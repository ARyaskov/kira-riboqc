use std::collections::BTreeMap;
use std::sync::OnceLock;

use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use tracing::info;

use crate::core::{Acc, PanelBit, PanelMask, PanelMaskBuilder, round6, round6_nan};
use crate::input::{InputBundle, MatrixSource, normalize_symbol};
use crate::metrics::translation_extension::aggregate::{
    TranslationExtensionCellSnapshot, TranslationExtensionSummary, aggregate_translation_extension,
};
use crate::metrics::translation_extension::panels::{
    MIN_GENES_PER_PANEL, PANEL_BIOGENESIS, PANEL_INITIATION, PANEL_ISR, PANEL_MTOR,
    PANEL_PROTEOSTASIS, PANEL_RIBOSOME_CORE,
};
use crate::metrics::translation_extension::scores::{
    PanelCore, TranslationExtensionCellScores, build_scores, robust_baseline, trimmed_mean,
};
use crate::model::axes::{OFFSET, SCALE, clamp01};
use crate::pipeline::stage2_axes::Stage2Output;
use crate::simd;

const LOW_COVERAGE_MIN_RATIO: f64 = 0.35;
const HIGH_SELECTIVE_THRESHOLD: f64 = 0.65;
const PER_CELL_MISSING_PANEL_THRESHOLD: u8 = 3;

#[derive(Debug, Clone)]
pub struct TranslationRegimeCell {
    pub cell_id: String,
    pub ribosome_loading_heterogeneity: f64,
    pub translation_selectivity_index: f64,
    pub isr_like_signature_score: f64,
    pub codon_bias_proxy: f64,
    pub translation_commitment_score: f64,
    pub translation_regime: &'static str,
    pub ribosome_core: f64,
    pub initiation_core: f64,
    pub bio_core: f64,
    pub mtor_core: f64,
    pub isr_core: f64,
    pub proteo_core: f64,
    pub tpi: f64,
    pub rbl: f64,
    pub mtor_p: f64,
    pub isr_a: f64,
    pub tpib: f64,
    pub tsm: f64,
    pub translation_high: bool,
    pub biogenesis_high: bool,
    pub isr_active: bool,
    pub proteotoxic_risk: bool,
    pub translational_stress_mode: bool,
    pub missing_panel_gene_count: u8,
    pub low_confidence: bool,
}

#[derive(Debug, Clone)]
pub struct StageTranslationRegimeOutput {
    pub cells: Vec<TranslationRegimeCell>,
    pub regime_fractions: BTreeMap<String, f64>,
    pub mean_translation_commitment: f64,
    pub high_selective_translation_fraction: f64,
    pub translation_extension_summary: TranslationExtensionSummary,
}

#[derive(Clone, Copy)]
struct SetCoverage {
    ratio: f64,
    low_confidence: bool,
}

struct AxesBits {
    ribosome: PanelBit,
    housekeeping: PanelBit,
    global: PanelBit,
    stress: PanelBit,
    isr: PanelBit,
    codon: PanelBit,
    bg: PanelBit,
}

struct ExtBits {
    ribosome: PanelBit,
    initiation: PanelBit,
    bio: PanelBit,
    mtor: PanelBit,
    isr: PanelBit,
    proteo: PanelBit,
}

#[derive(Clone)]
struct CellMetrics {
    cell_id: String,
    ribosome_loading_heterogeneity: f64,
    translation_selectivity_index: f64,
    isr_like_signature_score: f64,
    codon_bias_proxy: f64,
    translation_commitment_score: f64,
    translation_regime: &'static str,
    ext_core: PanelCore,
}

pub fn run_stage_translation_regime(
    input: &InputBundle,
    stage2: &Stage2Output,
) -> anyhow::Result<StageTranslationRegimeOutput> {
    run_stage_translation_regime_with_ln1p(input, stage2, simd::ln1p_f64)
}

pub fn run_stage_translation_regime_with_ln1p(
    input: &InputBundle,
    stage2: &Stage2Output,
    ln1p_fn: fn(f64) -> f64,
) -> anyhow::Result<StageTranslationRegimeOutput> {
    let _span = tracing::info_span!("stage_translation_regime").entered();
    info!("Translation regime extension stage start");

    let source = MatrixSource::from_input(input);
    let n_cells = source.n_cols();
    if n_cells != input.barcodes.len() || n_cells != stage2.libsize.len() {
        return Err(anyhow::anyhow!(
            "cell count mismatch in translation extension stage"
        ));
    }

    let n_genes = input.gene_index.genes.len();
    let map = &input.gene_index.map;

    let sets = ResolvedSets::resolve(map);
    let cov = sets.coverage();

    let low_coverage = cov.any_low();
    if low_coverage {
        info!(
            housekeeping_ratio = cov.housekeeping.ratio,
            global_ratio = cov.global.ratio,
            stress_ratio = cov.stress.ratio,
            isr_ratio = cov.isr.ratio,
            codon_ratio = cov.codon.ratio,
            "Translation extension is running with LOW_CONFIDENCE due to limited panel coverage"
        );
    }

    let mut axes_builder = PanelMaskBuilder::new(n_genes);
    let axes_bits = AxesBits {
        ribosome: axes_builder.add_panel(&sets.ribosome_core),
        housekeeping: axes_builder.add_panel(&sets.housekeeping),
        global: axes_builder.add_panel(&sets.global_translation),
        stress: axes_builder.add_panel(&sets.stress_responsive),
        isr: axes_builder.add_panel(&sets.isr_like),
        codon: axes_builder.add_panel(&sets.codon_suboptimal),
        bg: axes_builder.add_panel(&stage2.bg_gene_ids),
    };
    let axes_mask = axes_builder.build();

    let mut ext_builder = PanelMaskBuilder::new(n_genes);
    let ext_bits = ExtBits {
        ribosome: ext_builder.add_panel(&sets.ext_ribosome_core),
        initiation: ext_builder.add_panel(&sets.ext_initiation),
        bio: ext_builder.add_panel(&sets.ext_biogenesis),
        mtor: ext_builder.add_panel(&sets.ext_mtor),
        isr: ext_builder.add_panel(&sets.ext_isr),
        proteo: ext_builder.add_panel(&sets.ext_proteostasis),
    };
    let ext_mask = ext_builder.build();
    let row_to_gene = input.gene_index.row_to_gene.as_slice();

    let cell_metrics: Vec<CellMetrics> = (0..n_cells)
        .into_par_iter()
        .map(|col| {
            compute_cell(
                &source,
                col,
                stage2.libsize[col],
                &input.barcodes[col],
                row_to_gene,
                &axes_mask,
                &axes_bits,
                &ext_mask,
                &ext_bits,
                ln1p_fn,
            )
        })
        .collect();

    let baseline = robust_baseline(&cell_metrics.iter().map(|c| c.ext_core).collect::<Vec<_>>());

    let mut cells: Vec<TranslationRegimeCell> = cell_metrics
        .par_iter()
        .map(|cm| {
            let ext: TranslationExtensionCellScores = build_scores(cm.ext_core, &baseline);
            let missing_panel_gene_count = [
                ext.missing_ribosome_core,
                ext.missing_initiation_core,
                ext.missing_bio_core,
                ext.missing_mtor_core,
                ext.missing_isr_core,
                ext.missing_proteo_core,
            ]
            .iter()
            .filter(|x| **x)
            .count() as u8;

            let low_conf =
                low_coverage || missing_panel_gene_count >= PER_CELL_MISSING_PANEL_THRESHOLD;

            TranslationRegimeCell {
                cell_id: cm.cell_id.clone(),
                ribosome_loading_heterogeneity: round6(cm.ribosome_loading_heterogeneity),
                translation_selectivity_index: round6(cm.translation_selectivity_index),
                isr_like_signature_score: round6(cm.isr_like_signature_score),
                codon_bias_proxy: round6(cm.codon_bias_proxy),
                translation_commitment_score: round6(cm.translation_commitment_score),
                translation_regime: cm.translation_regime,
                ribosome_core: round6_nan(ext.ribosome_core),
                initiation_core: round6_nan(ext.initiation_core),
                bio_core: round6_nan(ext.bio_core),
                mtor_core: round6_nan(ext.mtor_core),
                isr_core: round6_nan(ext.isr_core),
                proteo_core: round6_nan(ext.proteo_core),
                tpi: round6_nan(ext.tpi),
                rbl: round6_nan(ext.rbl),
                mtor_p: round6_nan(ext.mtor_p),
                isr_a: round6_nan(ext.isr_a),
                tpib: round6_nan(ext.tpib),
                tsm: round6_nan(ext.tsm),
                translation_high: ext.translation_high,
                biogenesis_high: ext.biogenesis_high,
                isr_active: ext.isr_active,
                proteotoxic_risk: ext.proteotoxic_risk,
                translational_stress_mode: ext.translational_stress_mode,
                missing_panel_gene_count,
                low_confidence: low_conf,
            }
        })
        .collect();

    let mut regime_counts: BTreeMap<&'static str, u64> = BTreeMap::new();
    let mut commitment_sum = 0.0;
    let mut high_selective_count = 0u64;
    for cm in &cell_metrics {
        *regime_counts.entry(cm.translation_regime).or_insert(0) += 1;
        commitment_sum += cm.translation_commitment_score;
        if cm.translation_selectivity_index >= HIGH_SELECTIVE_THRESHOLD {
            high_selective_count += 1;
        }
    }

    let ext_snapshots: Vec<TranslationExtensionCellSnapshot> = cells
        .iter()
        .map(|cell| TranslationExtensionCellSnapshot {
            cell_id: cell.cell_id.clone(),
            tpi: cell.tpi,
            rbl: cell.rbl,
            isr_a: cell.isr_a,
            tpib: cell.tpib,
            tsm: cell.tsm,
            translation_high: cell.translation_high,
            biogenesis_high: cell.biogenesis_high,
            isr_active: cell.isr_active,
            proteotoxic_risk: cell.proteotoxic_risk,
            translational_stress_mode: cell.translational_stress_mode,
            missing_ribosome_core: cell.ribosome_core.is_nan(),
            missing_initiation_core: cell.initiation_core.is_nan(),
            missing_bio_core: cell.bio_core.is_nan(),
            missing_mtor_core: cell.mtor_core.is_nan(),
            missing_isr_core: cell.isr_core.is_nan(),
            missing_proteo_core: cell.proteo_core.is_nan(),
        })
        .collect();

    let translation_extension_summary =
        aggregate_translation_extension(&ext_snapshots, input.metadata.as_ref());

    cells.sort_by(|a, b| a.cell_id.cmp(&b.cell_id));

    let mut regime_fractions = BTreeMap::new();
    for regime in [
        "BalancedTranslation",
        "SelectiveTranslationStress",
        "ISR_DominatedTranslation",
        "TranslationFixation",
    ] {
        let n = *regime_counts.get(regime).unwrap_or(&0);
        let frac = if n_cells == 0 {
            0.0
        } else {
            n as f64 / n_cells as f64
        };
        regime_fractions.insert(regime.to_string(), round6(frac));
    }

    let mean_translation_commitment = if n_cells == 0 {
        0.0
    } else {
        round6(commitment_sum / n_cells as f64)
    };

    let high_selective_translation_fraction = if n_cells == 0 {
        0.0
    } else {
        round6(high_selective_count as f64 / n_cells as f64)
    };

    info!(
        mean_translation_commitment,
        high_selective_translation_fraction, "Translation regime extension stage end"
    );

    Ok(StageTranslationRegimeOutput {
        cells,
        regime_fractions,
        mean_translation_commitment,
        high_selective_translation_fraction,
        translation_extension_summary,
    })
}

#[allow(clippy::too_many_arguments)]
fn compute_cell(
    source: &MatrixSource<'_>,
    col: usize,
    libsize: u64,
    cell_id: &str,
    row_to_gene: &[u32],
    axes_mask: &PanelMask,
    axes_bits: &AxesBits,
    ext_mask: &PanelMask,
    ext_bits: &ExtBits,
    ln1p_fn: fn(f64) -> f64,
) -> CellMetrics {
    let col_view = source.column(col);
    let denom = if libsize == 0 { 1.0 } else { libsize as f64 };

    let mut axes_accs = [Acc::default(); 7];
    let mut ribosome_values = Vec::new();
    let mut initiation_values = Vec::new();
    let mut bio_values = Vec::new();
    let mut mtor_values = Vec::new();
    let mut isr_values = Vec::new();
    let mut proteo_values = Vec::new();

    for (gid, val) in col_view.iter_genes(row_to_gene) {
        let axes_m = axes_mask.get(gid);
        let ext_m = ext_mask.get(gid);
        if axes_m == 0 && ext_m == 0 {
            continue;
        }
        let cpm = 1_000_000.0 * (val as f64) / denom;
        let x = ln1p_fn(cpm);

        let mut bits_left = axes_m;
        while bits_left != 0 {
            let lo = bits_left.trailing_zeros() as usize;
            axes_accs[lo].add(x);
            bits_left &= bits_left - 1;
        }

        if ext_m != 0 {
            if ext_m & ext_bits.ribosome.bit != 0 {
                ribosome_values.push(x);
            }
            if ext_m & ext_bits.initiation.bit != 0 {
                initiation_values.push(x);
            }
            if ext_m & ext_bits.bio.bit != 0 {
                bio_values.push(x);
            }
            if ext_m & ext_bits.mtor.bit != 0 {
                mtor_values.push(x);
            }
            if ext_m & ext_bits.isr.bit != 0 {
                isr_values.push(x);
            }
            if ext_m & ext_bits.proteo.bit != 0 {
                proteo_values.push(x);
            }
        }
    }

    let bg_mean = axes_accs[axes_bits.bg.index].mean();
    let ribosome = enrich_zero(axes_accs[axes_bits.ribosome.index].mean(), bg_mean);
    let housekeeping = enrich_zero(axes_accs[axes_bits.housekeeping.index].mean(), bg_mean);
    let global = enrich_zero(axes_accs[axes_bits.global.index].mean(), bg_mean);
    let stress = enrich_zero(axes_accs[axes_bits.stress.index].mean(), bg_mean);
    let isr = enrich_zero(axes_accs[axes_bits.isr.index].mean(), bg_mean);
    let codon = enrich_zero(axes_accs[axes_bits.codon.index].mean(), bg_mean);

    let ribosome_loading_heterogeneity =
        clamp01(0.6 * (ribosome - housekeeping).abs() + 0.4 * (global - housekeeping).abs());

    let translation_selectivity_index = {
        let numerator = 0.6 * stress + 0.4 * isr;
        let denominator = 0.5 * global + 0.5 * ribosome + 0.2;
        clamp01(numerator / denominator)
    };

    let isr_like_signature_score = clamp01(0.7 * isr + 0.3 * stress - 0.15 * global + 0.15);
    let codon_bias_proxy = clamp01(0.7 * codon + 0.3 * stress - 0.25 * housekeeping + 0.1);
    let translation_commitment_score = clamp01(
        0.40 * ribosome_loading_heterogeneity
            + 0.35 * translation_selectivity_index
            + 0.25 * isr_like_signature_score,
    );

    let translation_regime = classify_translation_regime(
        translation_selectivity_index,
        isr_like_signature_score,
        ribosome_loading_heterogeneity,
        codon_bias_proxy,
        translation_commitment_score,
    );

    let ext_core = PanelCore {
        ribosome_core: trimmed_mean(&ribosome_values, MIN_GENES_PER_PANEL),
        initiation_core: trimmed_mean(&initiation_values, MIN_GENES_PER_PANEL),
        bio_core: trimmed_mean(&bio_values, MIN_GENES_PER_PANEL),
        mtor_core: trimmed_mean(&mtor_values, MIN_GENES_PER_PANEL),
        isr_core: trimmed_mean(&isr_values, MIN_GENES_PER_PANEL),
        proteo_core: trimmed_mean(&proteo_values, MIN_GENES_PER_PANEL),
    };

    CellMetrics {
        cell_id: cell_id.to_string(),
        ribosome_loading_heterogeneity,
        translation_selectivity_index,
        isr_like_signature_score,
        codon_bias_proxy,
        translation_commitment_score,
        translation_regime,
        ext_core,
    }
}

pub fn classify_translation_regime(
    selectivity: f64,
    isr: f64,
    heterogeneity: f64,
    codon_bias: f64,
    commitment: f64,
) -> &'static str {
    if commitment >= 0.78 && selectivity >= 0.62 && isr >= 0.58 && codon_bias >= 0.50 {
        return "TranslationFixation";
    }
    if isr >= 0.68 && selectivity >= 0.45 {
        return "ISR_DominatedTranslation";
    }
    if selectivity >= 0.58 && heterogeneity >= 0.42 {
        return "SelectiveTranslationStress";
    }
    "BalancedTranslation"
}

fn coverage(found: usize, total: usize) -> SetCoverage {
    let ratio = if total == 0 {
        0.0
    } else {
        found as f64 / total as f64
    };
    SetCoverage {
        ratio,
        low_confidence: found == 0 || ratio < LOW_COVERAGE_MIN_RATIO,
    }
}

#[inline]
fn enrich_zero(set_mean: f64, bg_mean: f64) -> f64 {
    if set_mean.is_nan() || bg_mean.is_nan() {
        return 0.0;
    }
    let raw = set_mean - bg_mean;
    let scaled = (raw + OFFSET) / SCALE;
    clamp01(scaled)
}

struct ResolvedSets {
    ribosome_core: Vec<u32>,
    housekeeping: Vec<u32>,
    global_translation: Vec<u32>,
    stress_responsive: Vec<u32>,
    isr_like: Vec<u32>,
    codon_suboptimal: Vec<u32>,
    ext_ribosome_core: Vec<u32>,
    ext_biogenesis: Vec<u32>,
    ext_mtor: Vec<u32>,
    ext_isr: Vec<u32>,
    ext_initiation: Vec<u32>,
    ext_proteostasis: Vec<u32>,
    sizes: SetSizes,
}

struct SetSizes {
    housekeeping_total: usize,
    global_total: usize,
    stress_total: usize,
    isr_total: usize,
    codon_total: usize,
}

struct Coverage {
    housekeeping: SetCoverage,
    global: SetCoverage,
    stress: SetCoverage,
    isr: SetCoverage,
    codon: SetCoverage,
}

impl Coverage {
    fn any_low(&self) -> bool {
        self.housekeeping.low_confidence
            || self.global.low_confidence
            || self.stress.low_confidence
            || self.isr.low_confidence
            || self.codon.low_confidence
    }
}

impl ResolvedSets {
    fn resolve(map: &FxHashMap<String, u32>) -> Self {
        let hk = housekeeping_pairs();
        let gl = global_translation_pairs();
        let sr = stress_responsive_pairs();
        let is = isr_like_pairs();
        let cd = codon_suboptimal_pairs();
        let rb = ribosome_core_pairs();

        let sizes = SetSizes {
            housekeeping_total: hk.len(),
            global_total: gl.len(),
            stress_total: sr.len(),
            isr_total: is.len(),
            codon_total: cd.len(),
        };

        Self {
            ribosome_core: resolve_from_pairs(map, rb),
            housekeeping: resolve_from_pairs(map, hk),
            global_translation: resolve_from_pairs(map, gl),
            stress_responsive: resolve_from_pairs(map, sr),
            isr_like: resolve_from_pairs(map, is),
            codon_suboptimal: resolve_from_pairs(map, cd),
            ext_ribosome_core: resolve_from_list(map, PANEL_RIBOSOME_CORE),
            ext_biogenesis: resolve_from_list(map, PANEL_BIOGENESIS),
            ext_mtor: resolve_from_list(map, PANEL_MTOR),
            ext_isr: resolve_from_list(map, PANEL_ISR),
            ext_initiation: resolve_from_list(map, PANEL_INITIATION),
            ext_proteostasis: resolve_from_list(map, PANEL_PROTEOSTASIS),
            sizes,
        }
    }

    fn coverage(&self) -> Coverage {
        Coverage {
            housekeeping: coverage(self.housekeeping.len(), self.sizes.housekeeping_total),
            global: coverage(self.global_translation.len(), self.sizes.global_total),
            stress: coverage(self.stress_responsive.len(), self.sizes.stress_total),
            isr: coverage(self.isr_like.len(), self.sizes.isr_total),
            codon: coverage(self.codon_suboptimal.len(), self.sizes.codon_total),
        }
    }
}

fn resolve_from_pairs(
    map: &FxHashMap<String, u32>,
    pairs: &'static [(&'static str, &'static str)],
) -> Vec<u32> {
    let mut ids = FxHashSet::default();
    for (human, mouse) in pairs {
        for sym in [human, mouse] {
            let normalized = normalize_symbol(sym);
            if normalized.is_empty() {
                continue;
            }
            if let Some(id) = map.get(&normalized) {
                ids.insert(*id);
            }
        }
    }
    let mut out: Vec<u32> = ids.into_iter().collect();
    out.sort_unstable();
    out
}

fn resolve_from_list(map: &FxHashMap<String, u32>, symbols: &[&str]) -> Vec<u32> {
    let mut ids = FxHashSet::default();
    for symbol in symbols {
        let normalized = normalize_symbol(symbol);
        if normalized.is_empty() {
            continue;
        }
        if let Some(id) = map.get(&normalized) {
            ids.insert(*id);
        }
    }
    let mut out: Vec<u32> = ids.into_iter().collect();
    out.sort_unstable();
    out
}

fn parse_two_column_gene_set(tsv: &'static str) -> Vec<(&'static str, &'static str)> {
    let mut rows = Vec::new();
    for (idx, line) in tsv.lines().enumerate() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        if idx == 0 && trimmed.to_ascii_lowercase().contains("human_symbol") {
            continue;
        }
        let mut parts = trimmed.split('\t');
        let human = parts.next().unwrap_or("").trim();
        let mouse = parts.next().unwrap_or("").trim();
        if human.is_empty() && mouse.is_empty() {
            continue;
        }
        rows.push((human, mouse));
    }
    rows
}

macro_rules! cached_pairs {
    ($fn_name:ident, $slot:ident, $tsv:expr) => {
        fn $fn_name() -> &'static [(&'static str, &'static str)] {
            static $slot: OnceLock<Vec<(&'static str, &'static str)>> = OnceLock::new();
            $slot
                .get_or_init(|| parse_two_column_gene_set($tsv))
                .as_slice()
        }
    };
}

cached_pairs!(housekeeping_pairs, HOUSEKEEPING_PAIRS, HOUSEKEEPING_SET);
cached_pairs!(
    global_translation_pairs,
    GLOBAL_TRANSLATION_PAIRS,
    GLOBAL_TRANSLATION_SET
);
cached_pairs!(
    stress_responsive_pairs,
    STRESS_RESPONSIVE_PAIRS,
    STRESS_RESPONSIVE_SET
);
cached_pairs!(isr_like_pairs, ISR_LIKE_PAIRS, ISR_LIKE_SET);
cached_pairs!(
    codon_suboptimal_pairs,
    CODON_SUBOPTIMAL_PAIRS,
    CODON_SUBOPTIMAL_SET
);
cached_pairs!(ribosome_core_pairs, RIBOSOME_CORE_PAIRS, RIBOSOME_CORE_SET);

const HOUSEKEEPING_SET: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/resources/ribosome/translation/housekeeping.tsv"
));
const GLOBAL_TRANSLATION_SET: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/resources/ribosome/translation/global_translation.tsv"
));
const STRESS_RESPONSIVE_SET: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/resources/ribosome/translation/stress_responsive.tsv"
));
const ISR_LIKE_SET: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/resources/ribosome/translation/isr_like.tsv"
));
const CODON_SUBOPTIMAL_SET: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/resources/ribosome/translation/codon_suboptimal.tsv"
));
const RIBOSOME_CORE_SET: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/resources/ribosome/translation/ribosome_core.tsv"
));
