use std::collections::{BTreeMap, BTreeSet};

use tracing::info;

use crate::core::membership::MembershipVec;
use crate::input::{InputBundle, SharedCacheData, normalize_symbol};
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

#[derive(Clone, Copy)]
struct ColumnView<'a> {
    row_idx: &'a [u32],
    values: &'a [u32],
}

impl<'a> ColumnView<'a> {
    fn iter(self, row_to_gene: &'a [u32]) -> impl Iterator<Item = (u32, u32)> + 'a {
        self.row_idx
            .iter()
            .zip(self.values.iter())
            .map(move |(row, val)| (row_to_gene[*row as usize], *val))
    }
}

enum MatrixSource<'a> {
    Matrix {
        col_ptr: &'a [u32],
        row_idx: &'a [u32],
        values: &'a [u32],
        n_cols: usize,
    },
    Cache(&'a SharedCacheData),
}

impl<'a> MatrixSource<'a> {
    fn from_input(input: &'a InputBundle) -> Self {
        if let Some(cache) = input.shared_cache.as_ref() {
            Self::Cache(cache)
        } else {
            Self::Matrix {
                col_ptr: &input.matrix.col_ptr,
                row_idx: &input.matrix.row_idx,
                values: &input.matrix.values,
                n_cols: input.matrix.n_cols as usize,
            }
        }
    }

    fn n_cols(&self) -> usize {
        match self {
            MatrixSource::Matrix { n_cols, .. } => *n_cols,
            MatrixSource::Cache(cache) => cache.n_cells as usize,
        }
    }

    fn column_view(&self, col: usize) -> ColumnView<'_> {
        match self {
            MatrixSource::Matrix {
                col_ptr,
                row_idx,
                values,
                ..
            } => {
                let start = col_ptr[col] as usize;
                let end = col_ptr[col + 1] as usize;
                ColumnView {
                    row_idx: &row_idx[start..end],
                    values: &values[start..end],
                }
            }
            MatrixSource::Cache(cache) => {
                let start = cache.col_ptr[col] as usize;
                let end = cache.col_ptr[col + 1] as usize;
                ColumnView {
                    row_idx: &cache.row_idx[start..end],
                    values: &cache.values_u32[start..end],
                }
            }
        }
    }
}

#[derive(Clone, Copy, Default)]
struct Acc {
    sum: f64,
    n: u32,
}

impl Acc {
    fn add(&mut self, x: f64) {
        self.sum += x;
        self.n += 1;
    }

    fn mean(&self) -> f64 {
        if self.n == 0 {
            f64::NAN
        } else {
            self.sum / self.n as f64
        }
    }
}

#[derive(Clone)]
struct GeneSets {
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
}

pub fn run_stage_translation_regime(
    input: &InputBundle,
    stage2: &Stage2Output,
) -> anyhow::Result<StageTranslationRegimeOutput> {
    run_stage_translation_regime_with_ln1p(input, stage2, simd::ln1p_f64)
}

fn run_stage_translation_regime_with_ln1p(
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

    let sets = resolve_gene_sets(&input.gene_index.map);

    let cov_housekeeping = coverage(
        sets.housekeeping.len(),
        parse_two_column_gene_set(HOUSEKEEPING_SET).len(),
    );
    let cov_global = coverage(
        sets.global_translation.len(),
        parse_two_column_gene_set(GLOBAL_TRANSLATION_SET).len(),
    );
    let cov_stress = coverage(
        sets.stress_responsive.len(),
        parse_two_column_gene_set(STRESS_RESPONSIVE_SET).len(),
    );
    let cov_isr = coverage(
        sets.isr_like.len(),
        parse_two_column_gene_set(ISR_LIKE_SET).len(),
    );
    let cov_codon = coverage(
        sets.codon_suboptimal.len(),
        parse_two_column_gene_set(CODON_SUBOPTIMAL_SET).len(),
    );

    let low_coverage = cov_housekeeping.low_confidence
        || cov_global.low_confidence
        || cov_stress.low_confidence
        || cov_isr.low_confidence
        || cov_codon.low_confidence;

    if low_coverage {
        info!(
            housekeeping_ratio = cov_housekeeping.ratio,
            global_ratio = cov_global.ratio,
            stress_ratio = cov_stress.ratio,
            isr_ratio = cov_isr.ratio,
            codon_ratio = cov_codon.ratio,
            "Translation extension is running with LOW_CONFIDENCE due to limited panel coverage"
        );
    }

    let ribosome_mem =
        MembershipVec::from_gene_ids(&sets.ribosome_core, input.gene_index.genes.len());
    let housekeeping_mem =
        MembershipVec::from_gene_ids(&sets.housekeeping, input.gene_index.genes.len());
    let global_mem =
        MembershipVec::from_gene_ids(&sets.global_translation, input.gene_index.genes.len());
    let stress_mem =
        MembershipVec::from_gene_ids(&sets.stress_responsive, input.gene_index.genes.len());
    let isr_mem = MembershipVec::from_gene_ids(&sets.isr_like, input.gene_index.genes.len());
    let codon_mem =
        MembershipVec::from_gene_ids(&sets.codon_suboptimal, input.gene_index.genes.len());
    let bg_mem = MembershipVec::from_gene_ids(&stage2.bg_gene_ids, input.gene_index.genes.len());
    let ext_ribosome_mem =
        MembershipVec::from_gene_ids(&sets.ext_ribosome_core, input.gene_index.genes.len());
    let ext_initiation_mem =
        MembershipVec::from_gene_ids(&sets.ext_initiation, input.gene_index.genes.len());
    let ext_bio_mem =
        MembershipVec::from_gene_ids(&sets.ext_biogenesis, input.gene_index.genes.len());
    let ext_mtor_mem = MembershipVec::from_gene_ids(&sets.ext_mtor, input.gene_index.genes.len());
    let ext_isr_mem = MembershipVec::from_gene_ids(&sets.ext_isr, input.gene_index.genes.len());
    let ext_proteo_mem =
        MembershipVec::from_gene_ids(&sets.ext_proteostasis, input.gene_index.genes.len());

    let mut cells = Vec::with_capacity(n_cells);
    let mut ext_panel_cores = Vec::with_capacity(n_cells);
    let mut regime_counts: BTreeMap<String, u64> = BTreeMap::new();
    let mut commitment_sum = 0.0;
    let mut high_selective_count = 0u64;
    let mut ribosome_values = Vec::with_capacity(sets.ext_ribosome_core.len());
    let mut initiation_values = Vec::with_capacity(sets.ext_initiation.len());
    let mut bio_values = Vec::with_capacity(sets.ext_biogenesis.len());
    let mut mtor_values = Vec::with_capacity(sets.ext_mtor.len());
    let mut isr_values = Vec::with_capacity(sets.ext_isr.len());
    let mut proteo_values = Vec::with_capacity(sets.ext_proteostasis.len());

    for col in 0..n_cells {
        let col_view = source.column_view(col);
        let ls = stage2.libsize[col];
        let denom = if ls == 0 { 1.0 } else { ls as f64 };

        let mut acc_ribosome = Acc::default();
        let mut acc_housekeeping = Acc::default();
        let mut acc_global = Acc::default();
        let mut acc_stress = Acc::default();
        let mut acc_isr = Acc::default();
        let mut acc_codon = Acc::default();
        let mut acc_bg = Acc::default();
        ribosome_values.clear();
        initiation_values.clear();
        bio_values.clear();
        mtor_values.clear();
        isr_values.clear();
        proteo_values.clear();

        for (gid, val) in col_view.iter(input.gene_index.row_to_gene.as_slice()) {
            let cpm = 1_000_000.0 * (val as f64) / denom;
            let x = ln1p_fn(cpm);

            if ribosome_mem.contains(gid) {
                acc_ribosome.add(x);
            }
            if housekeeping_mem.contains(gid) {
                acc_housekeeping.add(x);
            }
            if global_mem.contains(gid) {
                acc_global.add(x);
            }
            if stress_mem.contains(gid) {
                acc_stress.add(x);
            }
            if isr_mem.contains(gid) {
                acc_isr.add(x);
            }
            if codon_mem.contains(gid) {
                acc_codon.add(x);
            }
            if bg_mem.contains(gid) {
                acc_bg.add(x);
            }
            if ext_ribosome_mem.contains(gid) {
                ribosome_values.push(x);
            }
            if ext_initiation_mem.contains(gid) {
                initiation_values.push(x);
            }
            if ext_bio_mem.contains(gid) {
                bio_values.push(x);
            }
            if ext_mtor_mem.contains(gid) {
                mtor_values.push(x);
            }
            if ext_isr_mem.contains(gid) {
                isr_values.push(x);
            }
            if ext_proteo_mem.contains(gid) {
                proteo_values.push(x);
            }
        }

        let bg_mean = acc_bg.mean();
        let ribosome = enrich(acc_ribosome.mean(), bg_mean);
        let housekeeping = enrich(acc_housekeeping.mean(), bg_mean);
        let global = enrich(acc_global.mean(), bg_mean);
        let stress = enrich(acc_stress.mean(), bg_mean);
        let isr = enrich(acc_isr.mean(), bg_mean);
        let codon = enrich(acc_codon.mean(), bg_mean);

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

        if translation_selectivity_index >= HIGH_SELECTIVE_THRESHOLD {
            high_selective_count += 1;
        }
        commitment_sum += translation_commitment_score;

        *regime_counts
            .entry(translation_regime.to_string())
            .or_insert(0) += 1;

        ext_panel_cores.push(PanelCore {
            ribosome_core: trimmed_mean(&ribosome_values, MIN_GENES_PER_PANEL),
            initiation_core: trimmed_mean(&initiation_values, MIN_GENES_PER_PANEL),
            bio_core: trimmed_mean(&bio_values, MIN_GENES_PER_PANEL),
            mtor_core: trimmed_mean(&mtor_values, MIN_GENES_PER_PANEL),
            isr_core: trimmed_mean(&isr_values, MIN_GENES_PER_PANEL),
            proteo_core: trimmed_mean(&proteo_values, MIN_GENES_PER_PANEL),
        });

        cells.push(TranslationRegimeCell {
            cell_id: input.barcodes[col].clone(),
            ribosome_loading_heterogeneity: round6(ribosome_loading_heterogeneity),
            translation_selectivity_index: round6(translation_selectivity_index),
            isr_like_signature_score: round6(isr_like_signature_score),
            codon_bias_proxy: round6(codon_bias_proxy),
            translation_commitment_score: round6(translation_commitment_score),
            translation_regime,
            ribosome_core: f64::NAN,
            initiation_core: f64::NAN,
            bio_core: f64::NAN,
            mtor_core: f64::NAN,
            isr_core: f64::NAN,
            proteo_core: f64::NAN,
            tpi: f64::NAN,
            rbl: f64::NAN,
            mtor_p: f64::NAN,
            isr_a: f64::NAN,
            tpib: f64::NAN,
            tsm: f64::NAN,
            translation_high: false,
            biogenesis_high: false,
            isr_active: false,
            proteotoxic_risk: false,
            translational_stress_mode: false,
            missing_panel_gene_count: 0,
            low_confidence: low_coverage,
        });
    }

    let baseline = robust_baseline(&ext_panel_cores);
    for (cell, core) in cells.iter_mut().zip(ext_panel_cores.iter().copied()) {
        let ext: TranslationExtensionCellScores = build_scores(core, &baseline);
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

        cell.ribosome_core = round6_nan(ext.ribosome_core);
        cell.initiation_core = round6_nan(ext.initiation_core);
        cell.bio_core = round6_nan(ext.bio_core);
        cell.mtor_core = round6_nan(ext.mtor_core);
        cell.isr_core = round6_nan(ext.isr_core);
        cell.proteo_core = round6_nan(ext.proteo_core);
        cell.tpi = round6_nan(ext.tpi);
        cell.rbl = round6_nan(ext.rbl);
        cell.mtor_p = round6_nan(ext.mtor_p);
        cell.isr_a = round6_nan(ext.isr_a);
        cell.tpib = round6_nan(ext.tpib);
        cell.tsm = round6_nan(ext.tsm);
        cell.translation_high = ext.translation_high;
        cell.biogenesis_high = ext.biogenesis_high;
        cell.isr_active = ext.isr_active;
        cell.proteotoxic_risk = ext.proteotoxic_risk;
        cell.translational_stress_mode = ext.translational_stress_mode;
        cell.missing_panel_gene_count = missing_panel_gene_count;
    }

    let ext_snapshots = cells
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
        .collect::<Vec<_>>();

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

// Threshold-based classification for expression-only translation-regime proxies.
// Rule order is fixed and deterministic. Earlier rules have priority.
fn classify_translation_regime(
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

fn enrich(set_mean: f64, bg_mean: f64) -> f64 {
    if set_mean.is_nan() || bg_mean.is_nan() {
        return 0.0;
    }
    let raw = set_mean - bg_mean;
    let scaled = (raw + OFFSET) / SCALE;
    clamp01(scaled)
}

fn resolve_gene_sets(map: &BTreeMap<String, u32>) -> GeneSets {
    GeneSets {
        ribosome_core: resolve_gene_set_tsv(map, RIBOSOME_CORE_SET),
        housekeeping: resolve_gene_set_tsv(map, HOUSEKEEPING_SET),
        global_translation: resolve_gene_set_tsv(map, GLOBAL_TRANSLATION_SET),
        stress_responsive: resolve_gene_set_tsv(map, STRESS_RESPONSIVE_SET),
        isr_like: resolve_gene_set_tsv(map, ISR_LIKE_SET),
        codon_suboptimal: resolve_gene_set_tsv(map, CODON_SUBOPTIMAL_SET),
        ext_ribosome_core: resolve_gene_list(map, PANEL_RIBOSOME_CORE),
        ext_biogenesis: resolve_gene_list(map, PANEL_BIOGENESIS),
        ext_mtor: resolve_gene_list(map, PANEL_MTOR),
        ext_isr: resolve_gene_list(map, PANEL_ISR),
        ext_initiation: resolve_gene_list(map, PANEL_INITIATION),
        ext_proteostasis: resolve_gene_list(map, PANEL_PROTEOSTASIS),
    }
}

fn resolve_gene_set_tsv(map: &BTreeMap<String, u32>, tsv: &str) -> Vec<u32> {
    let mut ids = BTreeSet::new();
    for (human, mouse) in parse_two_column_gene_set(tsv) {
        for sym in [human, mouse] {
            let normalized = normalize_symbol(&sym);
            if normalized.is_empty() {
                continue;
            }
            if let Some(id) = map.get(&normalized) {
                ids.insert(*id);
            }
        }
    }
    ids.into_iter().collect()
}

fn parse_two_column_gene_set(tsv: &str) -> Vec<(String, String)> {
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
        let human = parts.next().unwrap_or("").trim().to_string();
        let mouse = parts.next().unwrap_or("").trim().to_string();
        if human.is_empty() && mouse.is_empty() {
            continue;
        }
        rows.push((human, mouse));
    }
    rows
}

fn resolve_gene_list(map: &BTreeMap<String, u32>, symbols: &[&str]) -> Vec<u32> {
    let mut ids = BTreeSet::new();
    for symbol in symbols {
        let normalized = normalize_symbol(symbol);
        if normalized.is_empty() {
            continue;
        }
        if let Some(id) = map.get(&normalized) {
            ids.insert(*id);
        }
    }
    ids.into_iter().collect()
}

fn round6(x: f64) -> f64 {
    (x * 1_000_000.0).round() / 1_000_000.0
}

fn round6_nan(x: f64) -> f64 {
    if x.is_nan() { f64::NAN } else { round6(x) }
}

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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::input::{CscMatrix, FeatureRow, InputFormat, build_gene_index};

    fn assert_f64_eq_or_nan(a: f64, b: f64) {
        if a.is_nan() && b.is_nan() {
            return;
        }
        assert_eq!(a, b);
    }

    fn synthetic_input() -> (InputBundle, Stage2Output) {
        let symbols = [
            "RPLP0", "RPS3", "EEF2", "ATF4", "DDIT3", "ASNS", "GAPDH", "ACTB", "SQSTM1", "XBP1",
        ];

        let features: Vec<FeatureRow> = symbols
            .iter()
            .enumerate()
            .map(|(i, s)| FeatureRow {
                raw_id: format!("G{}", i + 1),
                raw_name: s.to_string(),
                raw_type: String::new(),
                norm_symbol: s.to_string(),
            })
            .collect();

        let gene_index = build_gene_index(&features);

        let col_ptr = vec![0, 10, 17, 24];
        let row_idx = vec![
            0, 1, 2, 3, 4, 5, 6, 7, 8, 9, // cell A
            0, 1, 2, 6, 7, 8, 9, // cell B
            0, 1, 2, 3, 4, 5, 8, // cell C
        ];
        let values = vec![
            10, 8, 7, 2, 1, 1, 12, 11, 1, 2, // cell A
            4, 3, 2, 13, 12, 1, 1, // cell B
            9, 8, 7, 10, 8, 7, 5, // cell C
        ];

        let matrix = CscMatrix {
            n_rows: symbols.len() as u32,
            n_cols: 3,
            col_ptr,
            row_idx,
            values,
        };

        let input = InputBundle {
            format: InputFormat::TenXDir {
                matrix_path: "matrix.mtx".into(),
                barcodes_path: "barcodes.tsv".into(),
                features_path: "features.tsv".into(),
            },
            matrix,
            barcodes: vec!["A".to_string(), "B".to_string(), "C".to_string()],
            features,
            gene_index,
            metadata: None,
            shared_cache: None,
        };

        let stage2 = Stage2Output {
            bg_gene_ids: vec![6, 7, 8],
            libsize: vec![55, 36, 54],
            detected_genes: vec![10, 7, 7],
            components: Vec::new(),
            axes: Vec::new(),
        };

        (input, stage2)
    }

    #[test]
    fn deterministic_metrics() {
        let (input, stage2) = synthetic_input();
        let a = run_stage_translation_regime(&input, &stage2).unwrap();
        let b = run_stage_translation_regime(&input, &stage2).unwrap();

        assert_eq!(a.cells.len(), b.cells.len());
        for i in 0..a.cells.len() {
            assert_eq!(a.cells[i].cell_id, b.cells[i].cell_id);
            assert_f64_eq_or_nan(
                a.cells[i].ribosome_loading_heterogeneity,
                b.cells[i].ribosome_loading_heterogeneity,
            );
            assert_f64_eq_or_nan(
                a.cells[i].translation_selectivity_index,
                b.cells[i].translation_selectivity_index,
            );
            assert_f64_eq_or_nan(
                a.cells[i].isr_like_signature_score,
                b.cells[i].isr_like_signature_score,
            );
            assert_f64_eq_or_nan(a.cells[i].codon_bias_proxy, b.cells[i].codon_bias_proxy);
            assert_f64_eq_or_nan(
                a.cells[i].translation_commitment_score,
                b.cells[i].translation_commitment_score,
            );
            assert_eq!(a.cells[i].translation_regime, b.cells[i].translation_regime);
            assert_f64_eq_or_nan(a.cells[i].ribosome_core, b.cells[i].ribosome_core);
            assert_f64_eq_or_nan(a.cells[i].initiation_core, b.cells[i].initiation_core);
            assert_f64_eq_or_nan(a.cells[i].bio_core, b.cells[i].bio_core);
            assert_f64_eq_or_nan(a.cells[i].mtor_core, b.cells[i].mtor_core);
            assert_f64_eq_or_nan(a.cells[i].isr_core, b.cells[i].isr_core);
            assert_f64_eq_or_nan(a.cells[i].tpi, b.cells[i].tpi);
            assert_f64_eq_or_nan(a.cells[i].rbl, b.cells[i].rbl);
            assert_f64_eq_or_nan(a.cells[i].mtor_p, b.cells[i].mtor_p);
            assert_f64_eq_or_nan(a.cells[i].isr_a, b.cells[i].isr_a);
            assert_f64_eq_or_nan(a.cells[i].tpib, b.cells[i].tpib);
            assert_f64_eq_or_nan(a.cells[i].tsm, b.cells[i].tsm);
            assert_eq!(a.cells[i].translation_high, b.cells[i].translation_high);
            assert_eq!(a.cells[i].biogenesis_high, b.cells[i].biogenesis_high);
            assert_eq!(a.cells[i].isr_active, b.cells[i].isr_active);
            assert_eq!(a.cells[i].proteotoxic_risk, b.cells[i].proteotoxic_risk);
            assert_eq!(
                a.cells[i].translational_stress_mode,
                b.cells[i].translational_stress_mode
            );
            assert_eq!(
                a.cells[i].missing_panel_gene_count,
                b.cells[i].missing_panel_gene_count
            );
        }

        assert_eq!(a.regime_fractions, b.regime_fractions);
        assert_eq!(a.mean_translation_commitment, b.mean_translation_commitment);
        assert_eq!(
            a.high_selective_translation_fraction,
            b.high_selective_translation_fraction
        );
        assert_eq!(
            serde_json::to_string(&a.translation_extension_summary).unwrap(),
            serde_json::to_string(&b.translation_extension_summary).unwrap()
        );
    }

    #[test]
    fn regime_threshold_boundaries() {
        assert_eq!(
            classify_translation_regime(0.57, 0.67, 0.45, 0.49, 0.77),
            "BalancedTranslation"
        );
        assert_eq!(
            classify_translation_regime(0.58, 0.50, 0.42, 0.30, 0.60),
            "SelectiveTranslationStress"
        );
        assert_eq!(
            classify_translation_regime(0.50, 0.68, 0.20, 0.30, 0.60),
            "ISR_DominatedTranslation"
        );
        assert_eq!(
            classify_translation_regime(0.62, 0.58, 0.45, 0.50, 0.78),
            "TranslationFixation"
        );
    }

    #[test]
    fn simd_scalar_equivalence() {
        let (input, stage2) = synthetic_input();

        let simd_out =
            run_stage_translation_regime_with_ln1p(&input, &stage2, crate::simd::ln1p_f64).unwrap();
        let scalar_out =
            run_stage_translation_regime_with_ln1p(&input, &stage2, crate::simd::scalar::ln1p_f64)
                .unwrap();

        assert_eq!(simd_out.cells.len(), scalar_out.cells.len());
        for i in 0..simd_out.cells.len() {
            let a = &simd_out.cells[i];
            let b = &scalar_out.cells[i];
            assert_eq!(a.cell_id, b.cell_id);
            assert_f64_eq_or_nan(
                a.ribosome_loading_heterogeneity,
                b.ribosome_loading_heterogeneity,
            );
            assert_f64_eq_or_nan(
                a.translation_selectivity_index,
                b.translation_selectivity_index,
            );
            assert_f64_eq_or_nan(a.isr_like_signature_score, b.isr_like_signature_score);
            assert_f64_eq_or_nan(a.codon_bias_proxy, b.codon_bias_proxy);
            assert_f64_eq_or_nan(
                a.translation_commitment_score,
                b.translation_commitment_score,
            );
            assert_eq!(a.translation_regime, b.translation_regime);
            assert_f64_eq_or_nan(a.tpi, b.tpi);
            assert_f64_eq_or_nan(a.tsm, b.tsm);
        }
    }
}
