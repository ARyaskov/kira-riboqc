use std::collections::{BTreeMap, BTreeSet};

use tracing::info;

use crate::core::membership::MembershipVec;
use crate::input::{InputBundle, SharedCacheData, normalize_symbol};
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
    pub low_confidence: bool,
}

#[derive(Debug, Clone)]
pub struct StageTranslationRegimeOutput {
    pub cells: Vec<TranslationRegimeCell>,
    pub regime_fractions: BTreeMap<String, f64>,
    pub mean_translation_commitment: f64,
    pub high_selective_translation_fraction: f64,
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

    let cov_housekeeping = coverage(sets.housekeeping.len(), HOUSEKEEPING_SET.len());
    let cov_global = coverage(sets.global_translation.len(), GLOBAL_TRANSLATION_SET.len());
    let cov_stress = coverage(sets.stress_responsive.len(), STRESS_RESPONSIVE_SET.len());
    let cov_isr = coverage(sets.isr_like.len(), ISR_LIKE_SET.len());
    let cov_codon = coverage(sets.codon_suboptimal.len(), CODON_SUBOPTIMAL_SET.len());

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

    let mut cells = Vec::with_capacity(n_cells);
    let mut regime_counts: BTreeMap<String, u64> = BTreeMap::new();
    let mut commitment_sum = 0.0;
    let mut high_selective_count = 0u64;

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

        cells.push(TranslationRegimeCell {
            cell_id: input.barcodes[col].clone(),
            ribosome_loading_heterogeneity: round6(ribosome_loading_heterogeneity),
            translation_selectivity_index: round6(translation_selectivity_index),
            isr_like_signature_score: round6(isr_like_signature_score),
            codon_bias_proxy: round6(codon_bias_proxy),
            translation_commitment_score: round6(translation_commitment_score),
            translation_regime,
            low_confidence: low_coverage,
        });
    }

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

fn round6(x: f64) -> f64 {
    (x * 1_000_000.0).round() / 1_000_000.0
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
            assert_eq!(
                a.cells[i].ribosome_loading_heterogeneity,
                b.cells[i].ribosome_loading_heterogeneity
            );
            assert_eq!(
                a.cells[i].translation_selectivity_index,
                b.cells[i].translation_selectivity_index
            );
            assert_eq!(
                a.cells[i].isr_like_signature_score,
                b.cells[i].isr_like_signature_score
            );
            assert_eq!(a.cells[i].codon_bias_proxy, b.cells[i].codon_bias_proxy);
            assert_eq!(
                a.cells[i].translation_commitment_score,
                b.cells[i].translation_commitment_score
            );
            assert_eq!(a.cells[i].translation_regime, b.cells[i].translation_regime);
        }

        assert_eq!(a.regime_fractions, b.regime_fractions);
        assert_eq!(a.mean_translation_commitment, b.mean_translation_commitment);
        assert_eq!(
            a.high_selective_translation_fraction,
            b.high_selective_translation_fraction
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
            assert_eq!(
                a.ribosome_loading_heterogeneity,
                b.ribosome_loading_heterogeneity
            );
            assert_eq!(
                a.translation_selectivity_index,
                b.translation_selectivity_index
            );
            assert_eq!(a.isr_like_signature_score, b.isr_like_signature_score);
            assert_eq!(a.codon_bias_proxy, b.codon_bias_proxy);
            assert_eq!(
                a.translation_commitment_score,
                b.translation_commitment_score
            );
            assert_eq!(a.translation_regime, b.translation_regime);
        }
    }
}
