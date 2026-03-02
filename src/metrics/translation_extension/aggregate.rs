use std::collections::BTreeMap;

use serde::Serialize;

use crate::input::MetadataTable;

use super::panels::RIBO_EXTENSION_PANEL_V1;
use super::scores::{
    THRESHOLD_BIOGENESIS_HIGH, THRESHOLD_ISR_ACTIVE, THRESHOLD_PROTEOTOXIC_RISK,
    THRESHOLD_TRANSLATION_HIGH, THRESHOLD_TRANSLATIONAL_STRESS_MODE, median, percentile,
};

#[derive(Debug, Clone)]
pub struct TranslationExtensionCellSnapshot {
    pub cell_id: String,
    pub tpi: f64,
    pub rbl: f64,
    pub isr_a: f64,
    pub tpib: f64,
    pub tsm: f64,
    pub translation_high: bool,
    pub biogenesis_high: bool,
    pub isr_active: bool,
    pub proteotoxic_risk: bool,
    pub translational_stress_mode: bool,
    pub missing_ribosome_core: bool,
    pub missing_initiation_core: bool,
    pub missing_bio_core: bool,
    pub missing_mtor_core: bool,
    pub missing_isr_core: bool,
    pub missing_proteo_core: bool,
}

#[derive(Debug, Clone, Serialize)]
pub struct TranslationExtensionSummary {
    pub panel_version: &'static str,
    pub thresholds: TranslationExtensionThresholds,
    pub global_stats: TranslationGlobalStats,
    pub cluster_stats: BTreeMap<String, ClusterStats>,
    pub missingness: MissingnessStats,
}

#[derive(Debug, Clone, Serialize)]
pub struct TranslationExtensionThresholds {
    pub translation_high: f64,
    pub biogenesis_high: f64,
    pub isr_active: f64,
    pub proteotoxic_risk: f64,
    pub translational_stress_mode: f64,
}

#[derive(Debug, Clone, Serialize)]
pub struct CenterMad {
    pub median: f64,
    pub mad: f64,
}

#[derive(Debug, Clone, Serialize)]
pub struct TranslationGlobalStats {
    pub tpi: CenterMad,
    pub rbl: CenterMad,
    pub isr_a: CenterMad,
    pub tpib: CenterMad,
    pub tsm: CenterMad,
    pub top_clusters_by_translational_stress_mode: Vec<String>,
}

#[derive(Debug, Clone, Serialize)]
pub struct ClusterStats {
    pub n_cells: u64,
    pub tpi: Quantiles3,
    pub rbl: Quantiles3,
    pub isr_a: Quantiles3,
    pub tpib: Quantiles3,
    pub tsm: Quantiles3,
    pub fraction_translation_high: f64,
    pub fraction_biogenesis_high: f64,
    pub fraction_isr_active: f64,
    pub fraction_proteotoxic_risk: f64,
    pub fraction_translational_stress_mode: f64,
}

#[derive(Debug, Clone, Serialize)]
pub struct Quantiles3 {
    pub median: f64,
    pub p10: f64,
    pub p90: f64,
}

#[derive(Debug, Clone, Serialize)]
pub struct MissingnessStats {
    pub ribosome_core_fraction: f64,
    pub initiation_core_fraction: f64,
    pub bio_core_fraction: f64,
    pub mtor_core_fraction: f64,
    pub isr_core_fraction: f64,
    pub proteo_core_fraction: f64,
}

pub fn aggregate_translation_extension(
    cells: &[TranslationExtensionCellSnapshot],
    metadata: Option<&MetadataTable>,
) -> TranslationExtensionSummary {
    let tpi = collect_finite(cells.iter().map(|c| c.tpi));
    let rbl = collect_finite(cells.iter().map(|c| c.rbl));
    let isr_a = collect_finite(cells.iter().map(|c| c.isr_a));
    let tpib = collect_finite(cells.iter().map(|c| c.tpib));
    let tsm = collect_finite(cells.iter().map(|c| c.tsm));

    let cluster_stats = aggregate_cluster_stats(cells, metadata);
    let top_clusters_by_translational_stress_mode =
        top_clusters_by_stress_mode(&cluster_stats, 5usize);

    TranslationExtensionSummary {
        panel_version: RIBO_EXTENSION_PANEL_V1,
        thresholds: TranslationExtensionThresholds {
            translation_high: THRESHOLD_TRANSLATION_HIGH,
            biogenesis_high: THRESHOLD_BIOGENESIS_HIGH,
            isr_active: THRESHOLD_ISR_ACTIVE,
            proteotoxic_risk: THRESHOLD_PROTEOTOXIC_RISK,
            translational_stress_mode: THRESHOLD_TRANSLATIONAL_STRESS_MODE,
        },
        global_stats: TranslationGlobalStats {
            tpi: center_mad(&tpi),
            rbl: center_mad(&rbl),
            isr_a: center_mad(&isr_a),
            tpib: center_mad(&tpib),
            tsm: center_mad(&tsm),
            top_clusters_by_translational_stress_mode,
        },
        cluster_stats,
        missingness: missingness(cells),
    }
}

fn aggregate_cluster_stats(
    cells: &[TranslationExtensionCellSnapshot],
    metadata: Option<&MetadataTable>,
) -> BTreeMap<String, ClusterStats> {
    let mut by_cluster: BTreeMap<String, Vec<&TranslationExtensionCellSnapshot>> = BTreeMap::new();
    for cell in cells {
        if let Some(cluster) = infer_cluster(cell.cell_id.as_str(), metadata) {
            by_cluster.entry(cluster).or_default().push(cell);
        }
    }

    let mut out = BTreeMap::new();
    for (cluster, cluster_cells) in by_cluster {
        let n = cluster_cells.len() as u64;
        if n == 0 {
            continue;
        }
        let tpi = collect_finite(cluster_cells.iter().map(|c| c.tpi));
        let rbl = collect_finite(cluster_cells.iter().map(|c| c.rbl));
        let isr_a = collect_finite(cluster_cells.iter().map(|c| c.isr_a));
        let tpib = collect_finite(cluster_cells.iter().map(|c| c.tpib));
        let tsm = collect_finite(cluster_cells.iter().map(|c| c.tsm));

        out.insert(
            cluster,
            ClusterStats {
                n_cells: n,
                tpi: quantiles3(&tpi),
                rbl: quantiles3(&rbl),
                isr_a: quantiles3(&isr_a),
                tpib: quantiles3(&tpib),
                tsm: quantiles3(&tsm),
                fraction_translation_high: frac(
                    cluster_cells.iter().filter(|c| c.translation_high).count() as u64,
                    n,
                ),
                fraction_biogenesis_high: frac(
                    cluster_cells.iter().filter(|c| c.biogenesis_high).count() as u64,
                    n,
                ),
                fraction_isr_active: frac(
                    cluster_cells.iter().filter(|c| c.isr_active).count() as u64,
                    n,
                ),
                fraction_proteotoxic_risk: frac(
                    cluster_cells.iter().filter(|c| c.proteotoxic_risk).count() as u64,
                    n,
                ),
                fraction_translational_stress_mode: frac(
                    cluster_cells
                        .iter()
                        .filter(|c| c.translational_stress_mode)
                        .count() as u64,
                    n,
                ),
            },
        );
    }
    out
}

fn top_clusters_by_stress_mode(
    stats: &BTreeMap<String, ClusterStats>,
    top_n: usize,
) -> Vec<String> {
    let mut ranked: Vec<(&str, f64)> = stats
        .iter()
        .map(|(name, stat)| (name.as_str(), stat.tsm.median))
        .collect();
    ranked.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
    ranked
        .into_iter()
        .take(top_n)
        .map(|(name, _)| name.to_string())
        .collect()
}

fn infer_cluster(cell_id: &str, metadata: Option<&MetadataTable>) -> Option<String> {
    let meta = metadata.and_then(|m| m.rows.get(cell_id))?;
    for key in ["cluster", "seurat_clusters", "leiden", "louvain"] {
        for (field_key, value) in &meta.fields {
            if field_key.eq_ignore_ascii_case(key) {
                let trimmed = value.trim();
                if !trimmed.is_empty() {
                    return Some(trimmed.to_string());
                }
            }
        }
    }
    None
}

fn center_mad(values: &[f64]) -> CenterMad {
    let med = median(values);
    let mad = if med.is_nan() {
        f64::NAN
    } else {
        let deviations: Vec<f64> = values.iter().map(|v| (v - med).abs()).collect();
        median(&deviations)
    };
    CenterMad {
        median: round6_or_zero(med),
        mad: round6_or_zero(mad),
    }
}

fn quantiles3(values: &[f64]) -> Quantiles3 {
    Quantiles3 {
        median: round6_or_zero(percentile(values, 0.50)),
        p10: round6_or_zero(percentile(values, 0.10)),
        p90: round6_or_zero(percentile(values, 0.90)),
    }
}

fn missingness(cells: &[TranslationExtensionCellSnapshot]) -> MissingnessStats {
    let n = cells.len() as u64;
    MissingnessStats {
        ribosome_core_fraction: frac(
            cells.iter().filter(|c| c.missing_ribosome_core).count() as u64,
            n,
        ),
        initiation_core_fraction: frac(
            cells.iter().filter(|c| c.missing_initiation_core).count() as u64,
            n,
        ),
        bio_core_fraction: frac(
            cells.iter().filter(|c| c.missing_bio_core).count() as u64,
            n,
        ),
        mtor_core_fraction: frac(
            cells.iter().filter(|c| c.missing_mtor_core).count() as u64,
            n,
        ),
        isr_core_fraction: frac(
            cells.iter().filter(|c| c.missing_isr_core).count() as u64,
            n,
        ),
        proteo_core_fraction: frac(
            cells.iter().filter(|c| c.missing_proteo_core).count() as u64,
            n,
        ),
    }
}

fn frac(numerator: u64, denom: u64) -> f64 {
    if denom == 0 {
        0.0
    } else {
        round6_or_nan(numerator as f64 / denom as f64)
    }
}

fn collect_finite<I>(iter: I) -> Vec<f64>
where
    I: Iterator<Item = f64>,
{
    iter.filter(|v| v.is_finite()).collect()
}

fn round6_or_nan(value: f64) -> f64 {
    if value.is_nan() {
        f64::NAN
    } else {
        (value * 1_000_000.0).round() / 1_000_000.0
    }
}

fn round6_or_zero(value: f64) -> f64 {
    if value.is_nan() {
        0.0
    } else {
        (value * 1_000_000.0).round() / 1_000_000.0
    }
}
