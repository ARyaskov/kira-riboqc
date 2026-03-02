use serde::Serialize;
use std::collections::BTreeMap;

use crate::input::Stage1Stats;
use crate::metrics::translation_extension::aggregate::TranslationExtensionSummary;
use crate::pipeline::stage_translation_regime::StageTranslationRegimeOutput;
use crate::report::pipeline_contract::PipelineCellRow;
use crate::simd;

#[derive(Serialize)]
pub struct Summary {
    pub tool: ToolInfo,
    pub input: InputInfo,
    pub distributions: Distributions,
    pub regimes: Regimes,
    pub qc: QcSummary,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub translation: Option<TranslationSummary>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub translation_extension: Option<TranslationExtensionSummary>,
}

#[derive(Serialize)]
pub struct ToolInfo {
    pub name: &'static str,
    pub version: &'static str,
    pub simd: &'static str,
}

#[derive(Serialize)]
pub struct InputInfo {
    pub n_cells: u64,
    pub species: String,
}

#[derive(Serialize)]
pub struct Distributions {
    pub translation_load: Stat3,
    pub ribosome_density: Stat3,
    pub stress_translation_index: Stat3,
}

#[derive(Serialize)]
pub struct Stat3 {
    pub median: f64,
    pub p90: f64,
    pub p99: f64,
}

#[derive(Serialize)]
pub struct Regimes {
    pub counts: BTreeMap<String, u64>,
    pub fractions: BTreeMap<String, f64>,
}

#[derive(Serialize)]
pub struct QcSummary {
    pub low_confidence_fraction: f64,
    pub low_ribo_signal_fraction: f64,
}

#[derive(Serialize)]
pub struct TranslationSummary {
    pub regime_fractions: BTreeMap<String, f64>,
    pub mean_translation_commitment: f64,
    pub high_selective_translation_fraction: f64,
}

pub fn build_summary(
    _stage1: &Stage1Stats,
    rows: &[PipelineCellRow],
    stage_translation: Option<&StageTranslationRegimeOutput>,
) -> Summary {
    let n_cells = rows.len() as u64;
    let species = rows
        .first()
        .map(|r| r.species.clone())
        .unwrap_or_else(|| "unknown".to_string());

    let tl_vals: Vec<f64> = rows.iter().map(|r| r.translation_load).collect();
    let rd_vals: Vec<f64> = rows.iter().map(|r| r.ribosome_density).collect();
    let sti_vals: Vec<f64> = rows.iter().map(|r| r.stress_translation_index).collect();

    let mut counts = BTreeMap::new();
    for regime in [
        "HomeostaticTranslation",
        "GrowthDrivenTranslation",
        "StressAdaptiveTranslation",
        "TranslationalOverdrive",
        "TranslationalCollapse",
        "Unclassified",
    ] {
        counts.insert(regime.to_string(), 0);
    }
    for row in rows {
        if let Some(v) = counts.get_mut(row.regime) {
            *v += 1;
        }
    }

    let mut fractions = BTreeMap::new();
    for (k, v) in &counts {
        let frac = if n_cells == 0 {
            0.0
        } else {
            *v as f64 / n_cells as f64
        };
        fractions.insert(k.clone(), round6(frac));
    }

    let low_conf = rows.iter().filter(|r| r.has_low_confidence).count() as u64;
    let low_ribo = rows.iter().filter(|r| r.has_low_ribo_signal).count() as u64;

    Summary {
        tool: ToolInfo {
            name: "kira-riboqc",
            version: env!("CARGO_PKG_VERSION"),
            simd: simd::SIMD_KIND,
        },
        input: InputInfo { n_cells, species },
        distributions: Distributions {
            translation_load: stat3(&tl_vals),
            ribosome_density: stat3(&rd_vals),
            stress_translation_index: stat3(&sti_vals),
        },
        regimes: Regimes { counts, fractions },
        qc: QcSummary {
            low_confidence_fraction: round6(if n_cells == 0 {
                0.0
            } else {
                low_conf as f64 / n_cells as f64
            }),
            low_ribo_signal_fraction: round6(if n_cells == 0 {
                0.0
            } else {
                low_ribo as f64 / n_cells as f64
            }),
        },
        translation: stage_translation.map(|t| TranslationSummary {
            regime_fractions: t.regime_fractions.clone(),
            mean_translation_commitment: t.mean_translation_commitment,
            high_selective_translation_fraction: t.high_selective_translation_fraction,
        }),
        translation_extension: stage_translation.map(|t| t.translation_extension_summary.clone()),
    }
}

fn stat3(values: &[f64]) -> Stat3 {
    if values.is_empty() {
        return Stat3 {
            median: 0.0,
            p90: 0.0,
            p99: 0.0,
        };
    }

    let mut sorted = values
        .iter()
        .copied()
        .map(|v| if v.is_finite() { v } else { 0.0 })
        .collect::<Vec<_>>();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());

    Stat3 {
        median: round6(percentile_sorted(&sorted, 0.50)),
        p90: round6(percentile_sorted(&sorted, 0.90)),
        p99: round6(percentile_sorted(&sorted, 0.99)),
    }
}

fn percentile_sorted(sorted: &[f64], p: f64) -> f64 {
    let rank = (p * sorted.len() as f64).ceil() as usize;
    let idx = if rank == 0 { 0 } else { rank - 1 };
    sorted[idx.min(sorted.len() - 1)]
}

fn round6(x: f64) -> f64 {
    (x * 1_000_000.0).round() / 1_000_000.0
}
