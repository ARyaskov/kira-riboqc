use serde::Serialize;
use std::collections::BTreeMap;

use crate::core::math::{percentile_sorted, round6};
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

    let mut tl_vals = Vec::with_capacity(rows.len());
    let mut rd_vals = Vec::with_capacity(rows.len());
    let mut sti_vals = Vec::with_capacity(rows.len());
    for r in rows {
        if r.translation_load.is_finite() {
            tl_vals.push(r.translation_load);
        }
        if r.ribosome_density.is_finite() {
            rd_vals.push(r.ribosome_density);
        }
        if r.stress_translation_index.is_finite() {
            sti_vals.push(r.stress_translation_index);
        }
    }

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
            translation_load: stat3(&mut tl_vals),
            ribosome_density: stat3(&mut rd_vals),
            stress_translation_index: stat3(&mut sti_vals),
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

fn stat3(values: &mut [f64]) -> Stat3 {
    if values.is_empty() {
        return Stat3 {
            median: 0.0,
            p90: 0.0,
            p99: 0.0,
        };
    }
    values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    Stat3 {
        median: round6(percentile_sorted(values, 0.50)),
        p90: round6(percentile_sorted(values, 0.90)),
        p99: round6(percentile_sorted(values, 0.99)),
    }
}
