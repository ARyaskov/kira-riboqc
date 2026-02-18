use std::fs;
use std::path::Path;

use anyhow::Result;
use serde::Serialize;

#[derive(Serialize)]
struct PipelineStep {
    tool: Tool,
    artifacts: Artifacts,
    cell_metrics: CellMetrics,
    regimes: Vec<&'static str>,
}

#[derive(Serialize)]
struct Tool {
    name: &'static str,
    stage: &'static str,
    version: &'static str,
}

#[derive(Serialize)]
struct Artifacts {
    summary: &'static str,
    primary_metrics: &'static str,
    panels: &'static str,
    #[serde(skip_serializing_if = "Option::is_none")]
    translation_metrics: Option<&'static str>,
}

#[derive(Serialize)]
struct CellMetrics {
    file: &'static str,
    id_column: &'static str,
    regime_column: &'static str,
    confidence_column: &'static str,
    flag_column: &'static str,
    #[serde(skip_serializing_if = "Option::is_none")]
    translation_mode_column: Option<&'static str>,
}

pub fn write_pipeline_step(path: &Path, has_translation_extension: bool) -> Result<()> {
    let step = PipelineStep {
        tool: Tool {
            name: "kira-riboqc",
            stage: "translation",
            version: env!("CARGO_PKG_VERSION"),
        },
        artifacts: Artifacts {
            summary: "summary.json",
            primary_metrics: "riboqc.tsv",
            panels: "panels_report.tsv",
            translation_metrics: has_translation_extension
                .then_some("translation_regime_metrics.tsv"),
        },
        cell_metrics: CellMetrics {
            file: "riboqc.tsv",
            id_column: "barcode",
            regime_column: "regime",
            confidence_column: "confidence",
            flag_column: "flags",
            translation_mode_column: has_translation_extension.then_some("translation_regime"),
        },
        regimes: vec![
            "HomeostaticTranslation",
            "GrowthDrivenTranslation",
            "StressAdaptiveTranslation",
            "TranslationalOverdrive",
            "TranslationalCollapse",
            "Unclassified",
        ],
    };

    fs::write(path, serde_json::to_string_pretty(&step)?)?;
    Ok(())
}
