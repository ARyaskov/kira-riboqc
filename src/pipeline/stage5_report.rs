use std::fs;
use std::path::Path;

use tracing::info;

use crate::cli::RunMode;
use crate::input::{InputBundle, Stage1Stats};
use crate::pipeline::stage_translation_regime::StageTranslationRegimeOutput;
use crate::pipeline::stage2_axes::Stage2Output;
use crate::pipeline::stage3_scores::Stage3Output;
use crate::pipeline::stage4_classify::Stage4Output;
use crate::report::clinical::build_report;
use crate::report::panels::write_panels_report;
use crate::report::pipeline_contract::build_pipeline_rows;
use crate::report::pipeline_step::write_pipeline_step;
use crate::report::summary::build_summary;
use crate::report::tsv::{write_riboqc_tsv, write_translation_regime_tsv};

pub fn run_stage5(
    out_dir: &Path,
    run_mode: RunMode,
    input: &InputBundle,
    stage1_stats: &Stage1Stats,
    stage2: &Stage2Output,
    stage3: &Stage3Output,
    stage4: &Stage4Output,
    stage_translation: Option<&StageTranslationRegimeOutput>,
) -> anyhow::Result<()> {
    let _span = tracing::info_span!("stage5_report").entered();
    info!("Stage 5 start");

    let target_dir = if matches!(run_mode, RunMode::Pipeline)
        && out_dir
            .file_name()
            .and_then(|v| v.to_str())
            .is_some_and(|v| v != "kira-riboqc")
    {
        out_dir.join("kira-riboqc")
    } else {
        out_dir.to_path_buf()
    };

    fs::create_dir_all(&target_dir)?;

    let rows = build_pipeline_rows(input, stage2, stage3, stage4, stage_translation);

    let tsv_path = target_dir.join("riboqc.tsv");
    write_riboqc_tsv(&tsv_path, &rows)?;

    let summary = build_summary(stage1_stats, &rows, stage_translation);
    let summary_path = target_dir.join("summary.json");
    let summary_json = serde_json::to_string_pretty(&summary)?;
    fs::write(&summary_path, summary_json)?;

    let panels_path = target_dir.join("panels_report.tsv");
    write_panels_report(&panels_path, input)?;

    let translation_metrics_path = if let Some(translation) = stage_translation {
        let path = target_dir.join("translation_regime_metrics.tsv");
        write_translation_regime_tsv(&path, translation)?;
        Some(path)
    } else {
        None
    };

    let pipeline_step_path = if matches!(run_mode, RunMode::Pipeline) {
        let path = target_dir.join("pipeline_step.json");
        write_pipeline_step(&path, stage_translation.is_some())?;
        Some(path)
    } else {
        None
    };

    if matches!(run_mode, RunMode::Standalone) {
        let report_text = build_report(&summary);
        let report_path = target_dir.join("report.txt");
        fs::write(&report_path, report_text)?;
    }

    info!(
        riboqc_tsv = %tsv_path.display(),
        summary_json = %summary_path.display(),
        panels_report = %panels_path.display(),
        translation_metrics = translation_metrics_path
            .as_ref()
            .map(|p| p.display().to_string())
            .unwrap_or_else(|| "n/a".to_string()),
        pipeline_step = pipeline_step_path
            .as_ref()
            .map(|p| p.display().to_string())
            .unwrap_or_else(|| "n/a".to_string()),
        "Stage 5 outputs written"
    );
    info!("Stage 5 end");

    Ok(())
}
