use std::fs::{create_dir_all, read, read_to_string, write};

use kira_riboqc::cli::{Mode, RunArgs, RunMode};
use kira_riboqc::pipeline::stage_translation_regime::run_stage_translation_regime;
use kira_riboqc::pipeline::stage1_load::run_stage1;
use kira_riboqc::pipeline::stage2_axes::run_stage2;
use kira_riboqc::pipeline::stage3_scores::run_stage3;
use kira_riboqc::pipeline::stage4_classify::run_stage4;
use kira_riboqc::pipeline::stage5_report::run_stage5;
use kira_riboqc::report::tsv::PIPELINE_TSV_HEADER;
use serde_json::Value;
use tempfile::tempdir;

#[test]
fn stage5_pipeline_outputs_written() {
    let tmp = tempdir().unwrap();
    let input_dir = tmp.path().join("input");
    let out_dir = tmp.path().join("out");
    create_dir_all(&input_dir).unwrap();

    let matrix = "%%MatrixMarket matrix coordinate integer general\n3 2 2\n1 1 5\n2 2 3\n";
    let barcodes = "CELL1\nCELL2\n";
    let features = "G1\tRPLP0\nG2\tEEF2\nG3\tATF4\n";

    write(input_dir.join("matrix.mtx"), matrix).unwrap();
    write(input_dir.join("barcodes.tsv"), barcodes).unwrap();
    write(input_dir.join("features.tsv"), features).unwrap();

    let args = RunArgs {
        input: input_dir,
        out: out_dir.clone(),
        mode: Mode::Cell,
        metadata: None,
        run_mode: RunMode::Pipeline,
        disable_translation_extension: false,
    };

    let (input, stage1_stats) = run_stage1(&args).unwrap();
    let stage2 = run_stage2(&input).unwrap();
    let stage3 = run_stage3(&stage2).unwrap();
    let stage4 = run_stage4(&stage2, &stage3).unwrap();
    let stage_translation = run_stage_translation_regime(&input, &stage2).unwrap();
    run_stage5(
        &out_dir,
        RunMode::Pipeline,
        &input,
        &stage1_stats,
        &stage2,
        &stage3,
        &stage4,
        Some(&stage_translation),
    )
    .unwrap();

    let final_dir = out_dir.join("kira-riboqc");
    let tsv_path = final_dir.join("riboqc.tsv");
    let json_path = final_dir.join("summary.json");
    let panels_path = final_dir.join("panels_report.tsv");
    let step_path = final_dir.join("pipeline_step.json");
    let translation_metrics_path = final_dir.join("translation_regime_metrics.tsv");

    assert!(tsv_path.exists());
    assert!(json_path.exists());
    assert!(panels_path.exists());
    assert!(step_path.exists());
    assert!(translation_metrics_path.exists());

    let tsv = read_to_string(tsv_path).unwrap();
    let lines: Vec<&str> = tsv.lines().collect();
    assert_eq!(lines[0], PIPELINE_TSV_HEADER);
    assert_eq!(lines.len(), 3);

    let summary: Value = serde_json::from_str(&read_to_string(json_path).unwrap()).unwrap();
    assert!(summary.get("tool").is_some());
    assert!(summary.get("input").is_some());
    assert!(summary.get("distributions").is_some());
    assert!(summary.get("regimes").is_some());
    assert!(summary.get("qc").is_some());
    assert!(summary.get("translation").is_some());

    let step: Value = serde_json::from_str(&read_to_string(step_path).unwrap()).unwrap();
    assert_eq!(step["tool"]["name"], "kira-riboqc");
    assert_eq!(step["tool"]["stage"], "translation");
    assert_eq!(step["artifacts"]["primary_metrics"], "riboqc.tsv");
    assert_eq!(step["cell_metrics"]["file"], "riboqc.tsv");
    assert_eq!(step["cell_metrics"]["id_column"], "barcode");
    assert_eq!(step["cell_metrics"]["regime_column"], "regime");
    assert_eq!(step["cell_metrics"]["confidence_column"], "confidence");
    assert_eq!(step["cell_metrics"]["flag_column"], "flags");
    assert_eq!(
        step["artifacts"]["translation_metrics"],
        "translation_regime_metrics.tsv"
    );
    assert_eq!(
        step["cell_metrics"]["translation_mode_column"],
        "translation_regime"
    );

    let translation_metrics = read_to_string(translation_metrics_path).unwrap();
    assert!(
        translation_metrics
            .lines()
            .next()
            .unwrap()
            .starts_with("cell_id\tribosome_loading_heterogeneity")
    );

    let panels = read_to_string(panels_path).unwrap();
    assert!(
        panels
            .lines()
            .next()
            .unwrap()
            .starts_with("panel_id\tpanel_name\tpanel_group")
    );
}

#[test]
fn stage5_deterministic_outputs() {
    let tmp = tempdir().unwrap();
    let input_dir = tmp.path().join("input");
    let out_a = tmp.path().join("out_a");
    let out_b = tmp.path().join("out_b");
    create_dir_all(&input_dir).unwrap();

    write(
        input_dir.join("matrix.mtx"),
        "%%MatrixMarket matrix coordinate integer general\n3 2 3\n1 1 5\n2 2 3\n3 2 1\n",
    )
    .unwrap();
    write(input_dir.join("barcodes.tsv"), "CELL1\nCELL2\n").unwrap();
    write(
        input_dir.join("features.tsv"),
        "G1\tRPLP0\nG2\tEEF2\nG3\tATF4\n",
    )
    .unwrap();

    let run_once = |out_dir: &std::path::Path| {
        let args = RunArgs {
            input: input_dir.clone(),
            out: out_dir.to_path_buf(),
            mode: Mode::Cell,
            metadata: None,
            run_mode: RunMode::Pipeline,
            disable_translation_extension: false,
        };

        let (input, stage1_stats) = run_stage1(&args).unwrap();
        let stage2 = run_stage2(&input).unwrap();
        let stage3 = run_stage3(&stage2).unwrap();
        let stage4 = run_stage4(&stage2, &stage3).unwrap();
        let stage_translation = run_stage_translation_regime(&input, &stage2).unwrap();
        run_stage5(
            out_dir,
            RunMode::Pipeline,
            &input,
            &stage1_stats,
            &stage2,
            &stage3,
            &stage4,
            Some(&stage_translation),
        )
        .unwrap();
    };

    run_once(&out_a);
    run_once(&out_b);

    for name in [
        "riboqc.tsv",
        "summary.json",
        "panels_report.tsv",
        "pipeline_step.json",
        "translation_regime_metrics.tsv",
    ] {
        let a = read(out_a.join("kira-riboqc").join(name)).unwrap();
        let b = read(out_b.join("kira-riboqc").join(name)).unwrap();
        assert_eq!(a, b, "mismatch in {name}");
    }
}

#[test]
fn stage5_standalone_does_not_write_pipeline_step() {
    let tmp = tempdir().unwrap();
    let input_dir = tmp.path().join("input");
    let out_dir = tmp.path().join("out");
    create_dir_all(&input_dir).unwrap();

    let matrix = "%%MatrixMarket matrix coordinate integer general\n3 2 2\n1 1 5\n2 2 3\n";
    let barcodes = "CELL1\nCELL2\n";
    let features = "G1\tRPLP0\nG2\tEEF2\nG3\tATF4\n";

    write(input_dir.join("matrix.mtx"), matrix).unwrap();
    write(input_dir.join("barcodes.tsv"), barcodes).unwrap();
    write(input_dir.join("features.tsv"), features).unwrap();

    let args = RunArgs {
        input: input_dir,
        out: out_dir.clone(),
        mode: Mode::Cell,
        metadata: None,
        run_mode: RunMode::Standalone,
        disable_translation_extension: false,
    };

    let (input, stage1_stats) = run_stage1(&args).unwrap();
    let stage2 = run_stage2(&input).unwrap();
    let stage3 = run_stage3(&stage2).unwrap();
    let stage4 = run_stage4(&stage2, &stage3).unwrap();
    let stage_translation = run_stage_translation_regime(&input, &stage2).unwrap();
    run_stage5(
        &out_dir,
        RunMode::Standalone,
        &input,
        &stage1_stats,
        &stage2,
        &stage3,
        &stage4,
        Some(&stage_translation),
    )
    .unwrap();

    assert!(out_dir.join("riboqc.tsv").exists());
    assert!(out_dir.join("summary.json").exists());
    assert!(out_dir.join("panels_report.tsv").exists());
    assert!(out_dir.join("translation_regime_metrics.tsv").exists());
    assert!(out_dir.join("report.txt").exists());
    assert!(!out_dir.join("pipeline_step.json").exists());
}

#[test]
fn stage5_extension_disabled_keeps_legacy_artifacts_stable() {
    let tmp = tempdir().unwrap();
    let input_dir = tmp.path().join("input");
    let out_dir = tmp.path().join("out");
    create_dir_all(&input_dir).unwrap();

    write(
        input_dir.join("matrix.mtx"),
        "%%MatrixMarket matrix coordinate integer general\n3 2 2\n1 1 5\n2 2 3\n",
    )
    .unwrap();
    write(input_dir.join("barcodes.tsv"), "CELL1\nCELL2\n").unwrap();
    write(
        input_dir.join("features.tsv"),
        "G1\tRPLP0\nG2\tEEF2\nG3\tATF4\n",
    )
    .unwrap();

    let args = RunArgs {
        input: input_dir,
        out: out_dir.clone(),
        mode: Mode::Cell,
        metadata: None,
        run_mode: RunMode::Pipeline,
        disable_translation_extension: true,
    };

    let (input, stage1_stats) = run_stage1(&args).unwrap();
    let stage2 = run_stage2(&input).unwrap();
    let stage3 = run_stage3(&stage2).unwrap();
    let stage4 = run_stage4(&stage2, &stage3).unwrap();
    run_stage5(
        &out_dir,
        RunMode::Pipeline,
        &input,
        &stage1_stats,
        &stage2,
        &stage3,
        &stage4,
        None,
    )
    .unwrap();

    let final_dir = out_dir.join("kira-riboqc");
    assert!(final_dir.join("riboqc.tsv").exists());
    assert!(final_dir.join("summary.json").exists());
    assert!(final_dir.join("panels_report.tsv").exists());
    assert!(!final_dir.join("translation_regime_metrics.tsv").exists());

    let summary: Value =
        serde_json::from_str(&read_to_string(final_dir.join("summary.json")).unwrap()).unwrap();
    assert!(summary.get("translation").is_none());
}
