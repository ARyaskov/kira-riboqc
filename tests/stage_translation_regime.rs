use kira_riboqc::input::{CscMatrix, FeatureRow, InputBundle, InputFormat, build_gene_index};
use kira_riboqc::pipeline::stage_translation_regime::{
    classify_translation_regime, run_stage_translation_regime,
    run_stage_translation_regime_with_ln1p,
};
use kira_riboqc::pipeline::stage2_axes::Stage2Output;

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
        0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 1, 2, 6, 7, 8, 9, 0, 1, 2, 3, 4, 5, 8,
    ];
    let values = vec![
        10, 8, 7, 2, 1, 1, 12, 11, 1, 2, 4, 3, 2, 13, 12, 1, 1, 9, 8, 7, 10, 8, 7, 5,
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
        assert_f64_eq_or_nan(a.cells[i].tpi, b.cells[i].tpi);
        assert_f64_eq_or_nan(a.cells[i].tsm, b.cells[i].tsm);
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
        run_stage_translation_regime_with_ln1p(&input, &stage2, kira_riboqc::simd::ln1p_f64)
            .unwrap();
    let scalar_out = run_stage_translation_regime_with_ln1p(
        &input,
        &stage2,
        kira_riboqc::simd::scalar::ln1p_f64,
    )
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
        assert_eq!(a.translation_regime, b.translation_regime);
    }

    serde_json::to_string(&simd_out.translation_extension_summary).unwrap();
}

#[test]
fn requires_serde_runtime() {
    let json = serde_json::json!({});
    assert!(json.is_object());
}
