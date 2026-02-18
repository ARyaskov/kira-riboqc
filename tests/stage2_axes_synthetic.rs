use kira_riboqc::input::{CscMatrix, FeatureRow, InputBundle, InputFormat, build_gene_index};
use kira_riboqc::pipeline::stage2_axes::run_stage2;

fn build_features(symbols: &[&str]) -> Vec<FeatureRow> {
    symbols
        .iter()
        .enumerate()
        .map(|(i, s)| FeatureRow {
            raw_id: format!("G{}", i + 1),
            raw_name: s.to_string(),
            raw_type: String::new(),
            norm_symbol: s.to_string(),
        })
        .collect()
}

#[test]
fn tiny_matrix_axes_ordering() {
    let symbols = ["RPLP0", "EEF2", "ATF4", "ZNF598", "HSPA8", "UBB", "SLC2A1"];
    let features = build_features(&symbols);
    let gene_index = build_gene_index(&features);

    // 4 cells, 7 genes.
    // Cell0: RPLP0 + EEF2 (high TL)
    // Cell1: ATF4 + ZNF598 + SLC2A1 (high ST + RQC)
    // Cell2: HSPA8 + UBB (positive TPC)
    // Cell3: empty
    let col_ptr = vec![0, 2, 5, 7, 7];
    let row_idx = vec![0, 1, 2, 3, 6, 4, 5];
    let values = vec![100, 50, 40, 30, 20, 60, 60];
    let matrix = CscMatrix {
        n_rows: 7,
        n_cols: 4,
        col_ptr,
        row_idx,
        values,
    };

    let barcodes = vec![
        "C0".to_string(),
        "C1".to_string(),
        "C2".to_string(),
        "C3".to_string(),
    ];

    let input = InputBundle {
        format: InputFormat::TenXDir {
            matrix_path: "matrix.mtx".into(),
            barcodes_path: "barcodes.tsv".into(),
            features_path: "features.tsv".into(),
        },
        matrix,
        barcodes,
        features,
        gene_index,
        metadata: None,
        shared_cache: None,
    };

    let out = run_stage2(&input).unwrap();
    assert_eq!(out.axes.len(), 4);

    let tl = out.axes.iter().map(|a| a.tl).collect::<Vec<_>>();
    let st = out.axes.iter().map(|a| a.st).collect::<Vec<_>>();
    let rqc = out.axes.iter().map(|a| a.rqc).collect::<Vec<_>>();
    let tpc = out.axes.iter().map(|a| a.tpc).collect::<Vec<_>>();

    assert!(tl[0] > tl[1] && tl[0] > tl[2] && tl[0] > tl[3]);
    assert!(st[1].is_finite());
    assert!(rqc[1].is_finite());
    assert!(st[0].is_nan() && st[2].is_nan() && st[3].is_nan());
    assert!(rqc[0].is_nan() && rqc[2].is_nan() && rqc[3].is_nan());
    assert!(tpc[2].is_finite());
    assert!(tpc[0].is_nan() && tpc[1].is_nan() && tpc[3].is_nan());

    for ax in &out.axes {
        assert!(ax.tl >= 0.0 && ax.tl <= 1.0 || ax.tl.is_nan());
        assert!(ax.st >= 0.0 && ax.st <= 1.0 || ax.st.is_nan());
        assert!(ax.rqc >= 0.0 && ax.rqc <= 1.0 || ax.rqc.is_nan());
        assert!(ax.rqc_pressure >= 0.0 && ax.rqc_pressure <= 1.0 || ax.rqc_pressure.is_nan());
        assert!(ax.tpc >= -1.0 && ax.tpc <= 1.0 || ax.tpc.is_nan());
    }

    assert!(out.axes.iter().all(|a| a.st_low_confidence));
}

#[test]
fn background_excludes_rpl_and_mt() {
    let symbols = ["RPLP0", "MT-CO1", "GAPDH"];
    let features = build_features(&symbols);
    let gene_index = build_gene_index(&features);

    let col_ptr = vec![0, 3];
    let row_idx = vec![0, 1, 2];
    let values = vec![10, 10, 10];
    let matrix = CscMatrix {
        n_rows: 3,
        n_cols: 1,
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
        barcodes: vec!["C0".to_string()],
        features,
        gene_index,
        metadata: None,
        shared_cache: None,
    };

    let out = run_stage2(&input).unwrap();
    assert!(out.bg_gene_ids.contains(&2));
    assert!(!out.bg_gene_ids.contains(&0));
    assert!(!out.bg_gene_ids.contains(&1));
}
