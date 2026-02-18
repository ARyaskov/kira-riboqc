use kira_riboqc::core::membership::MembershipVec;
use kira_riboqc::input::{CscMatrix, FeatureRow, InputBundle, InputFormat, build_gene_index};
use kira_riboqc::pipeline::stage2_axes::run_stage2;
use kira_riboqc::simd::ln1p_f64;

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
fn membership_contains_and_bounds() {
    let mem = MembershipVec::from_gene_ids(&[0, 2, 5], 4);
    assert!(mem.contains(0));
    assert!(!mem.contains(1));
    assert!(mem.contains(2));
    assert!(!mem.contains(3));
    assert!(!mem.contains(5));
}

#[test]
fn stage2_axes_equivalence_reference() {
    let symbols = ["RPLP0", "EEF2", "ATF4", "SLC2A1"];
    let features = build_features(&symbols);
    let gene_index = build_gene_index(&features);

    let col_ptr = vec![0, 4, 7];
    let row_idx = vec![0, 1, 2, 3, 1, 2, 3];
    let values = vec![10, 5, 2, 4, 1, 9, 3];
    let matrix = CscMatrix {
        n_rows: 4,
        n_cols: 2,
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
        barcodes: vec!["C0".to_string(), "C1".to_string()],
        features,
        gene_index,
        metadata: None,
        shared_cache: None,
    };

    let out = run_stage2(&input).unwrap();

    let expected_cell0 = expected_components_for_cell(&[10.0, 5.0, 2.0, 4.0], 21.0);
    let expected_cell1 = expected_components_for_cell(&[0.0, 1.0, 9.0, 3.0], 13.0);

    assert_close(out.components[0].rp, expected_cell0.0);
    assert_close(out.components[0].mach, expected_cell0.1);
    assert_close(out.components[0].stress, expected_cell0.2);
    assert_close(out.components[0].metab, expected_cell0.3);
    assert_close(out.axes[0].tl, expected_cell0.4);
    assert_close(out.axes[0].st, expected_cell0.5);

    assert_close(out.components[1].rp, expected_cell1.0);
    assert_close(out.components[1].mach, expected_cell1.1);
    assert_close(out.components[1].stress, expected_cell1.2);
    assert_close(out.components[1].metab, expected_cell1.3);
    assert_close(out.axes[1].tl, expected_cell1.4);
    assert_close(out.axes[1].st, expected_cell1.5);
}

fn expected_components_for_cell(counts: &[f64; 4], libsize: f64) -> (f64, f64, f64, f64, f64, f64) {
    let denom = if libsize == 0.0 { 1.0 } else { libsize };
    let x: Vec<f64> = counts
        .iter()
        .map(|c| ln1p_f64(1_000_000.0 * *c / denom))
        .collect();

    let rp_mean = mean_for_present(counts[0], x[0]);
    let mach_mean = mean_for_present(counts[1], x[1]);
    let stress_mean = mean_for_present(counts[2], x[2]);
    let metab_mean = mean_for_present(counts[3], x[3]);
    let bg_mean = mean_for_bg(&[(counts[1], x[1]), (counts[2], x[2]), (counts[3], x[3])]);

    let rp = enrich(rp_mean, bg_mean);
    let mach = enrich(mach_mean, bg_mean);
    let stress = enrich(stress_mean, bg_mean);
    let metab = enrich(metab_mean, bg_mean);

    let tl = clamp01(0.5 * nan_to_zero(rp) + 0.3 * nan_to_zero(mach) + 0.2 * nan_to_zero(f64::NAN));
    let st = clamp01(0.55 * stress + 0.45 * metab);

    (rp, mach, stress, metab, tl, st)
}

fn mean_for_present(count: f64, x: f64) -> f64 {
    if count > 0.0 { x } else { f64::NAN }
}

fn mean_for_bg(entries: &[(f64, f64)]) -> f64 {
    let mut sum = 0.0;
    let mut n = 0usize;
    for (count, x) in entries {
        if *count > 0.0 {
            sum += *x;
            n += 1;
        }
    }
    if n == 0 { f64::NAN } else { sum / n as f64 }
}

fn enrich(set_mean: f64, bg_mean: f64) -> f64 {
    let raw = set_mean - bg_mean;
    let scaled = (raw + 0.5) / 1.5;
    clamp01(scaled)
}

fn clamp01(x: f64) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    if x < 0.0 {
        0.0
    } else if x > 1.0 {
        1.0
    } else {
        x
    }
}

fn nan_to_zero(x: f64) -> f64 {
    if x.is_nan() { 0.0 } else { x }
}

fn assert_close(a: f64, b: f64) {
    let eps = 1e-12;
    if a.is_nan() && b.is_nan() {
        return;
    }
    assert!((a - b).abs() < eps, "{a} vs {b}");
}
