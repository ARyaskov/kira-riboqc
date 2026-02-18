use kira_riboqc::pipeline::stage2_axes::{CellAxes, CellAxisComponents, Stage2Output};
use kira_riboqc::pipeline::stage3_scores::{explain, run_stage3};

fn synthetic_stage2() -> Stage2Output {
    let axes = vec![CellAxes {
        tl: 0.8,
        st: 0.6,
        rqc: 0.7,
        rqc_pressure: 0.0,
        tpc: -0.4,
        st_low_confidence: false,
    }];

    let components = vec![CellAxisComponents {
        rp: 0.0,
        mach: 0.0,
        rrna_proxy: 0.0,
        stress: 0.0,
        metab: 0.0,
        luxury: None,
        rqc_mach: 0.0,
        chaperone: 0.0,
        ubiquitin: 0.0,
    }];

    Stage2Output {
        bg_gene_ids: Vec::new(),
        libsize: vec![0],
        detected_genes: vec![0],
        components,
        axes,
    }
}

#[test]
fn deterministic_scores() {
    let stage2 = synthetic_stage2();
    let out = run_stage3(&stage2).unwrap();
    let cell = &out.scores[0];

    let fragility = 0.5 - (-0.4) / 2.0;
    let tss = 0.45 * 0.8 + 0.35 * 0.7 + 0.20 * 0.6;
    let lfti = 0.40 * 0.7 + 0.35 * 0.6 + 0.25 * fragility;
    let ras = 0.50 * 0.8 + 0.30 * 0.7 + 0.20 * fragility;

    let eps = 1e-6;
    assert!((cell.tss - tss).abs() < eps);
    assert!((cell.lfti - lfti).abs() < eps);
    assert!((cell.ras - ras).abs() < eps);
    assert!(!cell.ras_red_flag);
}

#[test]
fn explainability_format() {
    let out = explain(vec![("A", 0.2), ("B", 0.5), ("C", 0.5)]);
    assert_eq!(out, "+B(0.50) +C(0.50) +A(0.20)");
}
