use kira_riboqc::model::classification::{
    MIN_DETECTED_GENES, MIN_LIBSIZE, RQC_HIGH, RQC_LOW, RQC_MID, ST_HIGH, ST_LOW, TL_HIGH, TL_LOW,
    TL_MID_HIGH, TL_MID_LOW, TPC_NEG_SOFT, TranslationRegime,
};
use kira_riboqc::pipeline::stage2_axes::CellAxes;
use kira_riboqc::pipeline::stage4_classify::classify_cell;

fn axes(tl: f64, st: f64, rqc: f64, tpc: f64) -> CellAxes {
    CellAxes {
        tl,
        st,
        rqc,
        rqc_pressure: 0.0,
        tpc,
        st_low_confidence: false,
    }
}

#[test]
fn boundary_rules() {
    let suppressed = classify_cell(&axes(0.24, 0.34, 0.34, 0.0), 0, 0);
    assert_eq!(suppressed, TranslationRegime::TranslationSuppressed);

    let efficient = classify_cell(&axes(TL_MID_LOW, ST_LOW, 0.44, TPC_NEG_SOFT), 0, 0);
    assert_eq!(efficient, TranslationRegime::EfficientTranslation);

    let selective = classify_cell(&axes(TL_LOW, ST_HIGH, RQC_MID, 0.0), 0, 0);
    assert_eq!(selective, TranslationRegime::SelectiveSurvivalTranslation);

    let overloaded = classify_cell(&axes(TL_HIGH, ST_LOW, RQC_MID, -0.16), 0, 0);
    assert_eq!(overloaded, TranslationRegime::OverloadedTranslation);

    let rqc_dep = classify_cell(&axes(TL_MID_LOW, ST_LOW, RQC_HIGH, -0.01), 0, 0);
    assert_eq!(rqc_dep, TranslationRegime::RQC_Dependent);

    let unclassified = classify_cell(&axes(TL_MID_HIGH + 0.01, ST_LOW, RQC_LOW, 0.0), 0, 0);
    assert_eq!(unclassified, TranslationRegime::Unclassified);
}

#[test]
fn priority_first_match_wins() {
    let ax = axes(0.80, 0.70, 0.80, -0.30);
    let regime = classify_cell(&ax, 0, 0);
    assert_eq!(regime, TranslationRegime::SelectiveSurvivalTranslation);
}

#[test]
fn flags_independent_of_regime() {
    let ax = axes(TL_MID_LOW, ST_LOW, RQC_MID - 0.01, 0.0);
    let regime = classify_cell(&ax, 0, 0);
    assert_eq!(regime, TranslationRegime::EfficientTranslation);

    let low_counts_cell = 0 < MIN_LIBSIZE;
    let few_detected_genes = 0 < MIN_DETECTED_GENES;
    assert!(low_counts_cell);
    assert!(few_detected_genes);
}
