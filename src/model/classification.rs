use serde::Serialize;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[allow(non_camel_case_types)]
pub enum TranslationRegime {
    EfficientTranslation,
    SelectiveSurvivalTranslation,
    OverloadedTranslation,
    RQC_Dependent,
    TranslationSuppressed,
    Unclassified,
}

pub const TL_LOW: f64 = 0.25;
pub const TL_MID_LOW: f64 = 0.35;
pub const TL_MID_HIGH: f64 = 0.70;
pub const TL_HIGH: f64 = 0.75;

pub const ST_LOW: f64 = 0.35;
pub const ST_HIGH: f64 = 0.65;

pub const RQC_LOW: f64 = 0.35;
pub const RQC_MID: f64 = 0.45;
pub const RQC_HIGH: f64 = 0.70;

pub const TPC_NEG_SOFT: f64 = -0.05;
pub const TPC_NEG: f64 = -0.15;

pub const RAS_RED_TL: f64 = 0.70;
pub const RAS_RED_RQC: f64 = 0.70;
pub const RAS_RED_TPC: f64 = -0.20;

pub const MIN_LIBSIZE: u64 = 500;
pub const MIN_DETECTED_GENES: u32 = 200;

pub fn regime_label(regime: TranslationRegime) -> &'static str {
    match regime {
        TranslationRegime::EfficientTranslation => "EfficientTranslation",
        TranslationRegime::SelectiveSurvivalTranslation => "SelectiveSurvivalTranslation",
        TranslationRegime::OverloadedTranslation => "OverloadedTranslation",
        TranslationRegime::RQC_Dependent => "RQC_Dependent",
        TranslationRegime::TranslationSuppressed => "TranslationSuppressed",
        TranslationRegime::Unclassified => "Unclassified",
    }
}
