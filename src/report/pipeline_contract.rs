use rustc_hash::FxHashMap;

use crate::input::InputBundle;
use crate::model::axes::clamp01;
use crate::model::classification::TranslationRegime;
use crate::pipeline::stage_translation_regime::{
    StageTranslationRegimeOutput, TranslationRegimeCell,
};
use crate::pipeline::stage2_axes::Stage2Output;
use crate::pipeline::stage3_scores::Stage3Output;
use crate::pipeline::stage4_classify::Stage4Output;

#[derive(Debug, Clone)]
pub struct PipelineCellRow {
    pub barcode: String,
    pub sample: String,
    pub condition: String,
    pub species: String,
    pub libsize: u64,
    pub nnz: u32,
    pub expressed_genes: u32,
    pub translation_load: f64,
    pub ribosome_density: f64,
    pub elongation_pressure: f64,
    pub initiation_bias: f64,
    pub ribosomal_specialization: f64,
    pub stress_translation_index: f64,
    pub ribosome_core: f64,
    pub initiation_core: f64,
    pub bio_core: f64,
    pub mtor_core: f64,
    pub isr_core: f64,
    pub tpi: f64,
    pub rbl: f64,
    pub mtor_p: f64,
    pub isr_a: f64,
    pub tpib: f64,
    pub tsm: f64,
    pub translation_high: bool,
    pub biogenesis_high: bool,
    pub isr_active: bool,
    pub proteotoxic_risk: bool,
    pub translational_stress_mode: bool,
    pub regime: &'static str,
    pub flags: String,
    pub confidence: f64,
    pub has_low_confidence: bool,
    pub has_low_ribo_signal: bool,
}

const FLAG_LOW_COUNTS: u32 = 1 << 0;
const FLAG_FEW_DETECTED: u32 = 1 << 1;
const FLAG_LOW_CONFIDENCE: u32 = 1 << 2;
const FLAG_TRANSLATION_HIGH: u32 = 1 << 3;
const FLAG_BIOGENESIS_HIGH: u32 = 1 << 4;
const FLAG_ISR_ACTIVE: u32 = 1 << 5;
const FLAG_PROTEOTOXIC_RISK: u32 = 1 << 6;
const FLAG_TRANSLATIONAL_STRESS: u32 = 1 << 7;
const FLAG_LOW_RIBO_SIGNAL: u32 = 1 << 8;
const FLAG_RAS_RED_FLAG: u32 = 1 << 9;

const FLAGS: &[(u32, &str)] = &[
    (FLAG_BIOGENESIS_HIGH, "BIOGENESIS_HIGH"),
    (FLAG_FEW_DETECTED, "FEW_DETECTED_GENES"),
    (FLAG_ISR_ACTIVE, "ISR_ACTIVE"),
    (FLAG_LOW_CONFIDENCE, "LOW_CONFIDENCE"),
    (FLAG_LOW_COUNTS, "LOW_COUNTS_CELL"),
    (FLAG_LOW_RIBO_SIGNAL, "LOW_RIBO_SIGNAL"),
    (FLAG_PROTEOTOXIC_RISK, "PROTEOTOXIC_RISK"),
    (FLAG_RAS_RED_FLAG, "RAS_RED_FLAG"),
    (FLAG_TRANSLATIONAL_STRESS, "TRANSLATIONAL_STRESS_MODE"),
    (FLAG_TRANSLATION_HIGH, "TRANSLATION_HIGH"),
];

pub fn build_pipeline_rows(
    input: &InputBundle,
    stage2: &Stage2Output,
    stage3: &Stage3Output,
    stage4: &Stage4Output,
    stage_translation: Option<&StageTranslationRegimeOutput>,
) -> Vec<PipelineCellRow> {
    let species = infer_species(input);
    let mut rows = Vec::with_capacity(input.barcodes.len());

    let translation_cells: FxHashMap<&str, &TranslationRegimeCell> = stage_translation
        .map(|s| s.cells.iter().map(|c| (c.cell_id.as_str(), c)).collect())
        .unwrap_or_default();

    for i in 0..input.barcodes.len() {
        let barcode = &input.barcodes[i];
        let meta = input.metadata.as_ref().and_then(|m| m.rows.get(barcode));
        let sample = meta
            .and_then(|m| m.field("sample"))
            .unwrap_or("unknown")
            .to_string();
        let condition = meta
            .and_then(|m| m.field("condition"))
            .unwrap_or("unknown")
            .to_string();

        let translation_load = nan_to_zero(stage2.axes[i].tl);
        let ribosome_density = nan_to_zero(stage2.components[i].rp);
        let elongation_pressure = nan_to_zero(stage2.axes[i].rqc_pressure);
        let initiation_bias = nan_to_zero(stage2.components[i].mach);
        let ribosomal_specialization = nan_to_zero(stage2.axes[i].rqc);
        let stress_translation_index = nan_to_zero(stage2.axes[i].st);
        let tcell = translation_cells.get(barcode.as_str()).copied();

        let regime = map_regime(stage4.classification[i].regime, translation_load);

        let mut flag_bits = 0u32;
        if stage4.classification[i].low_counts_cell {
            flag_bits |= FLAG_LOW_COUNTS;
        }
        if stage4.classification[i].few_detected_genes {
            flag_bits |= FLAG_FEW_DETECTED;
        }
        if stage2.axes[i].st_low_confidence {
            flag_bits |= FLAG_LOW_CONFIDENCE;
        }
        if let Some(c) = tcell {
            if c.low_confidence {
                flag_bits |= FLAG_LOW_CONFIDENCE;
            }
            if c.translation_high {
                flag_bits |= FLAG_TRANSLATION_HIGH;
            }
            if c.biogenesis_high {
                flag_bits |= FLAG_BIOGENESIS_HIGH;
            }
            if c.isr_active {
                flag_bits |= FLAG_ISR_ACTIVE;
            }
            if c.proteotoxic_risk {
                flag_bits |= FLAG_PROTEOTOXIC_RISK;
            }
            if c.translational_stress_mode {
                flag_bits |= FLAG_TRANSLATIONAL_STRESS;
            }
        }
        let low_ribo_signal = ribosome_density < 0.2;
        if low_ribo_signal {
            flag_bits |= FLAG_LOW_RIBO_SIGNAL;
        }
        if stage3.scores[i].ras_red_flag {
            flag_bits |= FLAG_RAS_RED_FLAG;
        }

        let flags = format_flags(flag_bits);
        let has_low_confidence = flag_bits & FLAG_LOW_CONFIDENCE != 0;

        let confidence = calc_confidence(
            stage4.classification[i].low_counts_cell,
            stage4.classification[i].few_detected_genes,
            has_low_confidence,
            low_ribo_signal,
        );

        rows.push(PipelineCellRow {
            barcode: barcode.clone(),
            sample,
            condition,
            species: species.to_string(),
            libsize: stage2.libsize[i],
            nnz: stage2.detected_genes[i],
            expressed_genes: stage2.detected_genes[i],
            translation_load,
            ribosome_density,
            elongation_pressure,
            initiation_bias,
            ribosomal_specialization,
            stress_translation_index,
            ribosome_core: tcell.map(|c| c.ribosome_core).unwrap_or(f64::NAN),
            initiation_core: tcell.map(|c| c.initiation_core).unwrap_or(f64::NAN),
            bio_core: tcell.map(|c| c.bio_core).unwrap_or(f64::NAN),
            mtor_core: tcell.map(|c| c.mtor_core).unwrap_or(f64::NAN),
            isr_core: tcell.map(|c| c.isr_core).unwrap_or(f64::NAN),
            tpi: tcell.map(|c| c.tpi).unwrap_or(f64::NAN),
            rbl: tcell.map(|c| c.rbl).unwrap_or(f64::NAN),
            mtor_p: tcell.map(|c| c.mtor_p).unwrap_or(f64::NAN),
            isr_a: tcell.map(|c| c.isr_a).unwrap_or(f64::NAN),
            tpib: tcell.map(|c| c.tpib).unwrap_or(f64::NAN),
            tsm: tcell.map(|c| c.tsm).unwrap_or(f64::NAN),
            translation_high: tcell.map(|c| c.translation_high).unwrap_or(false),
            biogenesis_high: tcell.map(|c| c.biogenesis_high).unwrap_or(false),
            isr_active: tcell.map(|c| c.isr_active).unwrap_or(false),
            proteotoxic_risk: tcell.map(|c| c.proteotoxic_risk).unwrap_or(false),
            translational_stress_mode: tcell.map(|c| c.translational_stress_mode).unwrap_or(false),
            regime,
            flags,
            confidence,
            has_low_confidence,
            has_low_ribo_signal: low_ribo_signal,
        });
    }

    rows.sort_by(|a, b| a.barcode.cmp(&b.barcode));
    rows
}

fn format_flags(bits: u32) -> String {
    if bits == 0 {
        return String::new();
    }
    let mut out = String::with_capacity(64);
    let mut first = true;
    for (bit, name) in FLAGS {
        if bits & bit != 0 {
            if !first {
                out.push(',');
            }
            out.push_str(name);
            first = false;
        }
    }
    out
}

fn calc_confidence(
    low_counts_cell: bool,
    few_detected_genes: bool,
    low_conf: bool,
    low_ribo_signal: bool,
) -> f64 {
    let mut conf = 1.0;
    if low_counts_cell {
        conf -= 0.30;
    }
    if few_detected_genes {
        conf -= 0.25;
    }
    if low_conf {
        conf -= 0.20;
    }
    if low_ribo_signal {
        conf -= 0.15;
    }
    clamp01(conf).max(0.0)
}

fn map_regime(regime: TranslationRegime, translation_load: f64) -> &'static str {
    match regime {
        TranslationRegime::EfficientTranslation => {
            if translation_load >= 0.60 {
                "GrowthDrivenTranslation"
            } else {
                "HomeostaticTranslation"
            }
        }
        TranslationRegime::SelectiveSurvivalTranslation => "StressAdaptiveTranslation",
        TranslationRegime::OverloadedTranslation => "TranslationalOverdrive",
        TranslationRegime::RQC_Dependent => "StressAdaptiveTranslation",
        TranslationRegime::TranslationSuppressed => "TranslationalCollapse",
        TranslationRegime::Unclassified => "Unclassified",
    }
}

fn infer_species(input: &InputBundle) -> &'static str {
    let mut human_like = 0usize;
    let mut mouse_like = 0usize;

    for f in &input.features {
        let raw = f.raw_name.as_str();
        if raw.is_empty() {
            continue;
        }
        let bytes = raw.as_bytes();
        let mut is_all_upper = true;
        let mut first_upper = false;
        let mut any_lower_after_first = false;
        for (i, &b) in bytes.iter().enumerate() {
            if b.is_ascii_alphabetic() {
                if b.is_ascii_lowercase() {
                    is_all_upper = false;
                    if i > 0 {
                        any_lower_after_first = true;
                    }
                } else if i == 0 {
                    first_upper = true;
                }
            }
        }
        if is_all_upper {
            human_like += 1;
        } else if first_upper && any_lower_after_first {
            mouse_like += 1;
        }
    }

    if human_like > mouse_like * 2 {
        "human"
    } else if mouse_like > human_like * 2 {
        "mouse"
    } else {
        "unknown"
    }
}

#[inline]
fn nan_to_zero(v: f64) -> f64 {
    if v.is_nan() { 0.0 } else { v }
}
