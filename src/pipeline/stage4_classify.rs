use std::collections::BTreeMap;

use tracing::info;

use crate::model::classification::{
    MIN_DETECTED_GENES, MIN_LIBSIZE, RQC_HIGH, RQC_LOW, RQC_MID, ST_HIGH, ST_LOW, TL_HIGH, TL_LOW,
    TL_MID_HIGH, TL_MID_LOW, TPC_NEG, TPC_NEG_SOFT, TranslationRegime, regime_label,
};
use crate::pipeline::stage2_axes::{CellAxes, Stage2Output};
use crate::pipeline::stage3_scores::Stage3Output;

pub struct CellClassification {
    pub regime: TranslationRegime,
    pub low_counts_cell: bool,
    pub few_detected_genes: bool,
}

pub struct Stage4Output {
    pub classification: Vec<CellClassification>,
}

pub fn classify_cell(axes: &CellAxes, _libsize: u64, _detected_genes: u32) -> TranslationRegime {
    let tl = axes.tl;
    let st = axes.st;
    let rqc = axes.rqc;
    let tpc = axes.tpc;

    if tl < TL_LOW && st < ST_LOW && rqc < RQC_LOW {
        return TranslationRegime::TranslationSuppressed;
    }

    if tl >= TL_MID_LOW && tl <= TL_MID_HIGH && rqc < RQC_MID && tpc >= TPC_NEG_SOFT {
        return TranslationRegime::EfficientTranslation;
    }

    if st >= ST_HIGH && rqc >= RQC_MID {
        return TranslationRegime::SelectiveSurvivalTranslation;
    }

    if tl >= TL_HIGH && rqc >= RQC_MID && tpc < TPC_NEG {
        return TranslationRegime::OverloadedTranslation;
    }

    if rqc >= RQC_HIGH && tl >= TL_MID_LOW && tpc < 0.0 {
        return TranslationRegime::RQC_Dependent;
    }

    TranslationRegime::Unclassified
}

pub fn run_stage4(stage2: &Stage2Output, stage3: &Stage3Output) -> anyhow::Result<Stage4Output> {
    let _span = tracing::info_span!("stage4_classify").entered();
    info!("Stage 4 start");

    let n_cells = stage2.axes.len();
    if stage3.scores.len() != n_cells {
        return Err(anyhow::anyhow!(
            "stage3 length mismatch: {} vs {}",
            stage3.scores.len(),
            n_cells
        ));
    }

    let mut classification = Vec::with_capacity(n_cells);
    let mut regime_counts: BTreeMap<&'static str, u32> = BTreeMap::new();
    let mut low_counts = 0u32;
    let mut few_detected = 0u32;

    for idx in 0..n_cells {
        let axes = &stage2.axes[idx];
        let libsize = stage2.libsize[idx];
        let detected = stage2.detected_genes[idx];

        let regime = classify_cell(axes, libsize, detected);
        let low_counts_cell = libsize < MIN_LIBSIZE;
        let few_detected_genes = detected < MIN_DETECTED_GENES;

        if low_counts_cell {
            low_counts += 1;
        }
        if few_detected_genes {
            few_detected += 1;
        }

        *regime_counts.entry(regime_label(regime)).or_insert(0) += 1;

        classification.push(CellClassification {
            regime,
            low_counts_cell,
            few_detected_genes,
        });
    }

    info!(
        low_counts_frac = if n_cells == 0 {
            0.0
        } else {
            low_counts as f64 / n_cells as f64
        },
        few_detected_frac = if n_cells == 0 {
            0.0
        } else {
            few_detected as f64 / n_cells as f64
        },
        "Stage 4 summary"
    );

    for (regime, count) in regime_counts {
        info!(regime = regime, count, "Stage 4 regime count");
    }

    info!("Stage 4 end");

    Ok(Stage4Output { classification })
}
