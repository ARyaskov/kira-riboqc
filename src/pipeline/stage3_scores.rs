use std::fmt::Write;

use tracing::info;

use crate::core::math::median_non_nan;
use crate::model::scores::{fragility_from_tpc, lfti, ras, tss};
use crate::pipeline::stage2_axes::Stage2Output;

pub struct CellScores {
    pub tss: f64,
    pub lfti: f64,
    pub ras: f64,
    pub ras_red_flag: bool,
    pub tss_explain: String,
    pub lfti_explain: String,
    pub ras_explain: String,
}

pub struct Stage3Output {
    pub scores: Vec<CellScores>,
}

pub fn run_stage3(stage2: &Stage2Output) -> anyhow::Result<Stage3Output> {
    let _span = tracing::info_span!("stage3_scores").entered();
    info!("Stage 3 start");

    let mut scores = Vec::with_capacity(stage2.axes.len());
    let mut red_flags = 0u32;

    for ax in &stage2.axes {
        let tl = ax.tl;
        let st = ax.st;
        let rqc = ax.rqc;
        let tpc = ax.tpc;
        let fragility = fragility_from_tpc(tpc);

        let tss_val = tss(tl, rqc, st);
        let lfti_val = lfti(rqc, st, fragility);
        let ras_val = ras(tl, rqc, fragility);
        let ras_red_flag = tl > 0.70 && rqc > 0.70 && tpc < -0.20;

        if ras_red_flag {
            red_flags += 1;
        }

        let tss_explain = explain(&[("TL", 0.45 * tl), ("RQC", 0.35 * rqc), ("ST", 0.20 * st)]);
        let lfti_explain = explain(&[
            ("RQC", 0.40 * rqc),
            ("ST", 0.35 * st),
            ("fragility", 0.25 * fragility),
        ]);
        let ras_explain = explain(&[
            ("TL", 0.50 * tl),
            ("RQC", 0.30 * rqc),
            ("fragility", 0.20 * fragility),
        ]);

        scores.push(CellScores {
            tss: tss_val,
            lfti: lfti_val,
            ras: ras_val,
            ras_red_flag,
            tss_explain,
            lfti_explain,
            ras_explain,
        });
    }

    let tss_median = median_non_nan(scores.iter().map(|s| s.tss));
    let lfti_median = median_non_nan(scores.iter().map(|s| s.lfti));
    let ras_median = median_non_nan(scores.iter().map(|s| s.ras));
    let frac_red_flag = if scores.is_empty() {
        0.0
    } else {
        red_flags as f64 / scores.len() as f64
    };

    info!(
        tss_median = tss_median,
        lfti_median = lfti_median,
        ras_median = ras_median,
        ras_red_flag_frac = frac_red_flag,
        "Stage 3 summary"
    );
    info!("Stage 3 end");

    Ok(Stage3Output { scores })
}

pub fn explain(components: &[(&str, f64)]) -> String {
    let mut parts: [(&str, f64); 4] = [("", 0.0); 4];
    let n = components.len().min(parts.len());
    parts[..n].copy_from_slice(&components[..n]);
    parts[..n].sort_by(|a, b| {
        b.1.abs()
            .partial_cmp(&a.1.abs())
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| a.0.cmp(b.0))
    });

    let mut out = String::with_capacity(48);
    for (idx, (name, value)) in parts[..n].iter().enumerate() {
        if idx > 0 {
            out.push(' ');
        }
        out.push(if *value < 0.0 { '-' } else { '+' });
        out.push_str(name);
        let _ = write!(out, "({:.2})", value.abs());
    }
    out
}
