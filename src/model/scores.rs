use crate::model::axes::clamp01;

pub fn fragility_from_tpc(tpc: f64) -> f64 {
    clamp01(0.5 - tpc / 2.0)
}

pub fn tss(tl: f64, rqc: f64, st: f64) -> f64 {
    clamp01(0.45 * tl + 0.35 * rqc + 0.20 * st)
}

pub fn lfti(rqc: f64, st: f64, fragility: f64) -> f64 {
    clamp01(0.40 * rqc + 0.35 * st + 0.25 * fragility)
}

pub fn ras(tl: f64, rqc: f64, fragility: f64) -> f64 {
    clamp01(0.50 * tl + 0.30 * rqc + 0.20 * fragility)
}
