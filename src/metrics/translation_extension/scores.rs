use std::cmp::Ordering;

pub const TRIM_FRAC: f64 = 0.10;
pub const MAD_SCALE: f64 = 1.4826;
pub const EPS: f64 = 1e-9;

pub const THRESHOLD_TRANSLATION_HIGH: f64 = 2.0;
pub const THRESHOLD_BIOGENESIS_HIGH: f64 = 2.0;
pub const THRESHOLD_ISR_ACTIVE: f64 = 1.5;
pub const THRESHOLD_PROTEOTOXIC_RISK: f64 = 1.5;
pub const THRESHOLD_TRANSLATIONAL_STRESS_MODE: f64 = 2.0;

#[derive(Debug, Clone, Copy)]
pub struct TranslationExtensionCellScores {
    pub ribosome_core: f64,
    pub initiation_core: f64,
    pub bio_core: f64,
    pub mtor_core: f64,
    pub isr_core: f64,
    pub proteo_core: f64,
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
    pub missing_ribosome_core: bool,
    pub missing_initiation_core: bool,
    pub missing_bio_core: bool,
    pub missing_mtor_core: bool,
    pub missing_isr_core: bool,
    pub missing_proteo_core: bool,
}

#[derive(Debug, Clone, Copy)]
pub struct PanelCore {
    pub ribosome_core: f64,
    pub initiation_core: f64,
    pub bio_core: f64,
    pub mtor_core: f64,
    pub isr_core: f64,
    pub proteo_core: f64,
}

#[derive(Debug, Clone, Copy)]
pub struct RobustBaseline {
    pub ribosome_median: f64,
    pub ribosome_mad: f64,
    pub initiation_median: f64,
    pub initiation_mad: f64,
    pub bio_median: f64,
    pub bio_mad: f64,
    pub mtor_median: f64,
    pub mtor_mad: f64,
    pub isr_median: f64,
    pub isr_mad: f64,
    pub proteo_median: f64,
    pub proteo_mad: f64,
}

pub fn trimmed_mean(values: &[f64], min_genes: usize) -> f64 {
    if values.len() < min_genes {
        return f64::NAN;
    }
    let mut sorted = values.to_vec();
    sorted.sort_by(total_cmp);
    let trim = ((sorted.len() as f64) * TRIM_FRAC).floor() as usize;
    let start = trim.min(sorted.len());
    let end = sorted.len().saturating_sub(trim);
    if end <= start {
        return f64::NAN;
    }
    let sum: f64 = sorted[start..end].iter().sum();
    sum / (end - start) as f64
}

pub fn robust_baseline(panel_cores: &[PanelCore]) -> RobustBaseline {
    let ribosome = collect_finite(panel_cores.iter().map(|x| x.ribosome_core));
    let initiation = collect_finite(panel_cores.iter().map(|x| x.initiation_core));
    let bio = collect_finite(panel_cores.iter().map(|x| x.bio_core));
    let mtor = collect_finite(panel_cores.iter().map(|x| x.mtor_core));
    let isr = collect_finite(panel_cores.iter().map(|x| x.isr_core));
    let proteo = collect_finite(panel_cores.iter().map(|x| x.proteo_core));

    let (ribosome_median, ribosome_mad) = median_and_mad(&ribosome);
    let (initiation_median, initiation_mad) = median_and_mad(&initiation);
    let (bio_median, bio_mad) = median_and_mad(&bio);
    let (mtor_median, mtor_mad) = median_and_mad(&mtor);
    let (isr_median, isr_mad) = median_and_mad(&isr);
    let (proteo_median, proteo_mad) = median_and_mad(&proteo);

    RobustBaseline {
        ribosome_median,
        ribosome_mad,
        initiation_median,
        initiation_mad,
        bio_median,
        bio_mad,
        mtor_median,
        mtor_mad,
        isr_median,
        isr_mad,
        proteo_median,
        proteo_mad,
    }
}

pub fn robust_z(value: f64, median: f64, mad: f64) -> f64 {
    if value.is_nan() || median.is_nan() {
        return f64::NAN;
    }
    if mad == 0.0 || mad.is_nan() {
        return 0.0;
    }
    (value - median) / (MAD_SCALE * mad + EPS)
}

pub fn build_scores(core: PanelCore, baseline: &RobustBaseline) -> TranslationExtensionCellScores {
    let zr = robust_z(
        core.ribosome_core,
        baseline.ribosome_median,
        baseline.ribosome_mad,
    );
    let zi = robust_z(
        core.initiation_core,
        baseline.initiation_median,
        baseline.initiation_mad,
    );
    let zbio = robust_z(core.bio_core, baseline.bio_median, baseline.bio_mad);
    let zmtor = robust_z(core.mtor_core, baseline.mtor_median, baseline.mtor_mad);
    let zisr = robust_z(core.isr_core, baseline.isr_median, baseline.isr_mad);
    let zp = robust_z(
        core.proteo_core,
        baseline.proteo_median,
        baseline.proteo_mad,
    );

    let tpi = combine_weighted(&[(zr, 0.7), (zi, 0.3)]);
    let rbl = zbio;
    let mtor_p = combine_weighted(&[(zmtor, 0.6), (zr, 0.4)]);
    let isr_a = zisr;
    let tpib = if zp.is_nan() {
        if tpi.is_nan() { f64::NAN } else { tpi.max(0.0) }
    } else if tpi.is_nan() {
        f64::NAN
    } else {
        (tpi - zp).max(0.0)
    };
    let tsm = combine_weighted(&[(tpi, 0.5), (isr_a, 0.3), (tpib, 0.2)]);

    TranslationExtensionCellScores {
        ribosome_core: core.ribosome_core,
        initiation_core: core.initiation_core,
        bio_core: core.bio_core,
        mtor_core: core.mtor_core,
        isr_core: core.isr_core,
        proteo_core: core.proteo_core,
        tpi,
        rbl,
        mtor_p,
        isr_a,
        tpib,
        tsm,
        translation_high: tpi.is_finite() && tpi >= THRESHOLD_TRANSLATION_HIGH,
        biogenesis_high: rbl.is_finite() && rbl >= THRESHOLD_BIOGENESIS_HIGH,
        isr_active: isr_a.is_finite() && isr_a >= THRESHOLD_ISR_ACTIVE,
        proteotoxic_risk: tpib.is_finite() && tpib >= THRESHOLD_PROTEOTOXIC_RISK,
        translational_stress_mode: tsm.is_finite() && tsm >= THRESHOLD_TRANSLATIONAL_STRESS_MODE,
        missing_ribosome_core: core.ribosome_core.is_nan(),
        missing_initiation_core: core.initiation_core.is_nan(),
        missing_bio_core: core.bio_core.is_nan(),
        missing_mtor_core: core.mtor_core.is_nan(),
        missing_isr_core: core.isr_core.is_nan(),
        missing_proteo_core: core.proteo_core.is_nan(),
    }
}

pub fn median(values: &[f64]) -> f64 {
    let mut sorted = values.to_vec();
    if sorted.is_empty() {
        return f64::NAN;
    }
    sorted.sort_by(total_cmp);
    let mid = sorted.len() / 2;
    if sorted.len() % 2 == 0 {
        (sorted[mid - 1] + sorted[mid]) * 0.5
    } else {
        sorted[mid]
    }
}

pub fn percentile(values: &[f64], p: f64) -> f64 {
    if values.is_empty() {
        return f64::NAN;
    }
    let mut sorted = values.to_vec();
    sorted.sort_by(total_cmp);
    percentile_sorted(&sorted, p)
}

fn percentile_sorted(sorted: &[f64], p: f64) -> f64 {
    let rank = (p * sorted.len() as f64).ceil() as usize;
    let idx = if rank == 0 { 0 } else { rank - 1 };
    sorted[idx.min(sorted.len() - 1)]
}

fn median_and_mad(values: &[f64]) -> (f64, f64) {
    if values.is_empty() {
        return (f64::NAN, f64::NAN);
    }
    let med = median(values);
    if med.is_nan() {
        return (f64::NAN, f64::NAN);
    }
    let deviations: Vec<f64> = values.iter().map(|v| (v - med).abs()).collect();
    let mad = median(&deviations);
    (med, mad)
}

fn collect_finite<I>(iter: I) -> Vec<f64>
where
    I: Iterator<Item = f64>,
{
    iter.filter(|v| v.is_finite()).collect()
}

fn combine_weighted(parts: &[(f64, f64)]) -> f64 {
    let mut weighted_sum = 0.0;
    let mut weight_sum = 0.0;
    for (value, weight) in parts {
        if value.is_finite() {
            weighted_sum += value * weight;
            weight_sum += weight;
        }
    }
    if weight_sum == 0.0 {
        f64::NAN
    } else {
        weighted_sum / weight_sum
    }
}

fn total_cmp(a: &f64, b: &f64) -> Ordering {
    a.partial_cmp(b).unwrap_or(Ordering::Equal)
}
