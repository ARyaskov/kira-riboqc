use std::cmp::Ordering;

#[inline]
fn total_cmp(a: &f64, b: &f64) -> Ordering {
    a.partial_cmp(b).unwrap_or(Ordering::Equal)
}

pub fn median_sorted(sorted: &[f64]) -> f64 {
    if sorted.is_empty() {
        return f64::NAN;
    }
    let n = sorted.len();
    let mid = n / 2;
    if n.is_multiple_of(2) {
        (sorted[mid - 1] + sorted[mid]) * 0.5
    } else {
        sorted[mid]
    }
}

pub fn percentile_sorted(sorted: &[f64], p: f64) -> f64 {
    if sorted.is_empty() {
        return f64::NAN;
    }
    let rank = (p * sorted.len() as f64).ceil() as usize;
    let idx = rank.saturating_sub(1);
    sorted[idx.min(sorted.len() - 1)]
}

pub fn median(values: &[f64]) -> f64 {
    if values.is_empty() {
        return f64::NAN;
    }
    let mut buf = values.to_vec();
    buf.sort_by(total_cmp);
    median_sorted(&buf)
}

pub fn percentile(values: &[f64], p: f64) -> f64 {
    if values.is_empty() {
        return f64::NAN;
    }
    let mut buf = values.to_vec();
    buf.sort_by(total_cmp);
    percentile_sorted(&buf, p)
}

pub fn median_non_nan<I>(iter: I) -> f64
where
    I: Iterator<Item = f64>,
{
    let mut values: Vec<f64> = iter.filter(|v| !v.is_nan()).collect();
    if values.is_empty() {
        return f64::NAN;
    }
    values.sort_by(total_cmp);
    median_sorted(&values)
}

#[inline]
pub fn round6(x: f64) -> f64 {
    (x * 1_000_000.0).round() / 1_000_000.0
}

#[inline]
pub fn round6_nan(x: f64) -> f64 {
    if x.is_nan() { f64::NAN } else { round6(x) }
}

#[inline]
pub fn round6_or_zero(x: f64) -> f64 {
    if x.is_nan() { 0.0 } else { round6(x) }
}
