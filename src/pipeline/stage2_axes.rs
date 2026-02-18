use tracing::info;

use crate::core::membership::MembershipVec;
use crate::input::{InputBundle, SharedCacheData};
use crate::model::axes::{
    BG_TOP_K, CHAPERONE, LUXURY, MACHINERY, METAB, OFFSET, RQC_MACH, RRNA_PROXY, SCALE, STRESS,
    UBIQUITIN, clamp_to_minus1_plus1, clamp01, is_mito, is_ribosomal, resolve_gene_set,
};
use crate::simd::ln1p_f64;

pub struct CellAxisComponents {
    pub rp: f64,
    pub mach: f64,
    pub rrna_proxy: f64,
    pub stress: f64,
    pub metab: f64,
    pub luxury: Option<f64>,
    pub rqc_mach: f64,
    pub chaperone: f64,
    pub ubiquitin: f64,
}

pub struct CellAxes {
    pub tl: f64,
    pub st: f64,
    pub rqc: f64,
    pub rqc_pressure: f64,
    pub tpc: f64,
    pub st_low_confidence: bool,
}

pub struct Stage2Output {
    pub bg_gene_ids: Vec<u32>,
    pub libsize: Vec<u64>,
    pub detected_genes: Vec<u32>,
    pub components: Vec<CellAxisComponents>,
    pub axes: Vec<CellAxes>,
}

pub fn run_stage2(input: &InputBundle) -> anyhow::Result<Stage2Output> {
    let _span = tracing::info_span!("stage2_axes").entered();
    info!("Stage 2 start");

    let source = MatrixSource::from_input(input);
    let n_cells = source.n_cols();
    let n_genes_unique = input.gene_index.genes.len();

    let (libsize, detected_genes) = compute_cell_stats(&source);
    let gene_detect_counts = compute_gene_detection_counts(input, &source, n_genes_unique);
    let bg_gene_ids = select_background_genes(input, &gene_detect_counts, n_cells as u32);

    let rp_gene_ids = collect_ribosomal_gene_ids(input);
    let rrna_gene_ids = resolve_gene_set(&input.gene_index.map, RRNA_PROXY);
    let mach_gene_ids = resolve_gene_set(&input.gene_index.map, MACHINERY);
    let stress_gene_ids = resolve_gene_set(&input.gene_index.map, STRESS);
    let metab_gene_ids = resolve_gene_set(&input.gene_index.map, METAB);
    let luxury_gene_ids = resolve_gene_set(&input.gene_index.map, LUXURY);
    let rqc_gene_ids = resolve_gene_set(&input.gene_index.map, RQC_MACH);
    let chaperone_gene_ids = resolve_gene_set(&input.gene_index.map, CHAPERONE);
    let ubiquitin_gene_ids = resolve_gene_set(&input.gene_index.map, UBIQUITIN);

    let rp_mem = MembershipVec::from_gene_ids(&rp_gene_ids, n_genes_unique);
    let rrna_mem = MembershipVec::from_gene_ids(&rrna_gene_ids, n_genes_unique);
    let mach_mem = MembershipVec::from_gene_ids(&mach_gene_ids, n_genes_unique);
    let stress_mem = MembershipVec::from_gene_ids(&stress_gene_ids, n_genes_unique);
    let metab_mem = MembershipVec::from_gene_ids(&metab_gene_ids, n_genes_unique);
    let luxury_mem = MembershipVec::from_gene_ids(&luxury_gene_ids, n_genes_unique);
    let rqc_mem = MembershipVec::from_gene_ids(&rqc_gene_ids, n_genes_unique);
    let chaperone_mem = MembershipVec::from_gene_ids(&chaperone_gene_ids, n_genes_unique);
    let ubiquitin_mem = MembershipVec::from_gene_ids(&ubiquitin_gene_ids, n_genes_unique);
    let bg_mem = MembershipVec::from_gene_ids(&bg_gene_ids, n_genes_unique);

    let mut components = Vec::with_capacity(n_cells);
    let mut axes = Vec::with_capacity(n_cells);
    let mut st_low_conf_count: u32 = 0;

    for col in 0..n_cells {
        let col_view = source.column_view(col);
        let ls = libsize[col];

        let mut acc_rp = Acc::default();
        let mut acc_mach = Acc::default();
        let mut acc_rrna = Acc::default();
        let mut acc_stress = Acc::default();
        let mut acc_metab = Acc::default();
        let mut acc_luxury = Acc::default();
        let mut acc_rqc = Acc::default();
        let mut acc_chaperone = Acc::default();
        let mut acc_ubiquitin = Acc::default();
        let mut acc_bg = Acc::default();

        let denom = if ls == 0 { 1.0 } else { ls as f64 };
        for (gid, val) in col_view.iter(input.gene_index.row_to_gene.as_slice()) {
            let cpm = 1_000_000.0 * (val as f64) / denom;
            let x = ln1p_f64(cpm);

            if rp_mem.contains(gid) {
                acc_rp.add(x);
            }
            if mach_mem.contains(gid) {
                acc_mach.add(x);
            }
            if rrna_mem.contains(gid) {
                acc_rrna.add(x);
            }
            if stress_mem.contains(gid) {
                acc_stress.add(x);
            }
            if metab_mem.contains(gid) {
                acc_metab.add(x);
            }
            if luxury_mem.contains(gid) {
                acc_luxury.add(x);
            }
            if rqc_mem.contains(gid) {
                acc_rqc.add(x);
            }
            if chaperone_mem.contains(gid) {
                acc_chaperone.add(x);
            }
            if ubiquitin_mem.contains(gid) {
                acc_ubiquitin.add(x);
            }
            if bg_mem.contains(gid) {
                acc_bg.add(x);
            }
        }

        let bg_mean = acc_bg.mean();
        let rp = enrich(acc_rp.mean(), bg_mean);
        let mach = enrich(acc_mach.mean(), bg_mean);
        let rrna_proxy = enrich(acc_rrna.mean(), bg_mean);
        let stress = enrich(acc_stress.mean(), bg_mean);
        let metab = enrich(acc_metab.mean(), bg_mean);
        let luxury_raw = enrich(acc_luxury.mean(), bg_mean);
        let rqc_mach = enrich(acc_rqc.mean(), bg_mean);
        let chaperone = enrich(acc_chaperone.mean(), bg_mean);
        let ubiquitin = enrich(acc_ubiquitin.mean(), bg_mean);

        let luxury = if luxury_raw.is_nan() {
            None
        } else {
            Some(luxury_raw)
        };

        let (tl, _missing_component) = compute_tl(rp, mach, rrna_proxy);
        let (st, st_low_confidence) = compute_st(stress, metab, luxury);
        let (rqc, rqc_pressure) = compute_rqc(rqc_mach, stress, rp);
        let tpc = compute_tpc(chaperone, ubiquitin, tl);

        if st_low_confidence {
            st_low_conf_count += 1;
        }

        components.push(CellAxisComponents {
            rp,
            mach,
            rrna_proxy,
            stress,
            metab,
            luxury,
            rqc_mach,
            chaperone,
            ubiquitin,
        });

        axes.push(CellAxes {
            tl,
            st,
            rqc,
            rqc_pressure,
            tpc,
            st_low_confidence,
        });
    }

    let tl_median = median_non_nan(axes.iter().map(|a| a.tl));
    let st_median = median_non_nan(axes.iter().map(|a| a.st));
    let rqc_median = median_non_nan(axes.iter().map(|a| a.rqc));

    info!(
        bg_gene_count = bg_gene_ids.len(),
        st_low_conf_pct = if n_cells == 0 {
            0.0
        } else {
            (st_low_conf_count as f64) / (n_cells as f64)
        },
        tl_median = tl_median,
        st_median = st_median,
        rqc_median = rqc_median,
        "Stage 2 summary"
    );
    info!("Stage 2 end");

    Ok(Stage2Output {
        bg_gene_ids,
        libsize,
        detected_genes,
        components,
        axes,
    })
}

#[derive(Clone, Copy)]
struct ColumnView<'a> {
    row_idx: &'a [u32],
    values: &'a [u32],
}

impl<'a> ColumnView<'a> {
    fn iter(self, row_to_gene: &'a [u32]) -> impl Iterator<Item = (u32, u32)> + 'a {
        self.row_idx
            .iter()
            .zip(self.values.iter())
            .map(move |(row, val)| (row_to_gene[*row as usize], *val))
    }
}

enum MatrixSource<'a> {
    Matrix {
        col_ptr: &'a [u32],
        row_idx: &'a [u32],
        values: &'a [u32],
        n_cols: usize,
    },
    Cache(&'a SharedCacheData),
}

impl<'a> MatrixSource<'a> {
    fn from_input(input: &'a InputBundle) -> Self {
        if let Some(cache) = input.shared_cache.as_ref() {
            Self::Cache(cache)
        } else {
            Self::Matrix {
                col_ptr: &input.matrix.col_ptr,
                row_idx: &input.matrix.row_idx,
                values: &input.matrix.values,
                n_cols: input.matrix.n_cols as usize,
            }
        }
    }

    fn n_cols(&self) -> usize {
        match self {
            MatrixSource::Matrix { n_cols, .. } => *n_cols,
            MatrixSource::Cache(cache) => cache.n_cells as usize,
        }
    }

    fn column_view(&self, col: usize) -> ColumnView<'_> {
        match self {
            MatrixSource::Matrix {
                col_ptr,
                row_idx,
                values,
                ..
            } => {
                let start = col_ptr[col] as usize;
                let end = col_ptr[col + 1] as usize;
                ColumnView {
                    row_idx: &row_idx[start..end],
                    values: &values[start..end],
                }
            }
            MatrixSource::Cache(cache) => {
                let start = cache.col_ptr[col] as usize;
                let end = cache.col_ptr[col + 1] as usize;
                ColumnView {
                    row_idx: &cache.row_idx[start..end],
                    values: &cache.values_u32[start..end],
                }
            }
        }
    }
}

#[derive(Default)]
struct Acc {
    sum: f64,
    count: u32,
}

impl Acc {
    fn add(&mut self, x: f64) {
        self.sum += x;
        self.count += 1;
    }

    fn mean(&self) -> f64 {
        if self.count == 0 {
            f64::NAN
        } else {
            self.sum / (self.count as f64)
        }
    }
}

fn enrich(set_mean: f64, bg_mean: f64) -> f64 {
    if set_mean.is_nan() || bg_mean.is_nan() {
        return f64::NAN;
    }
    let raw = set_mean - bg_mean;
    let scaled = (raw + OFFSET) / SCALE;
    clamp01(scaled)
}

fn compute_tl(rp: f64, mach: f64, rrna_proxy: f64) -> (f64, bool) {
    let mut missing_component = false;
    let rp_val = if rp.is_nan() {
        missing_component = true;
        0.0
    } else {
        rp
    };
    let mach_val = if mach.is_nan() {
        missing_component = true;
        0.0
    } else {
        mach
    };
    let rrna_val = if rrna_proxy.is_nan() {
        missing_component = true;
        0.0
    } else {
        rrna_proxy
    };
    let tl = clamp01(0.5 * rp_val + 0.3 * mach_val + 0.2 * rrna_val);
    (tl, missing_component)
}

fn compute_st(stress: f64, metab: f64, luxury: Option<f64>) -> (f64, bool) {
    match luxury {
        Some(lux) => {
            let lux_supp = clamp01(1.0 - lux);
            let st = clamp01(0.45 * stress + 0.35 * metab + 0.20 * lux_supp);
            (st, false)
        }
        None => {
            let st = clamp01(0.55 * stress + 0.45 * metab);
            (st, true)
        }
    }
}

fn compute_rqc(rqc_mach: f64, stress: f64, rp: f64) -> (f64, f64) {
    let rqc = clamp01(0.7 * rqc_mach + 0.3 * stress);
    let rqc_pressure = clamp01(rqc_mach - rp);
    (rqc, rqc_pressure)
}

fn compute_tpc(chaperone: f64, ubiquitin: f64, tl: f64) -> f64 {
    let deg_cap = clamp01(0.6 * chaperone + 0.4 * ubiquitin);
    clamp_to_minus1_plus1(deg_cap - tl)
}

fn compute_cell_stats(source: &MatrixSource<'_>) -> (Vec<u64>, Vec<u32>) {
    let n_cells = source.n_cols();
    let mut libsize = vec![0u64; n_cells];
    let mut detected = vec![0u32; n_cells];

    for col in 0..n_cells {
        let view = source.column_view(col);
        let mut sum = 0u64;
        for val in view.values {
            sum += *val as u64;
        }
        libsize[col] = sum;
        detected[col] = view.values.len() as u32;
    }

    (libsize, detected)
}

fn compute_gene_detection_counts(
    input: &InputBundle,
    source: &MatrixSource<'_>,
    n_genes_unique: usize,
) -> Vec<u32> {
    let n_cells = source.n_cols();
    let mut counts = vec![0u32; n_genes_unique];
    let mut last_seen = vec![u32::MAX; n_genes_unique];

    for col in 0..n_cells {
        let view = source.column_view(col);
        for row in view.row_idx {
            let gene_id = input.gene_index.row_to_gene[*row as usize] as usize;
            if last_seen[gene_id] != col as u32 {
                last_seen[gene_id] = col as u32;
                counts[gene_id] += 1;
            }
        }
    }

    counts
}

fn select_background_genes(input: &InputBundle, counts: &[u32], n_cells: u32) -> Vec<u32> {
    let mut candidates: Vec<(u32, u32)> = Vec::new();

    for (gene_id, gene) in input.gene_index.genes.iter().enumerate() {
        let symbol = gene.symbol.as_str();
        if is_ribosomal(symbol) || is_mito(symbol) {
            continue;
        }
        let count = counts[gene_id];
        if n_cells == 0 {
            continue;
        }
        candidates.push((gene_id as u32, count));
    }

    candidates.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));

    let mut bg_gene_ids: Vec<u32> = candidates
        .iter()
        .take(BG_TOP_K)
        .map(|(gid, _)| *gid)
        .collect();

    bg_gene_ids.sort_unstable();
    bg_gene_ids
}

fn collect_ribosomal_gene_ids(input: &InputBundle) -> Vec<u32> {
    let mut ids = Vec::new();
    for entry in &input.gene_index.genes {
        if is_ribosomal(&entry.symbol) {
            ids.push(entry.gene_id);
        }
    }
    ids
}

fn median_non_nan<I>(iter: I) -> f64
where
    I: Iterator<Item = f64>,
{
    let mut values: Vec<f64> = iter.filter(|v| !v.is_nan()).collect();
    if values.is_empty() {
        return f64::NAN;
    }
    values.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mid = values.len() / 2;
    if values.len() % 2 == 1 {
        values[mid]
    } else {
        (values[mid - 1] + values[mid]) / 2.0
    }
}
