use rayon::prelude::*;
use tracing::info;

use crate::core::{Acc, PanelBit, PanelMask, PanelMaskBuilder, median_non_nan};
use crate::input::{InputBundle, MatrixSource};
use crate::model::axes::{
    BG_TOP_K, CHAPERONE, LUXURY, MACHINERY, METAB, OFFSET, RQC_MACH, RRNA_PROXY, SCALE, STRESS,
    UBIQUITIN, clamp_to_minus1_plus1, clamp01, is_mito, is_ribosomal,
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

struct StageBits {
    rp: PanelBit,
    mach: PanelBit,
    rrna: PanelBit,
    stress: PanelBit,
    metab: PanelBit,
    luxury: PanelBit,
    rqc: PanelBit,
    chaperone: PanelBit,
    ubiquitin: PanelBit,
    bg: PanelBit,
}

pub fn run_stage2(input: &InputBundle) -> anyhow::Result<Stage2Output> {
    let _span = tracing::info_span!("stage2_axes").entered();
    info!("Stage 2 start");

    let source = MatrixSource::from_input(input);
    let n_cells = source.n_cols();
    let n_genes_unique = input.gene_index.genes.len();

    let (libsize, detected_genes) = compute_cell_stats(&source, n_cells);
    let gene_detect_counts = compute_gene_detection_counts(input, &source, n_cells, n_genes_unique);
    let bg_gene_ids = select_background_genes(input, &gene_detect_counts);

    let (mask, bits) = build_stage_mask(input, &bg_gene_ids, n_genes_unique);
    let row_to_gene = input.gene_index.row_to_gene.as_slice();

    let results: Vec<(CellAxisComponents, CellAxes)> = (0..n_cells)
        .into_par_iter()
        .map(|col| compute_cell(&source, col, libsize[col], row_to_gene, &mask, &bits))
        .collect();

    let mut components = Vec::with_capacity(n_cells);
    let mut axes = Vec::with_capacity(n_cells);
    let mut st_low_conf_count: u32 = 0;
    for (comp, ax) in results {
        if ax.st_low_confidence {
            st_low_conf_count += 1;
        }
        components.push(comp);
        axes.push(ax);
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

fn build_stage_mask(
    input: &InputBundle,
    bg_gene_ids: &[u32],
    n_genes_unique: usize,
) -> (PanelMask, StageBits) {
    let rp_ids: Vec<u32> = input
        .gene_index
        .genes
        .iter()
        .filter(|g| is_ribosomal(&g.symbol))
        .map(|g| g.gene_id)
        .collect();

    let map = &input.gene_index.map;
    let resolve = |symbols: &[&str]| -> Vec<u32> {
        symbols
            .iter()
            .filter_map(|s| map.get(*s).copied())
            .collect()
    };

    let mut builder = PanelMaskBuilder::new(n_genes_unique);
    let bits = StageBits {
        rp: builder.add_panel(&rp_ids),
        mach: builder.add_panel(&resolve(MACHINERY)),
        rrna: builder.add_panel(&resolve(RRNA_PROXY)),
        stress: builder.add_panel(&resolve(STRESS)),
        metab: builder.add_panel(&resolve(METAB)),
        luxury: builder.add_panel(&resolve(LUXURY)),
        rqc: builder.add_panel(&resolve(RQC_MACH)),
        chaperone: builder.add_panel(&resolve(CHAPERONE)),
        ubiquitin: builder.add_panel(&resolve(UBIQUITIN)),
        bg: builder.add_panel(bg_gene_ids),
    };
    (builder.build(), bits)
}

fn compute_cell(
    source: &MatrixSource<'_>,
    col: usize,
    libsize: u64,
    row_to_gene: &[u32],
    mask: &PanelMask,
    bits: &StageBits,
) -> (CellAxisComponents, CellAxes) {
    let col_view = source.column(col);
    let denom = if libsize == 0 { 1.0 } else { libsize as f64 };

    let mut accs = [Acc::default(); 10];

    for (gid, val) in col_view.iter_genes(row_to_gene) {
        let m = mask.get(gid);
        if m == 0 {
            continue;
        }
        let cpm = 1_000_000.0 * (val as f64) / denom;
        let x = ln1p_f64(cpm);
        let mut bits_left = m;
        while bits_left != 0 {
            let lo = bits_left.trailing_zeros() as usize;
            accs[lo].add(x);
            bits_left &= bits_left - 1;
        }
    }

    let bg_mean = accs[bits.bg.index].mean();
    let rp = enrich(accs[bits.rp.index].mean(), bg_mean);
    let mach = enrich(accs[bits.mach.index].mean(), bg_mean);
    let rrna_proxy = enrich(accs[bits.rrna.index].mean(), bg_mean);
    let stress = enrich(accs[bits.stress.index].mean(), bg_mean);
    let metab = enrich(accs[bits.metab.index].mean(), bg_mean);
    let luxury_raw = enrich(accs[bits.luxury.index].mean(), bg_mean);
    let rqc_mach = enrich(accs[bits.rqc.index].mean(), bg_mean);
    let chaperone = enrich(accs[bits.chaperone.index].mean(), bg_mean);
    let ubiquitin = enrich(accs[bits.ubiquitin.index].mean(), bg_mean);

    let luxury = if luxury_raw.is_nan() {
        None
    } else {
        Some(luxury_raw)
    };

    let (tl, _missing_component) = compute_tl(rp, mach, rrna_proxy);
    let (st, st_low_confidence) = compute_st(stress, metab, luxury);
    let (rqc, rqc_pressure) = compute_rqc(rqc_mach, stress, rp);
    let tpc = compute_tpc(chaperone, ubiquitin, tl);

    (
        CellAxisComponents {
            rp,
            mach,
            rrna_proxy,
            stress,
            metab,
            luxury,
            rqc_mach,
            chaperone,
            ubiquitin,
        },
        CellAxes {
            tl,
            st,
            rqc,
            rqc_pressure,
            tpc,
            st_low_confidence,
        },
    )
}

#[inline]
fn enrich(set_mean: f64, bg_mean: f64) -> f64 {
    if set_mean.is_nan() || bg_mean.is_nan() {
        return f64::NAN;
    }
    let raw = set_mean - bg_mean;
    let scaled = (raw + OFFSET) / SCALE;
    clamp01(scaled)
}

fn compute_tl(rp: f64, mach: f64, rrna_proxy: f64) -> (f64, bool) {
    let mut missing = false;
    let rp_v = nan_to_zero(rp, &mut missing);
    let mach_v = nan_to_zero(mach, &mut missing);
    let rrna_v = nan_to_zero(rrna_proxy, &mut missing);
    let tl = clamp01(0.5 * rp_v + 0.3 * mach_v + 0.2 * rrna_v);
    (tl, missing)
}

#[inline]
fn nan_to_zero(x: f64, missing: &mut bool) -> f64 {
    if x.is_nan() {
        *missing = true;
        0.0
    } else {
        x
    }
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

fn compute_cell_stats(source: &MatrixSource<'_>, n_cells: usize) -> (Vec<u64>, Vec<u32>) {
    let stats: Vec<(u64, u32)> = (0..n_cells)
        .into_par_iter()
        .map(|col| {
            let view = source.column(col);
            let sum: u64 = view.values.iter().map(|v| *v as u64).sum();
            (sum, view.values.len() as u32)
        })
        .collect();

    let mut libsize = Vec::with_capacity(n_cells);
    let mut detected = Vec::with_capacity(n_cells);
    for (s, d) in stats {
        libsize.push(s);
        detected.push(d);
    }
    (libsize, detected)
}

fn compute_gene_detection_counts(
    input: &InputBundle,
    source: &MatrixSource<'_>,
    n_cells: usize,
    n_genes_unique: usize,
) -> Vec<u32> {
    let row_to_gene = input.gene_index.row_to_gene.as_slice();
    let has_dups = !input.gene_index.duplicates.is_empty();

    let mut counts = vec![0u32; n_genes_unique];

    if !has_dups {
        for col in 0..n_cells {
            let view = source.column(col);
            for row in view.row_idx {
                counts[*row as usize] += 1;
            }
        }
        return counts;
    }

    let mut last_seen = vec![u32::MAX; n_genes_unique];
    for col in 0..n_cells {
        let col_marker = col as u32;
        let view = source.column(col);
        for row in view.row_idx {
            let gene_id = row_to_gene[*row as usize] as usize;
            if last_seen[gene_id] != col_marker {
                last_seen[gene_id] = col_marker;
                counts[gene_id] += 1;
            }
        }
    }
    counts
}

fn select_background_genes(input: &InputBundle, counts: &[u32]) -> Vec<u32> {
    let mut candidates: Vec<(u32, u32)> = Vec::with_capacity(input.gene_index.genes.len());
    for (gene_id, gene) in input.gene_index.genes.iter().enumerate() {
        if is_ribosomal(&gene.symbol) || is_mito(&gene.symbol) {
            continue;
        }
        candidates.push((gene_id as u32, counts[gene_id]));
    }

    let take = BG_TOP_K.min(candidates.len());
    if take == 0 {
        return Vec::new();
    }
    if take < candidates.len() {
        candidates.select_nth_unstable_by(take - 1, |a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
    }
    let mut bg_gene_ids: Vec<u32> = candidates.iter().take(take).map(|(gid, _)| *gid).collect();
    bg_gene_ids.sort_unstable();
    bg_gene_ids
}
