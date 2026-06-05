use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::Result;
use rayon::prelude::*;

use crate::core::math::percentile_sorted;
use crate::input::{InputBundle, MatrixSource};
use crate::model::axes::{
    CHAPERONE, LUXURY, MACHINERY, METAB, RQC_MACH, RRNA_PROXY, STRESS, UBIQUITIN, is_ribosomal,
};

struct PanelDef {
    panel_id: &'static str,
    panel_name: &'static str,
    panel_group: &'static str,
    symbols_defined: &'static [&'static str],
    ribosomal_pattern: bool,
}

const PANEL_DEFS: &[PanelDef] = &[
    PanelDef {
        panel_id: "RP",
        panel_name: "Ribosomal proteins",
        panel_group: "translation_core",
        symbols_defined: &[],
        ribosomal_pattern: true,
    },
    PanelDef {
        panel_id: "MACH",
        panel_name: "Translation machinery",
        panel_group: "translation_core",
        symbols_defined: MACHINERY,
        ribosomal_pattern: false,
    },
    PanelDef {
        panel_id: "RRNA",
        panel_name: "rRNA processing",
        panel_group: "translation_core",
        symbols_defined: RRNA_PROXY,
        ribosomal_pattern: false,
    },
    PanelDef {
        panel_id: "STRESS",
        panel_name: "Stress translation",
        panel_group: "stress",
        symbols_defined: STRESS,
        ribosomal_pattern: false,
    },
    PanelDef {
        panel_id: "METAB",
        panel_name: "Metabolic survival",
        panel_group: "stress",
        symbols_defined: METAB,
        ribosomal_pattern: false,
    },
    PanelDef {
        panel_id: "LUXURY",
        panel_name: "Luxury anabolic",
        panel_group: "growth",
        symbols_defined: LUXURY,
        ribosomal_pattern: false,
    },
    PanelDef {
        panel_id: "RQC",
        panel_name: "RQC machinery",
        panel_group: "rqc",
        symbols_defined: RQC_MACH,
        ribosomal_pattern: false,
    },
    PanelDef {
        panel_id: "CHAPERONE",
        panel_name: "Chaperone",
        panel_group: "proteostasis",
        symbols_defined: CHAPERONE,
        ribosomal_pattern: false,
    },
    PanelDef {
        panel_id: "UBIQUITIN",
        panel_name: "Ubiquitin stress",
        panel_group: "proteostasis",
        symbols_defined: UBIQUITIN,
        ribosomal_pattern: false,
    },
];

struct PanelStats {
    panel_size_defined: usize,
    panel_size_mappable: usize,
    missing_genes: String,
    coverage_values: Vec<f64>,
    sum_values: Vec<f64>,
}

pub fn write_panels_report(path: &Path, input: &InputBundle) -> Result<()> {
    let file = File::create(path)?;
    let mut w = BufWriter::with_capacity(1 << 20, file);

    writeln!(
        w,
        "panel_id\tpanel_name\tpanel_group\tpanel_size_defined\tpanel_size_mappable\tmissing_genes\tcoverage_median\tcoverage_p10\tsum_median\tsum_p90\tsum_p99"
    )?;

    let source = MatrixSource::from_input(input);
    let n_cells = source.n_cols();
    let n_genes = input.gene_index.genes.len();
    let row_to_gene = input.gene_index.row_to_gene.as_slice();

    let n_panels = PANEL_DEFS.len();
    assert!(n_panels <= 16);

    let mut panel_stats: Vec<PanelStats> = Vec::with_capacity(n_panels);
    let mut gene_panel_mask = vec![0u16; n_genes];

    for (panel_idx, def) in PANEL_DEFS.iter().enumerate() {
        let (gene_ids, panel_size_defined, missing_genes) = resolve_panel(input, def);
        for gid in &gene_ids {
            let idx = *gid as usize;
            if idx < gene_panel_mask.len() {
                gene_panel_mask[idx] |= 1u16 << panel_idx;
            }
        }
        panel_stats.push(PanelStats {
            panel_size_defined,
            panel_size_mappable: gene_ids.len(),
            missing_genes,
            coverage_values: vec![0.0; n_cells],
            sum_values: vec![0.0; n_cells],
        });
    }

    let per_cell: Vec<Vec<(u32, f64)>> = (0..n_cells)
        .into_par_iter()
        .map(|col| {
            let view = source.column(col);
            let mut hits = vec![(0u32, 0.0f64); n_panels];
            for (k, &row) in view.row_idx.iter().enumerate() {
                let gid = row_to_gene[row as usize] as usize;
                let mask = *gene_panel_mask.get(gid).unwrap_or(&0);
                if mask == 0 {
                    continue;
                }
                let val = view.values[k] as f64;
                let mut bits = mask as u32;
                while bits != 0 {
                    let lo = bits.trailing_zeros() as usize;
                    hits[lo].0 += 1;
                    hits[lo].1 += val;
                    bits &= bits - 1;
                }
            }
            hits
        })
        .collect();

    for (col, hits) in per_cell.iter().enumerate() {
        for (panel_idx, (count, sum)) in hits.iter().enumerate() {
            let mappable = panel_stats[panel_idx].panel_size_mappable;
            panel_stats[panel_idx].coverage_values[col] = if mappable == 0 {
                0.0
            } else {
                *count as f64 / mappable as f64
            };
            panel_stats[panel_idx].sum_values[col] = *sum;
        }
    }

    let mut sort_buf = Vec::<f64>::new();
    let mut f_buf = ryu::Buffer::new();
    for (panel_idx, def) in PANEL_DEFS.iter().enumerate() {
        let stats = &panel_stats[panel_idx];

        let cov_med = sorted_percentile(&mut sort_buf, &stats.coverage_values, 0.50);
        let cov_p10 = sorted_percentile(&mut sort_buf, &stats.coverage_values, 0.10);
        let sum_med = sorted_percentile(&mut sort_buf, &stats.sum_values, 0.50);
        let sum_p90 = sorted_percentile(&mut sort_buf, &stats.sum_values, 0.90);
        let sum_p99 = sorted_percentile(&mut sort_buf, &stats.sum_values, 0.99);

        w.write_all(def.panel_id.as_bytes())?;
        w.write_all(b"\t")?;
        w.write_all(def.panel_name.as_bytes())?;
        w.write_all(b"\t")?;
        w.write_all(def.panel_group.as_bytes())?;
        w.write_all(b"\t")?;
        write_usize(&mut w, stats.panel_size_defined)?;
        w.write_all(b"\t")?;
        write_usize(&mut w, stats.panel_size_mappable)?;
        w.write_all(b"\t")?;
        w.write_all(stats.missing_genes.as_bytes())?;
        w.write_all(b"\t")?;
        write_f64_6(&mut w, &mut f_buf, cov_med)?;
        w.write_all(b"\t")?;
        write_f64_6(&mut w, &mut f_buf, cov_p10)?;
        w.write_all(b"\t")?;
        write_f64_6(&mut w, &mut f_buf, sum_med)?;
        w.write_all(b"\t")?;
        write_f64_6(&mut w, &mut f_buf, sum_p90)?;
        w.write_all(b"\t")?;
        write_f64_6(&mut w, &mut f_buf, sum_p99)?;
        w.write_all(b"\n")?;
    }

    Ok(())
}

fn resolve_panel(input: &InputBundle, def: &PanelDef) -> (Vec<u32>, usize, String) {
    if def.ribosomal_pattern {
        let mut ids = Vec::new();
        for g in &input.gene_index.genes {
            if is_ribosomal(&g.symbol) {
                ids.push(g.gene_id);
            }
        }
        let len = ids.len();
        return (ids, len, String::new());
    }

    let mut ids = Vec::new();
    let mut missing = Vec::new();
    for symbol in def.symbols_defined {
        if let Some(gid) = input.gene_index.map.get(*symbol) {
            ids.push(*gid);
        } else {
            missing.push(*symbol);
        }
    }
    let missing_str = if missing.is_empty() {
        String::new()
    } else {
        missing.join(",")
    };
    (ids, def.symbols_defined.len(), missing_str)
}

fn sorted_percentile(buf: &mut Vec<f64>, values: &[f64], p: f64) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    buf.clear();
    buf.extend_from_slice(values);
    buf.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    percentile_sorted(buf, p)
}

fn write_usize<W: Write>(w: &mut W, n: usize) -> std::io::Result<()> {
    let mut buf = itoa::Buffer::new();
    w.write_all(buf.format(n).as_bytes())
}

fn write_f64_6<W: Write>(w: &mut W, buf: &mut ryu::Buffer, v: f64) -> std::io::Result<()> {
    if v.is_nan() {
        return w.write_all(b"NaN");
    }
    let rounded = (v * 1_000_000.0).round() / 1_000_000.0;
    let s = buf.format(rounded);
    write_fixed6(w, s)
}

fn write_fixed6<W: Write>(w: &mut W, s: &str) -> std::io::Result<()> {
    if let Some(dot) = s.find('.') {
        let frac_len = s.len() - dot - 1;
        if frac_len >= 6 {
            return w.write_all(&s.as_bytes()[..dot + 7]);
        }
        w.write_all(s.as_bytes())?;
        for _ in frac_len..6 {
            w.write_all(b"0")?;
        }
        Ok(())
    } else {
        w.write_all(s.as_bytes())?;
        w.write_all(b".000000")
    }
}
