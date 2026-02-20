use std::collections::BTreeMap;
use std::path::Path;

use kira_scio::api::{Reader, ReaderOptions};
use kira_scio::detect::DetectedFormat;

use super::InputError;

#[derive(Debug, Clone)]
pub struct FeatureRow {
    pub raw_id: String,
    pub raw_name: String,
    pub raw_type: String,
    pub norm_symbol: String,
}

#[derive(Debug, Clone)]
pub struct GeneIndex {
    pub genes: Vec<GeneEntry>,
    pub map: BTreeMap<String, u32>,
    pub row_to_gene: Vec<u32>,
    pub duplicates: Vec<DuplicateGene>,
}

#[derive(Debug, Clone)]
pub struct GeneEntry {
    pub gene_id: u32,
    pub symbol: String,
    pub first_row: u32,
}

#[derive(Debug, Clone)]
pub struct DuplicateGene {
    pub symbol: String,
    pub first_row: u32,
    pub dup_row: u32,
}

pub fn normalize_symbol(s: &str) -> String {
    let collapsed = s.trim().split_whitespace().collect::<Vec<_>>().join(" ");
    collapsed.to_ascii_uppercase()
}

pub fn load_features(path: &Path) -> Result<Vec<FeatureRow>, InputError> {
    let md = Reader::with_options(
        path,
        ReaderOptions {
            force_format: Some(DetectedFormat::Mtx10x),
            strict: true,
        },
    )
    .read_metadata()
    .map_err(|e| InputError::Parse(e.message))?;

    let mut rows = Vec::with_capacity(md.gene_symbols.len());
    for (idx, symbol) in md.gene_symbols.iter().enumerate() {
        let raw_id = md
            .gene_ids
            .get(idx)
            .cloned()
            .unwrap_or_else(|| format!("GENE_{}", idx + 1));
        rows.push(FeatureRow {
            raw_id,
            raw_name: symbol.clone(),
            raw_type: String::new(),
            norm_symbol: normalize_symbol(symbol),
        });
    }
    Ok(rows)
}

pub fn build_gene_index(features: &[FeatureRow]) -> GeneIndex {
    let mut genes: Vec<GeneEntry> = Vec::new();
    let mut map: BTreeMap<String, u32> = BTreeMap::new();
    let mut row_to_gene = Vec::with_capacity(features.len());
    let mut duplicates = Vec::new();

    for (i, row) in features.iter().enumerate() {
        let row_idx = i as u32;
        if let Some(gene_id) = map.get(&row.norm_symbol) {
            let first_row = genes[*gene_id as usize].first_row;
            row_to_gene.push(*gene_id);
            duplicates.push(DuplicateGene {
                symbol: row.norm_symbol.clone(),
                first_row,
                dup_row: row_idx,
            });
        } else {
            let gene_id = genes.len() as u32;
            map.insert(row.norm_symbol.clone(), gene_id);
            genes.push(GeneEntry {
                gene_id,
                symbol: row.norm_symbol.clone(),
                first_row: row_idx,
            });
            row_to_gene.push(gene_id);
        }
    }

    GeneIndex {
        genes,
        map,
        row_to_gene,
        duplicates,
    }
}
