use std::collections::BTreeMap;
use std::io::{BufRead, BufReader};
use std::path::Path;

use super::InputError;
use super::mtx10x::open_maybe_gz;

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
    let reader = BufReader::new(open_maybe_gz(path)?);
    let mut rows = Vec::new();

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim_end();
        if trimmed.is_empty() {
            continue;
        }
        let parts: Vec<&str> = trimmed.split('\t').collect();
        if parts.len() < 2 {
            return Err(InputError::Parse(
                "features file must have at least 2 columns".to_string(),
            ));
        }
        let raw_id = parts[0].to_string();
        let raw_name = parts[1].to_string();
        let raw_type = if parts.len() >= 3 {
            parts[2].to_string()
        } else {
            String::new()
        };
        let norm_symbol = normalize_symbol(&raw_name);

        rows.push(FeatureRow {
            raw_id,
            raw_name,
            raw_type,
            norm_symbol,
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
