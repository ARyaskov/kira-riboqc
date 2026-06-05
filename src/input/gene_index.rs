use rustc_hash::FxHashMap;

#[derive(Debug, Clone)]
pub struct FeatureRow {
    pub raw_id: String,
    pub raw_name: String,
    pub raw_type: String,
    pub norm_symbol: String,
}

impl FeatureRow {
    pub fn from_raw(raw_id: String, raw_name: String) -> Self {
        let norm_symbol = normalize_symbol(&raw_name);
        Self {
            raw_id,
            raw_name,
            raw_type: String::new(),
            norm_symbol,
        }
    }
}

#[derive(Debug, Clone)]
pub struct GeneIndex {
    pub genes: Vec<GeneEntry>,
    pub map: FxHashMap<String, u32>,
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
    let trimmed = s.trim();
    let needs_collapse = trimmed
        .as_bytes()
        .windows(2)
        .any(|w| w[0] == b' ' && w[1] == b' ');
    let collapsed = if needs_collapse {
        trimmed.split_whitespace().collect::<Vec<_>>().join(" ")
    } else {
        trimmed.to_string()
    };
    collapsed.to_ascii_uppercase()
}

pub fn build_gene_index(features: &[FeatureRow]) -> GeneIndex {
    let n = features.len();
    let mut genes: Vec<GeneEntry> = Vec::with_capacity(n);
    let mut map: FxHashMap<String, u32> =
        FxHashMap::with_capacity_and_hasher(n, Default::default());
    let mut row_to_gene = Vec::with_capacity(n);
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

pub fn resolve_gene_ids(map: &FxHashMap<String, u32>, symbols: &[&str]) -> Vec<u32> {
    symbols
        .iter()
        .filter_map(|s| map.get(*s).copied())
        .collect()
}
