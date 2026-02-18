use std::path::PathBuf;

use super::cache::SharedCacheData;
use super::gene_index::{FeatureRow, GeneIndex};
use super::metadata::MetadataTable;
use super::mtx10x::CscMatrix;

#[derive(Debug, Clone)]
pub enum InputFormat {
    TenXDir {
        matrix_path: PathBuf,
        barcodes_path: PathBuf,
        features_path: PathBuf,
    },
    SharedCache {
        cache_path: PathBuf,
    },
}

#[derive(Debug, Clone)]
pub struct InputBundle {
    pub format: InputFormat,
    pub matrix: CscMatrix,
    pub barcodes: Vec<String>,
    pub features: Vec<FeatureRow>,
    pub gene_index: GeneIndex,
    pub metadata: Option<MetadataTable>,
    pub shared_cache: Option<SharedCacheData>,
}

#[derive(Debug, Clone, Copy)]
pub struct Stage1Stats {
    pub n_genes_raw: u32,
    pub n_genes_unique: u32,
    pub n_cells: u32,
    pub nnz: u32,
}
