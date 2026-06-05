use thiserror::Error;

pub mod cache;
pub mod detect;
pub mod gene_index;
pub mod matrix_source;
pub mod metadata;
pub mod mtx10x;
pub mod types;

pub use cache::{SharedCacheData, read_shared_cache};
pub use detect::{detect_input, detect_input_for_prefix, detect_prefix, resolve_cache_path};
pub use gene_index::{
    DuplicateGene, FeatureRow, GeneEntry, GeneIndex, build_gene_index, normalize_symbol,
    resolve_gene_ids,
};
pub use matrix_source::{ColumnView, MatrixSource};
pub use metadata::{CellMeta, MetadataTable, load_metadata};
pub use mtx10x::{CscMatrix, MtxDataset, load_mtx_dataset};
pub use types::{InputBundle, InputFormat, Stage1Stats};

#[derive(Error, Debug)]
pub enum InputError {
    #[error("io error: {0}")]
    Io(#[from] std::io::Error),
    #[error("invalid input format: {0}")]
    InvalidFormat(String),
    #[error("missing file: {0}")]
    MissingFile(String),
    #[error("parse error: {0}")]
    Parse(String),
    #[error("dimension mismatch: {0}")]
    Dimension(String),
    #[error("metadata error: {0}")]
    Metadata(String),
}
