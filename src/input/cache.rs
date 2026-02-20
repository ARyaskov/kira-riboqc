use std::fs::File;
use std::path::Path;
use std::sync::Arc;

use memmap2::Mmap;

use super::InputError;

#[derive(Debug, Clone)]
pub struct SharedCacheData {
    pub path: std::path::PathBuf,
    pub n_genes: u64,
    pub n_cells: u64,
    pub nnz: u64,
    pub genes: Vec<String>,
    pub barcodes: Vec<String>,
    pub col_ptr: Vec<u64>,
    pub row_idx: Vec<u32>,
    pub values_u32: Vec<u32>,
    pub mmap: Arc<Mmap>,
}

pub fn read_shared_cache(path: &Path) -> Result<SharedCacheData, InputError> {
    let shared = kira_shared_sc_cache::read_shared_cache_owned(path).map_err(map_err)?;

    let file = File::open(path)?;
    let mmap = unsafe { Mmap::map(&file) }
        .map_err(|e| InputError::Parse(format!("failed to mmap cache file: {e}")))?;

    Ok(SharedCacheData {
        path: shared.path,
        n_genes: shared.n_genes,
        n_cells: shared.n_cells,
        nnz: shared.nnz,
        genes: shared.genes,
        barcodes: shared.barcodes,
        col_ptr: shared.col_ptr,
        row_idx: shared.row_idx,
        values_u32: shared.values_u32,
        mmap: Arc::new(mmap),
    })
}

fn map_err(err: kira_shared_sc_cache::SharedCacheError) -> InputError {
    match err {
        kira_shared_sc_cache::SharedCacheError::Io { source, .. } => InputError::Io(source),
        kira_shared_sc_cache::SharedCacheError::Format { message, .. } => {
            InputError::Parse(message)
        }
    }
}
