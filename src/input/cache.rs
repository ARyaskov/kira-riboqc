use std::path::{Path, PathBuf};
use std::sync::Arc;

use kira_shared_sc_cache::SharedCacheMmap;

use super::InputError;

pub struct SharedCacheData {
    pub path: PathBuf,
    inner: Arc<SharedCacheMmap>,
}

impl SharedCacheData {
    #[inline]
    pub fn n_genes(&self) -> usize {
        self.inner.n_genes
    }

    #[inline]
    pub fn n_cells(&self) -> usize {
        self.inner.n_cells
    }

    #[inline]
    pub fn nnz(&self) -> usize {
        self.inner.nnz
    }

    #[inline]
    pub fn genes(&self) -> &[String] {
        &self.inner.genes
    }

    #[inline]
    pub fn barcodes(&self) -> &[String] {
        &self.inner.barcodes
    }

    #[inline]
    pub fn col_ptr(&self) -> &[u64] {
        self.inner.col_ptr()
    }

    #[inline]
    pub fn row_idx(&self) -> &[u32] {
        self.inner.row_idx()
    }

    #[inline]
    pub fn values_u32(&self) -> &[u32] {
        self.inner.values_u32()
    }
}

impl std::fmt::Debug for SharedCacheData {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("SharedCacheData")
            .field("path", &self.path)
            .field("n_genes", &self.n_genes())
            .field("n_cells", &self.n_cells())
            .field("nnz", &self.nnz())
            .finish()
    }
}

impl Clone for SharedCacheData {
    fn clone(&self) -> Self {
        Self {
            path: self.path.clone(),
            inner: Arc::clone(&self.inner),
        }
    }
}

pub fn read_shared_cache(path: &Path) -> Result<SharedCacheData, InputError> {
    let mapped = kira_shared_sc_cache::mmap_shared_cache(path).map_err(map_err)?;
    Ok(SharedCacheData {
        path: path.to_path_buf(),
        inner: Arc::new(mapped),
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
