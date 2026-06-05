use super::cache::SharedCacheData;
use super::mtx10x::CscMatrix;
use super::types::InputBundle;

pub enum MatrixSource<'a> {
    Owned {
        col_ptr: &'a [u32],
        row_idx: &'a [u32],
        values: &'a [u32],
        n_cols: usize,
    },
    Cache(&'a SharedCacheData),
}

impl<'a> MatrixSource<'a> {
    pub fn from_input(input: &'a InputBundle) -> Self {
        if let Some(cache) = input.shared_cache.as_ref() {
            Self::Cache(cache)
        } else {
            Self::from_matrix(&input.matrix)
        }
    }

    pub fn from_matrix(matrix: &'a CscMatrix) -> Self {
        Self::Owned {
            col_ptr: &matrix.col_ptr,
            row_idx: &matrix.row_idx,
            values: &matrix.values,
            n_cols: matrix.n_cols as usize,
        }
    }

    pub fn from_cache(cache: &'a SharedCacheData) -> Self {
        Self::Cache(cache)
    }

    #[inline]
    pub fn n_cols(&self) -> usize {
        match self {
            Self::Owned { n_cols, .. } => *n_cols,
            Self::Cache(cache) => cache.n_cells(),
        }
    }

    #[inline]
    pub fn column(&self, col: usize) -> ColumnView<'_> {
        match self {
            Self::Owned {
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
            Self::Cache(cache) => {
                let col_ptr = cache.col_ptr();
                let start = col_ptr[col] as usize;
                let end = col_ptr[col + 1] as usize;
                ColumnView {
                    row_idx: &cache.row_idx()[start..end],
                    values: &cache.values_u32()[start..end],
                }
            }
        }
    }
}

#[derive(Clone, Copy)]
pub struct ColumnView<'a> {
    pub row_idx: &'a [u32],
    pub values: &'a [u32],
}

impl<'a> ColumnView<'a> {
    #[inline]
    pub fn iter_genes(self, row_to_gene: &'a [u32]) -> impl Iterator<Item = (u32, u32)> + 'a {
        self.row_idx
            .iter()
            .zip(self.values.iter())
            .map(move |(row, val)| (row_to_gene[*row as usize], *val))
    }
}
