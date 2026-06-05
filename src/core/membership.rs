#[derive(Debug, Clone)]
pub struct MembershipVec {
    bits: Vec<u64>,
    len: usize,
}

impl MembershipVec {
    pub fn from_gene_ids(gene_ids: &[u32], n_genes: usize) -> Self {
        let words = n_genes.div_ceil(64);
        let mut bits = vec![0u64; words];
        for gid in gene_ids {
            let idx = *gid as usize;
            if idx < n_genes {
                bits[idx >> 6] |= 1u64 << (idx & 63);
            }
        }
        Self { bits, len: n_genes }
    }

    #[inline]
    pub fn contains(&self, gene_id: u32) -> bool {
        let idx = gene_id as usize;
        if idx >= self.len {
            return false;
        }
        (self.bits[idx >> 6] >> (idx & 63)) & 1 == 1
    }
}
