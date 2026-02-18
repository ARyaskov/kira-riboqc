#[derive(Debug, Clone)]
pub struct MembershipVec {
    data: Vec<u8>,
}

impl MembershipVec {
    pub fn from_gene_ids(gene_ids: &[u32], n_genes: usize) -> Self {
        let mut data = vec![0u8; n_genes];
        for gid in gene_ids {
            if let Some(slot) = data.get_mut(*gid as usize) {
                *slot = 1;
            }
        }
        Self { data }
    }

    pub fn contains(&self, gene_id: u32) -> bool {
        self.data
            .get(gene_id as usize)
            .map(|v| *v == 1)
            .unwrap_or(false)
    }
}
