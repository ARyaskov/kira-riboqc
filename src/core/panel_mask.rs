pub struct PanelMaskBuilder {
    mask: Vec<u32>,
    n_panels: u32,
}

impl PanelMaskBuilder {
    pub fn new(n_genes: usize) -> Self {
        Self {
            mask: vec![0u32; n_genes],
            n_panels: 0,
        }
    }

    pub fn add_panel(&mut self, gene_ids: &[u32]) -> PanelBit {
        assert!(self.n_panels < 32, "PanelMask supports up to 32 panels");
        let bit = 1u32 << self.n_panels;
        for &gid in gene_ids {
            let idx = gid as usize;
            if idx < self.mask.len() {
                self.mask[idx] |= bit;
            }
        }
        let panel = PanelBit {
            bit,
            index: self.n_panels as usize,
        };
        self.n_panels += 1;
        panel
    }

    pub fn build(self) -> PanelMask {
        PanelMask {
            mask: self.mask,
            n_panels: self.n_panels as usize,
        }
    }
}

#[derive(Clone, Copy)]
pub struct PanelBit {
    pub bit: u32,
    pub index: usize,
}

pub struct PanelMask {
    mask: Vec<u32>,
    n_panels: usize,
}

impl PanelMask {
    #[inline]
    pub fn n_panels(&self) -> usize {
        self.n_panels
    }

    #[inline]
    pub fn get(&self, gene_id: u32) -> u32 {
        *self.mask.get(gene_id as usize).unwrap_or(&0)
    }
}

#[derive(Clone, Copy, Default)]
pub struct Acc {
    pub sum: f64,
    pub count: u32,
}

impl Acc {
    #[inline]
    pub fn add(&mut self, x: f64) {
        self.sum += x;
        self.count += 1;
    }

    #[inline]
    pub fn mean(&self) -> f64 {
        if self.count == 0 {
            f64::NAN
        } else {
            self.sum / (self.count as f64)
        }
    }
}
