pub const RIBO_EXTENSION_PANEL_V1: &str = "RIBO_EXTENSION_PANEL_V1";

pub const PANEL_RIBOSOME_CORE: &[&str] = &[
    "RPL3", "RPL5", "RPL7", "RPL10", "RPL11", "RPL23", "RPS3", "RPS6", "RPS8", "RPS14", "RPS18",
];

pub const PANEL_BIOGENESIS: &[&str] = &[
    "NCL", "NPM1", "FBL", "DKC1", "UBTF", "POLR1A", "POLR1B", "WDR12", "PES1", "BOP1",
];

pub const PANEL_MTOR: &[&str] = &["MTOR", "RPTOR", "MLST8", "RPS6KB1", "EIF4EBP1", "AKT1"];

pub const PANEL_ISR: &[&str] = &["EIF2AK3", "EIF2AK4", "ATF4", "DDIT3", "PPP1R15A"];

pub const PANEL_INITIATION: &[&str] = &["EIF4E", "EIF4G1", "EIF3A", "EIF3B", "EEF1A1", "EEF2"];

pub const PANEL_PROTEOSTASIS: &[&str] = &["HSPA5", "HSP90AA1", "HSPB1", "PSMA1", "PSMB5"];

pub const MIN_GENES_PER_PANEL: usize = 3;
