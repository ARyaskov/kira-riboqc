use std::collections::BTreeMap;

// Constants for enrichment scaling.
pub const OFFSET: f64 = 0.5;
pub const SCALE: f64 = 1.5;

// Background selection size.
pub const BG_TOP_K: usize = 256;

// Translation load components.
// rRNA processing proxy genes.
pub const RRNA_PROXY: &[&str] = &["NOP56", "NOP58", "FBL", "DKC1"];
// Translation machinery genes.
pub const MACHINERY: &[&str] = &["EEF1A1", "EEF2", "EEF1B2", "EIF4A1", "EIF4G1"];

// Selective translation / stress response genes.
pub const STRESS: &[&str] = &["ATF4", "DDIT3", "XBP1", "HSPA5", "DNAJB9"];
pub const METAB: &[&str] = &["SLC2A1", "LDHA", "PFKP", "PDK1"];
pub const LUXURY: &[&str] = &["PCNA", "HMGB2", "TYMS"];

// RQC machinery genes.
pub const RQC_MACH: &[&str] = &["ZNF598", "ASCC3", "PELO", "HBS1L", "NEMF", "LTN1", "TCF25"];

// Proteostasis proxy genes for TPC.
pub const CHAPERONE: &[&str] = &["HSPA1A", "HSPA8", "HSP90AA1", "DNAJB1"];
pub const UBIQUITIN: &[&str] = &["UBB", "UBC", "SQSTM1"];

pub fn clamp01(x: f64) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    if x < 0.0 {
        0.0
    } else if x > 1.0 {
        1.0
    } else {
        x
    }
}

pub fn clamp_to_minus1_plus1(x: f64) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    if x < -1.0 {
        -1.0
    } else if x > 1.0 {
        1.0
    } else {
        x
    }
}

pub fn resolve_gene_set(map: &BTreeMap<String, u32>, symbols: &[&str]) -> Vec<u32> {
    let mut ids = Vec::new();
    for sym in symbols {
        if let Some(id) = map.get(&sym.to_string()) {
            ids.push(*id);
        }
    }
    ids
}

pub fn is_ribosomal(symbol: &str) -> bool {
    symbol.starts_with("RPL") || symbol.starts_with("RPS")
}

pub fn is_mito(symbol: &str) -> bool {
    symbol.starts_with("MT-") || symbol.starts_with("MT.")
}
