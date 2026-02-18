use kira_riboqc::simd::SIMD_KIND;

#[test]
fn simd_kind_is_valid() {
    assert!(matches!(SIMD_KIND, "avx2" | "neon" | "scalar"));
}
