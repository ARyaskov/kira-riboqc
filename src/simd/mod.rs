#[cfg(all(feature = "avx2", any(target_arch = "x86_64", target_arch = "x86")))]
pub use crate::simd::avx2::ln1p_f64;

#[cfg(all(
    not(all(feature = "avx2", any(target_arch = "x86_64", target_arch = "x86"))),
    feature = "neon",
    target_arch = "aarch64"
))]
pub use crate::simd::neon::ln1p_f64;

#[cfg(all(
    not(all(feature = "avx2", any(target_arch = "x86_64", target_arch = "x86"))),
    not(all(feature = "neon", target_arch = "aarch64"))
))]
pub use crate::simd::scalar::ln1p_f64;

#[cfg(all(feature = "avx2", any(target_arch = "x86_64", target_arch = "x86")))]
pub const SIMD_MODE: &str = "avx2";

#[cfg(all(
    not(all(feature = "avx2", any(target_arch = "x86_64", target_arch = "x86"))),
    feature = "neon",
    target_arch = "aarch64"
))]
pub const SIMD_MODE: &str = "neon";

#[cfg(all(
    not(all(feature = "avx2", any(target_arch = "x86_64", target_arch = "x86"))),
    not(all(feature = "neon", target_arch = "aarch64"))
))]
pub const SIMD_MODE: &str = "scalar";

#[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
pub const SIMD_KIND: &str = "avx2";

#[cfg(all(
    target_arch = "aarch64",
    not(all(target_arch = "x86_64", target_feature = "avx2"))
))]
pub const SIMD_KIND: &str = "neon";

#[cfg(not(any(
    all(target_arch = "x86_64", target_feature = "avx2"),
    target_arch = "aarch64"
)))]
pub const SIMD_KIND: &str = "scalar";

#[cfg(all(feature = "avx2", any(target_arch = "x86_64", target_arch = "x86")))]
pub mod avx2;
#[cfg(all(feature = "neon", target_arch = "aarch64"))]
pub mod neon;
pub mod scalar;
