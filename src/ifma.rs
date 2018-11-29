//! Definitions of (256-bit wide) IFMA intrinsics.
//!
//! Because Rust doesn't have support for AVX-512 yet, we manually
//! define our own intrinsics that lower to 256-bit IFMA instructions.
//! This doesn't really have a downside compared to using 512-bit
//! vectors, because on the only available Cannonlake processor, the
//! i3-8121U, executes IFMA at 256-bit wide anyways, and this nicely
//! matches the 4-way parallel Edwards formulas.

use packed_simd::u64x4;

// The `link_name`s below are pulled out of LLVM tablegen, have
// changed in the past, and might change again in the future.
#[allow(improper_ctypes)]
extern "C" {
    #[link_name = "llvm.x86.avx512.vpmadd52l.uq.256"]
    fn madd52lo_intrin(z: u64x4, x: u64x4, y: u64x4) -> u64x4;
    #[link_name = "llvm.x86.avx512.vpmadd52h.uq.256"]
    fn madd52hi_intrin(z: u64x4, x: u64x4, y: u64x4) -> u64x4;
}

/// A safe wrapper around `vpmadd52luq`.
///
/// The intrinsic itself is unsafe because it could generate SIGILL,
/// but this crate can't be compiled except for IFMA targets.
#[inline]
pub fn madd52lo(z: u64x4, x: u64x4, y: u64x4) -> u64x4 {
    unsafe { madd52lo_intrin(z, x, y) }
}

/// A safe wrapper around `vpmadd52huq`.
///
/// The intrinsic itself is unsafe because it could generate SIGILL,
/// but this crate can't be compiled except for IFMA targets.
#[inline]
pub fn madd52hi(z: u64x4, x: u64x4, y: u64x4) -> u64x4 {
    unsafe { madd52hi_intrin(z, x, y) }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_intrinsics() {
        let a = u64x4::new(1, 2, 3, 4);
        let b = u64x4::splat((1 << 52) + 3);
        let c = u64x4::new(5, 6, 7, 8);

        let x = madd52lo(a, b, c);
        let y = madd52hi(a, b, c);

        assert_eq!(x, u64x4::new(1 + 3 * 5, 2 + 3 * 6, 3 + 3 * 7, 4 + 3 * 8));
        assert_eq!(y, u64x4::new(1 + 1 * 5, 2 + 1 * 6, 3 + 1 * 7, 4 + 1 * 8));
    }
}
