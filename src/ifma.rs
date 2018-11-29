use packed_simd::{u64x4, IntoBits};

#[allow(improper_ctypes)]
extern "C" {
    #[link_name = "llvm.x86.avx512.vpmadd52l.uq.256"]
    fn madd52lo(z: u64x4, x: u64x4, y: u64x4) -> u64x4;
    #[link_name = "llvm.x86.avx512.vpmadd52h.uq.256"]
    fn madd52hi(z: u64x4, x: u64x4, y: u64x4) -> u64x4;
}

#[cfg(test)]
mod tests {
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
