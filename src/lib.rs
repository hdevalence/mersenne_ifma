//! An IFMA-based implementation of field operations over the Mersenne
//! prime `2**127 -1`.

//#![no_std]
#![feature(simd_ffi, link_llvm_intrinsics)]
#![deny(missing_docs)]

//#[cfg(not(target_feature = "avx512ifma"))]
//compile_error!("This crate requires AVX512-IFMA");

// The `packed_simd` crate contains what would have been the Rust SIMD
// code, except that it was decided to use untyped Intel __m256i
// bag-of-bits types instead of nice u64x4 types with arithmetic
// operations, so we use it instead of core::arch.
extern crate packed_simd;

#[cfg(target_feature = "avx512ifma")]
mod ifma;
pub mod serial;
#[cfg(target_feature = "avx512ifma")]
mod vector;
