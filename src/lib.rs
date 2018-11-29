//! An IFMA-based implementation of field operations over the Mersenne
//! prime `2**127 -1`.

#![feature(simd_ffi, link_llvm_intrinsics)]

extern crate packed_simd;

mod ifma;
mod serial;
mod vector;
