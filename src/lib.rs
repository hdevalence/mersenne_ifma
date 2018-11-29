//! An IFMA-based implementation of field operations over the Mersenne
//! prime `2**127 -1`.

extern crate packed_simd;

mod ifma;
mod serial;
mod vector;
