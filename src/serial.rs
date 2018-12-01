//! Serial arithmetic.
//!
//! This code implements prime-field and extension-field arithmetic
//! using `u128`s. Speed is not the highest priority, because the idea
//! is that the bulk of the work will be done using the vectorized
//! implementation.
//!
//!


pub struct F127(u128);


