//! A serial implementation of prime- and extension- field arithmetic.
//!
//! Speed is not the highest priority, because the idea is that the
//! bulk of the work will be done using the vectorized implementation.

mod ext_field;
mod prime_field;

pub use self::ext_field::ExtF127;
pub use self::prime_field::F127;
