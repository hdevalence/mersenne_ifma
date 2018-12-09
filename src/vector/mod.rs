//! Vectorized prime- and extension- field arithmetic.

mod ext_field;
mod prime_field;

pub use self::ext_field::ExtF127x4;
pub use self::prime_field::F127x4;
