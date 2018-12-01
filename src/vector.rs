//! Vectorized arithmetic.
//!
//! # Choice of radix
//!
//! We want to use a binary radix \\(2^r\\) for some \\(r\\).  For
//! IFMA we need the limbs to fit within 52 bits, so we could choose:
//!
//! * \\(r = 52\\) (saturated)
//! * \\(r = 51\\) (unsaturated, one carry bit)
//! * a much smaller radix.
//!
//! We need at least 127 bits, so we need at least 3 limbs. Choosing
//! \\(r = 52\\) is not ideal because saturated arithmetic forces
//! sequential dependencies between instructions, making ILP more
//! difficult.  Choosing \\(r = 51\\) means that the limb boundaries
//! are not aligned with the bitsize of the prime, which means that
//! reducing high product terms crosses the boundaries of the lower
//! limbs.  For instance, using \\( r = 51 \\), we would compute a
//! 52-bit product term \\(z_3 2^{153} = z_3 2^{127} 2^{26} = z_3
//! 2^{26} \\). This term needs to be placed at position \\(2^{26}\\),
//! crossing a limb boundary.  The problem here is that weight of the
//! value is not distributed evenly across the limbs when using \\(r =
//! 51\\).
//!
//! The third option is to choose a much smaller radix, \\( r = 43
//! \\), so that \\( 3r = 129 \\).  This means the weight is spread
//! evenly across the limbs and the limb boundaries are more closely
//! aligned with the bitsize of the prime.

use core::ops::Mul;

use packed_simd::u64x4;

use crate::ifma::{madd52hi, madd52lo};

pub struct F127x4(u64x4, u64x4, u64x4);

// XXX should this use u128 or F127 ?
impl From<(u128, u128, u128, u128)> for F127x4 {
    fn from(x: (u128, u128, u128, u128)) -> F127x4 {
        let low_43_bits = (1 << 43) - 1;

        F127x4(
            u64x4::new(
                (x.0 & low_43_bits) as u64,
                (x.1 & low_43_bits) as u64,
                (x.2 & low_43_bits) as u64,
                (x.3 & low_43_bits) as u64,
            ),
            u64x4::new(
                ((x.0 >> 43) & low_43_bits) as u64,
                ((x.1 >> 43) & low_43_bits) as u64,
                ((x.2 >> 43) & low_43_bits) as u64,
                ((x.3 >> 43) & low_43_bits) as u64,
            ),
            u64x4::new(
                ((x.0 >> 86) & low_43_bits) as u64,
                ((x.1 >> 86) & low_43_bits) as u64,
                ((x.2 >> 86) & low_43_bits) as u64,
                ((x.3 >> 86) & low_43_bits) as u64,
            ),
        )
    }
}

impl Into<(u128, u128, u128, u128)> for F127x4 {
    fn into(self) -> (u128, u128, u128, u128) {
        (
            (self.0.extract(0) as u128)
                + ((self.1.extract(0) as u128) << 43)
                + ((self.2.extract(0) as u128) << 86),
            (self.0.extract(1) as u128)
                + ((self.1.extract(1) as u128) << 43)
                + ((self.2.extract(1) as u128) << 86),
            (self.0.extract(2) as u128)
                + ((self.1.extract(2) as u128) << 43)
                + ((self.2.extract(2) as u128) << 86),
            (self.0.extract(3) as u128)
                + ((self.1.extract(3) as u128) << 43)
                + ((self.2.extract(3) as u128) << 86),
        )
    }
}

impl Mul<F127x4> for F127x4 {
    type Output = F127x4;
    #[inline]
    fn mul(self, other: F127x4) -> F127x4 {
        use super::ifma::{madd52hi, madd52lo};

        let (x0, y0) = (self.0, other.0);
        let (x1, y1) = (self.1, other.1);
        let (x2, y2) = (self.2, other.2);

        let mut z0_a = u64x4::splat(0);
        let mut z0_b = u64x4::splat(0);
        let mut z0_c = u64x4::splat(0);
        let mut z1_a = u64x4::splat(0);
        let mut z1_b = u64x4::splat(0);
        let mut z1_c = u64x4::splat(0);
        let mut z2_a = u64x4::splat(0);
        let mut z2_b = u64x4::splat(0);
        let mut z2_c = u64x4::splat(0);

        z0_a = madd52hi(z0_a, x2, y0);
        z0_a = z0_a << 11;
        z0_a = madd52lo(z0_a, x0, y0);

        z0_b = madd52lo(z0_b, x2, y1);
        z0_c = madd52hi(z0_c, x1, y1);

        unimplemented!();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn from_into_F127_round_trips() {
        let xs = (
            101054725971136791246222244709531340474u128,
            38188712660835962328561942614081743514u128,
            43654918112560223727172090912658261884u128,
            61331686004747624160469066397670963925u128,
        );

        let x_vec: F127x4 = xs.into();

        assert_eq!(xs, x_vec.into());
    }
}
