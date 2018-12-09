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

use core::ops::{Add, Mul, Neg};

use packed_simd::u64x4;

use crate::ifma::{madd52hi, madd52lo};

use crate::serial::F127;

#[derive(Copy, Clone)]
pub struct F127x4(u64x4, u64x4, u64x4);

impl From<(F127, F127, F127, F127)> for F127x4 {
    fn from(x: (F127, F127, F127, F127)) -> F127x4 {
        let low_43_bits = (1 << 43) - 1;

        F127x4(
            u64x4::new(
                ((x.0).0 & low_43_bits) as u64,
                ((x.1).0 & low_43_bits) as u64,
                ((x.2).0 & low_43_bits) as u64,
                ((x.3).0 & low_43_bits) as u64,
            ),
            u64x4::new(
                (((x.0).0 >> 43) & low_43_bits) as u64,
                (((x.1).0 >> 43) & low_43_bits) as u64,
                (((x.2).0 >> 43) & low_43_bits) as u64,
                (((x.3).0 >> 43) & low_43_bits) as u64,
            ),
            u64x4::new(
                (((x.0).0 >> 86) & low_43_bits) as u64,
                (((x.1).0 >> 86) & low_43_bits) as u64,
                (((x.2).0 >> 86) & low_43_bits) as u64,
                (((x.3).0 >> 86) & low_43_bits) as u64,
            ),
        )
    }
}

impl Into<(F127, F127, F127, F127)> for F127x4 {
    fn into(mut self) -> (F127, F127, F127, F127) {
        // We store self as a 129-bit (or more) integer,
        // so reduce it before conversion
        // XXX check bounds
        let mask = u64x4::splat((1 << 41) - 1);
        let c = self.2 >> 41;
        self.0 += c;
        self.2 &= mask;

        (
            F127::from(
                (self.0.extract(0) as u128)
                    + ((self.1.extract(0) as u128) << 43)
                    + ((self.2.extract(0) as u128) << 86),
            ),
            F127::from(
                (self.0.extract(1) as u128)
                    + ((self.1.extract(1) as u128) << 43)
                    + ((self.2.extract(1) as u128) << 86),
            ),
            F127::from(
                (self.0.extract(2) as u128)
                    + ((self.1.extract(2) as u128) << 43)
                    + ((self.2.extract(2) as u128) << 86),
            ),
            F127::from(
                (self.0.extract(3) as u128)
                    + ((self.1.extract(3) as u128) << 43)
                    + ((self.2.extract(3) as u128) << 86),
            ),
        )
    }
}

impl Neg for F127x4 {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self {
        F127x4(
            u64x4::splat(((1 << 43) - 1) << 3) - self.0,
            u64x4::splat(((1 << 43) - 1) << 3) - self.1,
            u64x4::splat(((1 << 41) - 1) << 3) - self.2,
        )
    }
}

impl Add<F127x4> for F127x4 {
    type Output = F127x4;
    #[inline]
    fn add(self, other: F127x4) -> F127x4 {
        F127x4(self.0 + other.0, self.1 + other.1, self.2 + other.2)
    }
}

impl Mul<F127x4> for F127x4 {
    type Output = F127x4;
    #[inline]
    fn mul(self, other: F127x4) -> F127x4 {
        let (x0, y0) = (self.0, other.0);
        let (x1, y1) = (self.1, other.1);
        let (x2, y2) = (self.2, other.2);

        // We have 18 multiplications, want 8 independent chains to
        // saturate the EUs, so split into 9 chains of length 2.

        let mut z0_a = u64x4::splat(0);
        let mut z0_b = u64x4::splat(0);
        let mut z0_c = u64x4::splat(0);
        let mut z1_a = u64x4::splat(0);
        let mut z1_b = u64x4::splat(0);
        let mut z1_c = u64x4::splat(0);
        let mut z2_a = u64x4::splat(0);
        let mut z2_b = u64x4::splat(0);
        let mut z2_c = u64x4::splat(0);

        z0_a = madd52hi(z0_a, x2, y0); // 2^11
        z0_b = madd52lo(z0_b, x2, y1); // 2^2
        z0_c = madd52hi(z0_c, x1, y1); // 2^11

        z1_a = madd52hi(z1_a, x0, y0); // 2^9
        z1_b = madd52hi(z1_b, x2, y1); // 2^11
        z1_c = madd52lo(z1_c, x1, y0); // 2^0

        z2_a = madd52hi(z2_a, x2, y2); // 2^11
        z2_b = madd52hi(z2_b, x0, y1); // 2^9
        z2_c = madd52lo(z2_c, x2, y0); // 2^0

        z0_a = z0_a << 11; // 2^11 -> 2^0
        z1_a = z1_a << 07; // 2^9  -> 2^2
        z2_a = z2_a << 11; // 2^11 -> 2^0

        z0_a = madd52lo(z0_a, x0, y0); // 2^0
        z0_b = madd52lo(z0_b, x1, y2); // 2^2
        z0_c = madd52hi(z0_c, x0, y2); // 2^11

        z1_a = madd52lo(z1_a, x2, y2); // 2^2
        z1_b = madd52hi(z1_b, x1, y2); // 2^11
        z1_c = madd52lo(z1_c, x0, y1); // 2^0

        z2_a = madd52lo(z2_a, x0, y2); // 2^0
        z2_b = madd52hi(z2_b, x1, y0); // 2^9
        z2_c = madd52lo(z2_c, x1, y1); // 2^0

        let z0 = z0_a + (z0_b << 2) + (z0_c << 11);
        let z1 = (z1_a << 2) + (z1_b << 11) + z1_c;
        let z2 = z2_a + (z2_b << 9) + z2_c;

        let c0 = z0 >> 43;
        let c1 = z1 >> 43;
        let c2 = z2 >> 43;

        let mask = u64x4::splat((1 << 43) - 1);

        F127x4((z0 & mask) + (c2 << 2), (z1 & mask) + c0, (z2 & mask) + c1)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[allow(non_snake_case)]
    #[test]
    fn from_into_F127_round_trips() {
        let xs: (F127, F127, F127, F127) = (
            101054725971136791246222244709531340474u128.into(),
            38188712660835962328561942614081743514u128.into(),
            43654918112560223727172090912658261884u128.into(),
            61331686004747624160469066397670963925u128.into(),
        );

        let x_vec: F127x4 = xs.into();

        assert_eq!(xs, x_vec.into());
    }

    #[test]
    fn mul_matches_serial() {
        let xs: (F127, F127, F127, F127) = (
            101054725971136791246222244709531340474u128.into(),
            38188712660835962328561942614081743514u128.into(),
            43654918112560223727172090912658261884u128.into(),
            61331686004747624160469066397670963925u128.into(),
        );

        let x_vec: F127x4 = xs.into();

        let z_vec = x_vec * x_vec;

        let zs: (F127, F127, F127, F127) = z_vec.into();

        assert_eq!(zs.0, xs.0 * xs.0);
        assert_eq!(zs.1, xs.1 * xs.1);
        assert_eq!(zs.2, xs.2 * xs.2);
        assert_eq!(zs.3, xs.3 * xs.3);
    }

    #[test]
    fn add_negation_is_zero() {
        let xs: (F127, F127, F127, F127) = (
            101054725971136791246222244709531340474u128.into(),
            38188712660835962328561942614081743514u128.into(),
            43654918112560223727172090912658261884u128.into(),
            61331686004747624160469066397670963925u128.into(),
        );

        let x_vec: F127x4 = xs.into();

        let zs: (F127, F127, F127, F127) = (x_vec + (-x_vec)).into();

        assert_eq!(zs.0, F127::zero());
        assert_eq!(zs.1, F127::zero());
        assert_eq!(zs.2, F127::zero());
        assert_eq!(zs.3, F127::zero());
    }

}
