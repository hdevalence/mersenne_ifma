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
//!
//!

pub struct F127x4(u64x4, u64x4, u64x4);

impl From<(u128, u128, u128, u128)> for F127x4 {
    fn from(x: (u128, u128, u128, u128)) -> F127x4 {
        let low_43_bits = (1 << 43) - 1;

        F127x4(
            u64x4::new(
                x.0 & low_43_bits,
                x.1 & low_43_bits,
                x.2 & low_43_bits,
                x.3 & low_43_bits,
            ),
            u64x4::new(
                (x.0 >> 43) & low_43_bits,
                (x.1 >> 43) & low_43_bits,
                (x.2 >> 43) & low_43_bits,
                (x.3 >> 43) & low_43_bits,
            ),
            u64x4::new(
                (x.0 >> 86) & low_43_bits,
                (x.1 >> 86) & low_43_bits,
                (x.2 >> 86) & low_43_bits,
                (x.3 >> 86) & low_43_bits,
            ),
        )
    }
}

impl Into<(u128, u128, u128, u128)> for F127x4 {
    fn into(x: F127x4) -> (u128, u128, u128, u128) {
        (
            x.0.extract(0) + (x.1.extract(0) << 43) + (x.2.extract(0) << 86),
            x.0.extract(1) + (x.1.extract(1) << 43) + (x.2.extract(1) << 86),
            x.0.extract(2) + (x.1.extract(2) << 43) + (x.2.extract(2) << 86),
            x.0.extract(3) + (x.1.extract(3) << 43) + (x.2.extract(3) << 86),
        )
    }
}

//impl Mul<F127x4> for F127x4(

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

        let x_vec: F127 = xs.into();

        assert_eq!(xs, x_vec.into());
    }
}
