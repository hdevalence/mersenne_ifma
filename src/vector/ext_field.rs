//! Vectorized arithmetic for the extension field

use super::F127x4;
use crate::serial::{ExtF127, F127};

/// A vector of four elements of the extension field.
pub struct ExtF127x4(F127x4, F127x4);

impl From<(ExtF127, ExtF127, ExtF127, ExtF127)> for ExtF127x4 {
    fn from(x: (ExtF127, ExtF127, ExtF127, ExtF127)) -> ExtF127x4 {
        ExtF127x4(
            ((x.0).0, (x.1).0, (x.2).0, (x.3).0).into(),
            ((x.0).1, (x.1).1, (x.2).1, (x.3).1).into(),
        )
    }
}

impl Into<(ExtF127, ExtF127, ExtF127, ExtF127)> for ExtF127x4 {
    fn into(self) -> (ExtF127, ExtF127, ExtF127, ExtF127) {
        let xs: (F127, F127, F127, F127) = self.0.into();
        let ys: (F127, F127, F127, F127) = self.1.into();

        (
            ExtF127(xs.0, ys.0),
            ExtF127(xs.1, ys.1),
            ExtF127(xs.2, ys.2),
            ExtF127(xs.3, ys.3),
        )
    }
}

use core::ops::{Add, Mul, Neg};

impl Add<ExtF127x4> for ExtF127x4 {
    type Output = ExtF127x4;
    #[inline]
    fn add(self, other: ExtF127x4) -> ExtF127x4 {
        ExtF127x4(self.0 + other.0, self.1 + other.1)
    }
}

impl Mul<ExtF127x4> for ExtF127x4 {
    type Output = ExtF127x4;
    #[inline]
    fn mul(self, other: ExtF127x4) -> ExtF127x4 {
        let (a, b) = (self.0, self.1);
        let (c, d) = (other.0, other.1);

        let ac = a * c;
        let bd = b * d;

        ExtF127x4(ac + (-bd), (b + (-a)) * (c + (-d)) + ac + bd)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mul_matches_serial() {
        let xs = (
            ExtF127::from((
                43654918112560223727172090912658261884u128,
                38188712660835962328561942614081743514u128,
            )),
            ExtF127::from((
                38188712660835962328561942614081743514u128,
                61331686004747624160469066397670963925u128,
            )),
            ExtF127::from((
                38188712660835962328561942614081743514u128,
                43654918112560223727172090912658261884u128,
            )),
            ExtF127::from((
                43654918112560223727172090912658261884u128,
                61331686004747624160469066397670963925u128,
            )),
        );

        let x_vec: ExtF127x4 = xs.into();

        let z_vec = x_vec * x_vec;

        let zs: (ExtF127, ExtF127, ExtF127, ExtF127) = z_vec.into();

        assert_eq!(zs.0, xs.0 * xs.0);
        assert_eq!(zs.1, xs.1 * xs.1);
        assert_eq!(zs.2, xs.2 * xs.2);
        assert_eq!(zs.3, xs.3 * xs.3);
    }
}
