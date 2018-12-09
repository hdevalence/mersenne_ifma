use super::F127;

/// An element of the quadratic extension field F127\[i\]
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub struct ExtF127(pub(crate) F127, pub(crate) F127);

use core::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

impl From<(u128, u128)> for ExtF127 {
    #[inline]
    fn from(x: (u128, u128)) -> ExtF127 {
        ExtF127(x.0.into(), x.1.into())
    }
}

impl Into<(u128, u128)> for ExtF127 {
    #[inline]
    fn into(self) -> (u128, u128) {
        (self.0.into(), self.1.into())
    }
}

impl Add<ExtF127> for ExtF127 {
    type Output = ExtF127;
    #[inline]
    fn add(self, other: ExtF127) -> ExtF127 {
        ExtF127(self.0 + other.0, self.1 + other.1)
    }
}

impl Sub<ExtF127> for ExtF127 {
    type Output = ExtF127;
    #[inline]
    fn sub(self, other: ExtF127) -> ExtF127 {
        ExtF127(self.0 - other.0, self.1 - other.1)
    }
}

impl Mul<ExtF127> for ExtF127 {
    type Output = ExtF127;
    #[inline]
    fn mul(self, other: ExtF127) -> ExtF127 {
        let (a, b) = (self.0, self.1);
        let (c, d) = (other.0, other.1);

        let ac = a * c;
        let bd = b * d;

        ExtF127(ac - bd, (b - a) * (c - d) + ac + bd)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mul_vs_sage() {
        let x = ExtF127::from((
            64602349736890547230188097686032968383u128,
            58401672467634577377614110902426170573u128,
        ));
        let y = ExtF127::from((
            36178516401130528447705023720593931265u128,
            57463319253223551344966612196770510351u128,
        ));
        let z = ExtF127::from((
            167087788139004297409615161698155907378u128,
            77896319433764489876703096387833153505u128,
        ));

        assert_eq!(x * y, z);
    }
}
