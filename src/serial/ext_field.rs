use super::F127;

/// An element of the quadratic extension field F127\[i\]
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub struct F127Ext(F127, F127);

use core::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

impl From<(u128, u128)> for F127Ext {
    #[inline]
    fn from(x: (u128, u128)) -> F127Ext {
        F127Ext(x.0.into(), x.1.into())
    }
}

impl Into<(u128, u128)> for F127Ext {
    #[inline]
    fn into(self) -> (u128, u128) {
        (self.0.into(), self.1.into())
    }
}

impl Add<F127Ext> for F127Ext {
    type Output = F127Ext;
    #[inline]
    fn add(self, other: F127Ext) -> F127Ext {
        F127Ext(self.0 + other.0, self.1 + other.1)
    }
}

impl Sub<F127Ext> for F127Ext {
    type Output = F127Ext;
    #[inline]
    fn sub(self, other: F127Ext) -> F127Ext {
        F127Ext(self.0 - other.0, self.1 - other.1)
    }
}

impl Mul<F127Ext> for F127Ext {
    type Output = F127Ext;
    #[inline]
    fn mul(self, other: F127Ext) -> F127Ext {
        let (a, b) = (self.0, self.1);
        let (c, d) = (other.0, other.1);

        let ac = a * c;
        let bd = b * d;

        F127Ext(ac - bd, (b-a)*(c-d) + ac + bd)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mul_vs_sage() {
        let x = F127Ext::from((
            64602349736890547230188097686032968383u128,
            58401672467634577377614110902426170573u128,
        ));
        let y = F127Ext::from((
            36178516401130528447705023720593931265u128,
            57463319253223551344966612196770510351u128,
        ));
        let z = F127Ext::from((
            167087788139004297409615161698155907378u128,
            77896319433764489876703096387833153505u128,
        ));

        assert_eq!(x * y, z);
    }
}
