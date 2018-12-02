//! Serial arithmetic.
//!
//! This code implements prime-field and extension-field arithmetic
//! using `u128`s. Speed is not the highest priority, because the idea
//! is that the bulk of the work will be done using the vectorized
//! implementation.
//!

use core::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

/// The Mersenne prime \\(2^{127} - 1\\).
const P: u128 = (1 << 127) - 1;

/// An element of the Mersenne field.
///
/// # Invariant
///
/// The inner `u128` always lies in the range \\([0, 2^{127} - 1]\\).
#[derive(Copy, Clone)]
pub struct F127(u128);

impl From<u128> for F127 {
    #[inline]
    fn from(x: u128) -> F127 {
        // Assume x in range [0, 2*P).
        debug_assert!(x < 2 * P);

        // Since x < 2*P, we have x - P < 2P - P = P < 2^127,
        // so the high bit of y is 1 if and only if x < P.
        let y = x.wrapping_sub(P);
        // Therefore sign_mask is 0b111...11 if x < P,
        //                        0b000...00 if x >= P
        let sign_mask = (0i128 - ((y >> 127) as i128)) as u128;

        // Conditionally add P back to y, to obtain
        // x_red = x      if 0 <= x < P
        // x_red = x - P  if P <= x < 2*P
        let x_red = y.wrapping_add(sign_mask & P);

        // In both cases x_red < P.
        F127(x_red)
    }
}

impl Into<u128> for F127 {
    #[inline]
    fn into(self) -> u128 {
        self.0
    }
}

impl Add<F127> for F127 {
    type Output = F127;
    #[inline]
    fn add(self, other: F127) -> F127 {
        let (x, y) = (self.0, other.0);

        // This cannot overflow:
        // x, y <= P, so z = x + y <= 2*P = 2^128 - 2 < 2^128.
        let z = x + y;

        // Write z = z0 + z1*2^127 with z0 < 2^127, z1 = 0,1
        let z0 = z & P;
        let z1 = z >> 127;

        // Since 2^127 = 1 mod P, z = z0 + z1 mod P.
        //
        // If z1 = 0 then z0 + z1 = z0 <= P.
        //
        // If z1 = 1 then z0 + z1 = (z - 2^127) + 1
        //                        < (2*(2^127 - 1) - 2^127) + 1
        //                        = (2^127 - 2) + 1 = P.
        //
        // In either case, z0 + z1 <= P.
        F127(z0 + z1)
    }
}

impl Sub<F127> for F127 {
    type Output = F127;
    #[inline]
    fn sub(self, other: F127) -> F127 {
        // XXX document
        let (x, y) = (self.0, other.0);
        let z = x.wrapping_sub(y);
        let z0 = z & P;
        let z1 = z >> 127;
        F127(z0 - z1)
    }
}

impl Neg for F127 {
    type Output = F127;
    #[inline]
    fn neg(self) -> F127 {
        // This does not underflow, since self.0 <= P
        F127(P - self.0)
    }
}

impl Mul<F127> for F127 {
    type Output = F127;
    #[inline]
    fn mul(self, other: F127) -> F127 {
        let (x, y) = (self.0, other.0);
        // x0 < 2^64, x1 < 2^63
        let (x0, x1) = (x as u64, (x >> 64) as u64);
        // y0 < 2^64, y1 < 2^63
        let (y0, y1) = (y as u64, (y >> 64) as u64);

        let m = |x: u64, y: u64| (x as u128) * (y as u128);

        // Write the product in mixed-radix
        //
        // z = z0 + z1*2^64 + z2*2^127
        //
        // with
        //
        // z0 = x0*y0
        // z1 = x1*y0 + x0*y1
        // z2 =         x1*y1*2

        let z0 = m(x0, y0);
        // The high 64 bits of z0 are accounted for in z1
        let z1 = m(x1, y0) + m(x0, y1) + (z0 >> 64);
        // The high 65 bits of z1 are accounted for in z2
        let z2 = m(x1, 2 * y1) + (z1 >> 63);

        // Now write z0, z1, z2 in radix 2^127 as w0 + w1*2^127:
        const MASK63: u64 = (1u64 << 63) - 1;
        // w0 is composed of the low 64 bits of z0 and the low 63 bits of z1
        let w0 = ((z0 as u64) as u128) | ((((z1 as u64) & MASK63) as u128) << 64);
        // w1 is just z2
        let w1 = z2;

        // Combine high and low halves, then reduce a carry bit
        let w = w0 + w1;
        F127((w & P) + (w >> 127))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn f127_from_u128_roundtrips_on_reduced_values() {
        let xs = [
            101054725971136791246222244709531340474u128,
            38188712660835962328561942614081743514u128,
            43654918112560223727172090912658261884u128,
            61331686004747624160469066397670963925u128,
        ];

        for x in &xs {
            let expected_x: u128 = F127::from(*x).into();
            assert_eq!(*x, expected_x);
        }
    }

    #[test]
    fn f127_from_u128_performs_reduction() {
        let unred_red_pairs = [
            (
                316359973995368844939217233962370990276u128,
                146218790534899613207529930246486884549u128,
            ),
            (
                323686882786572482033927420505583807954u128,
                153545699326103250302240116789699702227u128,
            ),
            (
                339808425234034235397003825875747537381u128,
                169667241773565003665316522159863431654u128,
            ),
            (
                78832676377202070809965755704023168971u128,
                78832676377202070809965755704023168971u128,
            ),
        ];

        for (x_unred, x_red) in &unred_red_pairs {
            let expected_x_red: u128 = F127::from(*x_unred).into();
            assert_eq!(*x_red, expected_x_red);
        }
    }

    #[test]
    fn iterated_add() {
        let x = F127::from(38188712660835962328561942614081743514u128);
        let mut z = F127::from(0);

        for _i in 0..1024 {
            z = z + x;
        }

        // XXX consider deriving Eq on F127
        let z_repr: u128 = z.into();
        assert_eq!(z_repr, 142910752248571357891036685882245146853u128);
    }

    #[test]
    fn iterated_sub() {
        let x = F127::from(38188712660835962328561942614081743514u128);
        let mut z = F127::from(0);

        for _i in 0..1024 {
            z = z - x;
        }

        // XXX consider deriving Eq on F127
        let z_repr: u128 = z.into();
        assert_eq!(z_repr, 27230431211897873840650617833638958874u128);
    }

    #[test]
    fn iterated_mul() {
        let x = F127::from(38188712660835962328561942614081743514u128);
        let mut z = F127::from(1);

        for _i in 0..1024 {
            z = z * x;
        }

        // XXX consider deriving Eq on F127
        let z_repr: u128 = z.into();
        assert_eq!(z_repr, 63115059284280959221284862234304285851u128);
    }
}
