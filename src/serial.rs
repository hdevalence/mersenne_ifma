//! Serial arithmetic.
//!
//! This code implements prime-field and extension-field arithmetic
//! using `u128`s. Speed is not the highest priority, because the idea
//! is that the bulk of the work will be done using the vectorized
//! implementation.
//!

/// The Mersenne prime \\(2^{127} - 1\\).
const P: u128 = (1 << 127) - 1;

/// An element of the Mersenne field.
pub struct F127(u128);

impl From<u128> for F127 {
    // XXX check behaviour in extremal case 0b111...11 (not fully reduced?)
    #[inline]
    fn from(x: u128) -> F127 {
        let y = x.wrapping_sub(P);
        // sign_mask is 0b111...11 if x < P, 0b000...00 if x >= P
        let sign_mask = (0i128 - ((y >> 127) as i128)) as u128;

        F127(y.wrapping_add(sign_mask & P))
    }
}

impl Into<u128> for F127 {
    #[inline]
    fn into(self) -> u128 {
        self.0
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
}
