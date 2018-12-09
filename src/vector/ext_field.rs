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
