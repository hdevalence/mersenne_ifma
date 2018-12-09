//! Vectorized arithmetic for the extension field

use super::F127x4;
use crate::serial::{F127, F127Ext};

/// A vector of four elements of the extension field.
pub struct ExtF127x4(F127x4, F127x4);

impl From<(F127Ext, F127Ext, F127Ext, F127Ext)> for ExtF127x4 {
    fn from(x: (F127Ext, F127Ext, F127Ext, F127Ext)) -> ExtF127x4 {
        ExtF127x4(
            ((x.0).0, (x.1).0, (x.2).0, (x.3).0).into(),
            ((x.0).1, (x.1).1, (x.2).1, (x.3).1).into(),
        )
    }
}

impl Into<(F127Ext, F127Ext, F127Ext, F127Ext)> for ExtF127x4 {
    fn into(self) -> (F127Ext, F127Ext, F127Ext, F127Ext) {
        let xs: (F127, F127, F127, F127) = self.0.into();
        let ys: (F127, F127, F127, F127) = self.1.into();

        (
            F127Ext(xs.0, ys.0),
            F127Ext(xs.1, ys.1),
            F127Ext(xs.2, ys.2),
            F127Ext(xs.3, ys.3),
        )
    }
}
