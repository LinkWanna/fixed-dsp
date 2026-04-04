use crate::common::error::{Error, Result};

/// Fixed-point division for Q15 values.
///
/// The returned pair `(quotient, shift)` represents:
/// `numerator / denominator = quotient * 2^shift`, where `quotient` is normalized.
///
/// On divide-by-zero, [`Error::NanInf`] is returned.
pub fn div_i16(numerator: i16, denominator: i16) -> Result<(i16, i16)> {
    let mut shift = 0_i16;

    let sign = (numerator < 0) ^ (denominator < 0);

    if denominator == 0 {
        return Err(Error::NanInf);
    }

    let numerator_abs = numerator.saturating_abs();
    let denominator_abs = denominator.saturating_abs();

    let mut temp = ((numerator_abs as i32) << 15) / (denominator_abs as i32);

    let shift_for_normalizing = 17 - (temp as u32).leading_zeros() as i16;
    if shift_for_normalizing > 0 {
        shift = shift_for_normalizing;
        temp >>= shift_for_normalizing;
    }

    if sign {
        temp = -temp;
    }

    let quotient = temp as i16;
    Ok((quotient, shift))
}

/// Fixed-point division for Q31 values.
///
/// The returned pair `(quotient, shift)` represents:
/// `numerator / denominator = quotient * 2^shift`, where `quotient` is normalized.
///
/// On divide-by-zero, [`Error::NanInf`] is returned.
pub fn div_i32(numerator: i32, denominator: i32) -> Result<(i32, i16)> {
    let mut shift = 0_i16;

    let sign = (numerator < 0) ^ (denominator < 0);

    if denominator == 0 {
        return Err(Error::NanInf);
    }

    let numerator_abs = numerator.saturating_abs();
    let denominator_abs = denominator.saturating_abs();

    let mut temp = ((numerator_abs as i64) << 31) / (denominator_abs as i64);

    let shift_for_normalizing = 32 - ((temp >> 31) as u32).leading_zeros() as i16;
    if shift_for_normalizing > 0 {
        shift = shift_for_normalizing;
        temp >>= shift_for_normalizing;
    }

    if sign {
        temp = -temp;
    }

    let quotient = temp as i32;
    Ok((quotient, shift))
}
