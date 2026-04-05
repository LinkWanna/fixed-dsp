use crate::common::error::{Error, Result};
use crate::common::tables::{SQRT_INITIAL_LUT_U16, SQRT_INITIAL_LUT_U32};

#[inline]
fn sqrt_newton_q15(number: i32, mut var1: i32) -> i32 {
    for _ in 0..3 {
        let mut temp = (var1 * var1) >> 12;
        temp = (number * temp) >> 15;
        temp = 0x3000 - temp;
        var1 = (var1 * temp) >> 13;
    }

    (number * var1) >> 12
}

#[inline]
fn sqrt_newton_q31(number: i64, mut var1: i64) -> i64 {
    for _ in 0..3 {
        let mut temp = (var1 * var1) >> 28;
        temp = (number * temp) >> 31;
        temp = 0x30000000 - temp;
        var1 = (var1 * temp) >> 29;
    }

    (number * var1) >> 28
}

/// Q15 square root using the CMSIS lookup-table + Newton iteration.
pub fn sqrt_i16(input: i16) -> Result<i16> {
    if input < 0 {
        return Err(Error::NanInf);
    }

    if input == 0 {
        return Ok(0);
    }

    let input_u = input as u16 as u32;
    let sign_bits = input_u.leading_zeros() as i32 - 17;
    let shift = if sign_bits & 1 == 0 {
        sign_bits
    } else {
        sign_bits - 1
    } as u32;

    let number = (input as i32).wrapping_shl(shift);
    let index = (((number >> 11) as usize) - (0x2000usize >> 11)) as usize;
    let var1 = SQRT_INITIAL_LUT_U16[index] as i32;
    let mut result = sqrt_newton_q15(number, var1);

    let final_shift = if sign_bits & 1 == 0 {
        (sign_bits / 2) as u32
    } else {
        ((sign_bits - 1) / 2) as u32
    };
    result >>= final_shift;

    Ok(result as i16)
}

/// Q31 square root using the CMSIS lookup-table + Newton iteration.
pub fn sqrt_i32(input: i32) -> Result<i32> {
    if input < 0 {
        return Err(Error::NanInf);
    }

    if input == 0 {
        return Ok(0);
    }

    let input_u = input as u32;
    let sign_bits = input_u.leading_zeros() as i32 - 1;
    let shift = if sign_bits & 1 == 0 {
        sign_bits
    } else {
        sign_bits - 1
    } as u32;

    let number = (input as i64).wrapping_shl(shift);
    let index = (((number >> 26) as usize) - (0x20000000usize >> 26)) as usize;
    let var1 = SQRT_INITIAL_LUT_U32[index] as i64;
    let mut result = sqrt_newton_q31(number, var1);

    let final_shift = if sign_bits & 1 == 0 {
        (sign_bits / 2) as u32
    } else {
        ((sign_bits - 1) / 2) as u32
    };
    result >>= final_shift;

    Ok(result as i32)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sqrt_q15_handles_sign_and_zero() {
        assert_eq!(sqrt_i16(-1), Err(Error::NanInf));
        assert_eq!(sqrt_i16(0), Ok(0));
    }

    #[test]
    fn sqrt_q31_handles_sign_and_zero() {
        assert_eq!(sqrt_i32(-1), Err(Error::NanInf));
        assert_eq!(sqrt_i32(0), Ok(0));
    }
}
