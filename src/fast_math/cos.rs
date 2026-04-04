use crate::common::tables::{SIN_TABLE_U16, SIN_TABLE_U32};

const FAST_MATH_Q15_SHIFT: i32 = 16 - 10;
const FAST_MATH_Q31_SHIFT: i32 = 32 - 10;

/// x is in Q15 format (1.15 fixed-point), representing angles in the range [0, 2pi).
pub fn cos_i16(x: i16) -> i16 {
    // Add 0.25 turn (pi/2) to map cosine onto sine lookup table.
    let mut x_bits = (x as u16).wrapping_add(0x2000);
    if (x_bits as i16) < 0 {
        x_bits = x_bits.wrapping_add(0x8000);
    }

    let x_mapped = x_bits as i32;
    let index = ((x_bits as u32) >> FAST_MATH_Q15_SHIFT) as usize;
    let fract = (x_mapped - ((index as i32) << FAST_MATH_Q15_SHIFT)) << 9;

    let a = SIN_TABLE_U16[index] as i16 as i32;
    let b = SIN_TABLE_U16[index + 1] as i16 as i32;

    let mut cos_val = ((0x8000_i32 - fract) * a) >> 16;
    cos_val = ((cos_val << 16) + (fract * b)) >> 16;

    (cos_val << 1) as i16
}

/// x is in Q31 format (1.31 fixed-point), representing angles in the range [0, 2pi).
pub fn cos_i32(x: i32) -> i32 {
    // Add 0.25 turn (pi/2) to map cosine onto sine lookup table.
    let mut x_bits = (x as u32).wrapping_add(0x2000_0000);
    if (x_bits as i32) < 0 {
        x_bits = x_bits.wrapping_add(0x8000_0000);
    }

    let x_mapped = x_bits as i64;
    let index = (x_bits >> FAST_MATH_Q31_SHIFT) as usize;
    let fract = (x_mapped - ((index as i64) << FAST_MATH_Q31_SHIFT)) << 9;

    let a = SIN_TABLE_U32[index] as i32 as i64;
    let b = SIN_TABLE_U32[index + 1] as i32 as i64;

    let mut cos_val = ((0x8000_0000_i64 - fract) * a) >> 32;
    cos_val = ((cos_val << 32) + (fract * b)) >> 32;

    (cos_val << 1) as i32
}
