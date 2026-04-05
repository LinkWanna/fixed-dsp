/// Shift Q15 vector elements by a specified number of bits using CMSIS `arm_shift_q15` semantics.
///
/// - Positive `shift_bits`: left shift with saturation to Q15 range.
/// - Negative `shift_bits`: right shift by `-shift_bits`.
pub fn shift_i16(val: i16, shift_bits: i8) -> i16 {
    let sign = (shift_bits as u8) & 0x80;

    if sign == 0 {
        // Left shift: saturate to Q15
        let shifted = (val as i32) << shift_bits;
        shifted.clamp(i16::MIN as i32, i16::MAX as i32) as i16
    } else {
        // Right shift by -shift_bits
        val >> (-shift_bits as u8)
    }
}

/// Shift Q31 vector elements by a specified number of bits using CMSIS `arm_shift_q31` semantics.
///
/// - Positive `shift_bits`: left shift with overflow detection and saturation.
/// - Negative `shift_bits`: right shift by `-shift_bits`.
pub fn shift_i32(val: i32, shift_bits: i8) -> i32 {
    let sign = (shift_bits as u8) & 0x80;

    if sign == 0 {
        // Left shift with overflow saturation
        let out = val << shift_bits;
        // Check for overflow: if val != (out >> shift_bits), saturation occurred
        if val != (out >> shift_bits) {
            // Return 0x7FFFFFFF if positive, 0x80000000 if negative
            0x7FFF_FFFF ^ (val >> 31)
        } else {
            out
        }
    } else {
        // Right shift by -shift_bits
        val >> (-shift_bits as u8)
    }
}
