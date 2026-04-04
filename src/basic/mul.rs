/// Multiply two Q15 fixed-point values and saturate result to Q15 range
///
/// a and b are in Q15 format (1.15 fixed-point), range [-32768, 32767] representing [-1.0, ~0.9999]
/// The multiplication produces a 2.30 format result, which must be rescaled back to Q15
/// by right-shifting 15 bits and saturating to the valid Q15 range.
///
/// Formula: result = saturate((a * b) >> 15)
pub fn mul_i16(a: i16, b: i16) -> i16 {
    let mul = (a as i32) * (b as i32);
    let result = mul >> 15;

    // Saturate to i16 range: [-2^15, 2^15-1] =  [-32768, 32767]
    result.clamp(i16::MIN as i32, i16::MAX as i32) as i16
}

/// Multiply two Q31 fixed-point values and saturate result to Q31 range
///
/// a and b are in Q31 format (1.31 fixed-point), range [-2147483648, 2147483647] representing [-1.0, ~0.9999]
/// The multiplication produces a 2.62 format result, which must be rescaled back to Q31
/// by right-shifting 31 bits and saturating to the valid Q31 range.
///
/// Formula: result = saturate((a * b) >> 31)
pub fn mul_i32(a: i32, b: i32) -> i32 {
    let mul = (a as i64) * (b as i64);
    let result = mul >> 31;

    // Saturate to i32 range: [-2^31, 2^31-1] = [-2147483648, 2147483647]
    result.clamp(i32::MIN as i64, i32::MAX as i64) as i32
}
