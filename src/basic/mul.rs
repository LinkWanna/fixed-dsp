/// Multiply two Q15 fixed-point values and saturate result to Q15 range
///
/// a and b are in Q15 format (1.15 fixed-point), range [-32768, 32767] representing [-1.0, ~0.9999]
/// The multiplication produces a 2.30 format result, which must be rescaled back to Q15
/// by right-shifting 15 bits and saturating to the valid Q15 range.
///
/// Formula: result = saturate((a * b) >> 15)
pub fn mul_i16(a: i16, b: i16) -> i16 {
    let mul = (a as i32) * (b as i32);
    let result = (mul >> 15) as i32;

    // Saturate to i16 range [-32768, 32767]
    result.clamp(i16::MIN as i32, i16::MAX as i32) as i16
}

/// Multiply two Q31 fixed-point values and saturate result to Q31 range
///
/// a and b are in Q31 format (1.31 fixed-point), range [-2147483648, 2147483647] representing [-1.0, ~0.9999]
/// The multiplication produces a 2.62 format result.
///
/// Algorithm (matching CMSIS arm_mult_q31):
/// 1. Compute 64-bit product: a * b
/// 2. Right-shift by 32 bits
/// 3. Saturate to 31-bit signed range [-2^30, 2^30-1]
/// 4. Left-shift by 1 bit to restore to Q31 format
///
/// This is different from Q15 which shifts by 15 directly; Q31 uses an intermediate
/// saturation step to handle the special case of i32::MIN * i32::MIN correctly.
pub fn mul_i32(a: i32, b: i32) -> i32 {
    let mul = (a as i64) * (b as i64);
    let intermediate = (mul >> 32) as i32;

    // Saturate to 31-bit signed range: [-2^30, 2^30-1] = [-1073741824, 1073741823]
    let saturated = intermediate.clamp(-1073741824, 1073741823);

    // Left-shift by 1 to restore to Q31 format
    saturated << 1
}
