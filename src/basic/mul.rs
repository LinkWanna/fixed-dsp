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

/// Multiply two Q31 fixed-point values using CMSIS `arm_mult_q31` scalar semantics.
///
/// CMSIS computes:
/// `out = ((a * b) >> 32); out = SSAT(out, 31); result = out << 1`.
pub fn mul_i32(a: i32, b: i32) -> i32 {
    let mul = (a as i64) * (b as i64);
    let out = (mul >> 32) as i32;
    let sat = out.clamp(-(1 << 30), (1 << 30) - 1);

    sat << 1
}
