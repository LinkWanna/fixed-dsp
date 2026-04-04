/// Dot product of two Q15 vectors.
///
/// Each multiplication is in 2.30 format and accumulated into a 64-bit sum.
/// The returned value matches CMSIS 34.30 accumulator semantics.
pub fn dot_i16(a: &[i16], b: &[i16]) -> i64 {
    assert_eq!(a.len(), b.len(), "input slices must have the same length");

    let mut sum: i64 = 0;
    for (&x, &y) in a.iter().zip(b.iter()) {
        sum = sum.wrapping_add((x as i32 as i64) * (y as i32 as i64));
    }

    sum
}

/// Dot product of two Q31 vectors.
///
/// CMSIS truncates each 2.62 product to 2.48 by discarding 14 LSBs,
/// then accumulates into a 64-bit 16.48 accumulator.
pub fn dot_i32(a: &[i32], b: &[i32]) -> i64 {
    assert_eq!(a.len(), b.len(), "input slices must have the same length");

    let mut sum: i64 = 0;
    for (&x, &y) in a.iter().zip(b.iter()) {
        let prod = ((x as i64) * (y as i64)) >> 14;
        sum = sum.wrapping_add(prod);
    }

    sum
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dot_i16_basic() {
        let a = [0x4000_i16, 0x4000_i16]; // 0.5, 0.5
        let b = [0x4000_i16, 0x2000_i16]; // 0.5, 0.25

        let out = dot_i16(&a, &b);
        assert_eq!(out, 402_653_184);
    }

    #[test]
    fn dot_i32_basic() {
        let a = [0x4000_0000_i32, 0x4000_0000_i32];
        let b = [0x4000_0000_i32, 0x2000_0000_i32];

        let out = dot_i32(&a, &b);
        let expected = (((0x4000_0000_i64 * 0x4000_0000_i64) >> 14)
            + ((0x4000_0000_i64 * 0x2000_0000_i64) >> 14)) as i64;
        assert_eq!(out, expected);
    }
}
