#[inline]
fn abs_sat_i16(x: i16) -> i16 {
    if x == i16::MIN { i16::MAX } else { x.abs() }
}

#[inline]
fn abs_sat_i32(x: i32) -> i32 {
    if x == i32::MIN { i32::MAX } else { x.abs() }
}

/// Maximum absolute value and its index for a Q15 slice prefix.
///
/// Only the first `size` samples are processed, matching CMSIS API style.
/// Returns the first index when there are ties.
pub fn absmax_i16(input: &[i16], size: usize) -> (i16, usize) {
    assert!(size <= input.len(), "size must be <= input length");

    let mut out = abs_sat_i16(input[0]);
    let mut out_index = 0usize;

    for (index, &value) in input[..size].iter().enumerate().skip(1) {
        let candidate = abs_sat_i16(value);
        if candidate > out {
            out = candidate;
            out_index = index;
        }
    }

    (out, out_index)
}

/// Maximum absolute value and its index for a Q31 slice prefix.
///
/// Only the first `size` samples are processed, matching CMSIS API style.
/// Returns the first index when there are ties.
pub fn absmax_i32(input: &[i32], size: usize) -> (i32, usize) {
    assert!(size <= input.len(), "size must be <= input length");

    let mut out = abs_sat_i32(input[0]);
    let mut out_index = 0usize;

    for (index, &value) in input[..size].iter().enumerate().skip(1) {
        let candidate = abs_sat_i32(value);
        if candidate > out {
            out = candidate;
            out_index = index;
        }
    }

    (out, out_index)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn absmax_i16_handles_q15_min_saturation() {
        let data = [0, -5, i16::MIN, 32766, -32767];
        let (max, idx) = absmax_i16(&data, data.len());
        assert_eq!(max, i16::MAX);
        assert_eq!(idx, 2);
    }

    #[test]
    fn absmax_i16_uses_first_index_on_tie() {
        let data = [-7, 7, -7, 3];
        let (max, idx) = absmax_i16(&data, data.len());
        assert_eq!(max, 7);
        assert_eq!(idx, 0);
    }

    #[test]
    fn absmax_i32_handles_q31_min_saturation() {
        let data = [0, -5, i32::MIN, i32::MAX - 1, -(i32::MAX - 1)];
        let (max, idx) = absmax_i32(&data, data.len());
        assert_eq!(max, i32::MAX);
        assert_eq!(idx, 2);
    }

    #[test]
    fn absmax_i32_respects_size_prefix() {
        let data = [10, -20, 30, -40, 1_000_000];
        let (max, idx) = absmax_i32(&data, 4);
        assert_eq!(max, 40);
        assert_eq!(idx, 3);
    }
}
