use crate::{sat_i16, sat_i32};

pub fn mul_i16(a: &[i16], b: &[i16], output: &mut [i16]) {
    output
        .iter_mut()
        .zip(a.iter().zip(b.iter()))
        .for_each(|(out, (&x, &y))| {
            *out = sat_i16((x as i32 * y as i32) >> 15);
        });
}

pub fn mul_i32(a: &[i32], b: &[i32], output: &mut [i32]) {
    output
        .iter_mut()
        .zip(a.iter().zip(b.iter()))
        .for_each(|(out, (&x, &y))| {
            *out = sat_i32((x as i64 * y as i64) >> 31);
        });
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mul_i16_basic_q15_scaling() {
        let a = [0x4000_i16, 0x4000_i16, -0x4000_i16]; // 0.5, 0.5, -0.5
        let b = [0x4000_i16, 0x2000_i16, 0x4000_i16]; // 0.5, 0.25, 0.5
        let mut out = [0_i16; 3];

        mul_i16(&a, &b, &mut out);

        assert_eq!(out, [0x2000, 0x1000, -0x2000]); // 0.25, 0.125, -0.25
    }

    #[test]
    fn mul_i16_saturates_at_positive_limit() {
        let a = [i16::MIN];
        let b = [i16::MIN];
        let mut out = [0_i16; 1];

        mul_i16(&a, &b, &mut out);

        assert_eq!(out[0], i16::MAX);
    }

    #[test]
    fn mul_i32_basic_q31_scaling() {
        let a = [0x4000_0000_i32, 0x4000_0000_i32, -0x4000_0000_i32]; // 0.5, 0.5, -0.5
        let b = [0x4000_0000_i32, 0x2000_0000_i32, 0x4000_0000_i32]; // 0.5, 0.25, 0.5
        let mut out = [0_i32; 3];

        mul_i32(&a, &b, &mut out);

        assert_eq!(out, [0x2000_0000, 0x1000_0000, -0x2000_0000]);
    }

    #[test]
    fn mul_i32_saturates_at_positive_limit() {
        let a = [i32::MIN];
        let b = [i32::MIN];
        let mut out = [0_i32; 1];

        mul_i32(&a, &b, &mut out);

        assert_eq!(out[0], i32::MAX);
    }
}
