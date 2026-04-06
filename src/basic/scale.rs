use crate::{sat_i16, sat_i32};

/// In-place Q15 scale, matching CMSIS `arm_scale_q15` scalar behavior.
pub fn scale_i16(data: &mut [i16], scale: i16, shift: i8) {
    let k_shift = 15 - shift;

    for sample in data.iter_mut() {
        *sample = sat_i16(((*sample as i32) * (scale as i32)) >> k_shift);
    }
}

/// In-place Q31 scale, matching CMSIS `arm_scale_q31` scalar behavior.
pub fn scale_i32(data: &mut [i32], scale: i32, shift: i8) {
    let k_shift = shift + 1;

    if k_shift >= 0 {
        for sample in data.iter_mut() {
            let in_val = ((*sample as i64) * (scale as i64)) >> 32;
            *sample = sat_i32(in_val << k_shift);
        }
    } else {
        for sample in data.iter_mut() {
            let in_val = ((*sample as i64) * (scale as i64)) >> 32;
            *sample = sat_i32(in_val >> -k_shift);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn scale_i16_in_place_changes_buffer() {
        let mut data = [0x4000_i16, -0x4000_i16, 0x2000_i16];
        scale_i16(&mut data, 0x4000, 0);
        assert_eq!(data, [0x2000, -0x2000, 0x1000]);
    }

    #[test]
    fn scale_i32_in_place_changes_buffer() {
        let mut data = [0x4000_0000_i32, -0x4000_0000_i32, 0x2000_0000_i32];
        scale_i32(&mut data, 0x4000_0000, 0);
        assert_eq!(data, [0x2000_0000, -0x2000_0000, 0x1000_0000]);
    }
}
