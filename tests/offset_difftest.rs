use fixed_dsp::basic::{offset_i16, offset_i32};

unsafe extern "C" {
    fn arm_offset_q15(pSrc: *const i16, offset: i16, pDst: *mut i16, blockSize: u32);
    fn arm_offset_q31(pSrc: *const i32, offset: i32, pDst: *mut i32, blockSize: u32);
}

fn sample_i16_data() -> Vec<i16> {
    (0..257)
        .map(|i| (i as i16).wrapping_mul(211).wrapping_sub(13_579))
        .collect()
}

fn sample_i32_data() -> Vec<i32> {
    (0..257)
        .map(|i| (i as i32).wrapping_mul(1_048_583).wrapping_sub(987_654_321))
        .collect()
}

#[test]
fn offset_q15_difftest_against_cmsis_in_place() {
    let src = sample_i16_data();
    let offsets = [0_i16, 1_i16, -1_i16, 0x4000_i16, -0x4000_i16, i16::MAX, i16::MIN];

    for &offset in &offsets {
        let mut rust_data = src.clone();
        let mut cmsis_data = src.clone();

        offset_i16(&mut rust_data, offset);

        unsafe {
            arm_offset_q15(
                cmsis_data.as_ptr(),
                offset,
                cmsis_data.as_mut_ptr(),
                cmsis_data.len() as u32,
            );
        }

        assert_eq!(
            rust_data, cmsis_data,
            "Q15 offset mismatch offset={}",
            offset
        );
    }
}

#[test]
fn offset_q31_difftest_against_cmsis_in_place() {
    let src = sample_i32_data();
    let offsets = [
        0_i32,
        1_i32,
        -1_i32,
        0x4000_0000_i32,
        -0x4000_0000_i32,
        i32::MAX,
        i32::MIN,
    ];

    for &offset in &offsets {
        let mut rust_data = src.clone();
        let mut cmsis_data = src.clone();

        offset_i32(&mut rust_data, offset);

        unsafe {
            arm_offset_q31(
                cmsis_data.as_ptr(),
                offset,
                cmsis_data.as_mut_ptr(),
                cmsis_data.len() as u32,
            );
        }

        assert_eq!(
            rust_data, cmsis_data,
            "Q31 offset mismatch offset={}",
            offset
        );
    }
}
