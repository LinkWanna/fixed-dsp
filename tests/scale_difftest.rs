use fixed_dsp::basic::{scale_i16, scale_i32};

unsafe extern "C" {
    fn arm_scale_q15(
        pSrc: *const i16,
        scaleFract: i16,
        shift: i8,
        pDst: *mut i16,
        blockSize: u32,
    );

    fn arm_scale_q31(
        pSrc: *const i32,
        scaleFract: i32,
        shift: i8,
        pDst: *mut i32,
        blockSize: u32,
    );
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
fn scale_q15_difftest_against_cmsis_in_place() {
    let src = sample_i16_data();
    let scales = [0x4000_i16, 0x7fff_i16, -0x4000_i16, -0x7fff_i16];
    let shifts = [-2_i8, -1_i8, 0_i8, 1_i8, 2_i8, 7_i8];

    for &scale in &scales {
        for &shift in &shifts {
            let mut rust_data = src.clone();
            let mut cmsis_data = src.clone();

            scale_i16(&mut rust_data, scale, shift);

            unsafe {
                arm_scale_q15(
                    cmsis_data.as_ptr(),
                    scale,
                    shift,
                    cmsis_data.as_mut_ptr(),
                    cmsis_data.len() as u32,
                );
            }

            assert_eq!(
                rust_data, cmsis_data,
                "Q15 scale mismatch scale={}, shift={}",
                scale, shift
            );
        }
    }
}

#[test]
fn scale_q31_difftest_against_cmsis_in_place() {
    let src = sample_i32_data();
    let scales = [
        0x4000_0000_i32,
        0x7fff_ffff_i32,
        -0x4000_0000_i32,
        i32::MIN + 1,
    ];
    let shifts = [-2_i8, -1_i8, 0_i8, 1_i8, 2_i8, 7_i8];

    for &scale in &scales {
        for &shift in &shifts {
            let mut rust_data = src.clone();
            let mut cmsis_data = src.clone();

            scale_i32(&mut rust_data, scale, shift);

            unsafe {
                arm_scale_q31(
                    cmsis_data.as_ptr(),
                    scale,
                    shift,
                    cmsis_data.as_mut_ptr(),
                    cmsis_data.len() as u32,
                );
            }

            assert_eq!(
                rust_data, cmsis_data,
                "Q31 scale mismatch scale={}, shift={}",
                scale, shift
            );
        }
    }
}
