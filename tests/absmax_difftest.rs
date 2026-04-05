use fixed_dsp::statistics::{absmax_i16, absmax_i32};

unsafe extern "C" {
    fn arm_absmax_q15(pSrc: *const i16, blockSize: u32, pResult: *mut i16, pIndex: *mut u32);
    fn arm_absmax_q31(pSrc: *const i32, blockSize: u32, pResult: *mut i32, pIndex: *mut u32);
}

fn sample_i16_data() -> Vec<i16> {
    let mut v: Vec<i16> = (0..257)
        .map(|i| (i as i16).wrapping_mul(101).wrapping_sub(11_000))
        .collect();
    v[13] = i16::MIN;
    v[89] = -32_767;
    v
}

fn sample_i32_data() -> Vec<i32> {
    let mut v: Vec<i32> = (0..257)
        .map(|i| (i as i32).wrapping_mul(1_048_583).wrapping_sub(123_456_789))
        .collect();
    v[17] = i32::MIN;
    v[144] = -(i32::MAX - 7);
    v
}

#[test]
fn absmax_q15_difftest_against_cmsis() {
    let input = sample_i16_data();

    for &block_size in &[1usize, 2, 3, 8, 31, 64, 128, 256, 257] {
        let (rust_val, rust_idx) = absmax_i16(&input, block_size);

        let mut cmsis_val = 0i16;
        let mut cmsis_idx = 0u32;
        unsafe {
            arm_absmax_q15(
                input.as_ptr(),
                block_size as u32,
                &mut cmsis_val,
                &mut cmsis_idx,
            );
        }

        assert_eq!(rust_val, cmsis_val, "Q15 absmax value mismatch block_size={}", block_size);
        assert_eq!(
            rust_idx,
            cmsis_idx as usize,
            "Q15 absmax index mismatch block_size={}",
            block_size
        );
    }
}

#[test]
fn absmax_q31_difftest_against_cmsis() {
    let input = sample_i32_data();

    for &block_size in &[1usize, 2, 3, 8, 31, 64, 128, 256, 257] {
        let (rust_val, rust_idx) = absmax_i32(&input, block_size);

        let mut cmsis_val = 0i32;
        let mut cmsis_idx = 0u32;
        unsafe {
            arm_absmax_q31(
                input.as_ptr(),
                block_size as u32,
                &mut cmsis_val,
                &mut cmsis_idx,
            );
        }

        assert_eq!(rust_val, cmsis_val, "Q31 absmax value mismatch block_size={}", block_size);
        assert_eq!(
            rust_idx,
            cmsis_idx as usize,
            "Q31 absmax index mismatch block_size={}",
            block_size
        );
    }
}
