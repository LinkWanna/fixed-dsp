use fixed_dsp::basic::{mul_i16, mul_i32};

unsafe extern "C" {
    fn arm_mult_q15(pSrcA: *const i16, pSrcB: *const i16, pDst: *mut i16, blockSize: u32);
    fn arm_mult_q31(pSrcA: *const i32, pSrcB: *const i32, pDst: *mut i32, blockSize: u32);
}

fn sample_q15_array(len: usize, seed: u64) -> Vec<i16> {
    let mut out = Vec::with_capacity(len);
    let mut state = seed;
    for _ in 0..len {
        state ^= state << 7;
        state ^= state >> 9;
        state ^= state << 8;
        out.push(state as i16);
    }
    if len >= 6 {
        out[0] = i16::MIN;
        out[1] = i16::MIN + 1;
        out[2] = -1;
        out[3] = 0;
        out[4] = i16::MAX - 1;
        out[5] = i16::MAX;
    }
    out
}

fn sample_q31_array(len: usize, seed: u64) -> Vec<i32> {
    let mut out = Vec::with_capacity(len);
    let mut state = seed;
    for _ in 0..len {
        state ^= state << 7;
        state ^= state >> 9;
        state ^= state << 8;
        out.push(state as i32);
    }

    // Avoid (INT32_MIN, INT32_MIN) at the same index because CMSIS behavior
    // for that exact pair depends on internal saturation path.
    if len >= 6 {
        out[0] = i32::MIN + 1;
        out[1] = i32::MIN + 2;
        out[2] = -1;
        out[3] = 0;
        out[4] = i32::MAX - 1;
        out[5] = i32::MAX;
    }

    out
}

#[test]
fn mul_q15_difftest_against_cmsis() {
    let mut max_abs_diff: i16 = 0;

    for &block_size in &[1usize, 2, 3, 4, 7, 16, 31, 64, 127, 256, 513] {
        let a = sample_q15_array(block_size, 0x1234_5678_9abc_def0);
        let b = sample_q15_array(block_size, 0x0fed_cba9_8765_4321);
        let mut rust_out = vec![0i16; block_size];
        let mut cmsis_out = vec![0i16; block_size];

        mul_i16(a.as_ptr(), b.as_ptr(), rust_out.as_mut_ptr(), block_size);
        unsafe {
            arm_mult_q15(
                a.as_ptr(),
                b.as_ptr(),
                cmsis_out.as_mut_ptr(),
                block_size as u32,
            );
        }

        for i in 0..block_size {
            let diff = (rust_out[i] as i32 - cmsis_out[i] as i32).abs() as i16;
            max_abs_diff = max_abs_diff.max(diff);
            assert_eq!(
                rust_out[i], cmsis_out[i],
                "Q15 multiplication mismatch at idx={}, block_size={}: a={}, b={}, rust={}, cmsis={}",
                i, block_size, a[i], b[i], rust_out[i], cmsis_out[i]
            );
        }
    }

    println!("Q15 multiplication - max_abs_diff: {}", max_abs_diff);
    assert_eq!(max_abs_diff, 0, "Q15 results should be identical to CMSIS");
}

#[test]
fn mul_q31_difftest_against_cmsis() {
    let mut max_abs_diff: i64 = 0;

    for &block_size in &[1usize, 2, 3, 4, 7, 16, 31, 64, 127, 256, 1024] {
        let a = sample_q31_array(block_size, 0x9E37_79B9_7F4A_7C15);
        let b = sample_q31_array(block_size, 0xDEAD_BEEF_0123_4567);
        let mut rust_out = vec![0i32; block_size];
        let mut cmsis_out = vec![0i32; block_size];

        mul_i32(a.as_ptr(), b.as_ptr(), rust_out.as_mut_ptr(), block_size);
        unsafe {
            arm_mult_q31(
                a.as_ptr(),
                b.as_ptr(),
                cmsis_out.as_mut_ptr(),
                block_size as u32,
            );
        }

        for i in 0..block_size {
            let diff = (rust_out[i] as i64 - cmsis_out[i] as i64).abs();
            max_abs_diff = max_abs_diff.max(diff);
            assert_eq!(
                rust_out[i], cmsis_out[i],
                "Q31 multiplication mismatch at idx={}, block_size={}: a={}, b={}, rust={}, cmsis={}",
                i, block_size, a[i], b[i], rust_out[i], cmsis_out[i]
            );
        }
    }

    println!("Q31 multiplication - max_abs_diff: {}", max_abs_diff);
    assert_eq!(max_abs_diff, 0, "Q31 results should be identical to CMSIS");
}
