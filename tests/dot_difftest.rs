use fixed_dsp::basic::{dot_i16, dot_i32};

unsafe extern "C" {
    fn arm_dot_prod_q15(pSrcA: *const i16, pSrcB: *const i16, blockSize: u32, result: *mut i64);
    fn arm_dot_prod_q31(pSrcA: *const i32, pSrcB: *const i32, blockSize: u32, result: *mut i64);
}

fn sample_q31_inputs() -> Vec<i32> {
    let mut data = Vec::with_capacity(20004);
    data.extend([i32::MIN, i32::MIN + 1, -1, 0, 1, i32::MAX - 1, i32::MAX]);

    let mut state: u64 = 0x9E37_79B9_7F4A_7C15;
    for _ in 0..20000 {
        state ^= state << 7;
        state ^= state >> 9;
        state ^= state << 8;
        data.push(state as i32);
    }

    data
}

#[test]
fn dot_q15_difftest_against_cmsis() {
    let data_a: Vec<i16> = (i16::MIN..=i16::MAX).step_by(97).collect();
    let data_b: Vec<i16> = (i16::MIN..=i16::MAX).step_by(193).collect();

    let len = data_a.len().min(data_b.len()).min(512);
    let a = &data_a[..len];
    let b = &data_b[..len];

    let rust_result = dot_i16(a, b);

    let mut cmsis_result = 0_i64;
    unsafe {
        arm_dot_prod_q15(a.as_ptr(), b.as_ptr(), len as u32, &mut cmsis_result);
    }

    assert_eq!(rust_result, cmsis_result, "Q15 dot product mismatch");
}

#[test]
fn dot_q31_difftest_against_cmsis() {
    let pool_a = sample_q31_inputs();
    let pool_b = sample_q31_inputs();

    // Keep length below 2^16 as documented by CMSIS to avoid accumulator overflow risk.
    let a: Vec<i32> = pool_a.into_iter().step_by(73).take(1024).collect();
    let b: Vec<i32> = pool_b.into_iter().step_by(151).take(1024).collect();

    let len = a.len().min(b.len());

    let rust_result = dot_i32(&a[..len], &b[..len]);

    let mut cmsis_result = 0_i64;
    unsafe {
        arm_dot_prod_q31(a.as_ptr(), b.as_ptr(), len as u32, &mut cmsis_result);
    }

    assert_eq!(rust_result, cmsis_result, "Q31 dot product mismatch");
}
