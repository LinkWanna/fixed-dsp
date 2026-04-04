use fixed_dsp::basic::{mul_i16, mul_i32};

unsafe extern "C" {
    fn arm_mult_q15(pSrcA: *const i16, pSrcB: *const i16, pDst: *mut i16, blockSize: u32);
    fn arm_mult_q31(pSrcA: *const i32, pSrcB: *const i32, pDst: *mut i32, blockSize: u32);
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
fn mul_q15_difftest_against_cmsis() {
    let inputs: Vec<i16> = (i16::MIN..=i16::MAX).step_by(257).collect();

    let mut max_abs_diff: i16 = 0;

    for &a in &inputs {
        for &b in &inputs {
            let rust_result = mul_i16(a, b);

            let mut cmsis_result = 0_i16;
            unsafe {
                arm_mult_q15(&a, &b, &mut cmsis_result, 1);
            }

            let diff = (rust_result as i32 - cmsis_result as i32).abs() as i16;
            max_abs_diff = max_abs_diff.max(diff);

            assert_eq!(
                rust_result, cmsis_result,
                "Q15 multiplication mismatch for a={}, b={}: rust={}, cmsis={}",
                a, b, rust_result, cmsis_result
            );
        }
    }

    println!("Q15 multiplication - max_abs_diff: {}", max_abs_diff);
    assert_eq!(max_abs_diff, 0, "Q15 results should be identical to CMSIS");
}

#[test]
fn mul_q31_difftest_against_cmsis() {
    let a_inputs: Vec<i32> = sample_q31_inputs()
        .into_iter()
        .step_by(97)
        .take(128)
        .collect();
    let b_inputs: Vec<i32> = sample_q31_inputs()
        .into_iter()
        .step_by(193)
        .take(128)
        .collect();

    let mut max_abs_diff: i64 = 0;

    for &a in &a_inputs {
        for &b in &b_inputs {
            let rust_result = mul_i32(a, b);

            let mut cmsis_result = 0_i32;
            unsafe {
                arm_mult_q31(&a, &b, &mut cmsis_result, 1);
            }

            let diff = (rust_result as i64 - cmsis_result as i64).abs();
            max_abs_diff = max_abs_diff.max(diff);

            assert_eq!(
                rust_result, cmsis_result,
                "Q31 multiplication mismatch for a={}, b={}: rust={}, cmsis={}",
                a, b, rust_result, cmsis_result
            );
        }
    }

    println!("Q31 multiplication - max_abs_diff: {}", max_abs_diff);
    assert_eq!(max_abs_diff, 0, "Q31 results should be identical to CMSIS");
}
