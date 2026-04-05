use fixed_dsp::basic::{shift_i16, shift_i32};

unsafe extern "C" {
    fn arm_shift_q15(pSrc: *const i16, shiftBits: i8, pDst: *mut i16, blockSize: u32);
    fn arm_shift_q31(pSrc: *const i32, shiftBits: i8, pDst: *mut i32, blockSize: u32);
}

#[test]
fn shift_q15_difftest_against_cmsis() {
    let test_values: Vec<i16> = vec![
        i16::MIN,
        i16::MIN + 1,
        -32767,
        -1024,
        -512,
        -256,
        -128,
        -64,
        -32,
        -16,
        -1,
        0,
        1,
        16,
        32,
        64,
        128,
        256,
        512,
        1024,
        32767,
        i16::MAX,
    ];

    let shift_values: Vec<i8> = vec![-8, -4, -2, -1, 0, 1, 2, 4, 8];

    let mut max_abs_diff: i16 = 0;

    for &val in &test_values {
        for &shift_bits in &shift_values {
            let rust_result = shift_i16(val, shift_bits);

            let mut cmsis_result = 0_i16;
            unsafe {
                arm_shift_q15(&val, shift_bits, &mut cmsis_result, 1);
            }

            let diff = (rust_result as i32 - cmsis_result as i32).abs() as i16;
            max_abs_diff = max_abs_diff.max(diff);

            assert_eq!(
                rust_result, cmsis_result,
                "Q15 shift mismatch for val={}, shift_bits={}: rust={}, cmsis={}",
                val, shift_bits, rust_result, cmsis_result
            );
        }
    }

    println!("Q15 shift - max_abs_diff: {}", max_abs_diff);
    assert_eq!(max_abs_diff, 0, "Q15 shift results should be identical to CMSIS");
}

#[test]
fn shift_q31_difftest_against_cmsis() {
    let test_values: Vec<i32> = vec![
        i32::MIN,
        i32::MIN + 1,
        -2147480000,
        -1048576,
        -65536,
        -256,
        -128,
        -64,
        -32,
        -16,
        -1,
        0,
        1,
        16,
        32,
        64,
        128,
        256,
        65536,
        1048576,
        2147480000,
        i32::MAX,
    ];

    let shift_values: Vec<i8> = vec![-16, -8, -4, -2, -1, 0, 1, 2, 4, 8, 16];

    let mut max_abs_diff: i64 = 0;

    for &val in &test_values {
        for &shift_bits in &shift_values {
            if val == i32::MIN && shift_bits > 0 {
                continue;
            }

            let rust_result = shift_i32(val, shift_bits);

            let mut cmsis_result = 0_i32;
            unsafe {
                arm_shift_q31(&val, shift_bits, &mut cmsis_result, 1);
            }

            let diff = (rust_result as i64 - cmsis_result as i64).abs();
            max_abs_diff = max_abs_diff.max(diff);

            assert_eq!(
                rust_result, cmsis_result,
                "Q31 shift mismatch for val={}, shift_bits={}: rust={}, cmsis={}",
                val, shift_bits, rust_result, cmsis_result
            );
        }
    }

    println!("Q31 shift - max_abs_diff: {}", max_abs_diff);
    assert_eq!(max_abs_diff, 0, "Q31 shift results should be identical to CMSIS");
}
