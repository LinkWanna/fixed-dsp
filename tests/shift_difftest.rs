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

    for &shift_bits in &shift_values {
        for &block_size in &[1usize, 2, 3, 4, 7, 16, 22] {
            let mut rust_data = test_values[..block_size].to_vec();
            let mut cmsis_data = rust_data.clone();

            shift_i16(&mut rust_data, shift_bits);
            unsafe {
                arm_shift_q15(
                    cmsis_data.as_ptr(),
                    shift_bits,
                    cmsis_data.as_mut_ptr(),
                    block_size as u32,
                );
            }

            for i in 0..block_size {
                let diff = (rust_data[i] as i32 - cmsis_data[i] as i32).abs() as i16;
                max_abs_diff = max_abs_diff.max(diff);

                assert_eq!(
                    rust_data[i], cmsis_data[i],
                    "Q15 shift mismatch idx={}, shift_bits={}: rust={}, cmsis={}",
                    i, shift_bits, rust_data[i], cmsis_data[i]
                );
            }
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

    for &shift_bits in &shift_values {
        for &block_size in &[1usize, 2, 3, 4, 7, 16, 22] {
            let mut rust_data = test_values[..block_size].to_vec();
            let mut cmsis_data = rust_data.clone();

            shift_i32(&mut rust_data, shift_bits);
            unsafe {
                arm_shift_q31(
                    cmsis_data.as_ptr(),
                    shift_bits,
                    cmsis_data.as_mut_ptr(),
                    block_size as u32,
                );
            }

            for i in 0..block_size {
                let diff = (rust_data[i] as i64 - cmsis_data[i] as i64).abs();
                max_abs_diff = max_abs_diff.max(diff);

                assert_eq!(
                    rust_data[i], cmsis_data[i],
                    "Q31 shift mismatch idx={}, shift_bits={}: rust={}, cmsis={}",
                    i, shift_bits, rust_data[i], cmsis_data[i]
                );
            }
        }
    }

    println!("Q31 shift - max_abs_diff: {}", max_abs_diff);
    assert_eq!(max_abs_diff, 0, "Q31 shift results should be identical to CMSIS");
}
