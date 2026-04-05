use fixed_dsp::basic::{sqrt_i16, sqrt_i32};
use fixed_dsp::common::error::Error;

unsafe extern "C" {
    fn arm_sqrt_q15(input: i16, pOut: *mut i16) -> i32;
    fn arm_sqrt_q31(input: i32, pOut: *mut i32) -> i32;
}

fn q15_samples() -> [i16; 8] {
    [0xFFFFu16 as i16, 0x0000, 0x0001, 0x0002, 0x0100, 0x1000, 0x4000, 0x7FFF]
}

fn q31_samples() -> [i32; 9] {
    [
        0xFFFFFFFFu32 as i32,
        0x00000000,
        0x00000001,
        0x00010000,
        0x10000000,
        0x40000000,
        0x60000000,
        0x7FFFFFFE,
        0x7FFFFFFF,
    ]
}

#[test]
fn sqrt_q15_difftest_against_cmsis() {
    for &input in &q15_samples() {
        let rust_result = sqrt_i16(input);
        let mut cmsis_out = 0_i16;
        let cmsis_status = unsafe { arm_sqrt_q15(input, &mut cmsis_out) };

        if input < 0 {
            assert_eq!(rust_result, Err(Error::NanInf));
            assert_eq!(cmsis_status, -1);
            assert_eq!(cmsis_out, 0);
        } else {
            assert_eq!(rust_result, Ok(cmsis_out), "Q15 sqrt mismatch for input={:#06X}", input as u16);
            assert_eq!(cmsis_status, 0);
        }
    }
}

#[test]
fn sqrt_q31_difftest_against_cmsis() {
    for &input in &q31_samples() {
        let rust_result = sqrt_i32(input);
        let mut cmsis_out = 0_i32;
        let cmsis_status = unsafe { arm_sqrt_q31(input, &mut cmsis_out) };

        if input < 0 {
            assert_eq!(rust_result, Err(Error::NanInf));
            assert_eq!(cmsis_status, -1);
            assert_eq!(cmsis_out, 0);
        } else {
            assert_eq!(rust_result, Ok(cmsis_out), "Q31 sqrt mismatch for input={:#010X}", input as u32);
            assert_eq!(cmsis_status, 0);
        }
    }
}
