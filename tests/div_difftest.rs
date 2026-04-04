use core::ffi::c_int;

use fixed_dsp::basic::{div_i16, div_i32};
use fixed_dsp::common::error::Error;

unsafe extern "C" {
    fn arm_divide_q15(
        numerator: i16,
        denominator: i16,
        quotient: *mut i16,
        shift: *mut i16,
    ) -> c_int;
    fn arm_divide_q31(
        numerator: i32,
        denominator: i32,
        quotient: *mut i32,
        shift: *mut i16,
    ) -> c_int;
}

const ARM_MATH_SUCCESS: c_int = 0;
const ARM_MATH_NANINF: c_int = -4;

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
fn div_q15_difftest_against_cmsis() {
    let numerators: Vec<i16> = (i16::MIN..=i16::MAX).step_by(257).collect();
    let denominators: Vec<i16> = (i16::MIN..=i16::MAX).step_by(257).collect();

    for &n in &numerators {
        for &d in &denominators {
            let rust_ret = div_i16(n, d);

            let mut q_ref = 0_i16;
            let mut s_ref = 0_i16;
            let ref_status = unsafe { arm_divide_q15(n, d, &mut q_ref, &mut s_ref) };

            match rust_ret {
                Ok((q_rust, s_rust)) => {
                    assert_eq!(
                        ref_status, ARM_MATH_SUCCESS,
                        "q15 status mismatch for n={}, d={}",
                        n, d
                    );
                    assert_eq!(q_rust, q_ref, "q15 quotient mismatch for n={}, d={}", n, d);
                    assert_eq!(s_rust, s_ref, "q15 shift mismatch for n={}, d={}", n, d);
                }
                Err(e) => {
                    assert_eq!(e, Error::NanInf, "q15 error mismatch for n={}, d={}", n, d);
                    assert_eq!(
                        ref_status, ARM_MATH_NANINF,
                        "q15 ref status mismatch for n={}, d={}",
                        n, d
                    );

                    let sign = (n < 0) ^ (d < 0);
                    let q_sat = if sign { i16::MIN } else { i16::MAX };
                    assert_eq!(
                        q_ref, q_sat,
                        "q15 saturated quotient mismatch for n={}, d={}",
                        n, d
                    );
                    assert_eq!(s_ref, 0, "q15 shift-on-error mismatch for n={}, d={}", n, d);
                }
            }
        }
    }
}

#[test]
fn div_q31_difftest_against_cmsis() {
    let numerators: Vec<i32> = sample_q31_inputs()
        .into_iter()
        .step_by(97)
        .take(256)
        .collect();
    let mut denominators: Vec<i32> = sample_q31_inputs()
        .into_iter()
        .step_by(193)
        .take(256)
        .collect();
    denominators.extend([0, 1, -1, i32::MIN, i32::MAX, i32::MIN + 1, i32::MAX - 1]);

    for &n in &numerators {
        for &d in &denominators {
            let rust_ret = div_i32(n, d);

            let mut q_ref = 0_i32;
            let mut s_ref = 0_i16;
            let ref_status = unsafe { arm_divide_q31(n, d, &mut q_ref, &mut s_ref) };

            match rust_ret {
                Ok((q_rust, s_rust)) => {
                    assert_eq!(
                        ref_status, ARM_MATH_SUCCESS,
                        "q31 status mismatch for n={}, d={}",
                        n, d
                    );
                    assert_eq!(q_rust, q_ref, "q31 quotient mismatch for n={}, d={}", n, d);
                    assert_eq!(s_rust, s_ref, "q31 shift mismatch for n={}, d={}", n, d);
                }
                Err(e) => {
                    assert_eq!(e, Error::NanInf, "q31 error mismatch for n={}, d={}", n, d);
                    assert_eq!(
                        ref_status, ARM_MATH_NANINF,
                        "q31 ref status mismatch for n={}, d={}",
                        n, d
                    );

                    let sign = (n < 0) ^ (d < 0);
                    let q_sat = if sign { i32::MIN } else { i32::MAX };
                    assert_eq!(
                        q_ref, q_sat,
                        "q31 saturated quotient mismatch for n={}, d={}",
                        n, d
                    );
                    assert_eq!(s_ref, 0, "q31 shift-on-error mismatch for n={}, d={}", n, d);
                }
            }
        }
    }
}
