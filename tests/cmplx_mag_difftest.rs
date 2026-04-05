use fixed_dsp::complex::{cmplx_mag_i16, cmplx_mag_i32};

unsafe extern "C" {
    fn arm_cmplx_mag_q15(pSrc: *const i16, pDst: *mut i16, numSamples: u32);
    fn arm_cmplx_mag_q31(pSrc: *const i32, pDst: *mut i32, numSamples: u32);
}

fn sample_i16_complex(num_samples: usize) -> Vec<i16> {
    let mut out = Vec::with_capacity(num_samples * 2);
    let mut state: u64 = 0x1234_5678_9ABC_DEF0;
    for _ in 0..num_samples {
        state ^= state << 7;
        state ^= state >> 9;
        state ^= state << 8;
        out.push(state as i16);
        state ^= state << 7;
        state ^= state >> 9;
        state ^= state << 8;
        out.push(state as i16);
    }

    if num_samples >= 4 {
        out[0] = i16::MIN;
        out[1] = i16::MIN;
        out[2] = i16::MAX;
        out[3] = 0;
        out[4] = 0;
        out[5] = i16::MAX;
    }

    out
}

fn sample_i32_complex(num_samples: usize) -> Vec<i32> {
    let mut out = Vec::with_capacity(num_samples * 2);
    let mut state: u64 = 0x9E37_79B9_7F4A_7C15;
    for _ in 0..num_samples {
        state ^= state << 7;
        state ^= state >> 9;
        state ^= state << 8;
        out.push(state as i32);
        state ^= state << 7;
        state ^= state >> 9;
        state ^= state << 8;
        out.push(state as i32);
    }

    if num_samples >= 4 {
        out[0] = i32::MIN;
        out[1] = i32::MIN;
        out[2] = i32::MAX;
        out[3] = 0;
        out[4] = 0;
        out[5] = i32::MAX;
    }

    out
}

#[test]
fn cmplx_mag_q15_difftest_against_cmsis() {
    for &num_samples in &[1usize, 2, 3, 4, 7, 16, 31, 64, 127, 256] {
        let input = sample_i16_complex(num_samples);
        let mut rust_out = vec![0i16; num_samples];
        let mut cmsis_out = vec![0i16; num_samples];

        cmplx_mag_i16(&input, &mut rust_out);
        unsafe {
            arm_cmplx_mag_q15(input.as_ptr(), cmsis_out.as_mut_ptr(), num_samples as u32);
        }

        assert_eq!(
            rust_out, cmsis_out,
            "Q15 complex magnitude mismatch for num_samples={}",
            num_samples
        );
    }
}

#[test]
fn cmplx_mag_q31_difftest_against_cmsis() {
    for &num_samples in &[1usize, 2, 3, 4, 7, 16, 31, 64, 127, 256] {
        let input = sample_i32_complex(num_samples);
        let mut rust_out = vec![0i32; num_samples];
        let mut cmsis_out = vec![0i32; num_samples];

        cmplx_mag_i32(&input, &mut rust_out);
        unsafe {
            arm_cmplx_mag_q31(input.as_ptr(), cmsis_out.as_mut_ptr(), num_samples as u32);
        }

        assert_eq!(
            rust_out, cmsis_out,
            "Q31 complex magnitude mismatch for num_samples={}",
            num_samples
        );
    }
}
