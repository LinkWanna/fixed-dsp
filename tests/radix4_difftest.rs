use fixed_dsp::common::tables::TWIDDLE_TABLE_4096_U16;
use fixed_dsp::transform::{radix4_butterfly_i16, radix4_butterfly_inverse_i16};

unsafe extern "C" {
    static twiddleCoef_4096_q15: [i16; 8192];

    fn arm_radix4_butterfly_q15(
        pSrc: *mut i16,
        fftLen: u32,
        pCoef: *const i16,
        twidCoefModifier: u32,
    );

    fn arm_radix4_butterfly_inverse_q15(
        pSrc: *mut i16,
        fftLen: u32,
        pCoef: *const i16,
        twidCoefModifier: u32,
    );
}

fn get_twiddle_modifier(fft_len: usize) -> u16 {
    match fft_len {
        16 => 256,
        64 => 64,
        256 => 16,
        1024 => 4,
        _ => panic!("Unsupported radix4 FFT length: {}", fft_len),
    }
}

fn sample_i16_buffer(fft_len: usize) -> Vec<i16> {
    (0..(fft_len * 2))
        .map(|index| (index as i16).wrapping_mul(37).wrapping_add(11))
        .collect()
}

#[test]
fn radix4_butterfly_q15_difftest_against_cmsis() {
    for &fft_len in &[16usize, 64, 256, 1024] {
        let mut rust_data = sample_i16_buffer(fft_len);
        let mut cmsis_data = rust_data.clone();

        let twiddle_modifier = get_twiddle_modifier(fft_len);

        radix4_butterfly_i16(
            &mut rust_data,
            fft_len,
            &TWIDDLE_TABLE_4096_U16,
            twiddle_modifier,
        );

        unsafe {
            arm_radix4_butterfly_q15(
                cmsis_data.as_mut_ptr(),
                fft_len as u32,
                twiddleCoef_4096_q15.as_ptr(),
                twiddle_modifier as u32,
            );
        }

        assert_eq!(
            rust_data, cmsis_data,
            "Q15 radix4 butterfly mismatch for fft_len={}",
            fft_len
        );
    }
}

#[test]
fn radix4_butterfly_inverse_q15_difftest_against_cmsis() {
    for &fft_len in &[16usize, 64, 256, 1024] {
        let mut rust_data = sample_i16_buffer(fft_len);
        let mut cmsis_data = rust_data.clone();

        let twiddle_modifier = get_twiddle_modifier(fft_len);

        radix4_butterfly_inverse_i16(
            &mut rust_data,
            fft_len,
            &TWIDDLE_TABLE_4096_U16,
            twiddle_modifier,
        );

        unsafe {
            arm_radix4_butterfly_inverse_q15(
                cmsis_data.as_mut_ptr(),
                fft_len as u32,
                twiddleCoef_4096_q15.as_ptr(),
                twiddle_modifier as u32,
            );
        }

        assert_eq!(
            rust_data, cmsis_data,
            "Q15 radix4 inverse butterfly mismatch for fft_len={}",
            fft_len
        );
    }
}
