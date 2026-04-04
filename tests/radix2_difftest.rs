use fixed_dsp::transform::{radix2_butterfly_i16, radix2_butterfly_inverse_i16};

unsafe extern "C" {
    // Twiddle coefficient table for 4096-point FFT (used for all sizes)
    static twiddleCoef_4096_q15: [i16; 8192];

    fn arm_radix2_butterfly_q15(
        pSrc: *mut i16,
        fftLen: u32,
        pCoef: *const i16,
        twidCoefModifier: u16,
    );

    fn arm_radix2_butterfly_inverse_q15(
        pSrc: *mut i16,
        fftLen: u32,
        pCoef: *const i16,
        twidCoefModifier: u16,
    );
}

/// Get the twiddle coefficient modifier for a given FFT size
/// Formula: twidCoefModifier = 4096 / fftLen
fn get_twiddle_modifier(fft_len: usize) -> u16 {
    match fft_len {
        16 => 256,
        32 => 128,
        64 => 64,
        128 => 32,
        256 => 16,
        512 => 8,
        1024 => 4,
        2048 => 2,
        4096 => 1,
        _ => panic!("Unsupported FFT length: {}", fft_len),
    }
}

fn sample_i16_buffer(fft_len: usize) -> Vec<i16> {
    (0..(fft_len * 2))
        .map(|index| (index as i16).wrapping_mul(37).wrapping_add(11))
        .collect()
}

#[test]
fn radix2_butterfly_q15_difftest_against_cmsis() {
    for &fft_len in &[16usize, 32, 64, 128, 256, 512, 1024] {
        let mut rust_data = sample_i16_buffer(fft_len);
        let mut cmsis_data = rust_data.clone();

        let twiddle_modifier = get_twiddle_modifier(fft_len);

        // Call Rust implementation
        radix2_butterfly_i16(&mut rust_data, fft_len, twiddle_modifier);

        // Call CMSIS reference implementation
        unsafe {
            arm_radix2_butterfly_q15(
                cmsis_data.as_mut_ptr(),
                fft_len as u32,
                twiddleCoef_4096_q15.as_ptr(),
                twiddle_modifier,
            );
        }

        // Compare results
        assert_eq!(
            rust_data, cmsis_data,
            "Q15 radix2 butterfly mismatch for fft_len={}",
            fft_len
        );

        println!("✓ Q15 radix2 butterfly test passed for fft_len={}", fft_len);
    }
}

#[test]
fn radix2_butterfly_inverse_q15_difftest_against_cmsis() {
    for &fft_len in &[16usize, 32, 64, 128, 256, 512, 1024] {
        let mut rust_data = sample_i16_buffer(fft_len);
        let mut cmsis_data = rust_data.clone();

        let twiddle_modifier = get_twiddle_modifier(fft_len);

        // Call Rust inverse implementation
        radix2_butterfly_inverse_i16(&mut rust_data, fft_len, twiddle_modifier);

        // Call CMSIS reference inverse implementation
        unsafe {
            arm_radix2_butterfly_inverse_q15(
                cmsis_data.as_mut_ptr(),
                fft_len as u32,
                twiddleCoef_4096_q15.as_ptr(),
                twiddle_modifier,
            );
        }

        // Compare results
        assert_eq!(
            rust_data, cmsis_data,
            "Q15 radix2 inverse butterfly mismatch for fft_len={}",
            fft_len
        );

        println!(
            "✓ Q15 radix2 inverse butterfly test passed for fft_len={}",
            fft_len
        );
    }
}
