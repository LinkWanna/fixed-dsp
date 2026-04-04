use fixed_dsp::transform::{bitreversal_i16, bitreversal_i32};

unsafe extern "C" {
    static armBitRevTable: [u16; 1024];

    fn arm_bitreversal_q15(
        pSrc16: *mut i16,
        fftLen: u32,
        bitRevFactor: u16,
        pBitRevTab: *const u16,
    );

    fn arm_bitreversal_q31(pSrc: *mut i32, fftLen: u32, bitRevFactor: u16, pBitRevTab: *const u16);
}

fn sample_i16_buffer(fft_len: usize) -> Vec<i16> {
    (0..(fft_len * 2))
        .map(|index| (index as i16).wrapping_mul(37).wrapping_add(11))
        .collect()
}

fn sample_i32_buffer(fft_len: usize) -> Vec<i32> {
    (0..(fft_len * 2))
        .map(|index| (index as i32).wrapping_mul(1_048_583).wrapping_add(29))
        .collect()
}

fn table_offset_and_factor(fft_len: usize) -> (usize, u16) {
    match fft_len {
        16 => (255, 256),
        32 => (127, 128),
        64 => (63, 64),
        128 => (31, 32),
        256 => (15, 16),
        512 => (7, 8),
        1024 => (3, 4),
        2048 => (1, 2),
        4096 => (0, 1),
        _ => panic!("unsupported FFT length for CMSIS bit reversal: {}", fft_len),
    }
}

#[test]
fn bitreversal_q15_difftest_against_cmsis() {
    for &fft_len in &[16usize, 32, 64, 128, 256, 512, 1024, 2048, 4096] {
        let mut rust_data = sample_i16_buffer(fft_len);
        let mut cmsis_data = rust_data.clone();
        let (offset, factor) = table_offset_and_factor(fft_len);

        bitreversal_i16(&mut rust_data, fft_len);

        unsafe {
            arm_bitreversal_q15(
                cmsis_data.as_mut_ptr(),
                fft_len as u32,
                factor,
                armBitRevTable.as_ptr().add(offset),
            );
        }

        assert_eq!(
            rust_data, cmsis_data,
            "Q15 bit reversal mismatch for fft_len={}",
            fft_len
        );
    }
}

#[test]
fn bitreversal_q31_difftest_against_cmsis() {
    for &fft_len in &[16usize, 32, 64, 128, 256, 512, 1024, 2048, 4096] {
        let mut rust_data = sample_i32_buffer(fft_len);
        let mut cmsis_data = rust_data.clone();
        let (offset, factor) = table_offset_and_factor(fft_len);

        bitreversal_i32(&mut rust_data, fft_len);

        unsafe {
            arm_bitreversal_q31(
                cmsis_data.as_mut_ptr(),
                fft_len as u32,
                factor,
                armBitRevTable.as_ptr().add(offset),
            );
        }

        assert_eq!(
            rust_data, cmsis_data,
            "Q31 bit reversal mismatch for fft_len={}",
            fft_len
        );
    }
}
