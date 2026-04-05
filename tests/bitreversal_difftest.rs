use fixed_dsp::transform::{bitreversal_i16, bitreversal_i32};

unsafe extern "C" {
    static armBitRevIndexTable_fixed_16: u16;
    static armBitRevIndexTable_fixed_32: u16;
    static armBitRevIndexTable_fixed_64: u16;
    static armBitRevIndexTable_fixed_128: u16;
    static armBitRevIndexTable_fixed_256: u16;
    static armBitRevIndexTable_fixed_512: u16;
    static armBitRevIndexTable_fixed_1024: u16;
    static armBitRevIndexTable_fixed_2048: u16;
    static armBitRevIndexTable_fixed_4096: u16;

    fn arm_bitreversal_16(pSrc: *mut u16, bitRevLen: u16, pBitRevTab: *const u16);

    fn arm_bitreversal_32(pSrc: *mut u32, bitRevLen: u16, pBitRevTab: *const u16);
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

fn table_ptr_and_len(fft_len: usize) -> (*const u16, u16) {
    match fft_len {
        16 => (&raw const armBitRevIndexTable_fixed_16, 12),
        32 => (&raw const armBitRevIndexTable_fixed_32, 24),
        64 => (&raw const armBitRevIndexTable_fixed_64, 56),
        128 => (&raw const armBitRevIndexTable_fixed_128, 112),
        256 => (&raw const armBitRevIndexTable_fixed_256, 240),
        512 => (&raw const armBitRevIndexTable_fixed_512, 480),
        1024 => (&raw const armBitRevIndexTable_fixed_1024, 992),
        2048 => (&raw const armBitRevIndexTable_fixed_2048, 1984),
        4096 => (&raw const armBitRevIndexTable_fixed_4096, 4032),
        _ => panic!(
            "unsupported FFT length for CMSIS new-entry bit reversal: {}",
            fft_len
        ),
    }
}

#[test]
fn bitreversal_q15_difftest_against_cmsis() {
    for &fft_len in &[16usize, 32, 64, 128, 256, 512, 1024, 2048, 4096] {
        let mut rust_data = sample_i16_buffer(fft_len);
        let mut cmsis_data = rust_data.clone();
        let (table_ptr, bit_rev_len) = table_ptr_and_len(fft_len);

        bitreversal_i16(&mut rust_data, fft_len);

        unsafe {
            arm_bitreversal_16(
                cmsis_data.as_mut_ptr().cast::<u16>(),
                bit_rev_len,
                table_ptr,
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
        let (table_ptr, bit_rev_len) = table_ptr_and_len(fft_len);

        bitreversal_i32(&mut rust_data, fft_len);

        unsafe {
            arm_bitreversal_32(
                cmsis_data.as_mut_ptr().cast::<u32>(),
                bit_rev_len,
                table_ptr,
            );
        }

        assert_eq!(
            rust_data, cmsis_data,
            "Q31 bit reversal mismatch for fft_len={}",
            fft_len
        );
    }
}
