use fixed_dsp::common::tables::{TWIDDLE_TABLE_4096_U16, TWIDDLE_TABLE_4096_U32};
use fixed_dsp::transform::{CfftConfigI16, CfftConfigI32};

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
    fn arm_cfft_radix4by2_q15(pSrc: *mut i16, fftLen: u32, pCoef: *const i16);
    fn arm_cfft_radix4by2_inverse_q15(pSrc: *mut i16, fftLen: u32, pCoef: *const i16);
    fn arm_bitreversal_16(pSrc: *mut u16, bitRevLen: u16, pBitRevTab: *const u16);

    fn arm_radix4_butterfly_q31(
        pSrc: *mut i32,
        fftLen: u32,
        pCoef: *const i32,
        twidCoefModifier: u32,
    );
    fn arm_radix4_butterfly_inverse_q31(
        pSrc: *mut i32,
        fftLen: u32,
        pCoef: *const i32,
        twidCoefModifier: u32,
    );
    fn arm_cfft_radix4by2_q31(pSrc: *mut i32, fftLen: u32, pCoef: *const i32);
    fn arm_cfft_radix4by2_inverse_q31(pSrc: *mut i32, fftLen: u32, pCoef: *const i32);
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
        _ => panic!("unsupported FFT length: {}", fft_len),
    }
}

fn twiddle_q15_for_fft_len(fft_len: usize) -> Vec<i16> {
    let modifier = 4096 / fft_len;
    let mut tw = Vec::with_capacity((fft_len * 3) / 2);
    for i in 0..((3 * fft_len) / 4) {
        let idx = i * modifier;
        tw.push(TWIDDLE_TABLE_4096_U16[2 * idx] as i16);
        tw.push(TWIDDLE_TABLE_4096_U16[2 * idx + 1] as i16);
    }
    tw
}

fn twiddle_q31_for_fft_len(fft_len: usize) -> Vec<i32> {
    let modifier = 4096 / fft_len;
    let mut tw = Vec::with_capacity((fft_len * 3) / 2);
    for i in 0..((3 * fft_len) / 4) {
        let idx = i * modifier;
        tw.push(TWIDDLE_TABLE_4096_U32[2 * idx] as i32);
        tw.push(TWIDDLE_TABLE_4096_U32[2 * idx + 1] as i32);
    }
    tw
}

fn cmsis_cfft_q15(data: &mut [i16], fft_len: usize, ifft_flag: bool, bit_reverse_flag: bool) {
    let twid_modifier = (4096 / fft_len) as u32;
    let twiddle = twiddle_q15_for_fft_len(fft_len);
    unsafe {
        if ifft_flag {
            match fft_len {
                16 | 64 | 256 | 1024 | 4096 => arm_radix4_butterfly_inverse_q15(
                    data.as_mut_ptr(),
                    fft_len as u32,
                    TWIDDLE_TABLE_4096_U16.as_ptr().cast::<i16>(),
                    twid_modifier,
                ),
                32 | 128 | 512 | 2048 => arm_cfft_radix4by2_inverse_q15(
                    data.as_mut_ptr(),
                    fft_len as u32,
                    twiddle.as_ptr(),
                ),
                _ => unreachable!(),
            }
        } else {
            match fft_len {
                16 | 64 | 256 | 1024 | 4096 => arm_radix4_butterfly_q15(
                    data.as_mut_ptr(),
                    fft_len as u32,
                    TWIDDLE_TABLE_4096_U16.as_ptr().cast::<i16>(),
                    twid_modifier,
                ),
                32 | 128 | 512 | 2048 => {
                    arm_cfft_radix4by2_q15(data.as_mut_ptr(), fft_len as u32, twiddle.as_ptr())
                }
                _ => unreachable!(),
            }
        }
    }

    if bit_reverse_flag {
        let (ptr, len) = table_ptr_and_len(fft_len);
        unsafe {
            arm_bitreversal_16(data.as_mut_ptr().cast::<u16>(), len, ptr);
        }
    }
}

fn cmsis_cfft_q31(data: &mut [i32], fft_len: usize, ifft_flag: bool, bit_reverse_flag: bool) {
    let twid_modifier = (4096 / fft_len) as u32;
    let twiddle = twiddle_q31_for_fft_len(fft_len);
    unsafe {
        if ifft_flag {
            match fft_len {
                16 | 64 | 256 | 1024 | 4096 => arm_radix4_butterfly_inverse_q31(
                    data.as_mut_ptr(),
                    fft_len as u32,
                    TWIDDLE_TABLE_4096_U32.as_ptr().cast::<i32>(),
                    twid_modifier,
                ),
                32 | 128 | 512 | 2048 => arm_cfft_radix4by2_inverse_q31(
                    data.as_mut_ptr(),
                    fft_len as u32,
                    twiddle.as_ptr(),
                ),
                _ => unreachable!(),
            }
        } else {
            match fft_len {
                16 | 64 | 256 | 1024 | 4096 => arm_radix4_butterfly_q31(
                    data.as_mut_ptr(),
                    fft_len as u32,
                    TWIDDLE_TABLE_4096_U32.as_ptr().cast::<i32>(),
                    twid_modifier,
                ),
                32 | 128 | 512 | 2048 => {
                    arm_cfft_radix4by2_q31(data.as_mut_ptr(), fft_len as u32, twiddle.as_ptr())
                }
                _ => unreachable!(),
            }
        }
    }

    if bit_reverse_flag {
        let (ptr, len) = table_ptr_and_len(fft_len);
        unsafe {
            arm_bitreversal_32(data.as_mut_ptr().cast::<u32>(), len, ptr);
        }
    }
}

#[test]
fn cfft_q15_difftest_against_cmsis() {
    for &fft_len in &[16usize, 32, 64, 128, 256, 512, 1024, 2048, 4096] {
        for &ifft_flag in &[false, true] {
            for &bit_reverse_flag in &[false, true] {
                let mut rust_data = sample_i16_buffer(fft_len);
                let mut cmsis_data = rust_data.clone();

                let cfg = CfftConfigI16::new(fft_len, ifft_flag, bit_reverse_flag);
                cfg.cfft_i16(&mut rust_data);
                cmsis_cfft_q15(&mut cmsis_data, fft_len, ifft_flag, bit_reverse_flag);

                assert_eq!(
                    rust_data, cmsis_data,
                    "Q15 CFFT mismatch fft_len={}, ifft={}, bitrev={}",
                    fft_len, ifft_flag, bit_reverse_flag
                );
            }
        }
    }
}

#[test]
fn cfft_q31_difftest_against_cmsis() {
    for &fft_len in &[16usize, 32, 64, 128, 256, 512, 1024, 2048, 4096] {
        for &ifft_flag in &[false, true] {
            for &bit_reverse_flag in &[false, true] {
                let mut rust_data = sample_i32_buffer(fft_len);
                let mut cmsis_data = rust_data.clone();

                let cfg = CfftConfigI32::new(fft_len, ifft_flag, bit_reverse_flag);
                cfg.cfft_i32(&mut rust_data);
                cmsis_cfft_q31(&mut cmsis_data, fft_len, ifft_flag, bit_reverse_flag);

                assert_eq!(
                    rust_data, cmsis_data,
                    "Q31 CFFT mismatch fft_len={}, ifft={}, bitrev={}",
                    fft_len, ifft_flag, bit_reverse_flag
                );
            }
        }
    }
}
