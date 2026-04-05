use fixed_dsp::common::tables::TWIDDLE_TABLE_4096_U16;
use fixed_dsp::transform::RfftConfigI16;

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

    static realCoefAQ15: i16;
    static realCoefBQ15: i16;

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

    fn arm_split_rfft_q15(
        pSrc: *mut i16,
        fftLen: u32,
        pATable: *const i16,
        pBTable: *const i16,
        pDst: *mut i16,
        modifier: u32,
    );

    fn arm_split_rifft_q15(
        pSrc: *mut i16,
        fftLen: u32,
        pATable: *const i16,
        pBTable: *const i16,
        pDst: *mut i16,
        modifier: u32,
    );

    fn arm_shift_q15(pSrc: *const i16, shiftBits: i8, pDst: *mut i16, blockSize: u32);
}

fn sample_real_i16_buffer(fft_len: usize) -> Vec<i16> {
    (0..fft_len)
        .map(|index| (index as i16).wrapping_mul(73).wrapping_add(17))
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

fn cmsis_rfft_q15_forward(input: &[i16], output: &mut [i16], bit_reverse_flag: bool) {
    let fft_len_real = input.len();
    let fft_len = fft_len_real >> 1;
    let modifier = (8192 / fft_len_real) as u32;

    let mut scratch = input.to_vec();
    cmsis_cfft_q15(&mut scratch, fft_len, false, bit_reverse_flag);

    unsafe {
        arm_split_rfft_q15(
            scratch.as_mut_ptr(),
            fft_len as u32,
            &raw const realCoefAQ15,
            &raw const realCoefBQ15,
            output.as_mut_ptr(),
            modifier,
        );
    }
}

fn cmsis_rfft_q15_inverse(input: &[i16], output: &mut [i16], bit_reverse_flag: bool) {
    let fft_len_real = output.len();
    let fft_len = fft_len_real >> 1;
    let modifier = (8192 / fft_len_real) as u32;

    let mut spectrum = input.to_vec();
    unsafe {
        arm_split_rifft_q15(
            spectrum.as_mut_ptr(),
            fft_len as u32,
            &raw const realCoefAQ15,
            &raw const realCoefBQ15,
            output.as_mut_ptr(),
            modifier,
        );
    }

    cmsis_cfft_q15(output, fft_len, true, bit_reverse_flag);

    unsafe {
        arm_shift_q15(output.as_ptr(), 1, output.as_mut_ptr(), fft_len_real as u32);
    }
}

#[test]
fn rfft_q15_forward_difftest_against_cmsis() {
    for &fft_len_real in &[32usize, 64, 128, 256, 512, 1024, 2048, 4096] {
        for &bit_reverse_flag in &[false, true] {
            let input = sample_real_i16_buffer(fft_len_real);

            let mut rust_output = vec![0i16; fft_len_real * 2];
            let mut cmsis_output = vec![0i16; fft_len_real * 2];

            let cfg = RfftConfigI16::new(fft_len_real, false, bit_reverse_flag);
            cfg.rfft_i16(&input, &mut rust_output);
            cmsis_rfft_q15_forward(&input, &mut cmsis_output, bit_reverse_flag);

            assert_eq!(
                rust_output, cmsis_output,
                "Q15 RFFT mismatch fft_len_real={}, bitrev={}",
                fft_len_real, bit_reverse_flag
            );
        }
    }
}

#[test]
fn rfft_q15_inverse_difftest_against_cmsis() {
    for &fft_len_real in &[32usize, 64, 128, 256, 512, 1024, 2048, 4096] {
        for &bit_reverse_flag in &[false, true] {
            let input_time = sample_real_i16_buffer(fft_len_real);
            let mut spectrum = vec![0i16; fft_len_real * 2];
            cmsis_rfft_q15_forward(&input_time, &mut spectrum, bit_reverse_flag);

            let mut rust_output = vec![0i16; fft_len_real];
            let mut cmsis_output = vec![0i16; fft_len_real];

            let cfg = RfftConfigI16::new(fft_len_real, true, bit_reverse_flag);
            cfg.rfft_i16(&spectrum, &mut rust_output);
            cmsis_rfft_q15_inverse(&spectrum, &mut cmsis_output, bit_reverse_flag);

            assert_eq!(
                rust_output, cmsis_output,
                "Q15 RIFFT mismatch fft_len_real={}, bitrev={}",
                fft_len_real, bit_reverse_flag
            );
        }
    }
}
