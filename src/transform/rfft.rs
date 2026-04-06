use crate::common::tables::{REAL_COEF_A_U16, REAL_COEF_A_U32, REAL_COEF_B_U16, REAL_COEF_B_U32};

use super::{CfftI16, CfftI32};

#[inline]
fn i32_mul(a: i32, b: i32) -> i32 {
    (((a as i64 * b as i64) + 0x8000_0000) >> 32) as i32
}

pub struct RfftI16 {
    pub n_fft_real: usize,
    pub ifft_flag: bool,
    pub bit_reverse_flag: bool,
    pub re_table: &'static [u16],
    pub im_table: &'static [u16],
}

impl RfftI16 {
    pub const fn new(n_fft_real: usize, ifft_flag: bool, bit_reverse_flag: bool) -> Self {
        assert!(n_fft_real.is_power_of_two(), "n_fft must be a power of two");
        assert!(
            32 <= n_fft_real && n_fft_real <= 4096,
            "n_fft must be in [32, 4096]",
        );

        let re_table: &'static [u16] = &REAL_COEF_A_U16;
        let im_table: &'static [u16] = &REAL_COEF_B_U16;

        Self {
            n_fft_real,
            ifft_flag,
            bit_reverse_flag,
            re_table,
            im_table,
        }
    }

    fn split_rfft_i16(
        input: &[i16],
        n_fft: usize,
        re_table: &[u16],
        im_table: &[u16],
        output: &mut [i16],
        modifier: u32,
    ) {
        let modifier = modifier as usize;

        for i in 1..n_fft {
            let src1 = 2 * i;
            let src2 = (2 * n_fft) - (2 * i);
            let coef = 2 * i * modifier;

            let a0 = re_table[coef] as i16 as i32;
            let a1 = re_table[coef + 1] as i16 as i32;
            let b0 = im_table[coef] as i16 as i32;
            let b1 = im_table[coef + 1] as i16 as i32;

            let s1r = input[src1] as i32;
            let s1i = input[src1 + 1] as i32;
            let s2r = input[src2] as i32;
            let s2i = input[src2 + 1] as i32;

            let out_r = ((s1r * a0) - (s1i * a1) + (s2r * b0) + (s2i * b1)) >> 16;
            let out_i = ((s2r * b1) - (s2i * b0) + (s1i * a0) + (s1r * a1)) >> 16;

            output[2 * i] = out_r as i16;
            output[2 * i + 1] = out_i as i16;

            let conj = (4 * n_fft) - (2 * i);
            output[conj] = out_r as i16;
            output[conj + 1] = (-(out_i as i16)) as i16;
        }

        output[2 * n_fft] = ((input[0] as i32 - input[1] as i32) >> 1) as i16;
        output[2 * n_fft + 1] = 0;

        output[0] = ((input[0] as i32 + input[1] as i32) >> 1) as i16;
        output[1] = 0;
    }

    fn split_rifft_i16(
        input: &[i16],
        n_fft: usize,
        re_table: &[u16],
        im_table: &[u16],
        output: &mut [i16],
        modifier: u32,
    ) {
        let modifier = modifier as usize;

        for i in 0..n_fft {
            let src1 = 2 * i;
            let src2 = (2 * n_fft) - (2 * i);
            let coef = 2 * i * modifier;

            let a0 = re_table[coef] as i16 as i32;
            let a1 = re_table[coef + 1] as i16 as i32;
            let b0 = im_table[coef] as i16 as i32;
            let b1 = im_table[coef + 1] as i16 as i32;

            let s1r = input[src1] as i32;
            let s1i = input[src1 + 1] as i32;
            let s2r = input[src2] as i32;
            let s2i = input[src2 + 1] as i32;

            let out_r = ((s2r * b0) - (s2i * b1) + (s1r * a0) + (s1i * a1)) >> 16;
            let out_i = ((s1i * a0) - (s1r * a1) - (s2r * b1) - (s2i * b0)) >> 16;

            output[src1] = out_r as i16;
            output[src1 + 1] = out_i as i16;
        }
    }

    pub fn run(&self, input: &[i16], output: &mut [i16]) {
        let n_fft = self.n_fft_real >> 1;
        let modifier = (8192 / self.n_fft_real) as u32;

        if self.ifft_flag {
            assert_eq!(
                input.len(),
                self.n_fft_real * 2,
                "RIFFT input length must be 2 * n_fft_real"
            );
            assert_eq!(
                output.len(),
                self.n_fft_real,
                "RIFFT output length must be n_fft_real"
            );

            Self::split_rifft_i16(input, n_fft, self.re_table, self.im_table, output, modifier);

            CfftI16::new(n_fft, true, self.bit_reverse_flag).run(output);

            for sample in output.iter_mut() {
                let v = (*sample as i32) << 1;
                *sample = v.clamp(i16::MIN as i32, i16::MAX as i32) as i16;
            }
        } else {
            assert_eq!(
                input.len(),
                self.n_fft_real,
                "RFFT input length must be n_fft_real"
            );
            assert_eq!(
                output.len(),
                self.n_fft_real * 2,
                "RFFT output length must be 2 * n_fft_real"
            );

            // Use a fixed-capacity local scratch to run the intermediate complex FFT.
            let mut scratch = [0i16; 4096];
            scratch[..self.n_fft_real].copy_from_slice(input);

            CfftI16::new(n_fft, false, self.bit_reverse_flag).run(&mut scratch[..self.n_fft_real]);

            Self::split_rfft_i16(
                &scratch[..self.n_fft_real],
                n_fft,
                self.re_table,
                self.im_table,
                output,
                modifier,
            );
        }
    }
}

pub struct RfftI32 {
    pub n_fft_real: usize,
    pub ifft_flag: bool,
    pub bit_reverse_flag: bool,
    pub re_table: &'static [u32],
    pub im_table: &'static [u32],
}

impl RfftI32 {
    pub const fn new(n_fft_real: usize, ifft_flag: bool, bit_reverse_flag: bool) -> Self {
        assert!(n_fft_real.is_power_of_two(), "n_fft must be a power of two");
        assert!(
            32 <= n_fft_real && n_fft_real <= 8192,
            "n_fft must be in [32, 8192]",
        );

        let re_table: &'static [u32] = &REAL_COEF_A_U32;
        let im_table: &'static [u32] = &REAL_COEF_B_U32;

        Self {
            n_fft_real,
            ifft_flag,
            bit_reverse_flag,
            re_table,
            im_table,
        }
    }

    fn split_rfft_i32(
        input: &[i32],
        n_fft: usize,
        re_table: &[u32],
        im_table: &[u32],
        output: &mut [i32],
        modifier: u32,
    ) {
        let modifier = modifier as usize;

        for i in 1..n_fft {
            let src1 = 2 * i;
            let src2 = (2 * n_fft) - (2 * i);
            let coef = 2 * i * modifier;

            let a0 = re_table[coef] as i32;
            let a1 = re_table[coef + 1] as i32;
            let b0 = im_table[coef] as i32;
            let b1 = im_table[coef + 1] as i32;

            let s1r = input[src1];
            let s1i = input[src1 + 1];
            let s2r = input[src2];
            let s2i = input[src2 + 1];

            let out_r = i32_mul(s1r, a0)
                .wrapping_sub(i32_mul(s1i, a1))
                .wrapping_add(i32_mul(s2r, b0))
                .wrapping_add(i32_mul(s2i, b1));
            let out_i = i32_mul(s2r, b1)
                .wrapping_sub(i32_mul(s2i, b0))
                .wrapping_add(i32_mul(s1i, a0))
                .wrapping_add(i32_mul(s1r, a1));

            output[2 * i] = out_r;
            output[2 * i + 1] = out_i;

            let conj = (4 * n_fft) - (2 * i);
            output[conj] = out_r;
            output[conj + 1] = out_i.wrapping_neg();
        }

        output[2 * n_fft] = (input[0] - input[1]) >> 1;
        output[2 * n_fft + 1] = 0;

        output[0] = (input[0] + input[1]) >> 1;
        output[1] = 0;
    }

    fn split_rifft_i32(
        input: &[i32],
        n_fft: usize,
        re_table: &[u32],
        im_table: &[u32],
        output: &mut [i32],
        modifier: u32,
    ) {
        let modifier = modifier as usize;

        for i in 0..n_fft {
            let src1 = 2 * i;
            let src2 = (2 * n_fft) - (2 * i);
            let coef = 2 * i * modifier;

            let a0 = re_table[coef] as i32;
            let a1 = re_table[coef + 1] as i32;
            let b0 = im_table[coef] as i32;
            let b1 = im_table[coef + 1] as i32;

            let s1r = input[src1];
            let s1i = input[src1 + 1];
            let s2r = input[src2];
            let s2i = input[src2 + 1];

            let out_r = i32_mul(s2r, b0)
                .wrapping_sub(i32_mul(s2i, b1))
                .wrapping_add(i32_mul(s1r, a0))
                .wrapping_add(i32_mul(s1i, a1));
            let out_i = i32_mul(s1i, a0)
                .wrapping_sub(i32_mul(s1r, a1))
                .wrapping_sub(i32_mul(s2r, b1))
                .wrapping_sub(i32_mul(s2i, b0));

            output[src1] = out_r;
            output[src1 + 1] = out_i;
        }
    }

    pub fn run(&self, input: &[i32], output: &mut [i32]) {
        let n_fft = self.n_fft_real >> 1;
        let modifier = (8192 / self.n_fft_real) as u32;

        if self.ifft_flag {
            assert_eq!(
                input.len(),
                self.n_fft_real * 2,
                "RIFFT input length must be 2 * n_fft_real"
            );
            assert_eq!(
                output.len(),
                self.n_fft_real,
                "RIFFT output length must be n_fft_real"
            );

            Self::split_rifft_i32(input, n_fft, self.re_table, self.im_table, output, modifier);

            CfftI32::new(n_fft, true, self.bit_reverse_flag).run(output);

            for sample in output.iter_mut() {
                *sample = sample.wrapping_shl(1);
            }
        } else {
            assert_eq!(
                input.len(),
                self.n_fft_real,
                "RFFT input length must be n_fft_real"
            );
            assert_eq!(
                output.len(),
                self.n_fft_real * 2,
                "RFFT output length must be 2 * n_fft_real"
            );

            // Use a fixed-capacity local scratch to run the intermediate complex FFT.
            let mut scratch = [0i32; 8192];
            scratch[..self.n_fft_real].copy_from_slice(input);

            CfftI32::new(n_fft, false, self.bit_reverse_flag).run(&mut scratch[..self.n_fft_real]);

            Self::split_rfft_i32(
                &scratch[..self.n_fft_real],
                n_fft,
                self.re_table,
                self.im_table,
                output,
                modifier,
            );
        }
    }
}
