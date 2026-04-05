use crate::common::tables::{
    BIT_REV_TABLE_16, BIT_REV_TABLE_32, BIT_REV_TABLE_64, BIT_REV_TABLE_128, BIT_REV_TABLE_256,
    BIT_REV_TABLE_512, BIT_REV_TABLE_1024, BIT_REV_TABLE_2048, BIT_REV_TABLE_4096,
};
use crate::common::tables::{
    TWIDDLE_TABLE_16_U16, TWIDDLE_TABLE_16_U32, TWIDDLE_TABLE_32_U16, TWIDDLE_TABLE_32_U32,
    TWIDDLE_TABLE_64_U16, TWIDDLE_TABLE_64_U32, TWIDDLE_TABLE_128_U16, TWIDDLE_TABLE_128_U32,
    TWIDDLE_TABLE_256_U16, TWIDDLE_TABLE_256_U32, TWIDDLE_TABLE_512_U16, TWIDDLE_TABLE_512_U32,
    TWIDDLE_TABLE_1024_U16, TWIDDLE_TABLE_1024_U32, TWIDDLE_TABLE_2048_U16, TWIDDLE_TABLE_2048_U32,
    TWIDDLE_TABLE_4096_U16, TWIDDLE_TABLE_4096_U32,
};

use super::{
    bitreversal_i16, bitreversal_i32, radix4_butterfly_i16, radix4_butterfly_i32,
    radix4_butterfly_inverse_i16, radix4_butterfly_inverse_i32,
};

#[inline]
fn i16_mul(a: i16, b: i16) -> i16 {
    (((a as i32) * (b as i32)) >> 16) as i16
}

#[inline]
fn i32_mul(a: i32, b: i32) -> i32 {
    (((a as i64 * b as i64) + 0x8000_0000) >> 32) as i32
}

pub struct CfftI16 {
    pub fft_len: usize,
    pub ifft_flag: bool,
    pub bit_reverse_flag: bool,
    pub twiddle: &'static [u16],
    pub bit_rev_table: &'static [u16],
}

impl CfftI16 {
    pub const fn new(fft_len: usize, ifft_flag: bool, bit_reverse_flag: bool) -> Self {
        assert!(fft_len.is_power_of_two(), "fft_len must be a power of two");
        assert!(
            16 <= fft_len && fft_len <= 4096,
            "fft_len must be in [16, 4096]",
        );

        let twiddle: &'static [u16] = match fft_len {
            16 => &TWIDDLE_TABLE_16_U16,
            32 => &TWIDDLE_TABLE_32_U16,
            64 => &TWIDDLE_TABLE_64_U16,
            128 => &TWIDDLE_TABLE_128_U16,
            256 => &TWIDDLE_TABLE_256_U16,
            512 => &TWIDDLE_TABLE_512_U16,
            1024 => &TWIDDLE_TABLE_1024_U16,
            2048 => &TWIDDLE_TABLE_2048_U16,
            4096 => &TWIDDLE_TABLE_4096_U16,
            _ => unreachable!(),
        };

        let bit_rev_table: &'static [u16] = match fft_len {
            16 => &BIT_REV_TABLE_16,
            32 => &BIT_REV_TABLE_32,
            64 => &BIT_REV_TABLE_64,
            128 => &BIT_REV_TABLE_128,
            256 => &BIT_REV_TABLE_256,
            512 => &BIT_REV_TABLE_512,
            1024 => &BIT_REV_TABLE_1024,
            2048 => &BIT_REV_TABLE_2048,
            4096 => &BIT_REV_TABLE_4096,
            _ => unreachable!(),
        };

        Self {
            fft_len,
            ifft_flag,
            bit_reverse_flag,
            twiddle,
            bit_rev_table,
        }
    }

    fn radix4by2_butterfly_i16(data: &mut [i16], fft_len: usize, twiddle: &[u16], modifier: u16) {
        let n2 = fft_len >> 1;

        for i in 0..n2 {
            let cos_val = twiddle[2 * i] as i16;
            let sin_val = twiddle[2 * i + 1] as i16;

            let l = i + n2;

            let t_re_half = (data[2 * i] as i32) >> 1;
            let t_im_half = (data[2 * i + 1] as i32) >> 1;
            let s_re_half = (data[2 * l] as i32) >> 1;
            let s_im_half = (data[2 * l + 1] as i32) >> 1;

            let xt = (t_re_half - s_re_half) as i16;
            let yt = (t_im_half - s_im_half) as i16;

            data[2 * i] = ((t_re_half + s_re_half) >> 1) as i16;
            data[2 * i + 1] = ((t_im_half + s_im_half) >> 1) as i16;

            let out_re = i16_mul(xt, cos_val).wrapping_add(i16_mul(yt, sin_val));
            let out_im = i16_mul(yt, cos_val).wrapping_sub(i16_mul(xt, sin_val));

            data[2 * l] = out_re;
            data[2 * l + 1] = out_im;
        }

        radix4_butterfly_i16(&mut data[..fft_len], n2, twiddle, modifier);
        radix4_butterfly_i16(&mut data[fft_len..], n2, twiddle, modifier);

        for i in 0..n2 {
            data[4 * i] = data[4 * i].wrapping_shl(1);
            data[4 * i + 1] = data[4 * i + 1].wrapping_shl(1);
            data[4 * i + 2] = data[4 * i + 2].wrapping_shl(1);
            data[4 * i + 3] = data[4 * i + 3].wrapping_shl(1);
        }
    }

    fn radix4by2_butterfly_inverse_i16(
        data: &mut [i16],
        fft_len: usize,
        twiddle: &[u16],
        modifier: u16,
    ) {
        let n2 = fft_len >> 1;

        for i in 0..n2 {
            let cos_val = twiddle[2 * i] as i16;
            let sin_val = twiddle[2 * i + 1] as i16;

            let l = i + n2;

            let t_re_half = (data[2 * i] as i32) >> 1;
            let t_im_half = (data[2 * i + 1] as i32) >> 1;
            let s_re_half = (data[2 * l] as i32) >> 1;
            let s_im_half = (data[2 * l + 1] as i32) >> 1;

            let xt = (t_re_half - s_re_half) as i16;
            let yt = (t_im_half - s_im_half) as i16;

            data[2 * i] = ((t_re_half + s_re_half) >> 1) as i16;
            data[2 * i + 1] = ((t_im_half + s_im_half) >> 1) as i16;

            let out_re = i16_mul(xt, cos_val).wrapping_sub(i16_mul(yt, sin_val));
            let out_im = i16_mul(yt, cos_val).wrapping_add(i16_mul(xt, sin_val));

            data[2 * l] = out_re;
            data[2 * l + 1] = out_im;
        }

        radix4_butterfly_inverse_i16(&mut data[..fft_len], n2, twiddle, modifier);
        radix4_butterfly_inverse_i16(&mut data[fft_len..], n2, twiddle, modifier);

        for i in 0..n2 {
            data[4 * i] = data[4 * i].wrapping_shl(1);
            data[4 * i + 1] = data[4 * i + 1].wrapping_shl(1);
            data[4 * i + 2] = data[4 * i + 2].wrapping_shl(1);
            data[4 * i + 3] = data[4 * i + 3].wrapping_shl(1);
        }
    }

    /// Q15 CFFT/CIFFT entry compatible with CMSIS `arm_cfft_q15` dispatch behavior.
    pub fn run(self, data: &mut [i16]) {
        assert_eq!(
            data.len(),
            self.fft_len * 2,
            "Q15 buffer length must be 2 * fft_len"
        );

        let radix4_modifier = 1u16;
        let radix4by2_modifier = 2u16;

        if self.ifft_flag {
            match self.fft_len {
                16 | 64 | 256 | 1024 | 4096 => {
                    radix4_butterfly_inverse_i16(data, self.fft_len, self.twiddle, radix4_modifier)
                }
                32 | 128 | 512 | 2048 => Self::radix4by2_butterfly_inverse_i16(
                    data,
                    self.fft_len,
                    self.twiddle,
                    radix4by2_modifier,
                ),
                _ => unreachable!(),
            }
        } else {
            match self.fft_len {
                16 | 64 | 256 | 1024 | 4096 => {
                    radix4_butterfly_i16(data, self.fft_len, self.twiddle, radix4_modifier)
                }
                32 | 128 | 512 | 2048 => Self::radix4by2_butterfly_i16(
                    data,
                    self.fft_len,
                    self.twiddle,
                    radix4by2_modifier,
                ),
                _ => unreachable!(),
            }
        }

        if self.bit_reverse_flag {
            bitreversal_i16(data, self.bit_rev_table);
        }
    }
}

pub struct CfftI32 {
    pub fft_len: usize,
    pub ifft_flag: bool,
    pub bit_reverse_flag: bool,
    pub twiddle: &'static [u32],
    pub bit_rev_table: &'static [u16],
}

impl CfftI32 {
    pub const fn new(fft_len: usize, ifft_flag: bool, bit_reverse_flag: bool) -> Self {
        assert!(fft_len.is_power_of_two(), "fft_len must be a power of two");
        assert!(
            16 <= fft_len && fft_len <= 4096,
            "fft_len must be in [16, 4096]",
        );

        let twiddle: &'static [u32] = match fft_len {
            16 => &TWIDDLE_TABLE_16_U32,
            32 => &TWIDDLE_TABLE_32_U32,
            64 => &TWIDDLE_TABLE_64_U32,
            128 => &TWIDDLE_TABLE_128_U32,
            256 => &TWIDDLE_TABLE_256_U32,
            512 => &TWIDDLE_TABLE_512_U32,
            1024 => &TWIDDLE_TABLE_1024_U32,
            2048 => &TWIDDLE_TABLE_2048_U32,
            4096 => &TWIDDLE_TABLE_4096_U32,
            _ => unreachable!(),
        };

        let bit_rev_table: &'static [u16] = match fft_len {
            16 => &BIT_REV_TABLE_16,
            32 => &BIT_REV_TABLE_32,
            64 => &BIT_REV_TABLE_64,
            128 => &BIT_REV_TABLE_128,
            256 => &BIT_REV_TABLE_256,
            512 => &BIT_REV_TABLE_512,
            1024 => &BIT_REV_TABLE_1024,
            2048 => &BIT_REV_TABLE_2048,
            4096 => &BIT_REV_TABLE_4096,
            _ => unreachable!(),
        };

        Self {
            fft_len,
            ifft_flag,
            bit_reverse_flag,
            twiddle,
            bit_rev_table,
        }
    }

    fn radix4by2_butterfly_i32(
        data: &mut [i32],
        fft_len: usize,
        twiddle: &[u32],
        radix4_modifier: u16,
    ) {
        let n2 = fft_len >> 1;

        for i in 0..n2 {
            let cos_val = twiddle[2 * i] as i32;
            let sin_val = twiddle[2 * i + 1] as i32;

            let l = i + n2;

            let xt = (data[2 * i] >> 2) - (data[2 * l] >> 2);
            data[2 * i] = (data[2 * i] >> 2) + (data[2 * l] >> 2);

            let yt = (data[2 * i + 1] >> 2) - (data[2 * l + 1] >> 2);
            data[2 * i + 1] = (data[2 * l + 1] >> 2) + (data[2 * i + 1] >> 2);

            let mut out_re = i32_mul(xt, cos_val);
            let mut out_im = i32_mul(yt, cos_val);
            out_re = out_re.wrapping_add(i32_mul(yt, sin_val));
            out_im = out_im.wrapping_sub(i32_mul(xt, sin_val));

            data[2 * l] = out_re.wrapping_shl(1);
            data[2 * l + 1] = out_im.wrapping_shl(1);
        }

        radix4_butterfly_i32(&mut data[..fft_len], n2, twiddle, radix4_modifier);
        radix4_butterfly_i32(&mut data[fft_len..], n2, twiddle, radix4_modifier);

        for i in 0..n2 {
            data[4 * i] = data[4 * i].wrapping_shl(1);
            data[4 * i + 1] = data[4 * i + 1].wrapping_shl(1);
            data[4 * i + 2] = data[4 * i + 2].wrapping_shl(1);
            data[4 * i + 3] = data[4 * i + 3].wrapping_shl(1);
        }
    }

    fn radix4by2_butterfly_inverse_i32(
        data: &mut [i32],
        fft_len: usize,
        twiddle: &[u32],
        radix4_modifier: u16,
    ) {
        let n2 = fft_len >> 1;

        for i in 0..n2 {
            let cos_val = twiddle[2 * i] as i32;
            let sin_val = twiddle[2 * i + 1] as i32;

            let l = i + n2;

            let xt = (data[2 * i] >> 2) - (data[2 * l] >> 2);
            data[2 * i] = (data[2 * i] >> 2) + (data[2 * l] >> 2);

            let yt = (data[2 * i + 1] >> 2) - (data[2 * l + 1] >> 2);
            data[2 * i + 1] = (data[2 * l + 1] >> 2) + (data[2 * i + 1] >> 2);

            let mut out_re = i32_mul(xt, cos_val);
            let mut out_im = i32_mul(yt, cos_val);
            out_re = out_re.wrapping_sub(i32_mul(yt, sin_val));
            out_im = out_im.wrapping_add(i32_mul(xt, sin_val));

            data[2 * l] = out_re.wrapping_shl(1);
            data[2 * l + 1] = out_im.wrapping_shl(1);
        }

        radix4_butterfly_inverse_i32(&mut data[..fft_len], n2, twiddle, radix4_modifier);
        radix4_butterfly_inverse_i32(&mut data[fft_len..], n2, twiddle, radix4_modifier);

        for i in 0..n2 {
            data[4 * i] = data[4 * i].wrapping_shl(1);
            data[4 * i + 1] = data[4 * i + 1].wrapping_shl(1);
            data[4 * i + 2] = data[4 * i + 2].wrapping_shl(1);
            data[4 * i + 3] = data[4 * i + 3].wrapping_shl(1);
        }
    }

    /// Q31 CFFT/CIFFT entry with the same dispatch pattern as CMSIS `arm_cfft_q31`.
    pub fn run(&self, data: &mut [i32]) {
        assert_eq!(
            data.len(),
            self.fft_len * 2,
            "Q31 buffer length must be 2 * fft_len"
        );

        let radix4_modifier = 1u16;
        let radix4by2_modifier = 2u16;

        if self.ifft_flag {
            match self.fft_len {
                16 | 64 | 256 | 1024 | 4096 => {
                    radix4_butterfly_inverse_i32(data, self.fft_len, self.twiddle, radix4_modifier)
                }
                32 | 128 | 512 | 2048 => Self::radix4by2_butterfly_inverse_i32(
                    data,
                    self.fft_len,
                    self.twiddle,
                    radix4by2_modifier,
                ),
                _ => unreachable!(),
            }
        } else {
            match self.fft_len {
                16 | 64 | 256 | 1024 | 4096 => {
                    radix4_butterfly_i32(data, self.fft_len, self.twiddle, radix4_modifier)
                }
                32 | 128 | 512 | 2048 => Self::radix4by2_butterfly_i32(
                    data,
                    self.fft_len,
                    self.twiddle,
                    radix4by2_modifier,
                ),
                _ => unreachable!(),
            }
        }

        if self.bit_reverse_flag {
            bitreversal_i32(data, self.bit_rev_table);
        }
    }
}
