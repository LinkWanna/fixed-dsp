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

const RADIX4_MODIFIER: u16 = 1;
const RADIX4BY2_MODIFIER: u16 = 2;

#[inline]
fn i16_mul(a: i16, b: i16) -> i16 {
    (((a as i32) * (b as i32)) >> 16) as i16
}

#[inline]
fn i32_mul(a: i32, b: i32) -> i32 {
    (((a as i64 * b as i64) + 0x8000_0000) >> 32) as i32
}

#[inline]
fn upshift_quads_i16(data: &mut [i16]) {
    for quad in data.chunks_exact_mut(4) {
        quad[0] = quad[0].wrapping_shl(1);
        quad[1] = quad[1].wrapping_shl(1);
        quad[2] = quad[2].wrapping_shl(1);
        quad[3] = quad[3].wrapping_shl(1);
    }
}

#[inline]
fn upshift_quads_i32(data: &mut [i32]) {
    for quad in data.chunks_exact_mut(4) {
        quad[0] = quad[0].wrapping_shl(1);
        quad[1] = quad[1].wrapping_shl(1);
        quad[2] = quad[2].wrapping_shl(1);
        quad[3] = quad[3].wrapping_shl(1);
    }
}

/// Input and Output formats for CFFT Q15
/// | CFFT Size  | Input Format  | Output Format  | Number of bits to upscale |
/// | ---------: | ------------: | -------------: | ------------------------: |
/// | 16         | 1.15          | 5.11           | 4                         |
/// | 64         | 1.15          | 7.9            | 6                         |
/// | 256        | 1.15          | 9.7            | 8                         |
/// | 1024       | 1.15          | 11.5           | 10                        |
///
/// Input and Output formats for CIFFT Q15
/// | CIFFT Size  | Input Format  | Output Format  | Number of bits to upscale |
/// | ----------: | ------------: | -------------: | ------------------------: |
/// | 16          | 1.15          | 5.11           | 0                         |
/// | 64          | 1.15          | 7.9            | 0                         |
/// | 256         | 1.15          | 9.7            | 0                         |
/// | 1024        | 1.15          | 11.5           | 0                         |
pub struct CfftI16 {
    pub n_fft: usize,
    pub ifft_flag: bool,
    pub bit_reverse_flag: bool,
    pub twiddle: &'static [u16],
    pub bit_rev_table: &'static [u16],
}

impl CfftI16 {
    pub const fn new(n_fft: usize, ifft_flag: bool, bit_reverse_flag: bool) -> Self {
        assert!(n_fft.is_power_of_two(), "n_fft must be a power of two");
        assert!(16 <= n_fft && n_fft <= 4096, "n_fft must be in [16, 4096]",);

        let twiddle: &'static [u16] = match n_fft {
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

        let bit_rev_table: &'static [u16] = match n_fft {
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
            n_fft,
            ifft_flag,
            bit_reverse_flag,
            twiddle,
            bit_rev_table,
        }
    }

    fn radix4by2_butterfly_i16(data: &mut [i16], n_fft: usize, twiddle: &[u16], modifier: u16) {
        let n2 = n_fft >> 1;

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

        radix4_butterfly_i16(&mut data[..n_fft], n2, twiddle, modifier);
        radix4_butterfly_i16(&mut data[n_fft..], n2, twiddle, modifier);

        upshift_quads_i16(data);
    }

    fn radix4by2_butterfly_inverse_i16(
        data: &mut [i16],
        n_fft: usize,
        twiddle: &[u16],
        modifier: u16,
    ) {
        let n2 = n_fft >> 1;

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

        radix4_butterfly_inverse_i16(&mut data[..n_fft], n2, twiddle, modifier);
        radix4_butterfly_inverse_i16(&mut data[n_fft..], n2, twiddle, modifier);

        upshift_quads_i16(data);
    }

    /// Q15 CFFT/CIFFT entry compatible with CMSIS `arm_cfft_q15` dispatch behavior.
    pub fn run(&self, data: &mut [i16]) {
        assert_eq!(
            data.len(),
            self.n_fft * 2,
            "Q15 buffer length must be 2 * n_fft"
        );

        let use_radix4by2 = matches!(self.n_fft, 32 | 128 | 512 | 2048);

        if self.ifft_flag {
            if use_radix4by2 {
                Self::radix4by2_butterfly_inverse_i16(
                    data,
                    self.n_fft,
                    self.twiddle,
                    RADIX4BY2_MODIFIER,
                );
            } else {
                radix4_butterfly_inverse_i16(data, self.n_fft, self.twiddle, RADIX4_MODIFIER);
            }
        } else {
            if use_radix4by2 {
                Self::radix4by2_butterfly_i16(data, self.n_fft, self.twiddle, RADIX4BY2_MODIFIER);
            } else {
                radix4_butterfly_i16(data, self.n_fft, self.twiddle, RADIX4_MODIFIER);
            }
        }

        if self.bit_reverse_flag {
            bitreversal_i16(data, self.bit_rev_table);
        }
    }
}

/// Input and Output formats for CFFT Q31
/// | CFFT Size  | Input Format  | Output Format  | Number of bits to upscale |
/// | ---------: | ------------: | -------------: | ------------------------: |
/// | 16         | 1.31          | 5.27           | 4                         |
/// | 64         | 1.31          | 7.25           | 6                         |
/// | 256        | 1.31          | 9.23           | 8                         |
/// | 1024       | 1.31          | 11.21          | 10                        |
///
/// Input and Output formats for CIFFT Q31
/// | CIFFT Size  | Input Format  | Output Format  | Number of bits to upscale |
/// | ----------: | ------------: | -------------: | ------------------------: |
/// | 16          | 1.31          | 5.27           | 0                         |
/// | 64          | 1.31          | 7.25           | 0                         |
/// | 256         | 1.31          | 9.23           | 0                         |
/// | 1024        | 1.31          | 11.21          | 0                         |
pub struct CfftI32 {
    pub n_fft: usize,
    pub ifft_flag: bool,
    pub bit_reverse_flag: bool,
    pub twiddle: &'static [u32],
    pub bit_rev_table: &'static [u16],
}

impl CfftI32 {
    pub const fn new(n_fft: usize, ifft_flag: bool, bit_reverse_flag: bool) -> Self {
        assert!(n_fft.is_power_of_two(), "n_fft must be a power of two");
        assert!(16 <= n_fft && n_fft <= 4096, "n_fft must be in [16, 4096]",);

        let twiddle: &'static [u32] = match n_fft {
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

        let bit_rev_table: &'static [u16] = match n_fft {
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
            n_fft,
            ifft_flag,
            bit_reverse_flag,
            twiddle,
            bit_rev_table,
        }
    }

    fn radix4by2_butterfly_i32(
        data: &mut [i32],
        n_fft: usize,
        twiddle: &[u32],
        radix4_modifier: u16,
    ) {
        let n2 = n_fft >> 1;

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

        radix4_butterfly_i32(&mut data[..n_fft], n2, twiddle, radix4_modifier);
        radix4_butterfly_i32(&mut data[n_fft..], n2, twiddle, radix4_modifier);

        upshift_quads_i32(data);
    }

    fn radix4by2_butterfly_inverse_i32(
        data: &mut [i32],
        n_fft: usize,
        twiddle: &[u32],
        radix4_modifier: u16,
    ) {
        let n2 = n_fft >> 1;

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

        radix4_butterfly_inverse_i32(&mut data[..n_fft], n2, twiddle, radix4_modifier);
        radix4_butterfly_inverse_i32(&mut data[n_fft..], n2, twiddle, radix4_modifier);

        upshift_quads_i32(data);
    }

    /// Q31 CFFT/CIFFT entry with the same dispatch pattern as CMSIS `arm_cfft_q31`.
    pub fn run(&self, data: &mut [i32]) {
        assert_eq!(
            data.len(),
            self.n_fft * 2,
            "Q31 buffer length must be 2 * n_fft"
        );

        let use_radix4by2 = matches!(self.n_fft, 32 | 128 | 512 | 2048);

        if self.ifft_flag {
            if use_radix4by2 {
                Self::radix4by2_butterfly_inverse_i32(
                    data,
                    self.n_fft,
                    self.twiddle,
                    RADIX4BY2_MODIFIER,
                );
            } else {
                radix4_butterfly_inverse_i32(data, self.n_fft, self.twiddle, RADIX4_MODIFIER);
            }
        } else {
            if use_radix4by2 {
                Self::radix4by2_butterfly_i32(data, self.n_fft, self.twiddle, RADIX4BY2_MODIFIER);
            } else {
                radix4_butterfly_i32(data, self.n_fft, self.twiddle, RADIX4_MODIFIER);
            }
        }

        if self.bit_reverse_flag {
            bitreversal_i32(data, self.bit_rev_table);
        }
    }
}
