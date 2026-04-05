use crate::common::tables::{TWIDDLE_TABLE_4096_U16, TWIDDLE_TABLE_4096_U32};

use super::{
    bitreversal_i16, bitreversal_i32, radix4_butterfly_i16, radix4_butterfly_i32,
    radix4_butterfly_inverse_i16, radix4_butterfly_inverse_i32,
};

#[inline]
fn q15_mul(a: i16, b: i16) -> i16 {
    (((a as i32) * (b as i32)) >> 16) as i16
}

#[inline]
fn q31_mul_r(a: i32, b: i32) -> i32 {
    (((a as i64 * b as i64) + 0x8000_0000) >> 32) as i32
}

fn radix4by2_butterfly_i16(data: &mut [i16], fft_len: usize, twiddle_modifier: usize) {
    let n2 = fft_len >> 1;

    for i in 0..n2 {
        let tw_idx = i * twiddle_modifier;
        let cos_val = TWIDDLE_TABLE_4096_U16[2 * tw_idx] as i16;
        let sin_val = TWIDDLE_TABLE_4096_U16[2 * tw_idx + 1] as i16;

        let l = i + n2;

        let t_re_half = (data[2 * i] as i32) >> 1;
        let t_im_half = (data[2 * i + 1] as i32) >> 1;
        let s_re_half = (data[2 * l] as i32) >> 1;
        let s_im_half = (data[2 * l + 1] as i32) >> 1;

        let xt = (t_re_half - s_re_half) as i16;
        let yt = (t_im_half - s_im_half) as i16;

        data[2 * i] = ((t_re_half + s_re_half) >> 1) as i16;
        data[2 * i + 1] = ((t_im_half + s_im_half) >> 1) as i16;

        let out_re = q15_mul(xt, cos_val).wrapping_add(q15_mul(yt, sin_val));
        let out_im = q15_mul(yt, cos_val).wrapping_sub(q15_mul(xt, sin_val));

        data[2 * l] = out_re;
        data[2 * l + 1] = out_im;
    }

    let radix4_modifier = (twiddle_modifier * 2) as u16;
    radix4_butterfly_i16(&mut data[..fft_len], n2, radix4_modifier);
    radix4_butterfly_i16(&mut data[fft_len..], n2, radix4_modifier);

    for i in 0..n2 {
        data[4 * i] = data[4 * i].wrapping_shl(1);
        data[4 * i + 1] = data[4 * i + 1].wrapping_shl(1);
        data[4 * i + 2] = data[4 * i + 2].wrapping_shl(1);
        data[4 * i + 3] = data[4 * i + 3].wrapping_shl(1);
    }
}

fn radix4by2_butterfly_inverse_i16(data: &mut [i16], fft_len: usize, twiddle_modifier: usize) {
    let n2 = fft_len >> 1;

    for i in 0..n2 {
        let tw_idx = i * twiddle_modifier;
        let cos_val = TWIDDLE_TABLE_4096_U16[2 * tw_idx] as i16;
        let sin_val = TWIDDLE_TABLE_4096_U16[2 * tw_idx + 1] as i16;

        let l = i + n2;

        let t_re_half = (data[2 * i] as i32) >> 1;
        let t_im_half = (data[2 * i + 1] as i32) >> 1;
        let s_re_half = (data[2 * l] as i32) >> 1;
        let s_im_half = (data[2 * l + 1] as i32) >> 1;

        let xt = (t_re_half - s_re_half) as i16;
        let yt = (t_im_half - s_im_half) as i16;

        data[2 * i] = ((t_re_half + s_re_half) >> 1) as i16;
        data[2 * i + 1] = ((t_im_half + s_im_half) >> 1) as i16;

        let out_re = q15_mul(xt, cos_val).wrapping_sub(q15_mul(yt, sin_val));
        let out_im = q15_mul(yt, cos_val).wrapping_add(q15_mul(xt, sin_val));

        data[2 * l] = out_re;
        data[2 * l + 1] = out_im;
    }

    let radix4_modifier = (twiddle_modifier * 2) as u16;
    radix4_butterfly_inverse_i16(&mut data[..fft_len], n2, radix4_modifier);
    radix4_butterfly_inverse_i16(&mut data[fft_len..], n2, radix4_modifier);

    for i in 0..n2 {
        data[4 * i] = data[4 * i].wrapping_shl(1);
        data[4 * i + 1] = data[4 * i + 1].wrapping_shl(1);
        data[4 * i + 2] = data[4 * i + 2].wrapping_shl(1);
        data[4 * i + 3] = data[4 * i + 3].wrapping_shl(1);
    }
}

/// Q15 CFFT/CIFFT entry compatible with CMSIS `arm_cfft_q15` dispatch behavior.
///
/// `data` is interleaved complex buffer with length `2 * fft_len`.
/// - `ifft_flag = false` => forward CFFT
/// - `ifft_flag = true` => inverse CIFFT
/// - `bit_reverse_flag` controls final in-place bit reversal.
pub fn cfft_i16(data: &mut [i16], fft_len: usize, ifft_flag: bool, bit_reverse_flag: bool) {
    assert_eq!(
        data.len(),
        fft_len * 2,
        "Q15 buffer length must be 2 * fft_len"
    );
    assert!(fft_len.is_power_of_two(), "fft_len must be a power of two");
    assert!(
        (16..=4096).contains(&fft_len),
        "fft_len must be in [16, 4096], got {}",
        fft_len
    );

    let twiddle_modifier = 4096 / fft_len;

    if ifft_flag {
        match fft_len {
            16 | 64 | 256 | 1024 | 4096 => {
                radix4_butterfly_inverse_i16(data, fft_len, twiddle_modifier as u16)
            }
            32 | 128 | 512 | 2048 => {
                radix4by2_butterfly_inverse_i16(data, fft_len, twiddle_modifier)
            }
            _ => unreachable!(),
        }
    } else {
        match fft_len {
            16 | 64 | 256 | 1024 | 4096 => {
                radix4_butterfly_i16(data, fft_len, twiddle_modifier as u16)
            }
            32 | 128 | 512 | 2048 => radix4by2_butterfly_i16(data, fft_len, twiddle_modifier),
            _ => unreachable!(),
        }
    }

    if bit_reverse_flag {
        bitreversal_i16(data, fft_len);
    }
}

fn radix4by2_butterfly_i32(data: &mut [i32], fft_len: usize, twiddle_modifier: usize) {
    let n2 = fft_len >> 1;

    for i in 0..n2 {
        let tw_idx = i * twiddle_modifier;
        let cos_val = TWIDDLE_TABLE_4096_U32[2 * tw_idx] as i32;
        let sin_val = TWIDDLE_TABLE_4096_U32[2 * tw_idx + 1] as i32;

        let l = i + n2;

        let xt = (data[2 * i] >> 2) - (data[2 * l] >> 2);
        data[2 * i] = (data[2 * i] >> 2) + (data[2 * l] >> 2);

        let yt = (data[2 * i + 1] >> 2) - (data[2 * l + 1] >> 2);
        data[2 * i + 1] = (data[2 * l + 1] >> 2) + (data[2 * i + 1] >> 2);

        let mut out_re = q31_mul_r(xt, cos_val);
        let mut out_im = q31_mul_r(yt, cos_val);
        out_re = out_re.wrapping_add(q31_mul_r(yt, sin_val));
        out_im = out_im.wrapping_sub(q31_mul_r(xt, sin_val));

        data[2 * l] = out_re.wrapping_shl(1);
        data[2 * l + 1] = out_im.wrapping_shl(1);
    }

    let radix4_modifier = (twiddle_modifier * 2) as u32;
    radix4_butterfly_i32(&mut data[..fft_len], n2, radix4_modifier as u16);
    radix4_butterfly_i32(&mut data[fft_len..], n2, radix4_modifier as u16);

    for i in 0..n2 {
        data[4 * i] = data[4 * i].wrapping_shl(1);
        data[4 * i + 1] = data[4 * i + 1].wrapping_shl(1);
        data[4 * i + 2] = data[4 * i + 2].wrapping_shl(1);
        data[4 * i + 3] = data[4 * i + 3].wrapping_shl(1);
    }
}

fn radix4by2_butterfly_inverse_i32(data: &mut [i32], fft_len: usize, twiddle_modifier: usize) {
    let n2 = fft_len >> 1;

    for i in 0..n2 {
        let tw_idx = i * twiddle_modifier;
        let cos_val = TWIDDLE_TABLE_4096_U32[2 * tw_idx] as i32;
        let sin_val = TWIDDLE_TABLE_4096_U32[2 * tw_idx + 1] as i32;

        let l = i + n2;

        let xt = (data[2 * i] >> 2) - (data[2 * l] >> 2);
        data[2 * i] = (data[2 * i] >> 2) + (data[2 * l] >> 2);

        let yt = (data[2 * i + 1] >> 2) - (data[2 * l + 1] >> 2);
        data[2 * i + 1] = (data[2 * l + 1] >> 2) + (data[2 * i + 1] >> 2);

        let mut out_re = q31_mul_r(xt, cos_val);
        let mut out_im = q31_mul_r(yt, cos_val);
        out_re = out_re.wrapping_sub(q31_mul_r(yt, sin_val));
        out_im = out_im.wrapping_add(q31_mul_r(xt, sin_val));

        data[2 * l] = out_re.wrapping_shl(1);
        data[2 * l + 1] = out_im.wrapping_shl(1);
    }

    let radix4_modifier = (twiddle_modifier * 2) as u32;
    radix4_butterfly_inverse_i32(&mut data[..fft_len], n2, radix4_modifier as u16);
    radix4_butterfly_inverse_i32(&mut data[fft_len..], n2, radix4_modifier as u16);

    for i in 0..n2 {
        data[4 * i] = data[4 * i].wrapping_shl(1);
        data[4 * i + 1] = data[4 * i + 1].wrapping_shl(1);
        data[4 * i + 2] = data[4 * i + 2].wrapping_shl(1);
        data[4 * i + 3] = data[4 * i + 3].wrapping_shl(1);
    }
}

/// Q31 CFFT/CIFFT entry with the same dispatch pattern as CMSIS `arm_cfft_q31`.
///
/// `data` is interleaved complex buffer with length `2 * fft_len`.
/// - `ifft_flag = false` => forward CFFT
/// - `ifft_flag = true` => inverse CIFFT
/// - `bit_reverse_flag` controls final in-place bit reversal.
pub fn cfft_i32(data: &mut [i32], fft_len: usize, ifft_flag: bool, bit_reverse_flag: bool) {
    assert_eq!(
        data.len(),
        fft_len * 2,
        "Q31 buffer length must be 2 * fft_len"
    );
    assert!(fft_len.is_power_of_two(), "fft_len must be a power of two");
    assert!(
        (16..=4096).contains(&fft_len),
        "fft_len must be in [16, 4096], got {}",
        fft_len
    );

    let twiddle_modifier = 4096 / fft_len;

    if ifft_flag {
        match fft_len {
            16 | 64 | 256 | 1024 | 4096 => {
                radix4_butterfly_inverse_i32(data, fft_len, twiddle_modifier as u16)
            }
            32 | 128 | 512 | 2048 => {
                radix4by2_butterfly_inverse_i32(data, fft_len, twiddle_modifier)
            }
            _ => unreachable!(),
        }
    } else {
        match fft_len {
            16 | 64 | 256 | 1024 | 4096 => {
                radix4_butterfly_i32(data, fft_len, twiddle_modifier as u16)
            }
            32 | 128 | 512 | 2048 => radix4by2_butterfly_i32(data, fft_len, twiddle_modifier),
            _ => unreachable!(),
        }
    }

    if bit_reverse_flag {
        bitreversal_i32(data, fft_len);
    }
}
