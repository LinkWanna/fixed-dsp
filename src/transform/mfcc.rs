use crate::basic::{
    div_i16, div_i32, dot_i16, dot_i32, offset_i32, scale_i16, scale_i32, shift_i32,
    vlog_i32_in_place,
};
use crate::complex::{cmplx_mag_i16, cmplx_mag_i32};
use crate::matrix::{Matrix, mat_vec_mul_i16, mat_vec_mul_i32};
use crate::statistics::{absmax_i16, absmax_i32};
use crate::transform::{RfftI16, RfftI32};
use crate::{sat_i16, sat_i32};

pub struct MfccI16 {
    pub n_fft: usize,
    pub n_mels: usize,
    pub n_mfcc: usize,
    pub dct: &'static [u16],
    pub filter: &'static [u16],
    pub filter_pos: &'static [u16],
    pub filter_len: &'static [u16],
    pub window: &'static [u16],
    pub rfft: RfftI16,
}

impl MfccI16 {
    pub const fn new(
        n_fft: usize,               // FFT size
        n_mels: usize,              // Number of Mel filters
        n_mfcc: usize,              // Number of MFCC coefficients
        dct: &'static [u16],        // DCT matrix for MFCC computation
        filter: &'static [u16],     // Mel filter bank coefficients
        filter_pos: &'static [u16], // Starting positions of each filter in the FFT output
        filter_len: &'static [u16], // Lengths of each filter
        window: &'static [u16],     // Window function coefficients
    ) -> Self {
        assert!(
            dct.len() == n_mfcc * n_mels,
            "DCT matrix size must be n_mfcc * n_mels"
        );
        assert!(
            filter_pos.len() == n_mels,
            "Filter position array size must be n_mels"
        );
        assert!(
            filter_len.len() == n_mels,
            "Filter length array size must be n_mels"
        );
        assert!(window.len() == n_fft, "Window size must be equal to n_fft");

        Self {
            n_fft,
            n_mels,
            n_mfcc,
            dct,
            filter,
            filter_pos,
            filter_len,
            window,
            rfft: RfftI16::new(n_fft, false, true),
        }
    }

    pub fn run(&self, input: &mut [i16], output: &mut [i16], tmp: &mut [i32]) {
        // CMSIS reuses the q31 tmp as q15 FFT scratch with a reinterpret cast.
        let tmp_q15 =
            unsafe { core::slice::from_raw_parts_mut(tmp.as_mut_ptr() as *mut i16, tmp.len() * 2) };
        let dct_i16 =
            unsafe { core::slice::from_raw_parts(self.dct.as_ptr() as *const i16, self.dct.len()) };
        let filter_i16 = unsafe {
            core::slice::from_raw_parts(self.filter.as_ptr() as *const i16, self.filter.len())
        };

        let (max_abs, _) = absmax_i16(input, self.n_fft);

        if max_abs != 0 && max_abs != i16::MAX {
            let (quotient, shift) = div_i16(i16::MAX, max_abs).expect("division must succeed");
            scale_i16(input, quotient, shift as i8);
        }

        for (sample, &window) in input.iter_mut().zip(self.window.iter()) {
            *sample = sat_i16(((*sample as i32) * (window as i32)) >> 15);
        }

        self.rfft.run(input, &mut tmp_q15[..self.n_fft * 2]);

        let filter_limit = 1 + (self.n_fft >> 1);
        cmplx_mag_i16(&tmp_q15[..filter_limit * 2], &mut input[..filter_limit]);

        let mut packed_pos = 0usize;
        for (index, (&filter_pos, &filter_len)) in self
            .filter_pos
            .iter()
            .zip(self.filter_len.iter())
            .enumerate()
        {
            let filter_pos = filter_pos as usize;
            let filter_len = filter_len as usize;
            let filter_end = packed_pos + filter_len;

            let acc = dot_i16(
                &input[filter_pos..filter_pos + filter_len],
                &filter_i16[packed_pos..filter_end],
            );
            packed_pos = filter_end;

            let acc = (acc + 0x219) >> 10;
            tmp[index] = sat_i32(acc);
        }

        if max_abs != 0 && max_abs != i16::MAX {
            scale_i32(&mut tmp[..self.n_mels], (max_abs as i32) << 16, 0);
        }

        vlog_i32_in_place(&mut tmp[..self.n_mels]);

        let fft_shift = self.n_fft.trailing_zeros() as i32;
        let log_exponent = (fft_shift + 2 + 10).wrapping_mul(0x02C5C860);
        offset_i32(&mut tmp[..self.n_mels], log_exponent);
        shift_i32(&mut tmp[..self.n_mels], -19);

        for (dst, &src) in input.iter_mut().zip(tmp.iter()).take(self.n_mels) {
            *dst = sat_i16(src);
        }

        let dct = Matrix {
            rows: self.n_mfcc,
            cols: self.n_mels,
            data: dct_i16.as_ptr() as *mut i16,
        };
        mat_vec_mul_i16(dct, &input[..self.n_mels], output);
    }
}

pub struct MfccI32 {
    pub n_fft: usize,
    pub n_mels: usize,
    pub n_mfcc: usize,
    pub dct: &'static [u32],
    pub filter: &'static [u32],
    pub filter_pos: &'static [u16],
    pub filter_len: &'static [u16],
    pub window: &'static [u32],
    pub rfft: RfftI32,
}

impl MfccI32 {
    pub const fn new(
        n_fft: usize,
        n_mels: usize,
        n_mfcc: usize,
        dct: &'static [u32],
        filter: &'static [u32],
        filter_pos: &'static [u16],
        filter_len: &'static [u16],
        window: &'static [u32],
    ) -> Self {
        assert!(
            dct.len() == n_mfcc * n_mels,
            "DCT matrix size must be n_mfcc * n_mels"
        );
        assert!(
            filter_pos.len() == n_mels,
            "Filter position array size must be n_mels"
        );
        assert!(
            filter_len.len() == n_mels,
            "Filter length array size must be n_mels"
        );
        assert!(window.len() == n_fft, "Window size must be equal to n_fft");

        Self {
            n_fft,
            n_mels,
            n_mfcc,
            dct,
            filter,
            filter_pos,
            filter_len,
            window,
            rfft: RfftI32::new(n_fft, false, true),
        }
    }

    pub fn run(&self, input: &mut [i32], output: &mut [i32], tmp: &mut [i32]) {
        const LOG2TOLOG_Q31: i32 = 0x02C5_C860;
        const MICRO_Q31: i64 = 0x0863_7BD0;
        const SHIFT_MELFILTER_SATURATION_Q31: i32 = 10;
        let dct =
            unsafe { core::slice::from_raw_parts(self.dct.as_ptr() as *const i32, self.dct.len()) };
        let filter = unsafe {
            core::slice::from_raw_parts(self.filter.as_ptr() as *const i32, self.filter.len())
        };

        assert!(
            tmp.len() >= self.n_fft * 2,
            "tmp length must be at least 2 * n_fft for q31 MFCC"
        );

        let (max_abs, _) = absmax_i32(input, self.n_fft);

        if max_abs != 0 && max_abs != i32::MAX {
            let (quotient, shift) = div_i32(i32::MAX, max_abs).expect("division must succeed");
            scale_i32(input, quotient, shift as i8);
        }

        for (sample, &window) in input.iter_mut().zip(self.window.iter()) {
            *sample = sat_i32((((*sample as i64) * (window as i64)) >> 32) << 1);
        }

        self.rfft.run(input, &mut tmp[..self.n_fft * 2]);

        let filter_limit = 1 + (self.n_fft >> 1);
        cmplx_mag_i32(&tmp[..filter_limit * 2], &mut input[..filter_limit]);

        let mut packed_pos = 0usize;
        for (index, (&filter_pos, &filter_len)) in self
            .filter_pos
            .iter()
            .zip(self.filter_len.iter())
            .enumerate()
        {
            let filter_pos = filter_pos as usize;
            let filter_len = filter_len as usize;
            let filter_end = packed_pos + filter_len;

            let mut acc = dot_i32(
                &input[filter_pos..filter_pos + filter_len],
                &filter[packed_pos..filter_end],
            );
            packed_pos = filter_end;

            acc = (acc + MICRO_Q31) >> (SHIFT_MELFILTER_SATURATION_Q31 + 18);
            tmp[index] = sat_i32(acc);
        }

        if max_abs != 0 && max_abs != i32::MAX {
            scale_i32(&mut tmp[..self.n_mels], max_abs, 0);
        }

        vlog_i32_in_place(&mut tmp[..self.n_mels]);

        let fft_shift = 31 - (self.n_fft as u32).leading_zeros() as i32;
        let log_exponent =
            (fft_shift + 2 + SHIFT_MELFILTER_SATURATION_Q31).wrapping_mul(LOG2TOLOG_Q31);
        offset_i32(&mut tmp[..self.n_mels], log_exponent);
        shift_i32(&mut tmp[..self.n_mels], -3);

        let dct = Matrix {
            rows: self.n_mfcc,
            cols: self.n_mels,
            data: dct.as_ptr() as *mut i32,
        };
        mat_vec_mul_i32(dct, &tmp[..self.n_mels], output);
    }
}
