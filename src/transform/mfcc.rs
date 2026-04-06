use crate::basic::{
    div_i16, dot_i16, offset_i32, scale_i16, scale_i32, shift_i32, vlog_i32_in_place,
};
use crate::complex::cmplx_mag_i16;
use crate::matrix::{Matrix, mat_vec_mul_i16};
use crate::statistics::absmax_i16;
use crate::transform::RfftI16;

#[inline]
fn sat_i16(x: i32) -> i16 {
    x.clamp(i16::MIN as i32, i16::MAX as i32) as i16
}

#[inline]
fn sat_i32(x: i64) -> i32 {
    x.clamp(i32::MIN as i64, i32::MAX as i64) as i32
}

#[inline]
fn mul_q15(a: i16, b: i16) -> i16 {
    (((a as i32) * (b as i32)) >> 15).clamp(i16::MIN as i32, i16::MAX as i32) as i16
}

pub struct MfccI16 {
    pub n_fft: usize,
    pub n_mels: usize,
    pub n_mfcc: usize,
    pub dct: &'static [i16],
    pub filter: &'static [i16],
    pub filter_pos: &'static [i16],
    pub filter_len: &'static [i16],
    pub window: &'static [i16],
    pub rfft: RfftI16,
}

impl MfccI16 {
    pub const fn new(
        n_fft: usize,               // FFT size
        n_mels: usize,              // Number of Mel filters
        n_mfcc: usize,              // Number of MFCC coefficients
        dct: &'static [i16],        // DCT matrix for MFCC computation
        filter: &'static [i16],     // Mel filter bank coefficients
        filter_pos: &'static [i16], // Starting positions of each filter in the FFT output
        filter_len: &'static [i16], // Lengths of each filter
        window: &'static [i16],     // Window function coefficients
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

        let (max_abs, _) = absmax_i16(input, self.n_fft);

        if max_abs != 0 && max_abs != i16::MAX {
            let (quotient, shift) = div_i16(i16::MAX, max_abs).expect("division must succeed");
            scale_i16(input, quotient, shift as i8);
        }

        for (sample, &window) in input.iter_mut().zip(self.window.iter()) {
            *sample = mul_q15(*sample, window);
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
                &self.filter[packed_pos..filter_end],
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
            data: self.dct.as_ptr() as *mut i16,
        };
        mat_vec_mul_i16(dct, &input[..self.n_mels], output);
    }
}
