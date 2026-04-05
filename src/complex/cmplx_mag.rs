use crate::basic::sqrt_i32;

/// Q15 complex magnitude.
///
/// Input is interleaved complex samples: `[re0, im0, re1, im1, ...]` in 1.15.
/// Output is magnitude in 2.14 format.
pub fn cmplx_mag_i16(input: &[i16], output: &mut [i16]) {
    assert!(
        input.len().is_multiple_of(2),
        "input length must be even (interleaved complex)"
    );
    let num_samples = input.len() / 2;
    assert_eq!(
        output.len(),
        num_samples,
        "output length must equal number of complex samples"
    );

    for (chunk, out) in input.chunks_exact(2).zip(output.iter_mut()) {
        let re = chunk[0] as i64;
        let im = chunk[1] as i64;

        // CMSIS Q15: sqrt_q31((re*re + im*im) >> 1), then downscale to 2.14.
        let sum = (re * re + im * im) as u64;
        let in_sqrt = (sum >> 1) as i32;
        let res = sqrt_i32(in_sqrt).unwrap_or(0);
        *out = (res >> 16) as i16;
    }
}

/// Q31 complex magnitude.
///
/// Input is interleaved complex samples: `[re0, im0, re1, im1, ...]` in 1.31.
/// Output is magnitude in 2.30 format.
pub fn cmplx_mag_i32(input: &[i32], output: &mut [i32]) {
    assert!(
        input.len().is_multiple_of(2),
        "input length must be even (interleaved complex)"
    );
    let num_samples = input.len() / 2;
    assert_eq!(
        output.len(),
        num_samples,
        "output length must equal number of complex samples"
    );

    for (chunk, out) in input.chunks_exact(2).zip(output.iter_mut()) {
        let re = chunk[0] as i64;
        let im = chunk[1] as i64;

        // CMSIS Q31: acc = (re*re >> 33) + (im*im >> 33), then sqrt_q31(acc).
        let acc0 = (re * re) >> 33;
        let acc1 = (im * im) >> 33;
        let in_sqrt = (acc0 + acc1) as i32;
        *out = sqrt_i32(in_sqrt).unwrap_or(0);
    }
}
