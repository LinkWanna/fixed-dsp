use crate::common::tables::TWIDDLE_TABLE_4096_U16;

/// Performs radix-2 DIF butterfly computation for Q15 fixed-point FFT.
///
/// # Arguments
/// - `data`: Mutable slice containing 2*n_fft i16 elements (interleaved complex)
/// - `n_fft`: FFT size (must be power of 2)
/// - `twiddle_modifier`: Stride for twiddle factor lookup
pub fn radix2_butterfly_i16(data: &mut [i16], n_fft: usize, mut twiddle_modifier: u16) {
    let twiddles = &TWIDDLE_TABLE_4096_U16;

    let mut n2 = n_fft >> 1;
    let mut ia;

    // First stage: groups loop
    ia = 0;
    let mut n1 = n_fft;
    for _j in 0..n2 {
        let cos_val = twiddles[ia * 2] as i16;
        let sin_val = twiddles[ia * 2 + 1] as i16;
        ia += twiddle_modifier as usize;

        let mut i = _j;
        while i < n_fft {
            let l = i + n2;

            // Read values with 1 bit shift for first stage
            let t_re_half = (data[2 * i] as i32) >> 1;
            let t_im_half = (data[2 * i + 1] as i32) >> 1;
            let s_re_half = (data[2 * l] as i32) >> 1;
            let s_im_half = (data[2 * l + 1] as i32) >> 1;

            // Compute differences
            let xt = (t_re_half - s_re_half) as i16;
            let yt = (t_im_half - s_im_half) as i16;

            // Write sums with additional right shift
            data[2 * i] = ((t_re_half + s_re_half) >> 1) as i16;
            data[2 * i + 1] = ((t_im_half + s_im_half) >> 1) as i16;

            // Multiply by twiddle factor
            // out1 = (xt * cos + yt * sin) >> 16
            // out2 = (yt * cos - xt * sin) >> 16
            let mul1 = (((xt as i32) * (cos_val as i32)) >> 16) as i16;
            let mul2 = (((yt as i32) * (sin_val as i32)) >> 16) as i16;
            let mul3 = (((yt as i32) * (cos_val as i32)) >> 16) as i16;
            let mul4 = (((xt as i32) * (sin_val as i32)) >> 16) as i16;

            data[2 * l] = mul1.wrapping_add(mul2);
            data[2 * l + 1] = mul3.wrapping_sub(mul4);

            i += n1;
        }
    }

    twiddle_modifier <<= 1;

    // Middle stages
    let mut k = n_fft >> 1;
    while k > 2 {
        n1 = n2;
        n2 >>= 1;
        ia = 0;

        // Loop for groups
        for _j in 0..n2 {
            let cos_val = twiddles[ia * 2] as i16;
            let sin_val = twiddles[ia * 2 + 1] as i16;
            ia += twiddle_modifier as usize;

            let mut i = _j;
            while i < n_fft {
                let l = i + n2;

                let t_re = data[2 * i] as i32;
                let t_im = data[2 * i + 1] as i32;
                let s_re = data[2 * l] as i32;
                let s_im = data[2 * l + 1] as i32;

                // Compute differences
                let xt = (t_re - s_re) as i16;
                let yt = (t_im - s_im) as i16;

                // Write sums with right shift
                data[2 * i] = ((t_re + s_re) >> 1) as i16;
                data[2 * i + 1] = ((t_im + s_im) >> 1) as i16;

                // Complex multiply: follow CMSIS exactly
                // Each multiply-shift result is cast to i16, then added/subtracted
                let mul1 = (((xt as i32) * (cos_val as i32)) >> 16) as i16;
                let mul2 = (((yt as i32) * (sin_val as i32)) >> 16) as i16;
                let mul3 = (((yt as i32) * (cos_val as i32)) >> 16) as i16;
                let mul4 = (((xt as i32) * (sin_val as i32)) >> 16) as i16;

                data[2 * l] = mul1.wrapping_add(mul2);
                data[2 * l + 1] = mul3.wrapping_sub(mul4);

                i += n1;
            }
        }

        twiddle_modifier <<= 1;
        k >>= 1;
    }

    // Final stage (no twiddle multiplication)
    n1 = n2;
    n2 >>= 1;

    let mut i = 0;
    while i < n_fft {
        let l = i + n2;

        let t_re = data[2 * i] as i32;
        let t_im = data[2 * i + 1] as i32;
        let s_re = data[2 * l] as i32;
        let s_im = data[2 * l + 1] as i32;

        let xt = (t_re - s_re) as i16;
        let yt = (t_im - s_im) as i16;

        // Final stage: no shift
        data[2 * i] = (t_re + s_re) as i16;
        data[2 * i + 1] = (t_im + s_im) as i16;

        data[2 * l] = xt;
        data[2 * l + 1] = yt;

        i += n1;
    }
}

/// Performs radix-2 DIF inverse butterfly computation for Q15 fixed-point FFT.
pub fn radix2_butterfly_inverse_i16(data: &mut [i16], n_fft: usize, mut twiddle_modifier: u16) {
    let twiddles = &TWIDDLE_TABLE_4096_U16;

    let mut n2 = n_fft >> 1;
    let mut ia;

    // First stage: groups loop
    ia = 0;
    let mut n1 = n_fft;
    for _j in 0..n2 {
        let cos_val = twiddles[ia * 2] as i16;
        let sin_val = twiddles[ia * 2 + 1] as i16;
        ia += twiddle_modifier as usize;

        let mut i = _j;
        while i < n_fft {
            let l = i + n2;

            // Read values with 1 bit shift for first stage
            let t_re_half = (data[2 * i] as i32) >> 1;
            let t_im_half = (data[2 * i + 1] as i32) >> 1;
            let s_re_half = (data[2 * l] as i32) >> 1;
            let s_im_half = (data[2 * l + 1] as i32) >> 1;

            // Compute differences
            let xt = (t_re_half - s_re_half) as i16;
            let yt = (t_im_half - s_im_half) as i16;

            // Write sums with additional right shift
            data[2 * i] = ((t_re_half + s_re_half) >> 1) as i16;
            data[2 * i + 1] = ((t_im_half + s_im_half) >> 1) as i16;

            // Complex multiply conjugate (for inverse): follow CMSIS exactly
            // Each multiply-shift result is cast to i16, then combined
            let mul1 = (((xt as i32) * (cos_val as i32)) >> 16) as i16;
            let mul2 = (((yt as i32) * (sin_val as i32)) >> 16) as i16;
            let mul3 = (((yt as i32) * (cos_val as i32)) >> 16) as i16;
            let mul4 = (((xt as i32) * (sin_val as i32)) >> 16) as i16;

            data[2 * l] = mul1.wrapping_sub(mul2);
            data[2 * l + 1] = mul3.wrapping_add(mul4);

            i += n1;
        }
    }

    twiddle_modifier <<= 1;

    // Middle stages
    let mut k = n_fft >> 1;
    while k > 2 {
        n1 = n2;
        n2 >>= 1;
        ia = 0;

        for _j in 0..n2 {
            let cos_val = twiddles[ia * 2] as i16;
            let sin_val = twiddles[ia * 2 + 1] as i16;
            ia += twiddle_modifier as usize;

            let mut i = _j;
            while i < n_fft {
                let l = i + n2;

                let t_re = data[2 * i] as i32;
                let t_im = data[2 * i + 1] as i32;
                let s_re = data[2 * l] as i32;
                let s_im = data[2 * l + 1] as i32;

                // Compute differences
                let xt = (t_re - s_re) as i16;
                let yt = (t_im - s_im) as i16;

                // Write sums with right shift
                data[2 * i] = ((t_re + s_re) >> 1) as i16;
                data[2 * i + 1] = ((t_im + s_im) >> 1) as i16;

                // Complex multiply conjugate (for inverse): follow CMSIS exactly
                // Individual multiply-shift results cast to i16, then combined
                let mul1 = (((xt as i32) * (cos_val as i32)) >> 16) as i16;
                let mul2 = (((yt as i32) * (sin_val as i32)) >> 16) as i16;
                let mul3 = (((yt as i32) * (cos_val as i32)) >> 16) as i16;
                let mul4 = (((xt as i32) * (sin_val as i32)) >> 16) as i16;

                data[2 * l] = mul1.wrapping_sub(mul2);
                data[2 * l + 1] = mul3.wrapping_add(mul4);

                i += n1;
            }
        }

        twiddle_modifier <<= 1;
        k >>= 1;
    }

    // Final stage (no twiddle multiplication)
    n1 = n2;
    n2 >>= 1;

    let mut i = 0;
    while i < n_fft {
        let l = i + n2;

        let t_re = data[2 * i] as i32;
        let t_im = data[2 * i + 1] as i32;
        let s_re = data[2 * l] as i32;
        let s_im = data[2 * l + 1] as i32;

        let xt = (t_re - s_re) as i16;
        let yt = (t_im - s_im) as i16;

        // Final stage: no shift
        data[2 * i] = (t_re + s_re) as i16;
        data[2 * i + 1] = (t_im + s_im) as i16;

        data[2 * l] = xt;
        data[2 * l + 1] = yt;

        i += n1;
    }
}
