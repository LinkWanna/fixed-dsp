#[inline]
fn i32_mul(a: i32, b: i32) -> i32 {
    ((a as i64 * b as i64) >> 32) as i32
}

#[inline]
fn twiddle_mul_fwd_i16(xr: i16, xi: i16, co: i16, si: i16) -> (i16, i16) {
    let re = (((co as i32) * (xr as i32) + (si as i32) * (xi as i32)) >> 16) as i16;
    let im = (((co as i32) * (xi as i32) - (si as i32) * (xr as i32)) >> 16) as i16;
    (re, im)
}

#[inline]
fn twiddle_mul_inv_i16(xr: i16, xi: i16, co: i16, si: i16) -> (i16, i16) {
    let re = (((co as i32) * (xr as i32) - (si as i32) * (xi as i32)) >> 16) as i16;
    let im = (((co as i32) * (xi as i32) + (si as i32) * (xr as i32)) >> 16) as i16;
    (re, im)
}

#[inline]
fn twiddle_mul_fwd_i32(xr: i32, xi: i32, co: i32, si: i32) -> (i32, i32) {
    let re = i32_mul(xr, co).wrapping_add(i32_mul(xi, si));
    let im = i32_mul(xi, co).wrapping_sub(i32_mul(xr, si));
    (re, im)
}

#[inline]
fn twiddle_mul_inv_i32(xr: i32, xi: i32, co: i32, si: i32) -> (i32, i32) {
    let re = i32_mul(xr, co).wrapping_sub(i32_mul(xi, si));
    let im = i32_mul(xi, co).wrapping_add(i32_mul(xr, si));
    (re, im)
}

/// Performs radix-4 DIF butterfly computation for Q15 fixed-point FFT.
pub fn radix4_butterfly_i16(
    data: &mut [i16],
    n_fft: usize,
    twiddles: &[u16],
    mut twiddle_modifier: u16,
) {
    let mut n2 = n_fft;
    let mut n1;
    n2 >>= 2;

    let mut i0 = 0usize;
    let mut ic = 0usize;
    let mut j = n2;

    // Stage 1: wide butterflies with initial downscaling and twiddle application.
    while j > 0 {
        let i1 = i0 + n2;
        let i2 = i1 + n2;
        let i3 = i2 + n2;

        let x0_re = data[2 * i0] >> 2;
        let x0_im = data[2 * i0 + 1] >> 2;
        let x2_re = data[2 * i2] >> 2;
        let x2_im = data[2 * i2 + 1] >> 2;

        let mut sum02_re = x0_re + x2_re;
        let mut sum02_im = x0_im + x2_im;
        let diff02_re = x0_re - x2_re;
        let diff02_im = x0_im - x2_im;

        let x1_re = data[i1 * 2] >> 2;
        let x1_im = data[i1 * 2 + 1] >> 2;
        let x3_re = data[i3 * 2] >> 2;
        let x3_im = data[i3 * 2 + 1] >> 2;

        let sum13_re = x1_re + x3_re;
        let sum13_im = x1_im + x3_im;

        data[i0 * 2] = (sum02_re >> 1) + (sum13_re >> 1);
        data[i0 * 2 + 1] = (sum02_im >> 1) + (sum13_im >> 1);

        sum02_re = sum02_re - sum13_re;
        sum02_im = sum02_im - sum13_im;

        let co1 = twiddles[2 * ic] as i16;
        let si1 = twiddles[2 * ic + 1] as i16;
        let co2 = twiddles[4 * ic] as i16;
        let si2 = twiddles[4 * ic + 1] as i16;
        let co3 = twiddles[6 * ic] as i16;
        let si3 = twiddles[6 * ic + 1] as i16;

        let (out1_re, out1_im) = twiddle_mul_fwd_i16(sum02_re, sum02_im, co2, si2);
        data[i1 * 2] = out1_re;
        data[i1 * 2 + 1] = out1_im;

        let diff13_re = x1_re - x3_re;
        let diff13_im = x1_im - x3_im;

        let acc2_re = diff02_re + diff13_im;
        let acc2_im = diff02_im - diff13_re;
        let acc3_re = diff02_re - diff13_im;
        let acc3_im = diff02_im + diff13_re;

        let (out2_re, out2_im) = twiddle_mul_fwd_i16(acc2_re, acc2_im, co1, si1);
        data[i2 * 2] = out2_re;
        data[i2 * 2 + 1] = out2_im;

        let (out3_re, out3_im) = twiddle_mul_fwd_i16(acc3_re, acc3_im, co3, si3);
        data[i3 * 2] = out3_re;
        data[i3 * 2 + 1] = out3_im;

        ic += twiddle_modifier as usize;
        i0 += 1;
        j -= 1;
    }

    // Stage 2: recursively process smaller radix-4 groups.
    twiddle_modifier <<= 2;

    let mut k = n_fft / 4;
    while k > 4 {
        n1 = n2;
        n2 >>= 2;
        ic = 0;

        for j in 0..n2 {
            let co1 = twiddles[2 * ic] as i16;
            let si1 = twiddles[2 * ic + 1] as i16;
            let co2 = twiddles[4 * ic] as i16;
            let si2 = twiddles[4 * ic + 1] as i16;
            let co3 = twiddles[6 * ic] as i16;
            let si3 = twiddles[6 * ic + 1] as i16;
            ic += twiddle_modifier as usize;

            let mut i0 = j;
            while i0 < n_fft {
                let i1 = i0 + n2;
                let i2 = i1 + n2;
                let i3 = i2 + n2;

                let mut t0 = data[i0 * 2];
                let mut t1 = data[i0 * 2 + 1];
                let s0_in = data[i2 * 2];
                let s1_in = data[i2 * 2 + 1];

                let mut r0 = t0 + s0_in;
                let mut r1 = t1 + s1_in;
                let s0 = t0 - s0_in;
                let s1 = t1 - s1_in;

                t0 = data[i1 * 2];
                t1 = data[i1 * 2 + 1];
                let u0 = data[i3 * 2];
                let u1 = data[i3 * 2 + 1];
                let t0_saved = t0;
                let t1_saved = t1;
                let u0_saved = u0;
                let u1_saved = u1;

                t0 = t0 + u0;
                t1 = t1 + u1;

                let out0_re = ((r0 >> 1) + (t0 >> 1)) >> 1;
                let out0_im = ((r1 >> 1) + (t1 >> 1)) >> 1;
                data[i0 * 2] = out0_re;
                data[i0 * 2 + 1] = out0_im;

                r0 = (r0 >> 1) - (t0 >> 1);
                r1 = (r1 >> 1) - (t1 >> 1);

                let (out1_re, out1_im) = twiddle_mul_fwd_i16(r0, r1, co2, si2);
                data[i1 * 2] = out1_re;
                data[i1 * 2 + 1] = out1_im;

                let td0 = t0_saved - u0_saved;
                let td1 = t1_saved - u1_saved;

                let r0_b = (s0 >> 1) - (td1 >> 1);
                let r1_b = (s1 >> 1) + (td0 >> 1);
                let s0_b = (s0 >> 1) + (td1 >> 1);
                let s1_b = (s1 >> 1) - (td0 >> 1);

                let (out2_re, out2_im) = twiddle_mul_fwd_i16(s0_b, s1_b, co1, si1);
                data[i2 * 2] = out2_re;
                data[i2 * 2 + 1] = out2_im;

                let (out3_re, out3_im) = twiddle_mul_fwd_i16(r0_b, r1_b, co3, si3);
                data[i3 * 2] = out3_re;
                data[i3 * 2 + 1] = out3_im;

                i0 += n1;
            }
        }

        twiddle_modifier <<= 2;
        k >>= 2;
    }

    // Stage 3: terminal butterflies without twiddle multiplications.
    n1 = n2;
    n2 >>= 2;

    let mut i0 = 0usize;
    while i0 <= n_fft - n1 {
        let i1 = i0 + n2;
        let i2 = i1 + n2;
        let i3 = i2 + n2;

        let t0 = data[i0 * 2] as i16;
        let t1 = data[i0 * 2 + 1] as i16;
        let s0_in = data[i2 * 2] as i16;
        let s1_in = data[i2 * 2 + 1] as i16;

        let mut r0 = t0 + s0_in;
        let mut r1 = t1 + s1_in;
        let s0 = t0 - s0_in;
        let s1 = t1 - s1_in;

        let t0b = data[i1 * 2] as i16;
        let t1b = data[i1 * 2 + 1] as i16;
        let u0 = data[i3 * 2] as i16;
        let u1 = data[i3 * 2 + 1] as i16;

        let tsum0 = t0b + u0;
        let tsum1 = t1b + u1;

        data[i0 * 2] = (r0 >> 1) + (tsum0 >> 1);
        data[i0 * 2 + 1] = (r1 >> 1) + (tsum1 >> 1);

        r0 = (r0 >> 1) - (tsum0 >> 1);
        r1 = (r1 >> 1) - (tsum1 >> 1);

        data[i1 * 2] = r0;
        data[i1 * 2 + 1] = r1;

        let tdiff0 = t0b - u0;
        let tdiff1 = t1b - u1;

        data[i2 * 2] = (s0 >> 1) + (tdiff1 >> 1);
        data[i2 * 2 + 1] = (s1 >> 1) - (tdiff0 >> 1);
        data[i3 * 2] = (s0 >> 1) - (tdiff1 >> 1);
        data[i3 * 2 + 1] = (s1 >> 1) + (tdiff0 >> 1);

        i0 += n1;
    }
}

/// Performs radix-4 DIF butterfly computation for inverse Q15 fixed-point FFT.
pub fn radix4_butterfly_inverse_i16(
    data: &mut [i16],
    n_fft: usize,
    twiddles: &[u16],
    mut twiddle_modifier: u16,
) {
    let mut n2 = n_fft;
    let mut n1;
    n2 >>= 2;

    let mut ic = 0usize;
    let mut i0 = 0usize;
    let mut j = n2;

    // Stage 1: inverse wide butterflies with initial downscaling and twiddle application.
    while j > 0 {
        let i1 = i0 + n2;
        let i2 = i1 + n2;
        let i3 = i2 + n2;

        let mut t0 = (data[i0 * 2] as i16) >> 2;
        let mut t1 = (data[i0 * 2 + 1] as i16) >> 2;
        let s0_in = (data[i2 * 2] as i16) >> 2;
        let s1_in = (data[i2 * 2 + 1] as i16) >> 2;

        let mut r0 = (t0 + s0_in) as i16;
        let mut r1 = (t1 + s1_in) as i16;
        let s0 = (t0 - s0_in) as i16;
        let s1 = (t1 - s1_in) as i16;

        t0 = (data[i1 * 2] as i16) >> 2;
        t1 = (data[i1 * 2 + 1] as i16) >> 2;
        let u0 = (data[i3 * 2] as i16) >> 2;
        let u1 = (data[i3 * 2 + 1] as i16) >> 2;
        let t0_saved = t0;
        let t1_saved = t1;
        let u0_saved = u0;
        let u1_saved = u1;

        t0 = t0 + u0;
        t1 = t1 + u1;

        data[i0 * 2] = (r0 >> 1) + (t0 >> 1);
        data[i0 * 2 + 1] = (r1 >> 1) + (t1 >> 1);

        r0 = r0 - t0;
        r1 = r1 - t1;

        let co1 = twiddles[2 * ic] as i16;
        let si1 = twiddles[2 * ic + 1] as i16;
        let co2 = twiddles[4 * ic] as i16;
        let si2 = twiddles[4 * ic + 1] as i16;
        let co3 = twiddles[6 * ic] as i16;
        let si3 = twiddles[6 * ic + 1] as i16;

        let (out1_re, out1_im) = twiddle_mul_inv_i16(r0, r1, co2, si2);
        data[i1 * 2] = out1_re;
        data[i1 * 2 + 1] = out1_im;

        let td0 = t0_saved - u0_saved;
        let td1 = t1_saved - u1_saved;

        let r0_b = s0 + td1;
        let r1_b = s1 - td0;
        let s0_b = s0 - td1;
        let s1_b = s1 + td0;

        let (out2_re, out2_im) = twiddle_mul_inv_i16(s0_b, s1_b, co1, si1);
        data[i2 * 2] = out2_re;
        data[i2 * 2 + 1] = out2_im;

        let (out3_re, out3_im) = twiddle_mul_inv_i16(r0_b, r1_b, co3, si3);
        data[i3 * 2] = out3_re;
        data[i3 * 2 + 1] = out3_im;

        ic += twiddle_modifier as usize;
        i0 += 1;
        j -= 1;
    }

    // Stage 2: recursively process smaller inverse radix-4 groups.
    twiddle_modifier <<= 2;
    let mut k = n_fft / 4;
    while k > 4 {
        n1 = n2;
        n2 >>= 2;
        ic = 0;

        for j in 0..n2 {
            let co1 = twiddles[2 * ic] as i16;
            let si1 = twiddles[2 * ic + 1] as i16;
            let co2 = twiddles[4 * ic] as i16;
            let si2 = twiddles[4 * ic + 1] as i16;
            let co3 = twiddles[6 * ic] as i16;
            let si3 = twiddles[6 * ic + 1] as i16;
            ic += twiddle_modifier as usize;

            let mut i0 = j;
            while i0 < n_fft {
                let i1 = i0 + n2;
                let i2 = i1 + n2;
                let i3 = i2 + n2;

                let mut t0 = data[i0 * 2] as i16;
                let mut t1 = data[i0 * 2 + 1] as i16;
                let s0_in = data[i2 * 2] as i16;
                let s1_in = data[i2 * 2 + 1] as i16;

                let mut r0 = t0 + s0_in;
                let mut r1 = t1 + s1_in;
                let s0 = t0 - s0_in;
                let s1 = t1 - s1_in;

                t0 = data[i1 * 2] as i16;
                t1 = data[i1 * 2 + 1] as i16;
                let u0 = data[i3 * 2] as i16;
                let u1 = data[i3 * 2 + 1] as i16;
                let t0_saved = t0;
                let t1_saved = t1;
                let u0_saved = u0;
                let u1_saved = u1;

                t0 = t0 + u0;
                t1 = t1 + u1;

                let out0_re = ((r0 >> 1) + (t0 >> 1)) >> 1;
                let out0_im = ((r1 >> 1) + (t1 >> 1)) >> 1;
                data[i0 * 2] = out0_re;
                data[i0 * 2 + 1] = out0_im;

                r0 = (r0 >> 1) - (t0 >> 1);
                r1 = (r1 >> 1) - (t1 >> 1);

                let (out1_re, out1_im) = twiddle_mul_inv_i16(r0, r1, co2, si2);
                data[i1 * 2] = out1_re;
                data[i1 * 2 + 1] = out1_im;

                let td0 = t0_saved - u0_saved;
                let td1 = t1_saved - u1_saved;

                let r0_b = (s0 >> 1) + (td1 >> 1);
                let r1_b = (s1 >> 1) - (td0 >> 1);
                let s0_b = (s0 >> 1) - (td1 >> 1);
                let s1_b = (s1 >> 1) + (td0 >> 1);

                let (out2_re, out2_im) = twiddle_mul_inv_i16(s0_b, s1_b, co1, si1);
                data[i2 * 2] = out2_re;
                data[i2 * 2 + 1] = out2_im;

                let (out3_re, out3_im) = twiddle_mul_inv_i16(r0_b, r1_b, co3, si3);
                data[i3 * 2] = out3_re;
                data[i3 * 2 + 1] = out3_im;

                i0 += n1;
            }
        }

        twiddle_modifier <<= 2;
        k >>= 2;
    }

    // Stage 3: terminal inverse butterflies without twiddle multiplications.
    n1 = n2;
    n2 >>= 2;

    let mut i0 = 0usize;
    while i0 <= n_fft - n1 {
        let i1 = i0 + n2;
        let i2 = i1 + n2;
        let i3 = i2 + n2;

        let t0 = data[i0 * 2] as i16;
        let t1 = data[i0 * 2 + 1] as i16;
        let s0_in = data[i2 * 2] as i16;
        let s1_in = data[i2 * 2 + 1] as i16;

        let mut r0 = (t0 + s0_in) as i16;
        let mut r1 = (t1 + s1_in) as i16;
        let s0 = (t0 - s0_in) as i16;
        let s1 = (t1 - s1_in) as i16;

        let t0b = data[i1 * 2] as i16;
        let t1b = data[i1 * 2 + 1] as i16;
        let u0 = data[i3 * 2] as i16;
        let u1 = data[i3 * 2 + 1] as i16;

        let tsum0 = t0b + u0;
        let tsum1 = t1b + u1;

        data[i0 * 2] = (r0 >> 1) + (tsum0 >> 1);
        data[i0 * 2 + 1] = (r1 >> 1) + (tsum1 >> 1);

        r0 = (r0 >> 1) - (tsum0 >> 1);
        r1 = (r1 >> 1) - (tsum1 >> 1);

        data[i1 * 2] = r0;
        data[i1 * 2 + 1] = r1;

        let tdiff0 = t0b - u0;
        let tdiff1 = t1b - u1;

        data[i2 * 2] = (s0 >> 1) - (tdiff1 >> 1);
        data[i2 * 2 + 1] = (s1 >> 1) + (tdiff0 >> 1);
        data[i3 * 2] = (s0 >> 1) + (tdiff1 >> 1);
        data[i3 * 2 + 1] = (s1 >> 1) - (tdiff0 >> 1);

        i0 += n1;
    }
}

/// Performs radix-4 DIF butterfly computation for Q31 fixed-point FFT.
pub fn radix4_butterfly_i32(
    data: &mut [i32],
    n_fft: usize,
    twiddles: &[u32],
    mut twiddle_modifier: u16,
) {
    let mut n2 = n_fft;
    let mut n1;
    n2 >>= 2;

    let mut i0 = 0usize;
    let mut ic = 0usize;
    let mut j = n2;

    // Stage 1: wide butterflies with coarse scaling and twiddle multiplication.
    while j > 0 {
        let i1 = i0 + n2;
        let i2 = i1 + n2;
        let i3 = i2 + n2;

        let x0_re = data[2 * i0] >> 4;
        let x0_im = data[2 * i0 + 1] >> 4;
        let x2_re = data[2 * i2] >> 4;
        let x2_im = data[2 * i2 + 1] >> 4;

        let mut sum02_re = x0_re + x2_re;
        let mut sum02_im = x0_im + x2_im;
        let diff02_re = x0_re - x2_re;
        let diff02_im = x0_im - x2_im;

        let x1_re = data[2 * i1] >> 4;
        let x1_im = data[2 * i1 + 1] >> 4;
        let x3_re = data[2 * i3] >> 4;
        let x3_im = data[2 * i3 + 1] >> 4;

        let sum13_re = x1_re + x3_re;
        let sum13_im = x1_im + x3_im;

        data[2 * i0] = sum02_re + sum13_re;
        data[2 * i0 + 1] = sum02_im + sum13_im;
        sum02_re -= sum13_re;
        sum02_im -= sum13_im;

        let diff13_re = x1_re - x3_re;
        let diff13_im = x1_im - x3_im;

        let co1 = twiddles[2 * ic] as i32;
        let si1 = twiddles[2 * ic + 1] as i32;
        let co2 = twiddles[4 * ic] as i32;
        let si2 = twiddles[4 * ic + 1] as i32;
        let co3 = twiddles[6 * ic] as i32;
        let si3 = twiddles[6 * ic + 1] as i32;

        let (out1_re, out1_im) = twiddle_mul_fwd_i32(sum02_re, sum02_im, co2, si2);
        data[2 * i1] = out1_re << 1;
        data[2 * i1 + 1] = out1_im << 1;

        let acc2_re = diff02_re + diff13_im;
        let acc2_im = diff02_im - diff13_re;
        let acc3_re = diff02_re - diff13_im;
        let acc3_im = diff02_im + diff13_re;

        let (out2_re, out2_im) = twiddle_mul_fwd_i32(acc2_re, acc2_im, co1, si1);
        data[2 * i2] = out2_re << 1;
        data[2 * i2 + 1] = out2_im << 1;

        let (out3_re, out3_im) = twiddle_mul_fwd_i32(acc3_re, acc3_im, co3, si3);
        data[2 * i3] = out3_re << 1;
        data[2 * i3 + 1] = out3_im << 1;

        ic += twiddle_modifier as usize;
        i0 += 1;
        j -= 1;
    }

    // Stage 2: descend through radix-4 groups and apply per-group twiddles.
    twiddle_modifier <<= 2;

    let mut k = n_fft / 4;
    while k > 4 {
        n1 = n2;
        n2 >>= 2;
        ic = 0;

        for j in 0..n2 {
            let co1 = twiddles[2 * ic] as i32;
            let si1 = twiddles[2 * ic + 1] as i32;
            let co2 = twiddles[4 * ic] as i32;
            let si2 = twiddles[4 * ic + 1] as i32;
            let co3 = twiddles[6 * ic] as i32;
            let si3 = twiddles[6 * ic + 1] as i32;
            ic += twiddle_modifier as usize;

            let mut i0 = j;
            while i0 < n_fft {
                let i1 = i0 + n2;
                let i2 = i1 + n2;
                let i3 = i2 + n2;

                let mut r0 = data[2 * i0] + data[2 * i2];
                let mut r1 = data[2 * i0] - data[2 * i2];
                let mut s0 = data[2 * i0 + 1] + data[2 * i2 + 1];
                let mut s1 = data[2 * i0 + 1] - data[2 * i2 + 1];

                let mut t0 = data[2 * i1] + data[2 * i3];
                data[2 * i0] = (r0 + t0) >> 2;
                r0 -= t0;

                let t1 = data[2 * i1 + 1] + data[2 * i3 + 1];
                data[2 * i0 + 1] = (s0 + t1) >> 2;
                s0 -= t1;

                t0 = data[2 * i1 + 1] - data[2 * i3 + 1];
                let t1 = data[2 * i1] - data[2 * i3];

                let (out1_re, out1_im) = twiddle_mul_fwd_i32(r0, s0, co2, si2);
                data[2 * i1] = out1_re >> 1;
                data[2 * i1 + 1] = out1_im >> 1;

                r0 = r1 + t0;
                r1 -= t0;
                s0 = s1 - t1;
                s1 += t1;

                let (out2_re, out2_im) = twiddle_mul_fwd_i32(r0, s0, co1, si1);
                data[2 * i2] = out2_re >> 1;
                data[2 * i2 + 1] = out2_im >> 1;

                let (out3_re, out3_im) = twiddle_mul_fwd_i32(r1, s1, co3, si3);
                data[2 * i3] = out3_re >> 1;
                data[2 * i3 + 1] = out3_im >> 1;

                i0 += n1;
            }
        }

        twiddle_modifier <<= 2;
        k >>= 2;
    }

    // Stage 3: final recombination of 4-point blocks.
    let mut ptr = 0usize;
    let mut j = n_fft >> 2;
    while j > 0 {
        let xa = data[ptr];
        let ya = data[ptr + 1];
        let xb = data[ptr + 2];
        let yb = data[ptr + 3];
        let xc = data[ptr + 4];
        let yc = data[ptr + 5];
        let xd = data[ptr + 6];
        let yd = data[ptr + 7];

        data[ptr] = xa + xb + xc + xd;
        data[ptr + 1] = ya + yb + yc + yd;
        data[ptr + 2] = xa - xb + xc - xd;
        data[ptr + 3] = ya - yb + yc - yd;
        data[ptr + 4] = xa + yb - xc - yd;
        data[ptr + 5] = ya - xb - yc + xd;
        data[ptr + 6] = xa - yb - xc + yd;
        data[ptr + 7] = ya + xb - yc - xd;

        ptr += 8;
        j -= 1;
    }
}

/// Performs radix-4 DIF butterfly computation for inverse Q31 fixed-point FFT.
pub fn radix4_butterfly_inverse_i32(
    data: &mut [i32],
    n_fft: usize,
    twiddles: &[u32],
    mut twiddle_modifier: u16,
) {
    let mut n2 = n_fft;
    let mut n1;
    n2 >>= 2;

    let mut i0 = 0usize;
    let mut ic = 0usize;
    let mut j = n2;

    // Stage 1: inverse wide butterflies with coarse scaling and twiddle multiplication.
    while j > 0 {
        let i1 = i0 + n2;
        let i2 = i1 + n2;
        let i3 = i2 + n2;

        let mut r0 = (data[2 * i0] >> 4) + (data[2 * i2] >> 4);
        let mut r1 = (data[2 * i0] >> 4) - (data[2 * i2] >> 4);
        let mut t0 = (data[2 * i1] >> 4) + (data[2 * i3] >> 4);
        let mut s0 = (data[2 * i0 + 1] >> 4) + (data[2 * i2 + 1] >> 4);
        let mut s1 = (data[2 * i0 + 1] >> 4) - (data[2 * i2 + 1] >> 4);

        data[2 * i0] = r0 + t0;
        r0 -= t0;
        let t1 = (data[2 * i1 + 1] >> 4) + (data[2 * i3 + 1] >> 4);
        data[2 * i0 + 1] = s0 + t1;
        s0 -= t1;

        t0 = (data[2 * i1 + 1] >> 4) - (data[2 * i3 + 1] >> 4);
        let t1 = (data[2 * i1] >> 4) - (data[2 * i3] >> 4);

        let co1 = twiddles[2 * ic] as i32;
        let si1 = twiddles[2 * ic + 1] as i32;
        let co2 = twiddles[4 * ic] as i32;
        let si2 = twiddles[4 * ic + 1] as i32;
        let co3 = twiddles[6 * ic] as i32;
        let si3 = twiddles[6 * ic + 1] as i32;

        let (out1_re, out1_im) = twiddle_mul_inv_i32(r0, s0, co2, si2);
        data[2 * i1] = out1_re << 1;
        data[2 * i1 + 1] = out1_im << 1;

        r0 = r1 - t0;
        r1 += t0;
        s0 = s1 + t1;
        s1 -= t1;

        let (out2_re, out2_im) = twiddle_mul_inv_i32(r0, s0, co1, si1);
        data[2 * i2] = out2_re << 1;
        data[2 * i2 + 1] = out2_im << 1;

        let (out3_re, out3_im) = twiddle_mul_inv_i32(r1, s1, co3, si3);
        data[2 * i3] = out3_re << 1;
        data[2 * i3 + 1] = out3_im << 1;

        ic += twiddle_modifier as usize;
        i0 += 1;
        j -= 1;
    }

    // Stage 2: descend through inverse radix-4 groups and apply per-group twiddles.
    twiddle_modifier <<= 2;

    let mut k = n_fft / 4;
    while k > 4 {
        n1 = n2;
        n2 >>= 2;
        ic = 0;

        for j in 0..n2 {
            let co1 = twiddles[2 * ic] as i32;
            let si1 = twiddles[2 * ic + 1] as i32;
            let co2 = twiddles[4 * ic] as i32;
            let si2 = twiddles[4 * ic + 1] as i32;
            let co3 = twiddles[6 * ic] as i32;
            let si3 = twiddles[6 * ic + 1] as i32;
            ic += twiddle_modifier as usize;

            let mut i0 = j;
            while i0 < n_fft {
                let i1 = i0 + n2;
                let i2 = i1 + n2;
                let i3 = i2 + n2;

                let mut r0 = data[2 * i0] + data[2 * i2];
                let mut r1 = data[2 * i0] - data[2 * i2];
                let mut s0 = data[2 * i0 + 1] + data[2 * i2 + 1];
                let mut s1 = data[2 * i0 + 1] - data[2 * i2 + 1];

                let mut t0 = data[2 * i1] + data[2 * i3];
                data[2 * i0] = (r0 + t0) >> 2;
                r0 -= t0;

                let t1 = data[2 * i1 + 1] + data[2 * i3 + 1];
                data[2 * i0 + 1] = (s0 + t1) >> 2;
                s0 -= t1;

                t0 = data[2 * i1 + 1] - data[2 * i3 + 1];
                let t1 = data[2 * i1] - data[2 * i3];

                let (out1_re, out1_im) = twiddle_mul_inv_i32(r0, s0, co2, si2);
                data[2 * i1] = out1_re >> 1;
                data[2 * i1 + 1] = out1_im >> 1;

                r0 = r1 - t0;
                r1 += t0;
                s0 = s1 + t1;
                s1 -= t1;

                let (out2_re, out2_im) = twiddle_mul_inv_i32(r0, s0, co1, si1);
                data[2 * i2] = out2_re >> 1;
                data[2 * i2 + 1] = out2_im >> 1;

                let (out3_re, out3_im) = twiddle_mul_inv_i32(r1, s1, co3, si3);
                data[2 * i3] = out3_re >> 1;
                data[2 * i3 + 1] = out3_im >> 1;

                i0 += n1;
            }
        }

        twiddle_modifier <<= 2;
        k >>= 2;
    }

    // Stage 3: final inverse recombination of 4-point blocks.
    let mut ptr = 0usize;
    let mut j = n_fft >> 2;
    while j > 0 {
        let xa = data[ptr];
        let ya = data[ptr + 1];
        let xb = data[ptr + 2];
        let yb = data[ptr + 3];
        let xc = data[ptr + 4];
        let yc = data[ptr + 5];
        let xd = data[ptr + 6];
        let yd = data[ptr + 7];

        data[ptr] = xa + xb + xc + xd;
        data[ptr + 1] = ya + yb + yc + yd;
        data[ptr + 2] = xa - xb + xc - xd;
        data[ptr + 3] = ya - yb + yc - yd;
        data[ptr + 4] = xa - yb - xc + yd;
        data[ptr + 5] = ya + xb - yc - xd;
        data[ptr + 6] = xa + yb - xc - yd;
        data[ptr + 7] = ya - xb - yc + xd;

        ptr += 8;
        j -= 1;
    }
}
