#[inline]
fn sat_i16(v: i32) -> i16 {
    v.clamp(i16::MIN as i32, i16::MAX as i32) as i16
}

#[inline]
fn i32_mul(a: i32, b: i32) -> i32 {
    ((a as i64 * b as i64) >> 32) as i32
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

    let mut ic = 0usize;
    let mut i0 = 0usize;
    let mut j = n2;

    // First stage
    while j > 0 {
        let i1 = i0 + n2;
        let i2 = i1 + n2;
        let i3 = i2 + n2;

        let mut t0 = (data[i0 * 2] as i32) >> 2;
        let mut t1 = (data[i0 * 2 + 1] as i32) >> 2;
        let s0_in = (data[i2 * 2] as i32) >> 2;
        let s1_in = (data[i2 * 2 + 1] as i32) >> 2;

        let mut r0 = sat_i16(t0 + s0_in);
        let mut r1 = sat_i16(t1 + s1_in);
        let s0 = sat_i16(t0 - s0_in);
        let s1 = sat_i16(t1 - s1_in);

        t0 = (data[i1 * 2] as i32) >> 2;
        t1 = (data[i1 * 2 + 1] as i32) >> 2;
        let u0 = (data[i3 * 2] as i32) >> 2;
        let u1 = (data[i3 * 2 + 1] as i32) >> 2;
        let t0_saved = t0;
        let t1_saved = t1;
        let u0_saved = u0;
        let u1_saved = u1;

        t0 = sat_i16(t0 + u0) as i32;
        t1 = sat_i16(t1 + u1) as i32;

        data[i0 * 2] = ((r0 as i32 >> 1) + (t0 >> 1)) as i16;
        data[i0 * 2 + 1] = ((r1 as i32 >> 1) + (t1 >> 1)) as i16;

        r0 = sat_i16(r0 as i32 - t0);
        r1 = sat_i16(r1 as i32 - t1);

        let co1 = twiddles[2 * ic] as i16;
        let si1 = twiddles[2 * ic + 1] as i16;
        let co2 = twiddles[4 * ic] as i16;
        let si2 = twiddles[4 * ic + 1] as i16;
        let co3 = twiddles[6 * ic] as i16;
        let si3 = twiddles[6 * ic + 1] as i16;

        let out1_i1 = (((co2 as i32) * (r0 as i32) + (si2 as i32) * (r1 as i32)) >> 16) as i16;
        let out2_i1 = (((-(si2 as i32)) * (r0 as i32) + (co2 as i32) * (r1 as i32)) >> 16) as i16;
        data[i1 * 2] = out1_i1;
        data[i1 * 2 + 1] = out2_i1;

        let td0 = sat_i16(t0_saved - u0_saved);
        let td1 = sat_i16(t1_saved - u1_saved);

        let r0_b = sat_i16(s0 as i32 - td1 as i32);
        let r1_b = sat_i16(s1 as i32 + td0 as i32);
        let s0_b = sat_i16(s0 as i32 + td1 as i32);
        let s1_b = sat_i16(s1 as i32 - td0 as i32);

        let out1_i2 = (((si1 as i32) * (s1_b as i32) + (co1 as i32) * (s0_b as i32)) >> 16) as i16;
        let out2_i2 =
            (((-(si1 as i32)) * (s0_b as i32) + (co1 as i32) * (s1_b as i32)) >> 16) as i16;
        data[i2 * 2] = out1_i2;
        data[i2 * 2 + 1] = out2_i2;

        let out1_i3 = (((si3 as i32) * (r1_b as i32) + (co3 as i32) * (r0_b as i32)) >> 16) as i16;
        let out2_i3 =
            (((-(si3 as i32)) * (r0_b as i32) + (co3 as i32) * (r1_b as i32)) >> 16) as i16;
        data[i3 * 2] = out1_i3;
        data[i3 * 2 + 1] = out2_i3;

        ic += twiddle_modifier as usize;
        i0 += 1;
        j -= 1;
    }

    // Middle stages
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

                let mut t0 = data[i0 * 2] as i32;
                let mut t1 = data[i0 * 2 + 1] as i32;
                let s0_in = data[i2 * 2] as i32;
                let s1_in = data[i2 * 2 + 1] as i32;

                let mut r0 = sat_i16(t0 + s0_in);
                let mut r1 = sat_i16(t1 + s1_in);
                let s0 = sat_i16(t0 - s0_in);
                let s1 = sat_i16(t1 - s1_in);

                t0 = data[i1 * 2] as i32;
                t1 = data[i1 * 2 + 1] as i32;
                let u0 = data[i3 * 2] as i32;
                let u1 = data[i3 * 2 + 1] as i32;
                let t0_saved = t0;
                let t1_saved = t1;
                let u0_saved = u0;
                let u1_saved = u1;

                t0 = sat_i16(t0 + u0) as i32;
                t1 = sat_i16(t1 + u1) as i32;

                let out0_re = (((r0 as i32 >> 1) + (t0 >> 1)) >> 1) as i16;
                let out0_im = (((r1 as i32 >> 1) + (t1 >> 1)) >> 1) as i16;
                data[i0 * 2] = out0_re;
                data[i0 * 2 + 1] = out0_im;

                r0 = ((r0 as i32 >> 1) - (t0 >> 1)) as i16;
                r1 = ((r1 as i32 >> 1) - (t1 >> 1)) as i16;

                let out1_re =
                    (((co2 as i32) * (r0 as i32) + (si2 as i32) * (r1 as i32)) >> 16) as i16;
                let out1_im =
                    (((-(si2 as i32)) * (r0 as i32) + (co2 as i32) * (r1 as i32)) >> 16) as i16;
                data[i1 * 2] = out1_re;
                data[i1 * 2 + 1] = out1_im;

                let td0 = sat_i16(t0_saved - u0_saved);
                let td1 = sat_i16(t1_saved - u1_saved);

                let r0_b = ((s0 as i32 >> 1) - (td1 as i32 >> 1)) as i16;
                let r1_b = ((s1 as i32 >> 1) + (td0 as i32 >> 1)) as i16;
                let s0_b = ((s0 as i32 >> 1) + (td1 as i32 >> 1)) as i16;
                let s1_b = ((s1 as i32 >> 1) - (td0 as i32 >> 1)) as i16;

                let out2_re =
                    (((co1 as i32) * (s0_b as i32) + (si1 as i32) * (s1_b as i32)) >> 16) as i16;
                let out2_im =
                    (((-(si1 as i32)) * (s0_b as i32) + (co1 as i32) * (s1_b as i32)) >> 16) as i16;
                data[i2 * 2] = out2_re;
                data[i2 * 2 + 1] = out2_im;

                let out3_re =
                    (((co3 as i32) * (r0_b as i32) + (si3 as i32) * (r1_b as i32)) >> 16) as i16;
                let out3_im =
                    (((-(si3 as i32)) * (r0_b as i32) + (co3 as i32) * (r1_b as i32)) >> 16) as i16;
                data[i3 * 2] = out3_re;
                data[i3 * 2 + 1] = out3_im;

                i0 += n1;
            }
        }

        twiddle_modifier <<= 2;
        k >>= 2;
    }

    // Final stage
    n1 = n2;
    n2 >>= 2;

    let mut i0 = 0usize;
    while i0 <= n_fft - n1 {
        let i1 = i0 + n2;
        let i2 = i1 + n2;
        let i3 = i2 + n2;

        let t0 = data[i0 * 2] as i32;
        let t1 = data[i0 * 2 + 1] as i32;
        let s0_in = data[i2 * 2] as i32;
        let s1_in = data[i2 * 2 + 1] as i32;

        let mut r0 = sat_i16(t0 + s0_in);
        let mut r1 = sat_i16(t1 + s1_in);
        let s0 = sat_i16(t0 - s0_in);
        let s1 = sat_i16(t1 - s1_in);

        let t0b = data[i1 * 2] as i32;
        let t1b = data[i1 * 2 + 1] as i32;
        let u0 = data[i3 * 2] as i32;
        let u1 = data[i3 * 2 + 1] as i32;

        let tsum0 = sat_i16(t0b + u0) as i32;
        let tsum1 = sat_i16(t1b + u1) as i32;

        data[i0 * 2] = ((r0 as i32 >> 1) + (tsum0 >> 1)) as i16;
        data[i0 * 2 + 1] = ((r1 as i32 >> 1) + (tsum1 >> 1)) as i16;

        r0 = ((r0 as i32 >> 1) - (tsum0 >> 1)) as i16;
        r1 = ((r1 as i32 >> 1) - (tsum1 >> 1)) as i16;

        data[i1 * 2] = r0;
        data[i1 * 2 + 1] = r1;

        let tdiff0 = sat_i16(t0b - u0) as i32;
        let tdiff1 = sat_i16(t1b - u1) as i32;

        data[i2 * 2] = ((s0 as i32 >> 1) + (tdiff1 >> 1)) as i16;
        data[i2 * 2 + 1] = ((s1 as i32 >> 1) - (tdiff0 >> 1)) as i16;
        data[i3 * 2] = ((s0 as i32 >> 1) - (tdiff1 >> 1)) as i16;
        data[i3 * 2 + 1] = ((s1 as i32 >> 1) + (tdiff0 >> 1)) as i16;

        i0 += n1;
    }
}

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

    // First stage
    while j > 0 {
        let i1 = i0 + n2;
        let i2 = i1 + n2;
        let i3 = i2 + n2;

        let mut t0 = (data[i0 * 2] as i32) >> 2;
        let mut t1 = (data[i0 * 2 + 1] as i32) >> 2;
        let s0_in = (data[i2 * 2] as i32) >> 2;
        let s1_in = (data[i2 * 2 + 1] as i32) >> 2;

        let mut r0 = sat_i16(t0 + s0_in);
        let mut r1 = sat_i16(t1 + s1_in);
        let s0 = sat_i16(t0 - s0_in);
        let s1 = sat_i16(t1 - s1_in);

        t0 = (data[i1 * 2] as i32) >> 2;
        t1 = (data[i1 * 2 + 1] as i32) >> 2;
        let u0 = (data[i3 * 2] as i32) >> 2;
        let u1 = (data[i3 * 2 + 1] as i32) >> 2;
        let t0_saved = t0;
        let t1_saved = t1;
        let u0_saved = u0;
        let u1_saved = u1;

        t0 = sat_i16(t0 + u0) as i32;
        t1 = sat_i16(t1 + u1) as i32;

        data[i0 * 2] = ((r0 as i32 >> 1) + (t0 >> 1)) as i16;
        data[i0 * 2 + 1] = ((r1 as i32 >> 1) + (t1 >> 1)) as i16;

        r0 = sat_i16(r0 as i32 - t0);
        r1 = sat_i16(r1 as i32 - t1);

        let co1 = twiddles[2 * ic] as i16;
        let si1 = twiddles[2 * ic + 1] as i16;
        let co2 = twiddles[4 * ic] as i16;
        let si2 = twiddles[4 * ic + 1] as i16;
        let co3 = twiddles[6 * ic] as i16;
        let si3 = twiddles[6 * ic + 1] as i16;

        let out1_i1 = (((co2 as i32) * (r0 as i32) - (si2 as i32) * (r1 as i32)) >> 16) as i16;
        let out2_i1 = (((si2 as i32) * (r0 as i32) + (co2 as i32) * (r1 as i32)) >> 16) as i16;
        data[i1 * 2] = out1_i1;
        data[i1 * 2 + 1] = out2_i1;

        let td0 = sat_i16(t0_saved - u0_saved);
        let td1 = sat_i16(t1_saved - u1_saved);

        let r0_b = sat_i16(s0 as i32 + td1 as i32);
        let r1_b = sat_i16(s1 as i32 - td0 as i32);
        let s0_b = sat_i16(s0 as i32 - td1 as i32);
        let s1_b = sat_i16(s1 as i32 + td0 as i32);

        let out1_i2 = (((co1 as i32) * (s0_b as i32) - (si1 as i32) * (s1_b as i32)) >> 16) as i16;
        let out2_i2 = (((si1 as i32) * (s0_b as i32) + (co1 as i32) * (s1_b as i32)) >> 16) as i16;
        data[i2 * 2] = out1_i2;
        data[i2 * 2 + 1] = out2_i2;

        let out1_i3 = (((co3 as i32) * (r0_b as i32) - (si3 as i32) * (r1_b as i32)) >> 16) as i16;
        let out2_i3 = (((si3 as i32) * (r0_b as i32) + (co3 as i32) * (r1_b as i32)) >> 16) as i16;
        data[i3 * 2] = out1_i3;
        data[i3 * 2 + 1] = out2_i3;

        ic += twiddle_modifier as usize;
        i0 += 1;
        j -= 1;
    }

    // Middle stages
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

                let mut t0 = data[i0 * 2] as i32;
                let mut t1 = data[i0 * 2 + 1] as i32;
                let s0_in = data[i2 * 2] as i32;
                let s1_in = data[i2 * 2 + 1] as i32;

                let mut r0 = sat_i16(t0 + s0_in);
                let mut r1 = sat_i16(t1 + s1_in);
                let s0 = sat_i16(t0 - s0_in);
                let s1 = sat_i16(t1 - s1_in);

                t0 = data[i1 * 2] as i32;
                t1 = data[i1 * 2 + 1] as i32;
                let u0 = data[i3 * 2] as i32;
                let u1 = data[i3 * 2 + 1] as i32;
                let t0_saved = t0;
                let t1_saved = t1;
                let u0_saved = u0;
                let u1_saved = u1;

                t0 = sat_i16(t0 + u0) as i32;
                t1 = sat_i16(t1 + u1) as i32;

                let out0_re = (((r0 as i32 >> 1) + (t0 >> 1)) >> 1) as i16;
                let out0_im = (((r1 as i32 >> 1) + (t1 >> 1)) >> 1) as i16;
                data[i0 * 2] = out0_re;
                data[i0 * 2 + 1] = out0_im;

                r0 = ((r0 as i32 >> 1) - (t0 >> 1)) as i16;
                r1 = ((r1 as i32 >> 1) - (t1 >> 1)) as i16;

                let out1_re =
                    (((co2 as i32) * (r0 as i32) - (si2 as i32) * (r1 as i32)) >> 16) as i16;
                let out1_im =
                    (((si2 as i32) * (r0 as i32) + (co2 as i32) * (r1 as i32)) >> 16) as i16;
                data[i1 * 2] = out1_re;
                data[i1 * 2 + 1] = out1_im;

                let td0 = sat_i16(t0_saved - u0_saved);
                let td1 = sat_i16(t1_saved - u1_saved);

                let r0_b = ((s0 as i32 >> 1) + (td1 as i32 >> 1)) as i16;
                let r1_b = ((s1 as i32 >> 1) - (td0 as i32 >> 1)) as i16;
                let s0_b = ((s0 as i32 >> 1) - (td1 as i32 >> 1)) as i16;
                let s1_b = ((s1 as i32 >> 1) + (td0 as i32 >> 1)) as i16;

                let out2_re =
                    (((co1 as i32) * (s0_b as i32) - (si1 as i32) * (s1_b as i32)) >> 16) as i16;
                let out2_im =
                    (((si1 as i32) * (s0_b as i32) + (co1 as i32) * (s1_b as i32)) >> 16) as i16;
                data[i2 * 2] = out2_re;
                data[i2 * 2 + 1] = out2_im;

                let out3_re =
                    (((co3 as i32) * (r0_b as i32) - (si3 as i32) * (r1_b as i32)) >> 16) as i16;
                let out3_im =
                    (((si3 as i32) * (r0_b as i32) + (co3 as i32) * (r1_b as i32)) >> 16) as i16;
                data[i3 * 2] = out3_re;
                data[i3 * 2 + 1] = out3_im;

                i0 += n1;
            }
        }

        twiddle_modifier <<= 2;
        k >>= 2;
    }

    // Final stage
    n1 = n2;
    n2 >>= 2;

    let mut i0 = 0usize;
    while i0 <= n_fft - n1 {
        let i1 = i0 + n2;
        let i2 = i1 + n2;
        let i3 = i2 + n2;

        let t0 = data[i0 * 2] as i32;
        let t1 = data[i0 * 2 + 1] as i32;
        let s0_in = data[i2 * 2] as i32;
        let s1_in = data[i2 * 2 + 1] as i32;

        let mut r0 = sat_i16(t0 + s0_in);
        let mut r1 = sat_i16(t1 + s1_in);
        let s0 = sat_i16(t0 - s0_in);
        let s1 = sat_i16(t1 - s1_in);

        let t0b = data[i1 * 2] as i32;
        let t1b = data[i1 * 2 + 1] as i32;
        let u0 = data[i3 * 2] as i32;
        let u1 = data[i3 * 2 + 1] as i32;

        let tsum0 = sat_i16(t0b + u0) as i32;
        let tsum1 = sat_i16(t1b + u1) as i32;

        data[i0 * 2] = ((r0 as i32 >> 1) + (tsum0 >> 1)) as i16;
        data[i0 * 2 + 1] = ((r1 as i32 >> 1) + (tsum1 >> 1)) as i16;

        r0 = ((r0 as i32 >> 1) - (tsum0 >> 1)) as i16;
        r1 = ((r1 as i32 >> 1) - (tsum1 >> 1)) as i16;

        data[i1 * 2] = r0;
        data[i1 * 2 + 1] = r1;

        let tdiff0 = sat_i16(t0b - u0) as i32;
        let tdiff1 = sat_i16(t1b - u1) as i32;

        data[i2 * 2] = ((s0 as i32 >> 1) - (tdiff1 >> 1)) as i16;
        data[i2 * 2 + 1] = ((s1 as i32 >> 1) + (tdiff0 >> 1)) as i16;
        data[i3 * 2] = ((s0 as i32 >> 1) + (tdiff1 >> 1)) as i16;
        data[i3 * 2 + 1] = ((s1 as i32 >> 1) - (tdiff0 >> 1)) as i16;

        i0 += n1;
    }
}

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
    let mut ia1 = 0usize;
    let mut j = n2;

    while j > 0 {
        let i1 = i0 + n2;
        let i2 = i1 + n2;
        let i3 = i2 + n2;

        let mut r1 = (data[2 * i0] >> 4) + (data[2 * i2] >> 4);
        let mut r2 = (data[2 * i0] >> 4) - (data[2 * i2] >> 4);
        let mut t1 = (data[2 * i1] >> 4) + (data[2 * i3] >> 4);
        let mut s1 = (data[2 * i0 + 1] >> 4) + (data[2 * i2 + 1] >> 4);
        let mut s2 = (data[2 * i0 + 1] >> 4) - (data[2 * i2 + 1] >> 4);

        data[2 * i0] = r1 + t1;
        r1 -= t1;
        let t2 = (data[2 * i1 + 1] >> 4) + (data[2 * i3 + 1] >> 4);
        data[2 * i0 + 1] = s1 + t2;
        s1 -= t2;

        t1 = (data[2 * i1 + 1] >> 4) - (data[2 * i3 + 1] >> 4);
        let t2 = (data[2 * i1] >> 4) - (data[2 * i3] >> 4);

        let ia2 = 2 * ia1;
        let co2 = twiddles[2 * ia2] as i32;
        let si2 = twiddles[2 * ia2 + 1] as i32;

        data[2 * i1] = (i32_mul(r1, co2) + i32_mul(s1, si2)) << 1;
        data[2 * i1 + 1] = (i32_mul(s1, co2) - i32_mul(r1, si2)) << 1;

        r1 = r2 + t1;
        r2 -= t1;
        s1 = s2 - t2;
        s2 += t2;

        let co1 = twiddles[2 * ia1] as i32;
        let si1 = twiddles[2 * ia1 + 1] as i32;
        data[2 * i2] = (i32_mul(r1, co1) + i32_mul(s1, si1)) << 1;
        data[2 * i2 + 1] = (i32_mul(s1, co1) - i32_mul(r1, si1)) << 1;

        let ia3 = 3 * ia1;
        let co3 = twiddles[2 * ia3] as i32;
        let si3 = twiddles[2 * ia3 + 1] as i32;
        data[2 * i3] = (i32_mul(r2, co3) + i32_mul(s2, si3)) << 1;
        data[2 * i3 + 1] = (i32_mul(s2, co3) - i32_mul(r2, si3)) << 1;

        ia1 += twiddle_modifier as usize;
        i0 += 1;
        j -= 1;
    }

    twiddle_modifier <<= 2;

    let mut k = n_fft / 4;
    while k > 4 {
        n1 = n2;
        n2 >>= 2;
        ia1 = 0;

        for j in 0..n2 {
            let ia2 = ia1 + ia1;
            let ia3 = ia2 + ia1;

            let co1 = twiddles[2 * ia1] as i32;
            let si1 = twiddles[2 * ia1 + 1] as i32;
            let co2 = twiddles[2 * ia2] as i32;
            let si2 = twiddles[2 * ia2 + 1] as i32;
            let co3 = twiddles[2 * ia3] as i32;
            let si3 = twiddles[2 * ia3 + 1] as i32;
            ia1 += twiddle_modifier as usize;

            let mut i0 = j;
            while i0 < n_fft {
                let i1 = i0 + n2;
                let i2 = i1 + n2;
                let i3 = i2 + n2;

                let mut r1 = data[2 * i0] + data[2 * i2];
                let mut r2 = data[2 * i0] - data[2 * i2];
                let mut s1 = data[2 * i0 + 1] + data[2 * i2 + 1];
                let mut s2 = data[2 * i0 + 1] - data[2 * i2 + 1];

                let mut t1 = data[2 * i1] + data[2 * i3];
                data[2 * i0] = (r1 + t1) >> 2;
                r1 -= t1;

                let t2 = data[2 * i1 + 1] + data[2 * i3 + 1];
                data[2 * i0 + 1] = (s1 + t2) >> 2;
                s1 -= t2;

                t1 = data[2 * i1 + 1] - data[2 * i3 + 1];
                let t2 = data[2 * i1] - data[2 * i3];

                data[2 * i1] = (i32_mul(r1, co2) + i32_mul(s1, si2)) >> 1;
                data[2 * i1 + 1] = (i32_mul(s1, co2) - i32_mul(r1, si2)) >> 1;

                r1 = r2 + t1;
                r2 -= t1;
                s1 = s2 - t2;
                s2 += t2;

                data[2 * i2] = (i32_mul(r1, co1) + i32_mul(s1, si1)) >> 1;
                data[2 * i2 + 1] = (i32_mul(s1, co1) - i32_mul(r1, si1)) >> 1;

                data[2 * i3] = (i32_mul(r2, co3) + i32_mul(s2, si3)) >> 1;
                data[2 * i3 + 1] = (i32_mul(s2, co3) - i32_mul(r2, si3)) >> 1;

                i0 += n1;
            }
        }

        twiddle_modifier <<= 2;
        k >>= 2;
    }

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
    let mut ia1 = 0usize;
    let mut j = n2;

    while j > 0 {
        let i1 = i0 + n2;
        let i2 = i1 + n2;
        let i3 = i2 + n2;

        let mut r1 = (data[2 * i0] >> 4) + (data[2 * i2] >> 4);
        let mut r2 = (data[2 * i0] >> 4) - (data[2 * i2] >> 4);
        let mut t1 = (data[2 * i1] >> 4) + (data[2 * i3] >> 4);
        let mut s1 = (data[2 * i0 + 1] >> 4) + (data[2 * i2 + 1] >> 4);
        let mut s2 = (data[2 * i0 + 1] >> 4) - (data[2 * i2 + 1] >> 4);

        data[2 * i0] = r1 + t1;
        r1 -= t1;
        let t2 = (data[2 * i1 + 1] >> 4) + (data[2 * i3 + 1] >> 4);
        data[2 * i0 + 1] = s1 + t2;
        s1 -= t2;

        t1 = (data[2 * i1 + 1] >> 4) - (data[2 * i3 + 1] >> 4);
        let t2 = (data[2 * i1] >> 4) - (data[2 * i3] >> 4);

        let ia2 = 2 * ia1;
        let co2 = twiddles[2 * ia2] as i32;
        let si2 = twiddles[2 * ia2 + 1] as i32;

        data[2 * i1] = (i32_mul(r1, co2) - i32_mul(s1, si2)) << 1;
        data[2 * i1 + 1] = (i32_mul(s1, co2) + i32_mul(r1, si2)) << 1;

        r1 = r2 - t1;
        r2 += t1;
        s1 = s2 + t2;
        s2 -= t2;

        let co1 = twiddles[2 * ia1] as i32;
        let si1 = twiddles[2 * ia1 + 1] as i32;
        data[2 * i2] = (i32_mul(r1, co1) - i32_mul(s1, si1)) << 1;
        data[2 * i2 + 1] = (i32_mul(s1, co1) + i32_mul(r1, si1)) << 1;

        let ia3 = 3 * ia1;
        let co3 = twiddles[2 * ia3] as i32;
        let si3 = twiddles[2 * ia3 + 1] as i32;
        data[2 * i3] = (i32_mul(r2, co3) - i32_mul(s2, si3)) << 1;
        data[2 * i3 + 1] = (i32_mul(s2, co3) + i32_mul(r2, si3)) << 1;

        ia1 += twiddle_modifier as usize;
        i0 += 1;
        j -= 1;
    }

    twiddle_modifier <<= 2;

    let mut k = n_fft / 4;
    while k > 4 {
        n1 = n2;
        n2 >>= 2;
        ia1 = 0;

        for j in 0..n2 {
            let ia2 = ia1 + ia1;
            let ia3 = ia2 + ia1;

            let co1 = twiddles[2 * ia1] as i32;
            let si1 = twiddles[2 * ia1 + 1] as i32;
            let co2 = twiddles[2 * ia2] as i32;
            let si2 = twiddles[2 * ia2 + 1] as i32;
            let co3 = twiddles[2 * ia3] as i32;
            let si3 = twiddles[2 * ia3 + 1] as i32;
            ia1 += twiddle_modifier as usize;

            let mut i0 = j;
            while i0 < n_fft {
                let i1 = i0 + n2;
                let i2 = i1 + n2;
                let i3 = i2 + n2;

                let mut r1 = data[2 * i0] + data[2 * i2];
                let mut r2 = data[2 * i0] - data[2 * i2];
                let mut s1 = data[2 * i0 + 1] + data[2 * i2 + 1];
                let mut s2 = data[2 * i0 + 1] - data[2 * i2 + 1];

                let mut t1 = data[2 * i1] + data[2 * i3];
                data[2 * i0] = (r1 + t1) >> 2;
                r1 -= t1;

                let t2 = data[2 * i1 + 1] + data[2 * i3 + 1];
                data[2 * i0 + 1] = (s1 + t2) >> 2;
                s1 -= t2;

                t1 = data[2 * i1 + 1] - data[2 * i3 + 1];
                let t2 = data[2 * i1] - data[2 * i3];

                data[2 * i1] = (i32_mul(r1, co2) - i32_mul(s1, si2)) >> 1;
                data[2 * i1 + 1] = (i32_mul(s1, co2) + i32_mul(r1, si2)) >> 1;

                r1 = r2 - t1;
                r2 += t1;
                s1 = s2 + t2;
                s2 -= t2;

                data[2 * i2] = (i32_mul(r1, co1) - i32_mul(s1, si1)) >> 1;
                data[2 * i2 + 1] = (i32_mul(s1, co1) + i32_mul(r1, si1)) >> 1;

                data[2 * i3] = (i32_mul(r2, co3) - i32_mul(s2, si3)) >> 1;
                data[2 * i3 + 1] = (i32_mul(s2, co3) + i32_mul(r2, si3)) >> 1;

                i0 += n1;
            }
        }

        twiddle_modifier <<= 2;
        k >>= 2;
    }

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
