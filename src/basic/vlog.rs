const LOG_Q15_ACCURACY: i32 = 15;
const LOG_Q15_INTEGER_PART: i32 = 4;
const LOQ_Q15_THRESHOLD: u32 = 1u32 << LOG_Q15_ACCURACY;
const LOQ_Q15_Q16_HALF: u32 = LOQ_Q15_THRESHOLD;
const LOG_Q15_INVLOG2EXP: i32 = 0x58b9; // 1.0 / ln(2) in Q15 format

const LOG_Q31_ACCURACY: i32 = 31;
const LOG_Q31_INTEGER_PART: i32 = 5;
const LOQ_Q31_THRESHOLD: u64 = 1u64 << LOG_Q31_ACCURACY;
const LOQ_Q31_Q32_HALF: u64 = LOQ_Q31_THRESHOLD;
const LOG_Q31_INVLOG2EXP: i64 = 0x58b90bfb; // 1.0 / ln(2) in Q31 format

fn log_i16(x: i16) -> i16 {
    let src = x as u16 as u32;

    let c = src.leading_zeros() as i32 - 16;
    let normalization = c;

    let mut inc: u32 = LOQ_Q15_Q16_HALF >> (LOG_Q15_INTEGER_PART + 1);

    let mut xn = src;
    if (c - 1) < 0 {
        xn >>= 1 - c;
    } else {
        xn <<= c - 1;
    }

    let mut y: u32 = 0;

    for _ in 0..LOG_Q15_ACCURACY {
        xn = ((xn as i32 * xn as i32) >> (LOG_Q15_ACCURACY - 1)) as u32;

        if xn >= LOQ_Q15_THRESHOLD {
            y = y.wrapping_add(inc);
            xn >>= 1;
        }
        inc >>= 1;
    }

    let tmp = y as i32 - (normalization << (LOG_Q15_ACCURACY - LOG_Q15_INTEGER_PART));
    (((tmp as i64 * LOG_Q15_INVLOG2EXP as i64) >> 15) as i32) as i16
}

pub fn vlog_i16(input: &[i16], output: &mut [i16]) {
    assert_eq!(
        input.len(),
        output.len(),
        "input and output lengths must match"
    );
    for (x, y) in input.iter().zip(output.iter_mut()) {
        *y = log_i16(*x);
    }
}

fn log_i32(x: i32) -> i32 {
    let src = x as u32;

    let c = src.leading_zeros() as i32;
    let normalization = c;

    let mut inc: u64 = LOQ_Q31_Q32_HALF >> (LOG_Q31_INTEGER_PART + 1);

    let mut xn = src as u64;
    if (c - 1) < 0 {
        xn >>= 1 - c;
    } else {
        xn <<= c - 1;
    }

    let mut y: u64 = 0;

    for _ in 0..LOG_Q31_ACCURACY {
        xn = ((xn * xn) >> (LOG_Q31_ACCURACY - 1)) as u64;

        if xn >= LOQ_Q31_THRESHOLD {
            y = y.wrapping_add(inc);
            xn >>= 1;
        }
        inc >>= 1;
    }

    let tmp = y as i64 - ((normalization as i64) << (LOG_Q31_ACCURACY - LOG_Q31_INTEGER_PART));
    ((tmp * LOG_Q31_INVLOG2EXP) >> 31) as i32
}

pub fn vlog_i32(input: &[i32], output: &mut [i32]) {
    assert_eq!(
        input.len(),
        output.len(),
        "input and output lengths must match"
    );
    for (x, y) in input.iter().zip(output.iter_mut()) {
        *y = log_i32(*x);
    }
}
