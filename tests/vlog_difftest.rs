use fixed_dsp::basic::{vlog_i16, vlog_i32};

unsafe extern "C" {
    fn arm_vlog_q15(pSrc: *const i16, pDst: *mut i16, blockSize: u32);
    fn arm_vlog_q31(pSrc: *const i32, pDst: *mut i32, blockSize: u32);
}

fn sample_q15(block_size: usize) -> Vec<i16> {
    let mut v = Vec::with_capacity(block_size);
    let mut state: u32 = 0x1234_5678;
    for _ in 0..block_size {
        state ^= state << 7;
        state ^= state >> 9;
        state ^= state << 8;
        v.push(state as u16 as i16);
    }
    if block_size >= 8 {
        v[0] = 0x0000;
        v[1] = 0x0001;
        v[2] = 0x7FFF;
        v[3] = 0x4000;
        v[4] = 0x2000;
        v[5] = 0x1000;
        v[6] = 0xFFFFu16 as i16;
        v[7] = 0x8000u16 as i16;
    }
    v
}

fn sample_q31(block_size: usize) -> Vec<i32> {
    let mut v = Vec::with_capacity(block_size);
    let mut state: u64 = 0x9E37_79B9_7F4A_7C15;
    for _ in 0..block_size {
        state ^= state << 7;
        state ^= state >> 9;
        state ^= state << 8;
        v.push(state as u32 as i32);
    }
    if block_size >= 8 {
        v[0] = 0x00000000;
        v[1] = 0x00000001;
        v[2] = 0x7FFFFFFF;
        v[3] = 0x40000000;
        v[4] = 0x20000000;
        v[5] = 0x10000000;
        v[6] = 0xFFFFFFFFu32 as i32;
        v[7] = 0x80000000u32 as i32;
    }
    v
}

#[test]
fn vlog_q15_difftest_against_cmsis() {
    for &block_size in &[1usize, 2, 3, 4, 7, 16, 31, 64, 127, 256] {
        let input = sample_q15(block_size);
        let mut rust_out = vec![0i16; block_size];
        let mut cmsis_out = vec![0i16; block_size];

        vlog_i16(&input, &mut rust_out);
        unsafe {
            arm_vlog_q15(input.as_ptr(), cmsis_out.as_mut_ptr(), block_size as u32);
        }

        assert_eq!(
            rust_out, cmsis_out,
            "Q15 vlog mismatch for block_size={}",
            block_size
        );
    }
}

#[test]
fn vlog_q31_difftest_against_cmsis() {
    for &block_size in &[1usize, 2, 3, 4, 7, 16, 31, 64, 127, 256] {
        let input = sample_q31(block_size);
        let mut rust_out = vec![0i32; block_size];
        let mut cmsis_out = vec![0i32; block_size];

        vlog_i32(&input, &mut rust_out);
        unsafe {
            arm_vlog_q31(input.as_ptr(), cmsis_out.as_mut_ptr(), block_size as u32);
        }

        assert_eq!(
            rust_out, cmsis_out,
            "Q31 vlog mismatch for block_size={}",
            block_size
        );
    }
}
