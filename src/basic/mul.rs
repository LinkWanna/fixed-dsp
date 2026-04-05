pub fn mul_i16(a: *const i16, b: *const i16, output: *mut i16, size: usize) {
    for i in 0..size {
        unsafe {
            let a_val = *a.add(i) as i32;
            let b_val = *b.add(i) as i32;
            let mul_val = (a_val * b_val) >> 15;
            *output.add(i) = mul_val.clamp(i16::MIN as i32, i16::MAX as i32) as i16; // Saturate to i16 range
        }
    }
}

pub fn mul_i32(a: *const i32, b: *const i32, output: *mut i32, size: usize) {
    for i in 0..size {
        unsafe {
            let a_val = *a.add(i) as i64;
            let b_val = *b.add(i) as i64;
            let mul_val = (a_val * b_val) >> 32;
            *output.add(i) = (mul_val.clamp(i32::MIN as i64, i32::MAX as i64) << 1) as i32; // Saturate to i32 range
        }
    }
}
