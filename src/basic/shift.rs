pub fn shift_i16(data: &mut [i16], shift: i8) {
    if shift >= 0 {
        let l = shift as u32;
        for x in data.iter_mut() {
            let shifted = if l >= 31 {
                if *x >= 0 { i32::MAX } else { i32::MIN }
            } else {
                (*x as i32) << l
            };
            *x = shifted.clamp(i16::MIN as i32, i16::MAX as i32) as i16;
        }
    } else {
        let r = (-shift) as u32;
        for x in data.iter_mut() {
            *x = if r >= 15 {
                if *x < 0 { -1 } else { 0 }
            } else {
                ((*x as i32) >> r) as i16
            };
        }
    }
}

pub fn shift_i32(data: &mut [i32], shift: i8) {
    if shift >= 0 {
        let l = shift as u32;
        for x in data.iter_mut() {
            let in_val = *x as i64;
            let shifted = if l >= 63 {
                if in_val >= 0 { i64::MAX } else { i64::MIN }
            } else {
                in_val << l
            };
            *x = shifted.clamp(i32::MIN as i64, i32::MAX as i64) as i32;
        }
    } else {
        let r = (-shift) as u32;
        for x in data.iter_mut() {
            *x = if r >= 31 {
                if *x < 0 { -1 } else { 0 }
            } else {
                *x >> r
            };
        }
    }
}
