pub fn offset_i16(data: &mut [i16], offset: i16) {
    for sample in data.iter_mut() {
        *sample = sample.saturating_add(offset);
    }
}

pub fn offset_i32(data: &mut [i32], offset: i32) {
    for sample in data.iter_mut() {
        *sample = sample.saturating_add(offset);
    }
}
