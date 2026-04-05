pub fn bitreversal_i16(data: &mut [i16], bit_rev_table: &'static [u16]) {
    let table = bit_rev_table;

    for i in (0..table.len()).step_by(2) {
        let a = (table[i] as usize) >> 2;
        let b = (table[i + 1] as usize) >> 2;

        assert!(a + 1 < data.len(), "bit reversal index out of range: {}", a);
        assert!(b + 1 < data.len(), "bit reversal index out of range: {}", b);

        data.swap(a, b);
        data.swap(a + 1, b + 1);
    }
}

pub fn bitreversal_i32(data: &mut [i32], bit_rev_table: &'static [u16]) {
    let table = bit_rev_table;

    for i in (0..table.len()).step_by(2) {
        let a = (table[i] as usize) >> 2;
        let b = (table[i + 1] as usize) >> 2;

        assert!(a + 1 < data.len(), "bit reversal index out of range: {}", a);
        assert!(b + 1 < data.len(), "bit reversal index out of range: {}", b);

        data.swap(a, b);
        data.swap(a + 1, b + 1);
    }
}
