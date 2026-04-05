use crate::common::tables::{
    BIT_REV_TABLE_16, BIT_REV_TABLE_32, BIT_REV_TABLE_64, BIT_REV_TABLE_128, BIT_REV_TABLE_256,
    BIT_REV_TABLE_512, BIT_REV_TABLE_1024, BIT_REV_TABLE_2048, BIT_REV_TABLE_4096,
};

fn bit_rev_table(fft_len: usize) -> &'static [u16] {
    match fft_len {
        16 => &BIT_REV_TABLE_16,
        32 => &BIT_REV_TABLE_32,
        64 => &BIT_REV_TABLE_64,
        128 => &BIT_REV_TABLE_128,
        256 => &BIT_REV_TABLE_256,
        512 => &BIT_REV_TABLE_512,
        1024 => &BIT_REV_TABLE_1024,
        2048 => &BIT_REV_TABLE_2048,
        4096 => &BIT_REV_TABLE_4096,
        _ => panic!(
            "unsupported fft_len: {} (expected power-of-two in [16, 4096])",
            fft_len
        ),
    }
}

pub fn bitreversal_i16(data: &mut [i16], fft_len: usize) {
    assert_eq!(
        data.len(),
        fft_len * 2,
        "Q15 buffer length must be 2 * fft_len"
    );
    assert!(fft_len.is_power_of_two(), "fft_len must be a power of two");
    let table = bit_rev_table(fft_len);
    assert!(
        table.len().is_multiple_of(2),
        "bit reversal table length must be even"
    );

    for i in (0..table.len()).step_by(2) {
        let a = (table[i] as usize) >> 2;
        let b = (table[i + 1] as usize) >> 2;

        assert!(a + 1 < data.len(), "bit reversal index out of range: {}", a);
        assert!(b + 1 < data.len(), "bit reversal index out of range: {}", b);

        data.swap(a, b);
        data.swap(a + 1, b + 1);
    }
}

pub fn bitreversal_i32(data: &mut [i32], fft_len: usize) {
    assert_eq!(
        data.len(),
        fft_len * 2,
        "Q31 buffer length must be 2 * fft_len"
    );
    assert!(fft_len.is_power_of_two(), "fft_len must be a power of two");
    let table = bit_rev_table(fft_len);
    assert!(
        table.len().is_multiple_of(2),
        "bit reversal table length must be even"
    );

    for i in (0..table.len()).step_by(2) {
        let a = (table[i] as usize) >> 2;
        let b = (table[i + 1] as usize) >> 2;

        assert!(a + 1 < data.len(), "bit reversal index out of range: {}", a);
        assert!(b + 1 < data.len(), "bit reversal index out of range: {}", b);

        data.swap(a, b);
        data.swap(a + 1, b + 1);
    }
}
