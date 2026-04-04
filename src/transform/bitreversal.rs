use crate::common::tables::BIT_REV_TABLE;

fn bit_reversal_table(fft_len: usize) -> (&'static [u16], u16) {
    match fft_len {
        16 => (&BIT_REV_TABLE[255..], 256),
        32 => (&BIT_REV_TABLE[127..], 128),
        64 => (&BIT_REV_TABLE[63..], 64),
        128 => (&BIT_REV_TABLE[31..], 32),
        256 => (&BIT_REV_TABLE[15..], 16),
        512 => (&BIT_REV_TABLE[7..], 8),
        1024 => (&BIT_REV_TABLE[3..], 4),
        2048 => (&BIT_REV_TABLE[1..], 2),
        4096 => (&BIT_REV_TABLE[0..], 1),
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

    let (bit_rev_tab, bit_rev_factor) = bit_reversal_table(fft_len);
    let fft_len_by2 = fft_len / 2;
    let fft_len_by2p1 = fft_len_by2 + 1;

    let mut j: usize = 0;
    let mut bit_rev_index: usize = 0;

    for i in (0..=(fft_len_by2 - 2)).step_by(2) {
        if i < j {
            data.swap(2 * i, 2 * j);
            data.swap(2 * i + 1, 2 * j + 1);

            data.swap(2 * (i + fft_len_by2p1), 2 * (j + fft_len_by2p1));
            data.swap(2 * (i + fft_len_by2p1) + 1, 2 * (j + fft_len_by2p1) + 1);
        }

        data.swap(2 * (i + 1), 2 * (j + fft_len_by2));
        data.swap(2 * (i + 1) + 1, 2 * (j + fft_len_by2) + 1);

        j = bit_rev_tab[bit_rev_index] as usize;
        bit_rev_index += bit_rev_factor as usize;
    }
}

pub fn bitreversal_i32(data: &mut [i32], fft_len: usize) {
    assert_eq!(
        data.len(),
        fft_len * 2,
        "Q31 buffer length must be 2 * fft_len"
    );
    assert!(fft_len.is_power_of_two(), "fft_len must be a power of two");

    let (bit_rev_tab, bit_rev_factor) = bit_reversal_table(fft_len);
    let fft_len_by2 = fft_len / 2;
    let fft_len_by2p1 = fft_len_by2 + 1;

    let mut j: usize = 0;
    let mut bit_rev_index: usize = 0;

    for i in (0..=(fft_len_by2 - 2)).step_by(2) {
        if i < j {
            data.swap(2 * i, 2 * j);
            data.swap(2 * i + 1, 2 * j + 1);

            data.swap(2 * (i + fft_len_by2p1), 2 * (j + fft_len_by2p1));
            data.swap(2 * (i + fft_len_by2p1) + 1, 2 * (j + fft_len_by2p1) + 1);
        }

        data.swap(2 * (i + 1), 2 * (j + fft_len_by2));
        data.swap(2 * (i + 1) + 1, 2 * (j + fft_len_by2) + 1);

        j = bit_rev_tab[bit_rev_index] as usize;
        bit_rev_index += bit_rev_factor as usize;
    }
}
