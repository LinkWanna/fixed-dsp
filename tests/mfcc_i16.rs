mod mfcc_data;

use core::ffi::c_void;
use core::mem::MaybeUninit;

use fixed_dsp::transform::MfccI16;

#[repr(C)]
struct ArmRfftInstanceQ15 {
    fft_len_real: u32,
    ifft_flag_r: u8,
    bit_reverse_flag_r: u8,
    twid_coef_rmodifier: u32,
    p_twiddle_a_real: *const i16,
    p_twiddle_b_real: *const i16,
    p_cfft: *const c_void,
}

#[repr(C)]
struct ArmMfccInstanceQ15 {
    dct_coefs: *const i16,
    filter_coefs: *const i16,
    window_coefs: *const i16,
    filter_pos: *const u32,
    filter_lengths: *const u32,
    fft_len: u32,
    nb_mel_filters: u32,
    nb_dct_outputs: u32,
    rfft: ArmRfftInstanceQ15,
}

unsafe extern "C" {
    fn arm_mfcc_init_q15(
        s: *mut ArmMfccInstanceQ15,
        fft_len: u32,
        nb_mel_filters: u32,
        nb_dct_outputs: u32,
        dct_coefs: *const i16,
        filter_pos: *const u32,
        filter_lengths: *const u32,
        filter_coefs: *const i16,
        window_coefs: *const i16,
    ) -> i32;

    fn arm_mfcc_q15(
        s: *const ArmMfccInstanceQ15,
        p_src: *mut i16,
        p_dst: *mut i16,
        p_tmp: *mut i32,
    ) -> i32;
}

fn leak_q15(data: &[u16]) -> &'static [i16] {
    let converted: Vec<i16> = data.iter().map(|&value| value as i16).collect();
    Box::leak(converted.into_boxed_slice())
}

fn leak_u32(data: &[u16]) -> &'static [u32] {
    let converted: Vec<u32> = data.iter().map(|&value| value as u32).collect();
    Box::leak(converted.into_boxed_slice())
}

fn parse_array(source: &str, start_marker: &str, end_marker: &str) -> Vec<f32> {
    let start = source.find(start_marker).expect("start marker missing") + start_marker.len();
    let end = source[start..]
        .find(end_marker)
        .expect("end marker missing")
        + start;

    source[start..end]
        .split(',')
        .filter_map(|token| token.trim().parse::<f32>().ok())
        .collect()
}

fn q15_to_f32(sample: i16) -> f32 {
    sample as f32 / 32768.0
}

fn run_cmsis_mfcc(
    input: &[i16],
    dct: &[i16],
    filter: &[i16],
    filter_pos: &[u32],
    filter_len: &[u32],
    window: &[i16],
) -> (i32, Vec<i16>) {
    let mut instance = MaybeUninit::<ArmMfccInstanceQ15>::uninit();
    let init_status = unsafe {
        arm_mfcc_init_q15(
            instance.as_mut_ptr(),
            256,
            20,
            13,
            dct.as_ptr(),
            filter_pos.as_ptr(),
            filter_len.as_ptr(),
            filter.as_ptr(),
            window.as_ptr(),
        )
    };

    let mut src = input.to_vec();
    let mut dst = vec![0_i16; 13];
    let mut tmp = vec![0_i32; 512];

    let run_status = unsafe {
        arm_mfcc_q15(
            instance.as_ptr(),
            src.as_mut_ptr(),
            dst.as_mut_ptr(),
            tmp.as_mut_ptr(),
        )
    };

    (init_status | run_status, dst)
}

#[test]
fn mfcc_i16_difftest_against_cmsis() {
    let source = include_str!("../CMSIS-DSP/PythonWrapper/examples/mfccdebugdata.py");
    let debug = parse_array(source, "debug=np.array([", "])");

    assert_eq!(debug.len(), 256);

    let input: Vec<i16> = debug
        .iter()
        .map(|&value| {
            (value * 32768.0)
                .round()
                .clamp(i16::MIN as f32, i16::MAX as f32) as i16
        })
        .collect();

    let dct = &mfcc_data::MCFF_DCT_U16;
    let filter = &mfcc_data::MEL_FILTER_U16;
    let filter_pos_q15 = &mfcc_data::MEL_FILTER_POS_U16;
    let filter_len_q15 = &mfcc_data::MEL_FILTER_LEN_U16;
    let filter_pos = leak_u32(&mfcc_data::MEL_FILTER_POS_U16);
    let filter_len = leak_u32(&mfcc_data::MEL_FILTER_LEN_U16);
    let window = &mfcc_data::MCFF_WINDOW_U16;

    let dct_cmsis = leak_q15(dct);
    let filter_cmsis = leak_q15(filter);
    let window_cmsis = leak_q15(window);

    let mfcc = MfccI16::new(
        256,
        20,
        13,
        dct,
        filter,
        filter_pos_q15,
        filter_len_q15,
        window,
    );
    let mut rust_input = input.clone();
    let mut rust_output = [0_i16; 13];
    let mut rust_tmp = [0_i32; 256];
    mfcc.run(&mut rust_input, &mut rust_output, &mut rust_tmp);

    let (status, cmsis_output) = run_cmsis_mfcc(
        &input,
        dct_cmsis,
        filter_cmsis,
        &filter_pos,
        &filter_len,
        window_cmsis,
    );
    assert_eq!(status, 0, "CMSIS MFCC init/run failed");
    assert_eq!(
        rust_output.as_slice(),
        cmsis_output.as_slice(),
        "Rust MFCC output diverged from CMSIS"
    );
}

#[test]
fn mfcc_i16_matches_cmsis_reference_scale() {
    let source = include_str!("../CMSIS-DSP/PythonWrapper/examples/mfccdebugdata.py");
    let debug = parse_array(source, "debug=np.array([", "])");
    let reference = parse_array(source, "ref = [", "]");

    assert_eq!(debug.len(), 256);
    assert_eq!(reference.len(), 13);

    let input: Vec<i16> = debug
        .iter()
        .map(|&value| {
            (value * 32768.0)
                .round()
                .clamp(i16::MIN as f32, i16::MAX as f32) as i16
        })
        .collect();

    let dct = &mfcc_data::MCFF_DCT_U16;
    let filter = &mfcc_data::MEL_FILTER_U16;
    let filter_pos = leak_u32(&mfcc_data::MEL_FILTER_POS_U16);
    let filter_len = leak_u32(&mfcc_data::MEL_FILTER_LEN_U16);
    let window = &mfcc_data::MCFF_WINDOW_U16;

    let dct_cmsis = leak_q15(dct);
    let filter_cmsis = leak_q15(filter);
    let window_cmsis = leak_q15(window);

    let (status, cmsis_output) = run_cmsis_mfcc(
        &input,
        dct_cmsis,
        filter_cmsis,
        &filter_pos,
        &filter_len,
        window_cmsis,
    );
    assert_eq!(status, 0, "CMSIS MFCC init/run failed");

    for (index, (&actual, &expected)) in cmsis_output.iter().zip(reference.iter()).enumerate() {
        let observed = q15_to_f32(actual) * 256.0;
        assert!(
            (observed - expected).abs() < 0.25,
            "mfcc[{index}] expected {expected}, got {observed}"
        );
    }
}
