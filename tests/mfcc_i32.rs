mod mfcc_data;

use core::ffi::c_void;

use fixed_dsp::transform::MfccI32;

unsafe extern "C" {
    fn arm_mfcc_init_q31(
        s: *mut c_void,
        fft_len: u32,
        nb_mel_filters: u32,
        nb_dct_outputs: u32,
        dct_coefs: *const i32,
        filter_pos: *const u32,
        filter_lengths: *const u32,
        filter_coefs: *const i32,
        window_coefs: *const i32,
    ) -> i32;

    fn arm_mfcc_q31(s: *const c_void, p_src: *mut i32, p_dst: *mut i32, p_tmp: *mut i32) -> i32;
}

fn leak_q31(data: &[u32]) -> &'static [i32] {
    let converted: Vec<i32> = data.iter().map(|&value| value as i32).collect();
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

fn q23_to_f32(sample: i32) -> f32 {
    sample as f32 / ((1u32 << 23) as f32)
}

fn run_cmsis_mfcc(
    input: &[i32],
    dct: &[i32],
    filter: &[i32],
    filter_pos: &[u16],
    filter_len: &[u16],
    window: &[i32],
) -> (i32, Vec<i32>) {
    let mut instance_words = [0u64; 128];
    // 构造 i32 版本的 filter_pos 和 filter_len
    let filter_pos_u32: Vec<u32> = filter_pos.iter().map(|&value| value as u32).collect();
    let filter_len_u32: Vec<u32> = filter_len.iter().map(|&value| value as u32).collect();

    let init_status = unsafe {
        arm_mfcc_init_q31(
            instance_words.as_mut_ptr() as *mut c_void,
            256,
            20,
            13,
            dct.as_ptr(),
            filter_pos_u32.as_ptr(),
            filter_len_u32.as_ptr(),
            filter.as_ptr(),
            window.as_ptr(),
        )
    };

    let mut src = input.to_vec();
    let mut dst = vec![0_i32; 13];
    let mut tmp = vec![0_i32; 512];

    let run_status = unsafe {
        arm_mfcc_q31(
            instance_words.as_ptr() as *const c_void,
            src.as_mut_ptr(),
            dst.as_mut_ptr(),
            tmp.as_mut_ptr(),
        )
    };

    (init_status | run_status, dst)
}

#[test]
fn mfcc_i32_difftest_against_cmsis() {
    let source = include_str!("../CMSIS-DSP/PythonWrapper/examples/mfccdebugdata.py");
    let debug = parse_array(source, "debug=np.array([", "])");

    assert_eq!(debug.len(), 256);

    let input: Vec<i32> = debug
        .iter()
        .map(|&value| {
            (value as f64 * 2147483648.0)
                .round()
                .clamp(i32::MIN as f64, i32::MAX as f64) as i32
        })
        .collect();

    let dct = &mfcc_data::MCFF_DCT_U32;
    let filter = &mfcc_data::MEL_FILTER_U32;
    let filter_pos = &mfcc_data::MEL_FILTER_POS_U16; // Changed to U32
    let filter_len = &mfcc_data::MEL_FILTER_LEN_U16; // Changed to U32
    let window = &mfcc_data::MCFF_WINDOW_U32;

    let dct_cmsis = leak_q31(dct);
    let filter_cmsis = leak_q31(filter);
    let window_cmsis = leak_q31(window);

    let mfcc = MfccI32::new(256, 20, 13, dct, filter, filter_pos, filter_len, window);
    let mut rust_input = input.clone();
    let mut rust_output = [0_i32; 13];
    let mut rust_tmp = [0_i32; 512];
    mfcc.run(&mut rust_input, &mut rust_output, &mut rust_tmp);

    let (status, cmsis_output) = run_cmsis_mfcc(
        &input,
        dct_cmsis,
        filter_cmsis,
        filter_pos,
        filter_len,
        window_cmsis,
    );
    assert_eq!(status, 0, "CMSIS MFCC init/run failed");
    assert_eq!(
        rust_output.as_slice(),
        cmsis_output.as_slice(),
        "Rust MFCC i32 output diverged from CMSIS"
    );
}

#[test]
fn mfcc_i32_matches_cmsis_reference_scale() {
    let source = include_str!("../CMSIS-DSP/PythonWrapper/examples/mfccdebugdata.py");
    let debug = parse_array(source, "debug=np.array([", "])");
    let reference = parse_array(source, "ref = [", "]");

    assert_eq!(debug.len(), 256);
    assert_eq!(reference.len(), 13);

    let input: Vec<i32> = debug
        .iter()
        .map(|&value| {
            (value as f64 * 2147483648.0)
                .round()
                .clamp(i32::MIN as f64, i32::MAX as f64) as i32
        })
        .collect();

    let dct = &mfcc_data::MCFF_DCT_U32;
    let filter = &mfcc_data::MEL_FILTER_U32;
    let filter_pos = &mfcc_data::MEL_FILTER_POS_U16;
    let filter_len = &mfcc_data::MEL_FILTER_LEN_U16;
    let window = &mfcc_data::MCFF_WINDOW_U32;

    let dct_cmsis = leak_q31(dct);
    let filter_cmsis = leak_q31(filter);
    let window_cmsis = leak_q31(window);

    let (status, cmsis_output) = run_cmsis_mfcc(
        &input,
        dct_cmsis,
        filter_cmsis,
        filter_pos,
        filter_len,
        window_cmsis,
    );
    assert_eq!(status, 0, "CMSIS MFCC init/run failed");

    for (index, (&actual, &expected)) in cmsis_output.iter().zip(reference.iter()).enumerate() {
        let observed = q23_to_f32(actual);
        assert!(
            (observed - expected).abs() < 0.02,
            "mfcc[{index}] expected {expected}, got {observed}"
        );
    }
}
