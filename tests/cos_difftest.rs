mod difftest;

use difftest::{assert_max_abs_diff, run_i16_difftest, run_i32_difftest};
use fixed_dsp::basic::{cos_i16, cos_i32};

unsafe extern "C" {
    fn arm_cos_q15(x: i16) -> i16;
    fn arm_cos_q31(x: i32) -> i32;
}

fn sample_q31_inputs() -> Vec<i32> {
    let mut data = Vec::with_capacity(20004);
    data.extend([i32::MIN, i32::MIN + 1, -1, 0, 1, i32::MAX - 1, i32::MAX]);

    let mut state: u64 = 0x9E37_79B9_7F4A_7C15;
    for _ in 0..20000 {
        state ^= state << 7;
        state ^= state >> 9;
        state ^= state << 8;
        data.push(state as i32);
    }

    data
}

#[test]
fn cos_q15_difftest_against_cmsis() {
    let inputs: Vec<i16> = (i16::MIN..=i16::MAX).step_by(17).collect();
    let stats = run_i16_difftest(&inputs, cos_i16, |x| unsafe { arm_cos_q15(x) });
    println!("q15 stats: {:?}", stats);

    assert_max_abs_diff(stats, 6, "cos_i16 vs arm_cos_q15");
}

#[test]
fn cos_q31_difftest_against_cmsis() {
    let inputs = sample_q31_inputs();
    let stats = run_i32_difftest(&inputs, cos_i32, |x| unsafe { arm_cos_q31(x) });
    println!("q31 stats: {:?}", stats);

    assert_max_abs_diff(stats, 50_000, "cos_i32 vs arm_cos_q31");
}
