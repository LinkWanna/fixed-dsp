#[derive(Debug, Clone, Copy)]
pub struct DiffStats {
    pub total: usize,
    pub max_abs_diff: i64,
    pub worst_input: i64,
    pub rust_output: i64,
    pub reference_output: i64,
}

pub fn run_i16_difftest<FDut, FRef>(inputs: &[i16], dut: FDut, reference: FRef) -> DiffStats
where
    FDut: Fn(i16) -> i16,
    FRef: Fn(i16) -> i16,
{
    let mut stats = DiffStats {
        total: inputs.len(),
        max_abs_diff: -1,
        worst_input: 0,
        rust_output: 0,
        reference_output: 0,
    };

    for &x in inputs {
        let y_rust = dut(x);
        let y_ref = reference(x);
        let diff = (y_rust as i64 - y_ref as i64).abs();

        if diff > stats.max_abs_diff {
            stats.max_abs_diff = diff;
            stats.worst_input = x as i64;
            stats.rust_output = y_rust as i64;
            stats.reference_output = y_ref as i64;
        }
    }

    stats
}

pub fn run_i32_difftest<FDut, FRef>(inputs: &[i32], dut: FDut, reference: FRef) -> DiffStats
where
    FDut: Fn(i32) -> i32,
    FRef: Fn(i32) -> i32,
{
    let mut stats = DiffStats {
        total: inputs.len(),
        max_abs_diff: -1,
        worst_input: 0,
        rust_output: 0,
        reference_output: 0,
    };

    for &x in inputs {
        let y_rust = dut(x);
        let y_ref = reference(x);
        let diff = (y_rust as i64 - y_ref as i64).abs();

        if diff > stats.max_abs_diff {
            stats.max_abs_diff = diff;
            stats.worst_input = x as i64;
            stats.rust_output = y_rust as i64;
            stats.reference_output = y_ref as i64;
        }
    }

    stats
}

pub fn assert_max_abs_diff(stats: DiffStats, limit: i64, label: &str) {
    assert!(
        stats.max_abs_diff <= limit,
        "{} diff failed: max_abs_diff={} (limit={}), worst_input={}, rust_output={}, reference_output={}, total={}",
        label,
        stats.max_abs_diff,
        limit,
        stats.worst_input,
        stats.rust_output,
        stats.reference_output,
        stats.total
    );
}
