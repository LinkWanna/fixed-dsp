use fixed_dsp::matrix::{mat_mul_i16, mat_mul_i32, Matrix};

#[repr(C)]
#[allow(non_snake_case)]
struct ArmMatrixQ15 {
    numRows: u16,
    numCols: u16,
    pData: *mut i16,
}

#[repr(C)]
#[allow(non_snake_case)]
struct ArmMatrixQ31 {
    numRows: u16,
    numCols: u16,
    pData: *mut i32,
}

unsafe extern "C" {
    fn arm_mat_mult_q15(
        pSrcA: *const ArmMatrixQ15,
        pSrcB: *const ArmMatrixQ15,
        pDst: *mut ArmMatrixQ15,
        pState: *mut i16,
    ) -> i32;

    fn arm_mat_mult_q31(
        pSrcA: *const ArmMatrixQ31,
        pSrcB: *const ArmMatrixQ31,
        pDst: *mut ArmMatrixQ31,
    ) -> i32;
}

fn make_matrix_i16(rows: usize, cols: usize, data: &mut [i16]) -> Matrix<i16> {
    Matrix {
        rows,
        cols,
        data: data.as_mut_ptr(),
    }
}

fn make_matrix_i32(rows: usize, cols: usize, data: &mut [i32]) -> Matrix<i32> {
    Matrix {
        rows,
        cols,
        data: data.as_mut_ptr(),
    }
}

fn sample_i16_data(len: usize, seed: i16) -> Vec<i16> {
    (0..len)
        .map(|i| (i as i16).wrapping_mul(173).wrapping_add(seed))
        .collect()
}

fn sample_i32_data(len: usize, seed: i32) -> Vec<i32> {
    (0..len)
        .map(|i| (i as i32).wrapping_mul(1_048_583).wrapping_add(seed))
        .collect()
}

#[test]
fn mat_mul_q15_difftest_against_cmsis() {
    let cases = [(1usize, 1usize, 1usize), (1, 3, 2), (2, 3, 4), (3, 4, 2), (4, 4, 4)];

    for &(rows_a, cols_a, cols_b) in &cases {
        let len_a = rows_a * cols_a;
        let len_b = cols_a * cols_b;
        let len_out = rows_a * cols_b;

        let mut a = sample_i16_data(len_a, -12_345);
        let mut b = sample_i16_data(len_b, 22_222);
        let mut out_rust = vec![0_i16; len_out];
        let mut out_cmsis = vec![0_i16; len_out];
        let mut state = vec![0_i16; 1];

        mat_mul_i16(
            make_matrix_i16(rows_a, cols_a, &mut a),
            make_matrix_i16(cols_a, cols_b, &mut b),
            make_matrix_i16(rows_a, cols_b, &mut out_rust),
        )
        .expect("rust mat_mul_i16 should succeed");

        let src_a = ArmMatrixQ15 {
            numRows: rows_a as u16,
            numCols: cols_a as u16,
            pData: a.as_mut_ptr(),
        };
        let src_b = ArmMatrixQ15 {
            numRows: cols_a as u16,
            numCols: cols_b as u16,
            pData: b.as_mut_ptr(),
        };
        let mut dst = ArmMatrixQ15 {
            numRows: rows_a as u16,
            numCols: cols_b as u16,
            pData: out_cmsis.as_mut_ptr(),
        };

        let status = unsafe { arm_mat_mult_q15(&src_a, &src_b, &mut dst, state.as_mut_ptr()) };
        assert_eq!(status, 0, "CMSIS q15 mat_mul returned non-success status");

        assert_eq!(
            out_rust, out_cmsis,
            "Q15 mat_mul mismatch at shape {}x{} * {}x{}",
            rows_a, cols_a, cols_a, cols_b
        );
    }
}

#[test]
fn mat_mul_q31_difftest_against_cmsis() {
    let cases = [(1usize, 1usize, 1usize), (1, 3, 2), (2, 3, 4), (3, 4, 2), (4, 4, 4)];

    for &(rows_a, cols_a, cols_b) in &cases {
        let len_a = rows_a * cols_a;
        let len_b = cols_a * cols_b;
        let len_out = rows_a * cols_b;

        let mut a = sample_i32_data(len_a, -987_654_321);
        let mut b = sample_i32_data(len_b, 1_234_567_890);
        let mut out_rust = vec![0_i32; len_out];
        let mut out_cmsis = vec![0_i32; len_out];

        mat_mul_i32(
            make_matrix_i32(rows_a, cols_a, &mut a),
            make_matrix_i32(cols_a, cols_b, &mut b),
            make_matrix_i32(rows_a, cols_b, &mut out_rust),
        )
        .expect("rust mat_mul_i32 should succeed");

        let src_a = ArmMatrixQ31 {
            numRows: rows_a as u16,
            numCols: cols_a as u16,
            pData: a.as_mut_ptr(),
        };
        let src_b = ArmMatrixQ31 {
            numRows: cols_a as u16,
            numCols: cols_b as u16,
            pData: b.as_mut_ptr(),
        };
        let mut dst = ArmMatrixQ31 {
            numRows: rows_a as u16,
            numCols: cols_b as u16,
            pData: out_cmsis.as_mut_ptr(),
        };

        let status = unsafe { arm_mat_mult_q31(&src_a, &src_b, &mut dst) };
        assert_eq!(status, 0, "CMSIS q31 mat_mul returned non-success status");

        assert_eq!(
            out_rust, out_cmsis,
            "Q31 mat_mul mismatch at shape {}x{} * {}x{}",
            rows_a, cols_a, cols_a, cols_b
        );
    }
}
