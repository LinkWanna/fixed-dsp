use fixed_dsp::matrix::{mat_add_i16, mat_add_i32, Matrix};

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
    fn arm_mat_add_q15(pSrcA: *const ArmMatrixQ15, pSrcB: *const ArmMatrixQ15, pDst: *mut ArmMatrixQ15) -> i32;
    fn arm_mat_add_q31(pSrcA: *const ArmMatrixQ31, pSrcB: *const ArmMatrixQ31, pDst: *mut ArmMatrixQ31) -> i32;
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

#[test]
fn mat_add_q15_difftest_against_cmsis() {
    let shapes = [(1usize, 1usize), (1, 7), (2, 3), (3, 5), (8, 4)];

    for &(rows, cols) in &shapes {
        let len = rows * cols;

        let mut a: Vec<i16> = (0..len)
            .map(|i| (i as i16).wrapping_mul(257).wrapping_sub(12_345))
            .collect();
        let mut b: Vec<i16> = (0..len)
            .map(|i| (i as i16).wrapping_mul(-193).wrapping_add(22_222))
            .collect();
        let mut out_rust = vec![0_i16; len];
        let mut out_cmsis = vec![0_i16; len];

        mat_add_i16(
            make_matrix_i16(rows, cols, &mut a),
            make_matrix_i16(rows, cols, &mut b),
            make_matrix_i16(rows, cols, &mut out_rust),
        );

        let src_a = ArmMatrixQ15 {
            numRows: rows as u16,
            numCols: cols as u16,
            pData: a.as_mut_ptr(),
        };
        let src_b = ArmMatrixQ15 {
            numRows: rows as u16,
            numCols: cols as u16,
            pData: b.as_mut_ptr(),
        };
        let mut dst = ArmMatrixQ15 {
            numRows: rows as u16,
            numCols: cols as u16,
            pData: out_cmsis.as_mut_ptr(),
        };

        let status = unsafe { arm_mat_add_q15(&src_a, &src_b, &mut dst) };
        assert_eq!(status, 0, "CMSIS q15 mat_add returned non-success status");

        assert_eq!(out_rust, out_cmsis, "Q15 mat_add mismatch at shape {}x{}", rows, cols);
    }
}

#[test]
fn mat_add_q31_difftest_against_cmsis() {
    let shapes = [(1usize, 1usize), (1, 7), (2, 3), (3, 5), (8, 4)];

    for &(rows, cols) in &shapes {
        let len = rows * cols;

        let mut a: Vec<i32> = (0..len)
            .map(|i| (i as i32).wrapping_mul(1_048_583).wrapping_sub(987_654_321))
            .collect();
        let mut b: Vec<i32> = (0..len)
            .map(|i| (i as i32).wrapping_mul(-1_333_777).wrapping_add(1_234_567_890))
            .collect();
        let mut out_rust = vec![0_i32; len];
        let mut out_cmsis = vec![0_i32; len];

        mat_add_i32(
            make_matrix_i32(rows, cols, &mut a),
            make_matrix_i32(rows, cols, &mut b),
            make_matrix_i32(rows, cols, &mut out_rust),
        );

        let src_a = ArmMatrixQ31 {
            numRows: rows as u16,
            numCols: cols as u16,
            pData: a.as_mut_ptr(),
        };
        let src_b = ArmMatrixQ31 {
            numRows: rows as u16,
            numCols: cols as u16,
            pData: b.as_mut_ptr(),
        };
        let mut dst = ArmMatrixQ31 {
            numRows: rows as u16,
            numCols: cols as u16,
            pData: out_cmsis.as_mut_ptr(),
        };

        let status = unsafe { arm_mat_add_q31(&src_a, &src_b, &mut dst) };
        assert_eq!(status, 0, "CMSIS q31 mat_add returned non-success status");

        assert_eq!(out_rust, out_cmsis, "Q31 mat_add mismatch at shape {}x{}", rows, cols);
    }
}
