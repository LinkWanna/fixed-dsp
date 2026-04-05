use fixed_dsp::matrix::{mat_vec_mul_i16, mat_vec_mul_i32, Matrix};

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
    fn arm_mat_vec_mult_q15(pSrcMat: *const ArmMatrixQ15, pSrcVec: *const i16, pDstVec: *mut i16);
    fn arm_mat_vec_mult_q31(pSrcMat: *const ArmMatrixQ31, pSrcVec: *const i32, pDstVec: *mut i32);
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

fn sample_i16_matrix(rows: usize, cols: usize) -> Vec<i16> {
    (0..rows * cols)
        .map(|i| (i as i16).wrapping_mul(127).wrapping_sub(4_321))
        .collect()
}

fn sample_i32_matrix(rows: usize, cols: usize) -> Vec<i32> {
    (0..rows * cols)
        .map(|i| (i as i32).wrapping_mul(65_537).wrapping_sub(123_456_789))
        .collect()
}

fn sample_i16_vec(cols: usize) -> Vec<i16> {
    (0..cols)
        .map(|i| (i as i16).wrapping_mul(251).wrapping_add(3_333))
        .collect()
}

fn sample_i32_vec(cols: usize) -> Vec<i32> {
    (0..cols)
        .map(|i| (i as i32).wrapping_mul(1_048_573).wrapping_add(7_777_777))
        .collect()
}

#[test]
fn mat_vec_mul_q15_difftest_against_cmsis() {
    let cases = [(1usize, 1usize), (2, 3), (3, 4), (5, 2), (7, 5)];

    for &(rows, cols) in &cases {
        let mut mat_data = sample_i16_matrix(rows, cols);
        let vec_data = sample_i16_vec(cols);
        let mut out_rust = vec![0_i16; rows];
        let mut out_cmsis = vec![0_i16; rows];

        mat_vec_mul_i16(
            make_matrix_i16(rows, cols, &mut mat_data),
            &vec_data,
            &mut out_rust,
        );

        let src_mat = ArmMatrixQ15 {
            numRows: rows as u16,
            numCols: cols as u16,
            pData: mat_data.as_mut_ptr(),
        };

        unsafe {
            arm_mat_vec_mult_q15(&src_mat, vec_data.as_ptr(), out_cmsis.as_mut_ptr());
        }

        assert_eq!(
            out_rust, out_cmsis,
            "Q15 mat_vec_mul mismatch at shape {}x{}",
            rows, cols
        );
    }
}

#[test]
fn mat_vec_mul_q31_difftest_against_cmsis() {
    let cases = [(1usize, 1usize), (2, 3), (3, 4), (5, 2), (7, 5)];

    for &(rows, cols) in &cases {
        let mut mat_data = sample_i32_matrix(rows, cols);
        let vec_data = sample_i32_vec(cols);
        let mut out_rust = vec![0_i32; rows];
        let mut out_cmsis = vec![0_i32; rows];

        mat_vec_mul_i32(
            make_matrix_i32(rows, cols, &mut mat_data),
            &vec_data,
            &mut out_rust,
        );

        let src_mat = ArmMatrixQ31 {
            numRows: rows as u16,
            numCols: cols as u16,
            pData: mat_data.as_mut_ptr(),
        };

        unsafe {
            arm_mat_vec_mult_q31(&src_mat, vec_data.as_ptr(), out_cmsis.as_mut_ptr());
        }

        assert_eq!(
            out_rust, out_cmsis,
            "Q31 mat_vec_mul mismatch at shape {}x{}",
            rows, cols
        );
    }
}
