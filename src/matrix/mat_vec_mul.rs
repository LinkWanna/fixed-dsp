use crate::matrix::Matrix;

#[inline]
fn sat_i16(x: i64) -> i16 {
    x.clamp(i16::MIN as i64, i16::MAX as i64) as i16
}

pub fn mat_vec_mul_i16(mat: Matrix<i16>, vec: &[i16], output: &mut [i16]) {
    let len = mat.rows * mat.cols;

    if vec.len() < mat.cols || output.len() < mat.rows {
        return;
    }

    unsafe {
        let mat_data = core::slice::from_raw_parts(mat.data as *const i16, len);

        for row in 0..mat.rows {
            let mut acc = 0_i64;

            for col in 0..mat.cols {
                let a = mat_data[row * mat.cols + col] as i64;
                let b = vec[col] as i64;
                acc = acc.wrapping_add(a * b);
            }

            output[row] = sat_i16(acc >> 15);
        }
    }
}

pub fn mat_vec_mul_i32(mat: Matrix<i32>, vec: &[i32], output: &mut [i32]) {
    let len = mat.rows * mat.cols;

    if vec.len() < mat.cols || output.len() < mat.rows {
        return;
    }

    unsafe {
        let mat_data = core::slice::from_raw_parts(mat.data as *const i32, len);

        for row in 0..mat.rows {
            let mut acc = 0_i64;

            for col in 0..mat.cols {
                let a = mat_data[row * mat.cols + col] as i64;
                let b = vec[col] as i64;
                acc = acc.wrapping_add(a * b);
            }

            output[row] = (acc >> 31) as i32;
        }
    }
}
