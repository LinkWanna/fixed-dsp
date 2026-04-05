use crate::common::error::Error;
use crate::common::error::Result;
use crate::matrix::Matrix;

pub fn mat_mul_i16(a: Matrix<i16>, b: Matrix<i16>, output: Matrix<i16>) -> Result<()> {
    if a.cols != b.rows || a.rows != output.rows || b.cols != output.cols {
        return Err(Error::SizeMismatch);
    }

    let rows = a.rows;
    let cols = a.cols;
    let out_cols = output.cols;

    unsafe {
        for row in 0..rows {
            for col in 0..b.cols {
                let mut acc = 0_i64;

                for k in 0..cols {
                    let a_val = *a.data.add(row * cols + k) as i64;
                    let b_val = *b.data.add(k * b.cols + col) as i64;
                    acc = acc.wrapping_add(a_val * b_val);
                }

                output
                    .data
                    .add(row * out_cols + col)
                    .write(((acc >> 15) as i32).clamp(i16::MIN as i32, i16::MAX as i32) as i16);
            }
        }
    }

    Ok(())
}

pub fn mat_mul_i32(a: Matrix<i32>, b: Matrix<i32>, output: Matrix<i32>) -> Result<()> {
    if a.cols != b.rows || a.rows != output.rows || b.cols != output.cols {
        return Err(Error::SizeMismatch);
    }

    let rows = a.rows;
    let cols = a.cols;
    let out_cols = output.cols;

    unsafe {
        for row in 0..rows {
            for col in 0..b.cols {
                let mut acc = 0_i64;

                for k in 0..cols {
                    let a_val = *a.data.add(row * cols + k) as i64;
                    let b_val = *b.data.add(k * b.cols + col) as i64;
                    acc = acc.wrapping_add(a_val * b_val);
                }

                output
                    .data
                    .add(row * out_cols + col)
                    .write((acc >> 31) as i32);
            }
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mat_mul_reports_size_mismatch() {
        let mut a = [0_i16; 4];
        let mut b = [0_i16; 4];
        let mut out = [0_i16; 4];

        let a = Matrix {
            rows: 2,
            cols: 2,
            data: a.as_mut_ptr(),
        };
        let b = Matrix {
            rows: 1,
            cols: 4,
            data: b.as_mut_ptr(),
        };
        let out = Matrix {
            rows: 2,
            cols: 2,
            data: out.as_mut_ptr(),
        };

        assert_eq!(mat_mul_i16(a, b, out), Err(Error::SizeMismatch));
    }
}
