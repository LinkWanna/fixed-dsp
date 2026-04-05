use crate::matrix::Matrix;

pub fn mat_add_i16(a: Matrix<i16>, b: Matrix<i16>, output: Matrix<i16>) {
    let len = a.rows * a.cols;

    unsafe {
        for i in 0..len {
            let a_val = *a.data.add(i);
            let b_val = *b.data.add(i);
            let sum = a_val.saturating_add(b_val);
            output.data.add(i).write(sum);
        }
    }
}

pub fn mat_add_i32(a: Matrix<i32>, b: Matrix<i32>, output: Matrix<i32>) {
    let len = a.rows * a.cols;

    unsafe {
        for i in 0..len {
            let a_val = *a.data.add(i);
            let b_val = *b.data.add(i);
            let sum = a_val.saturating_add(b_val);
            output.data.add(i).write(sum);
        }
    }
}
