use crate::matrix::Matrix;

pub fn mat_add_i16(a: Matrix<i16>, b: Matrix<i16>, output: Matrix<i16>) {
    let len = a.rows * a.cols;

    unsafe {
        let in_a = core::slice::from_raw_parts(a.data as *const i16, len);
        let in_b = core::slice::from_raw_parts(b.data as *const i16, len);
        let out = core::slice::from_raw_parts_mut(output.data, len);

        for i in 0..len {
            out[i] = in_a[i].saturating_add(in_b[i]);
        }
    }
}

pub fn mat_add_i32(a: Matrix<i32>, b: Matrix<i32>, output: Matrix<i32>) {
    let len = a.rows * a.cols;

    unsafe {
        let in_a = core::slice::from_raw_parts(a.data as *const i32, len);
        let in_b = core::slice::from_raw_parts(b.data as *const i32, len);
        let out = core::slice::from_raw_parts_mut(output.data, len);

        for i in 0..len {
            out[i] = in_a[i].saturating_add(in_b[i]);
        }
    }
}
