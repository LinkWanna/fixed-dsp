mod mat_add;

pub use mat_add::*;

pub struct Matrix<T> {
    pub rows: usize,
    pub cols: usize,
    pub data: *mut T,
}
