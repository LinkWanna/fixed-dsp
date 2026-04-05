mod mat_add;
mod mat_mul;

pub use mat_add::*;
pub use mat_mul::*;

pub struct Matrix<T> {
    pub rows: usize,
    pub cols: usize,
    pub data: *mut T,
}
