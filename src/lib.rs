#![no_std]

pub mod basic;
pub mod common;
pub mod complex;
pub mod matrix;
pub mod statistics;
pub mod transform;

#[inline]
pub(crate) fn sat_i16(x: i32) -> i16 {
    x.clamp(i16::MIN as i32, i16::MAX as i32) as i16
}

#[inline]
pub(crate) fn sat_i32(x: i64) -> i32 {
    x.clamp(i32::MIN as i64, i32::MAX as i64) as i32
}
