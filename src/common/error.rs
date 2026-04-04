pub type Result<T> = core::result::Result<T, Error>;

/// Error type for fixed-dsp.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Error {
	/// Not-a-number / infinity equivalent in fixed-point APIs.
	NanInf,
}
