# fixed-dsp

`fixed-dsp` is a Rust fixed-point DSP library focused on Q15/Q31 algorithms with behavior close to CMSIS-DSP. It also supports differential testing (difftest) against CMSIS C reference implementations.

[中文](https://github.com/LinkWanna/fixed-dsp/blob/main/README_zh.md)

## Highlights

- `#![no_std]` friendly for embedded targets
- Q15/Q31 implementations (`i16` / `i32`)
- Modular structure: basic ops, matrix, statistics, complex math, transforms, MFCC
- Optional `difftest` feature: when enabled, compiles selected CMSIS-DSP C sources for differential validation

## Implemented Modules

### Basic

Located in [src/basic](src/basic).

- Trigonometric functions: `sin_i16` / `sin_i32`, `cos_i16` / `cos_i32`
- Multiply and vector ops: `mul_i16` / `mul_i32`, `dot_i16` / `dot_i32`
- Scale and shift: `scale_i16` / `scale_i32`, `shift_i16` / `shift_i32`
- Offset and log: `offset_i16` / `offset_i32`, `vlog_i16` / `vlog_i32`
- Sqrt and division: `sqrt_i16` / `sqrt_i32`, `div_i16` / `div_i32`

### Matrix

Located in [src/matrix](src/matrix).

- `Matrix<T>` view struct
- `mat_add_i16` / `mat_add_i32`
- `mat_mul_i16` / `mat_mul_i32`
- `mat_vec_mul_i16` / `mat_vec_mul_i32`

### Statistics

Located in [src/statistics](src/statistics).

- `absmax_i16` / `absmax_i32`

### Complex

Located in [src/complex](src/complex).

- `cmplx_mag_i16` / `cmplx_mag_i32`

### Transform

Located in [src/transform](src/transform).

- Bit-reversal: `bitreversal_i16` / `bitreversal_i32`
- CFFT: `CfftI16` / `CfftI32`
- RFFT: `RfftI16` / `RfftI32`
- Radix-4 butterfly: forward/inverse Q15/Q31 kernels
- MFCC: `MfccI16` / `MfccI32`

## Quick Start

### 1. Add dependency

For local workspace development, use a path dependency:

```toml
[dependencies]
fixed-dsp = "0.2"
```

### 2. Basic example (Q15 sine)

```rust
use fixed_dsp::basic::sin_i16;

fn main() {
	// In Q15, 0x4000 is approximately pi/2.
	let x: i16 = 0x4000;
	let y = sin_i16(x);
	println!("sin_q15 = {}", y);
}
```

### 3. CFFT example (Q15)

```rust
use fixed_dsp::transform::CfftI16;

fn main() {
	let fft = CfftI16::new(16, false, true);

	// Interleaved complex format: [re0, im0, re1, im1, ...]
	let mut data = [0i16; 32];
	data[0] = 0x4000; // simple impulse

	fft.run(&mut data);
}
```

## Build and Test

With the `difftest` feature enabled, [build.rs](build.rs) compiles referenced C files under [CMSIS-DSP/Source](CMSIS-DSP/Source) and links them into tests.

```bash
cargo test --lib --tests --features difftest
```

Notes:

- Requires a working host C compiler (`cc`/`clang`)
- Requires a complete [CMSIS-DSP](CMSIS-DSP) directory
- Without `difftest`, integration tests that rely on CMSIS extern symbols are expected to fail at link time

## Design Goals

- Provide readable and maintainable fixed-point DSP implementations in Rust
- Continuously validate numerical behavior against CMSIS-DSP through difftest
- Serve as a stable base for embedded deployment and algorithm iteration
