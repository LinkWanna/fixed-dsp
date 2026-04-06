# fixed-dsp

Rust 固定点 DSP 库，目标是提供一组与 CMSIS-DSP 行为接近的 Q15/Q31 算法实现，并支持与 CMSIS C 参考实现进行差分测试（difftest）。

当前 crate 名称是 `fixed-dsp`，代码入口在 [src/lib.rs](src/lib.rs)。

## 特性概览

- `#![no_std]`：适合嵌入式环境，无需标准库
- Q15/Q31 双精度族实现（`i16`/`i32`）
- 模块化结构：基础算子、矩阵、统计、复数运算、变换、MFCC
- 可选 `difftest` 特性：启用后通过 [build.rs](build.rs) 构建 CMSIS-DSP C 参考函数，用于差分验证

## 已实现模块

### Basic

位于 [src/basic](src/basic)。

- 三角函数：`sin_i16` / `sin_i32`，`cos_i16` / `cos_i32`
- 乘法与向量运算：`mul_i16` / `mul_i32`，`dot_i16` / `dot_i32`
- 缩放与位移：`scale_i16` / `scale_i32`，`shift_i16` / `shift_i32`
- 偏置与对数：`offset_i16` / `offset_i32`，`vlog_i16` / `vlog_i32`
- 开方与除法：`sqrt_i16` / `sqrt_i32`，`div_i16` / `div_i32`

### Matrix

位于 [src/matrix](src/matrix)。

- `Matrix<T>` 视图结构
- `mat_add_i16` / `mat_add_i32`
- `mat_mul_i16` / `mat_mul_i32`
- `mat_vec_mul_i16` / `mat_vec_mul_i32`

### Statistics

位于 [src/statistics](src/statistics)。

- `absmax_i16` / `absmax_i32`

### Complex

位于 [src/complex](src/complex)。

- `cmplx_mag_i16` / `cmplx_mag_i32`

### Transform

位于 [src/transform](src/transform)。

- 位反转：`bitreversal_i16` / `bitreversal_i32`
- CFFT：`CfftI16` / `CfftI32`
- RFFT：`RfftI16` / `RfftI32`
- Radix-4 Butterfly：Q15/Q31 正反变换实现
- MFCC：`MfccI16` / `MfccI32`

### Common

位于 [src/common](src/common)。

- `Error` 与 `Result` 错误类型
- 各类查表常量（sin/twiddle/bit-reversal/real-FFT 系数）

## 快速开始

### 1. 添加依赖

如果你在同一工作区开发，可先用路径依赖：

```toml
[dependencies]
fixed-dsp = { path = "../fixed-dsp" }
```

### 2. 基础示例（Q15 正弦）

```rust
use fixed_dsp::basic::sin_i16;

fn main() {
	// Q15 中 0x4000 近似对应 pi/2
	let x: i16 = 0x4000;
	let y = sin_i16(x);
	println!("sin_q15 = {}", y);
}
```

### 3. CFFT 示例（Q15）

```rust
use fixed_dsp::transform::CfftI16;

fn main() {
	let fft = CfftI16::new(16, false, true);

	// 复数交错格式: [re0, im0, re1, im1, ...]
	let mut data = [0i16; 32];
	data[0] = 0x4000; // 一个简单脉冲

	fft.run(&mut data);
}
```

## 构建与测试

### 普通构建

```bash
cargo check
```

### 不依赖 CMSIS C 参考库的测试

```bash
cargo test --lib
```

### 运行 difftest（与 CMSIS-DSP 对比）

启用 `difftest` 特性后，`build.rs` 会编译 `CMSIS-DSP/Source` 下被引用的 C 文件并链接测试。

```bash
cargo test --lib --tests --features difftest
```

说明：

- 该模式需要可用的本地 C 编译器（如 `cc`/`clang`）
- 需要仓库中的 `CMSIS-DSP` 子目录保持完整
- 如果不启用 `difftest`，依赖 CMSIS extern 符号的集成测试会链接失败，这是预期行为

## 仓库结构

- 核心实现： [src](src)
- 差分测试： [tests](tests)
- CMSIS 参考实现与数据： [CMSIS-DSP](CMSIS-DSP)
- 构建脚本： [build.rs](build.rs)
- 包配置： [Cargo.toml](Cargo.toml)

## 设计目标

- 在 Rust 中提供可读、可维护的固定点 DSP 实现
- 用 difftest 持续校验与 CMSIS-DSP 的数值行为一致性
- 为后续嵌入式部署和算法迭代提供稳定基础
