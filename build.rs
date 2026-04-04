fn main() {
    println!("cargo:rerun-if-changed=CMSIS-DSP/Source/FastMathFunctions/arm_sin_q15.c");
    println!("cargo:rerun-if-changed=CMSIS-DSP/Source/FastMathFunctions/arm_sin_q31.c");
    println!("cargo:rerun-if-changed=CMSIS-DSP/Source/FastMathFunctions/arm_cos_q15.c");
    println!("cargo:rerun-if-changed=CMSIS-DSP/Source/FastMathFunctions/arm_cos_q31.c");
    println!("cargo:rerun-if-changed=CMSIS-DSP/Source/CommonTables/arm_common_tables.c");
    println!("cargo:rerun-if-changed=CMSIS-DSP/Include");
    println!("cargo:rerun-if-changed=CMSIS-DSP/PrivateInclude");

    cc::Build::new()
        .define("__GNUC_PYTHON__", Some("1"))
        .include("CMSIS-DSP/Include")
        .include("CMSIS-DSP/PrivateInclude")
        .file("CMSIS-DSP/Source/FastMathFunctions/arm_sin_q15.c")
        .file("CMSIS-DSP/Source/FastMathFunctions/arm_sin_q31.c")
        .file("CMSIS-DSP/Source/FastMathFunctions/arm_cos_q15.c")
        .file("CMSIS-DSP/Source/FastMathFunctions/arm_cos_q31.c")
        .file("CMSIS-DSP/Source/CommonTables/arm_common_tables.c")
        .flag_if_supported("-Wno-unused-parameter")
        .compile("cmsis_dsp_ref");
}
