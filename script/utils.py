import numpy as np


def to_i16(value, scale: bool) -> str:
    if scale:
        scaled = int(np.round(value * 2**15))
        if scaled > 0x7FFF:
            scaled = 0x7FFF
        if scaled < -0x8000:
            scaled = -0x8000
        value = scaled

    # convert to unsigned
    if value < 0:
        value = value + 0x10000
    return f"0x{value:04X}"


def to_i32(value, scale: bool) -> str:
    if scale:
        scaled = int(np.round(value * 2**31))
        if scaled > 0x7FFFFFFF:
            scaled = 0x7FFFFFFF
        if scaled < -0x80000000:
            scaled = -0x80000000
        value = scaled

    # convert to unsigned
    if value < 0:
        value = value + 0x100000000
    return f"0x{value:08X}"


def format_array(values, ty: str, scale: bool) -> str:
    def cvt(x) -> str:
        if ty == "i16":
            return to_i16(x, scale)
        elif ty == "i32":
            return to_i32(x, scale)
        else:
            raise ValueError(f"Unsupported source type: {ty}")

    if ty == "i16":
        step = 12
    elif ty == "i32":
        step = 8

    num = 0
    content = ""
    for value in values:
        content += cvt(value)
        content += ", "
        num = num + 1
        if num == step:
            content += "\n\t"
            num = 0
    return content
