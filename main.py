import numpy as np
import scipy.signal.windows as sw

from script.mfcc import dct, mel_filter, render_mfcc_template


def main():
    # Mel滤波器组
    filt_len, filt_pos, total_len, packed_filters, filters = mel_filter(
        sample_rate=16000,
        n_fft=256,
        f_min=64.0,
        f_max=8000.0,
        n_mels=20,
        htk=True,
    )

    # 离散余弦变换矩阵
    dct_matrix = dct(13, 20)

    # 窗函数
    window = sw.hamming(256, sym=False)

    rust_code = render_mfcc_template(
        ty="i16",
        window=window,
        dct_matrix=dct_matrix,
        filt_pos=filt_pos,
        filt_len=filt_len,
        packed_filters=packed_filters,
    )
    print(rust_code)


if __name__ == "__main__":
    main()
