import numpy as np

from script.utils import format_array

template = """
pub const MCFF_WINDOW_{TY}: [{ty}; {window_len}] = [
    {window}
];

pub const MCFF_DCT_{TY}: [{ty}; {dct_len}] = [
    {dct}
];

pub const MEL_FILTER_POS_U16: [u16; {n_mels}] = [
    {filter_pos}
];

pub const MEL_FILTER_LEN_U16: [u16; {n_mels}] = [
    {filter_len_array}
];

pub const MEL_FILTER_{TY}: [{ty}; {total_filter_len}] = [
    {filter}
];
"""


def freq_to_mel(freq, htk: bool = True):
    if htk:
        return 1127.0 * np.log(1.0 + freq / 700.0)
    else:
        freq_step = 200.0 / 3
        min_hz = 1000.0
        log_step = 1.8562979903656 / 27.0  # log(6.4)
        min_mel = min_hz / freq_step
        mel = freq / freq_step
        if freq >= min_hz:
            mel = min_mel + np.log(freq / min_hz) / log_step
        return mel


def mel_to_freq(mels, htk: bool = True):
    if htk:
        return 700.0 * (np.exp(mels / 1127.0) - 1.0)
    else:
        freq_step = 200.0 / 3
        min_hz = 1000.0
        log_step = 1.8562979903656 / 27.0  # log(6.4)
        min_mel = min_hz / freq_step
        freqs = freq_step * mels
        if mels >= min_mel:
            freqs = min_hz * np.exp(log_step * (mels - min_mel))
        return freqs


def dct(n_mfcc: int, n_mels: int) -> np.ndarray:
    """Compute the DCT matrix for MFCC computation."""
    result = np.zeros((n_mfcc, n_mels))
    s = (np.linspace(1, n_mels, n_mels) - 0.5) / n_mels

    for i in range(0, n_mfcc):
        result[i, :] = np.cos(i * np.pi * s) * np.sqrt(2.0 / n_mels)

    return result


def mel_filter(
    sample_rate: int,
    n_fft: int,
    f_min: float,
    f_max: float,
    n_mels: int,
    htk: bool = True,
) -> tuple[list[int], list[int], int, list[float], np.ndarray]:
    filters = np.zeros((n_mels, int(n_fft / 2 + 1)))
    zeros = np.zeros(int(n_fft // 2))

    fmin_mel = freq_to_mel(f_min, htk)
    fmax_mel = freq_to_mel(f_max, htk)
    mels = np.linspace(fmin_mel, fmax_mel, num=n_mels + 2)

    linearfreqs = np.linspace(0, sample_rate / 2.0, int(n_fft // 2 + 1))
    spectrogrammels = freq_to_mel(linearfreqs, htk)[1:]

    filt_pos = []
    filt_len = []
    packed_filters = []
    total_len = 0
    for n in range(n_mels):
        upper = (spectrogrammels - mels[n]) / (mels[n + 1] - mels[n])
        lower = (mels[n + 2] - spectrogrammels) / (mels[n + 2] - mels[n + 1])

        filters[n, :] = np.hstack([0, np.maximum(zeros, np.minimum(upper, lower))])
        num = 0
        first = True
        start_pos = 0
        end_pos = 0
        for sample in filters[n, :]:
            if first and sample != 0.0:
                first = False
                start_pos = num

            if not first and sample == 0.0:
                end_pos = num - 1
                break
            num = num + 1
        filt_len.append(end_pos - start_pos + 1)
        total_len += end_pos - start_pos + 1
        filt_pos.append(start_pos)
        packed_filters += list(filters[n, start_pos : end_pos + 1])

    return filt_len, filt_pos, total_len, packed_filters, filters


def render_mfcc_template(
    ty: str,
    window: np.ndarray,
    dct_matrix: np.ndarray,
    filt_pos: list[int],
    filt_len: list[int],
    packed_filters: list[float],
) -> str:
    dct_flat = np.asarray(dct_matrix).reshape(-1)
    window_arr = np.asarray(window)

    return template.format(
        TY=ty.replace("i", "u").upper(),
        ty=ty.replace("i", "u"),
        window_len=window_arr.size,
        window=format_array(window_arr, ty, scale=True),
        dct_len=dct_flat.size,
        dct=format_array(dct_flat, ty, scale=True),
        n_mels=len(filt_pos),
        filter_pos=format_array(filt_pos, "i16", scale=False),
        filter_len_array=format_array(filt_len, "i16", scale=False),
        total_filter_len=len(packed_filters),
        filter=format_array(packed_filters, ty, scale=True),
    )
