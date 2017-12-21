import numpy as np
from mtspec import mtspec, mt_coherence


def shift_bit_length(x):
    return 1 << (x - 1).bit_length()


def padwithzeros(vector, pad_width, iaxis, kwargs):
    vector[:pad_width[0]] = 0
    vector[-pad_width[1]:] = 0
    return vector


def get_periodicity(values):
    tbp = 4
    kspec = 3
    nf = 30
    p = 0.90
    length = len(values)
    next_pow2_length = shift_bit_length(length)
    values = np.lib.pad(values,  (0, next_pow2_length -
                                  len(values) % next_pow2_length), padwithzeros)
    mean_centered_values = values - np.nanmean(values)
    normalized_values = mean_centered_values / \
        np.max(np.abs(mean_centered_values))
    uniform_signal = [1, -0.6, -0.4] * (next_pow2_length // 3)
    uniform_signal = np.lib.pad(
        uniform_signal, (0, next_pow2_length - len(uniform_signal) % next_pow2_length), padwithzeros)

    out = mt_coherence(1,
                       normalized_values,
                       uniform_signal,
                       tbp,
                       kspec,
                       nf,
                       p,
                       freq=True,
                       phase=True,
                       cohe=True,
                       iadapt=1)

    return out['cohe'][np.argmin(np.abs(out['freq'] - 1 / 3.0))]
