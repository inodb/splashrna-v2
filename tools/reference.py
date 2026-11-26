"""
Reference implementation of the SplashRNA SVM scoring function, without Shogun.

This module is the *authority* on what the TypeScript engine must reproduce. It is
deliberately written for clarity over speed; `tools/validate.py` checks it against
ground truth captured from the original Shogun 2.1.0 models.

Kernel specification (recovered from the original serialized models -- see
docs/MODEL-ARCHAEOLOGY.md for how these constants were determined):

  Classifier: SVMLight, combined kernel of three sub-kernels, equal weights,
  identity normalizer on the combination, sqrt-diagonal normalizer on each part.

  1. Weighted-degree-position string kernel over the full 22-mer
       degree d = 10
       weights w_k = 2(d - k + 1) / (d(d + 1))      for k = 1..10
       shift s = 1 at every position, beta_s = 1/(s + 1)
  2. Weighted spectrum kernel over nucleotides 1..15
       order 3, SUFFIX masking, weights [1/2, 1/3, 1/6]
  3. Weighted spectrum kernel over nucleotides 16..22, same parameters

  score(x) = sum_i alpha_i * sum_j K_j(sv_i, x) / (||sv_i||_j * ||x||_j) + bias
"""

from __future__ import annotations

import numpy as np

DEGREE = 10
WD_WEIGHTS = np.array(
    [2.0 * (DEGREE - k + 1) / (DEGREE * (DEGREE + 1)) for k in range(1, DEGREE + 1)]
)
# WD_CUMULATIVE[k] = sum of the first k degree weights; a match of run length k
# contributes exactly this much.
WD_CUMULATIVE = np.concatenate([[0.0], np.cumsum(WD_WEIGHTS)])

SPECTRUM_ORDER = 3
SPECTRUM_WEIGHTS = np.array([1 / 2.0, 1 / 3.0, 1 / 6.0])  # suffix lengths 1, 2, 3

BASES = "ACGT"
BASE_INDEX = {c: i for i, c in enumerate(BASES)}

SPLIT = 15  # boundary between the two spectrum sub-kernels


def encode(seqs: list[str]) -> np.ndarray:
    """Encode ACGT strings as an (n, L) array of base indices."""
    return np.array([[BASE_INDEX[c] for c in s] for s in seqs], dtype=np.int8)


def _run_lengths(eq: np.ndarray, cap: int) -> np.ndarray:
    """For a boolean (n, L) match array, the run of Trues starting at each column."""
    n, length = eq.shape
    out = np.zeros((n, length), dtype=np.int32)
    out[:, length - 1] = eq[:, length - 1]
    for j in range(length - 2, -1, -1):
        out[:, j] = np.where(eq[:, j], np.minimum(out[:, j + 1] + 1, cap), 0)
    return out


def wd_kernel(mat: np.ndarray, vec: np.ndarray, shift: int = 1) -> np.ndarray:
    """Unnormalized weighted-degree-position kernel of every row of `mat` with `vec`."""
    _, length = mat.shape
    total = WD_CUMULATIVE[_run_lengths(mat == vec[None, :], DEGREE)].sum(axis=1)
    for s in range(1, shift + 1):
        beta = 1.0 / (s + 1.0)
        total += beta * WD_CUMULATIVE[
            _run_lengths(mat[:, s:] == vec[None, : length - s], DEGREE)
        ].sum(axis=1)
        total += beta * WD_CUMULATIVE[
            _run_lengths(mat[:, : length - s] == vec[None, s:], DEGREE)
        ].sum(axis=1)
    return total


def spectrum_counts(mat: np.ndarray) -> list[np.ndarray]:
    """Suffix-masked k-mer count vectors for suffix lengths 1, 2 and 3.

    Shogun builds order-3 words with `reverse=True`, which makes the degree
    masking keep the *last* k characters of each word rather than the first.
    """
    n, length = mat.shape
    n_words = length - SPECTRUM_ORDER + 1
    word3 = mat[:, 0:n_words] * 16 + mat[:, 1 : n_words + 1] * 4 + mat[:, 2 : n_words + 2]
    word2 = mat[:, 1 : n_words + 1] * 4 + mat[:, 2 : n_words + 2]
    word1 = mat[:, 2 : n_words + 2]

    counts = []
    for words, size in ((word1, 4), (word2, 16), (word3, 64)):
        table = np.zeros((n, size))
        for row in range(n):
            table[row] = np.bincount(words[row], minlength=size)
        counts.append(table)
    return counts


def spectrum_kernel(counts_a: list[np.ndarray], counts_b: list[np.ndarray]) -> np.ndarray:
    """Unnormalized weighted spectrum kernel between count matrices and a count vector."""
    return sum(SPECTRUM_WEIGHTS[i] * (counts_a[i] @ counts_b[i]) for i in range(3))


class SplashSVM:
    """One of the two SplashRNA sub-classifiers (miR-30 or miR-E)."""

    def __init__(self, sv_strings: list[str], alphas, bias: float):
        self.sv_strings = list(sv_strings)
        self.alphas = np.asarray(alphas, dtype=np.float64)
        self.bias = float(bias)

        self.enc = encode(self.sv_strings)
        self.enc_left = self.enc[:, :SPLIT]
        self.enc_right = self.enc[:, SPLIT:]
        self.counts_left = spectrum_counts(self.enc_left)
        self.counts_right = spectrum_counts(self.enc_right)

        n = len(self.enc)
        self.norm_wd = np.sqrt(
            np.array([wd_kernel(self.enc[i : i + 1], self.enc[i])[0] for i in range(n)])
        )
        self.norm_left = np.sqrt(self._self_spectrum(self.counts_left, n))
        self.norm_right = np.sqrt(self._self_spectrum(self.counts_right, n))

    @staticmethod
    def _self_spectrum(counts: list[np.ndarray], n: int) -> np.ndarray:
        return np.array(
            [
                sum(SPECTRUM_WEIGHTS[k] * counts[k][i] @ counts[k][i] for k in range(3))
                for i in range(n)
            ]
        )

    def score(self, seqs: list[str]) -> np.ndarray:
        out = np.empty(len(seqs))
        for t, seq in enumerate(seqs):
            vec = encode([seq])[0]
            left, right = vec[:SPLIT], vec[SPLIT:]
            c_left = [c[0] for c in spectrum_counts(left[None, :])]
            c_right = [c[0] for c in spectrum_counts(right[None, :])]

            n_wd = np.sqrt(wd_kernel(vec[None, :], vec)[0])
            n_l = np.sqrt(sum(SPECTRUM_WEIGHTS[k] * c_left[k] @ c_left[k] for k in range(3)))
            n_r = np.sqrt(sum(SPECTRUM_WEIGHTS[k] * c_right[k] @ c_right[k] for k in range(3)))

            k_wd = wd_kernel(self.enc, vec) / (self.norm_wd * n_wd)
            k_l = spectrum_kernel(self.counts_left, c_left) / (self.norm_left * n_l)
            k_r = spectrum_kernel(self.counts_right, c_right) / (self.norm_right * n_r)

            out[t] = self.alphas @ (k_wd + k_l + k_r) + self.bias
        return out


# --- SplashRNA sequential combination -------------------------------------------------
# From the original server (server.py:277-285): the miR-30 classifier gates, and the
# miR-E classifier only contributes for candidates that clear the threshold.
COMBINE_THRESHOLD = 1.1
COMBINE_WEIGHT = 0.6


def combine(mir30: np.ndarray, mire: np.ndarray) -> np.ndarray:
    mir30 = np.asarray(mir30, dtype=np.float64)
    mire = np.asarray(mire, dtype=np.float64)
    return np.where(
        mir30 > COMBINE_THRESHOLD,
        mir30 * COMBINE_WEIGHT + mire * (1 - COMBINE_WEIGHT),
        mir30 * COMBINE_WEIGHT,
    )
