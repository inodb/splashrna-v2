"""
Build the browser-ready model bundle from the extracted Shogun parameters.

    python3 tools/build_model.py

Reads  models/raw/{mir30,mirE}.npz   (committed; produced by tools/extract_shogun.py)
Writes public/splashrna-model.json   (loaded by the static site at runtime)

Two pre-computations happen here so the browser does not have to:

  * The sqrt-diagonal normalizers are folded into the support-vector coefficients,
    giving  coef_i = alpha_i / ||sv_i||  for the weighted-degree kernel.
  * The two spectrum kernels are *collapsed* into a single weight vector each.
    Because their feature space is only 4 + 16 + 64 = 84 dimensions, the entire
    sum over support vectors reduces to one dot product per sequence:

        sum_i alpha_i K(sv_i,x) / (||sv_i|| ||x||)
          = (1/||x||) * sum_k w_k <  sum_i (alpha_i/||sv_i||) c_i[k],  c_x[k] >

    This removes roughly two thirds of the per-tile work.
"""

from __future__ import annotations

import base64
import json
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from reference import (  # noqa: E402
    COMBINE_THRESHOLD,
    COMBINE_WEIGHT,
    SPECTRUM_WEIGHTS,
    SplashSVM,
)

ROOT = pathlib.Path(__file__).resolve().parent.parent
RAW = ROOT / "models" / "raw"
OUT = ROOT / "public" / "splashrna-model.json"


def load_raw(name: str) -> SplashSVM:
    data = np.load(RAW / f"{name}.npz", allow_pickle=True)
    strings = [
        s.decode() if isinstance(s, bytes) else str(s) for s in data["sv_strings"]
    ]
    return SplashSVM(strings, data["alphas"], float(data["bias"][0]))


def collapse_spectrum(counts: list[np.ndarray], coef: np.ndarray) -> list[float]:
    """Fold the support vectors into one 84-dimensional weight vector."""
    return np.concatenate([coef @ counts[k] for k in range(3)]).tolist()


def pack(model: SplashSVM) -> dict:
    coef_wd = model.alphas / model.norm_wd
    coef_l = model.alphas / model.norm_left
    coef_r = model.alphas / model.norm_right

    return {
        "bias": model.bias,
        "numSupportVectors": len(model.sv_strings),
        # Concatenated ACGT; 22 characters per support vector. Plain text because
        # it compresses to ~2 bits per base over HTTP and stays human-inspectable.
        "supportVectors": "".join(model.sv_strings),
        # float64 little-endian, base64 encoded. float32 would halve the size but
        # costs an order of magnitude of accuracy against the original scores.
        "coefWd": base64.b64encode(coef_wd.astype("<f8").tobytes()).decode(),
        "spectrumLeft": collapse_spectrum(model.counts_left, coef_l),
        "spectrumRight": collapse_spectrum(model.counts_right, coef_r),
    }


def main() -> None:
    bundle = {
        "formatVersion": 1,
        "source": "Shogun 2.1.0 SVMLight models from splashrna.mskcc.org (2017)",
        "citation": (
            "Pelossof R, Fairchild L, et al. Prediction of potent shRNAs with a "
            "sequential classification algorithm. Nat Biotechnol. 2017;35(4):350-353."
        ),
        "kernel": {
            "sequenceLength": 22,
            "split": 15,
            "wdDegree": 10,
            "wdShift": 1,
            "spectrumOrder": 3,
            "spectrumWeights": SPECTRUM_WEIGHTS.tolist(),
        },
        "combine": {"threshold": COMBINE_THRESHOLD, "weight": COMBINE_WEIGHT},
        "models": {},
    }

    for name in ("mir30", "mirE"):
        model = load_raw(name)
        bundle["models"][name] = pack(model)
        print(f"{name}: {len(model.sv_strings)} SVs, bias {model.bias:.12f}")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(bundle, separators=(",", ":")))
    print(f"wrote {OUT} ({OUT.stat().st_size / 1024:.0f} KB)")


if __name__ == "__main__":
    main()
