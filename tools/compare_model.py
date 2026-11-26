"""
Compare two model bundles numerically.

    python3 tools/compare_model.py <a.json> <b.json> [tolerance]

Why this exists instead of a byte comparison: the bundle's float fields are computed
through NumPy, and different NumPy/BLAS builds differ in the last bit or two of a dot
product or square root. A fresh build on another machine is therefore *not* byte
identical to the committed one, even though both are correct. CI used to assert byte
equality and failed the moment the runner's NumPy (2.5.5) differed from the one the
bundle was built with (1.26.4).

What actually matters is that the two bundles are numerically the same model, so that
is what this checks. Structural fields must match exactly; floats must agree to the
given tolerance, which defaults to far tighter than anything that could move a score.

Exits non-zero if they differ.
"""

from __future__ import annotations

import base64
import json
import pathlib
import sys

import numpy as np

DEFAULT_TOLERANCE = 1e-12


def decode_coef(encoded: str) -> np.ndarray:
    return np.frombuffer(base64.b64decode(encoded), dtype="<f8")


def compare(a_path: pathlib.Path, b_path: pathlib.Path, tolerance: float) -> int:
    a = json.loads(a_path.read_text())
    b = json.loads(b_path.read_text())
    problems: list[str] = []

    for field in ("formatVersion", "kernel", "combine"):
        if a.get(field) != b.get(field):
            problems.append(f"{field} differs: {a.get(field)!r} vs {b.get(field)!r}")

    if set(a["models"]) != set(b["models"]):
        problems.append(f"model sets differ: {set(a['models'])} vs {set(b['models'])}")
        print("\n".join(problems))
        return 1

    worst = 0.0
    for name in sorted(a["models"]):
        ma, mb = a["models"][name], b["models"][name]

        if ma["supportVectors"] != mb["supportVectors"]:
            problems.append(f"{name}: support vector sequences differ")
        if ma["numSupportVectors"] != mb["numSupportVectors"]:
            problems.append(
                f"{name}: support vector count {ma['numSupportVectors']} vs {mb['numSupportVectors']}"
            )

        bias = abs(ma["bias"] - mb["bias"])
        worst = max(worst, bias)
        if bias > tolerance:
            problems.append(f"{name}: bias differs by {bias:.3e}")

        for field in ("coefWd",):
            va, vb = decode_coef(ma[field]), decode_coef(mb[field])
            if va.shape != vb.shape:
                problems.append(f"{name}.{field}: shape {va.shape} vs {vb.shape}")
                continue
            delta = float(np.abs(va - vb).max()) if va.size else 0.0
            worst = max(worst, delta)
            if delta > tolerance:
                problems.append(f"{name}.{field}: max difference {delta:.3e}")

        for field in ("spectrumLeft", "spectrumRight"):
            va, vb = np.array(ma[field]), np.array(mb[field])
            if va.shape != vb.shape:
                problems.append(f"{name}.{field}: shape {va.shape} vs {vb.shape}")
                continue
            delta = float(np.abs(va - vb).max())
            worst = max(worst, delta)
            if delta > tolerance:
                problems.append(f"{name}.{field}: max difference {delta:.3e}")

        print(f"  {name}: {ma['numSupportVectors']:,} support vectors")

    print(f"\nlargest numeric difference: {worst:.3e}  (tolerance {tolerance:.0e})")
    if problems:
        print("\nDIFFERENT:")
        for problem in problems:
            print(f"  - {problem}")
        return 1
    print("Bundles are numerically identical.")
    return 0


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(__doc__)
        raise SystemExit(2)
    tol = float(sys.argv[3]) if len(sys.argv) > 3 else DEFAULT_TOLERANCE
    raise SystemExit(compare(pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2]), tol))
