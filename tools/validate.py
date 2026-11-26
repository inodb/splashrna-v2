"""
Check the Python reference implementation against ground truth from the original
Shogun models.

    python3 tools/validate.py

The TypeScript engine has its own tests over the same fixtures. This script keeps the
Python reference honest, so it stays usable as the specification and as a second
opinion when the two ever disagree.
"""

from __future__ import annotations

import json
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from reference import SplashSVM, combine  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parent.parent
FIXTURES = ROOT / "tests" / "fixtures"
RAW = ROOT / "models" / "raw"

# Shogun caches kernel values as float32, so the recorded scores carry about four
# float32 ulps of noise. Anything below this is agreement.
TOLERANCE = 1e-6


def load(name: str) -> SplashSVM:
    data = np.load(RAW / f"{name}.npz", allow_pickle=True)
    strings = [s.decode() if isinstance(s, bytes) else str(s) for s in data["sv_strings"]]
    return SplashSVM(strings, data["alphas"], float(data["bias"][0]))


def check(label: str, worst: float) -> bool:
    ok = worst < TOLERANCE
    print(f"  {'PASS' if ok else 'FAIL'}  {label:<44} max error {worst:.3e}")
    return ok


def main() -> int:
    failures = 0
    models = {name: load(name) for name in ("mir30", "mirE")}

    print("Individual guide scores")
    scores = json.loads((FIXTURES / "scores.json").read_text())
    for name, key in (("mir30", "mir30_scores"), ("mirE", "mirE_scores")):
        got = models[name].score(scores["tests"])
        worst = float(np.abs(got - np.array(scores[key])).max())
        failures += not check(f"{name} classifier, {len(scores['tests'])} guides", worst)

    print("\nWhole transcript")
    tiles = json.loads((FIXTURES / "tiles.json").read_text())
    guides = tiles["tiles"]
    new30 = models["mir30"].score(guides)
    new_e = models["mirE"].score(guides)
    old30 = np.array(tiles["mir30"])
    old_e = np.array(tiles["mirE"])

    failures += not check(f"miR-30 over {len(guides)} tiles", float(np.abs(new30 - old30).max()))
    failures += not check(f"miR-E over {len(guides)} tiles", float(np.abs(new_e - old_e).max()))

    new_combined = combine(new30, new_e)
    old_combined = combine(old30, old_e)
    failures += not check("combined SplashRNA score", float(np.abs(new_combined - old_combined).max()))

    print("\nRanking")
    order_new = np.argsort(-new_combined)
    order_old = np.argsort(-old_combined)
    for n in (5, 10, 25, 100):
        identical = bool((order_new[:n] == order_old[:n]).all())
        print(f"  {'PASS' if identical else 'FAIL'}  top-{n} ranking identical")
        failures += not identical

    print()
    if failures:
        print(f"{failures} check(s) FAILED")
        return 1
    print("All checks passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
