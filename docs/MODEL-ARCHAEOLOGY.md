# Model archaeology

How the SplashRNA classifiers were recovered from the 2017 server, and why the
constants in `tools/reference.py` and `src/engine/kernels.ts` are what they are.

Read this before changing anything in the scoring path.

## The problem

The published classifiers existed only as bzip2'd Python 2 pickles of **Shogun 2.1.0**
objects:

| File | Size | Model |
|---|---|---|
| `svm_VS3M1pos_2015_02_12.bz2` | 1.0 MB | miR-30 classifier |
| `mirE_plus_ultramir_c15_2015_04_27.bz2` | 48 KB | miR-E classifier |

They could only be unpickled by the exact Shogun build on the original host. Shogun is
effectively unmaintained, that build cannot be reproduced, and nothing about the models
was recoverable from the paper alone. Until this was solved, no rewrite was possible.

## Step 1 — reading the pickles

The pickles reference `__main__.ShRNAtor` and `modshogun.SVMLight`. Loading them needs
the legacy environment, but once loaded, the SVM exposes its parameters directly:

```
miR-30:  28,042 training examples, 11,585 support vectors, bias -0.840949519872
miR-E:    1,177 training examples,  1,016 support vectors, bias -0.084442271170
```

`tools/extract_shogun.py` does this and writes `models/raw/*.npz`. **Those files are
committed**, so the legacy environment is never needed again.

## Step 2 — recovering the kernel

Shogun's Python bindings would not downcast the kernel objects, so the hyper-parameters
came from a full ASCII serialization (`svm.save_serializable(SerializableAsciiFile(...))`),
which dumps every field:

```
CombinedKernel, IdentityKernelNormalizer, subkernel weights [1, 1, 1]
├── WeightedDegreePositionStringKernel   SqrtDiagKernelNormalizer
│     degree 10, shift [1]*22, max_shift 1
│     weights [0.18181…, 0.16363…, …]       = 2(d-k+1)/(d(d+1))
│     input: the full 22-mer
├── WeightedCommWordStringKernel          SqrtDiagKernelNormalizer
│     order 3, use_sign false, weights [0.7071, 0.5774, 0.4082]
│     input: nucleotides 1–15  (13 order-3 words)
└── WeightedCommWordStringKernel          same parameters
      input: nucleotides 16–22 (5 order-3 words)
```

So the decision function is

```
score(x) = Σᵢ αᵢ · Σⱼ Kⱼ(svᵢ, x) / (‖svᵢ‖ⱼ · ‖x‖ⱼ) + bias
```

## Step 3 — the two constants the dump did not settle

A serialized field list is not an implementation. Two details had to be pinned down
empirically, by probing the original Shogun with controlled inputs and matching numbers.

### The shift weight is 1/(s+1), not 1/(2(s+1))

Shogun's documentation gives the shift weighting as `β_s = 1/(2(s+1))`. Implementing
that reproduced the *unshifted* kernel exactly — a self-comparison of a non-repetitive
22-mer gives `Σ_k w_k (L-k+1) = 19.0` on the nose — but every shifted value came out
short by exactly a factor of two in the shift term. The installed build uses
`β_s = 1/(s+1)`, so `β_1 = 1/2`.

Verified against a recorded 3×3 kernel matrix; see `tests/kernels.test.ts`.

### The spectrum kernel masks suffixes, with weights 1/2, 1/3, 1/6

`get_spectrum_features` builds order-3 words with `reverse=True`. The consequence is
subtle and is the easiest way to get this wrong: **degree masking retains the last k
characters of each word, not the first.**

It was settled by probing Shogun with controlled 7-mers:

| Probe | Shogun | Suffix masking | Prefix masking |
|---|---|---|---|
| `AAAAAAA` vs `AAAAAAA` | 25.000 | 25.000 ✓ | 25.000 ✓ |
| `AAAAAAA` vs `AAAAAAC` | 20.000 | 20.000 ✓ | 24.167 ✗ |
| `TTCAGTT` vs itself | 6.000 | 6.000 ✓ | 9.878 ✗ |
| `AAAAACC` vs `CCCCCCC` | 6.667 | 6.667 ✓ | — |

The weights are the *squares* of the values in the dump: the stored
`[0.7071, 0.5774, 0.4082]` are `sqrt(3/6), sqrt(2/6), sqrt(1/6)`, and the kernel applies
them squared, giving `[1/2, 1/3, 1/6]` — which sum to exactly 1. The `AAAAAAA`
self-comparison of 25.000 is what proves it: all three masked levels contribute 25, so
the weights must sum to one.

These probes are preserved as `tests/fixtures/spectrum-probe.json` and asserted in the
test suite, so the orientation can never silently flip.

## Step 4 — validation

| Check | Result |
|---|---|
| 44 individual guides, both models | max error 3.3 × 10⁻⁷ |
| 1,979 tiles of a real transcript (Entrez 100129520) | max error 3.1 × 10⁻⁷ |
| Combined SplashRNA score | max error 1.9 × 10⁻⁷ |
| Top 5 / 10 / 25 / 100 ranking | **identical, including order** |

The residual is about four float32 ulps. Shogun passes kernel values through a float32
cache (`KERNELCACHE_ELEM`), so the *original* scores carry that noise and this
implementation is marginally more precise than the one it replaces. The errors are
mixed-sign rounding, not systematic bias.

Ground truth lives in `tests/fixtures/`. Regenerating it requires the legacy server;
the scripts that produced it are in `tools/`.

## Step 5 — collapsing the support-vector sum

Scoring needs only `S(x) = Σᵢ coefᵢ · K(svᵢ, x)`, never individual kernel values. Since

```
K(sv, x) = Σⱼ Σₖ wₖ · [ sv[j..j+k] == x[j..j+k] ]
```

the sum over support vectors can move inside, giving position-indexed tables
`T[j][k][code]`. Both shifted alignments reuse the same tables, because the table for
"support vector k-mer at position j+1" is just `T[j+1]`.

Degrees 1–7 are dense arrays; degrees 8–10 are rare enough to resolve through an 8-mer
index over the handful of candidate support vectors. This turned 11,585 kernel
evaluations per guide into a few hundred lookups — **an 18× end-to-end speedup**
(6.2 s → 0.35 s for a 2 kb transcript) with bit-identical output, asserted by
`tests/collapsed.test.ts`.

## Step 6 — a bug found in the off-target tables

While building the static off-target index, the legacy loader turned out to assign the
two match tables to the wrong species:

```python
server.py:662   self.mouseMatch = pickle.load(open("dat/hg38_22mer_matches.pickle"))
server.py:663   self.humanMatch = pickle.load(open("dat/mm10_22mer_matches.pickle"))
```

`hg38` is the human assembly and `mm10` the mouse one, so these are crossed. Confirmed
by checking marker genes in each file's Entrez ID set:

| File | TP53/ACTB/EGFR/KRAS (human) | Trp53/Myc/Actb/Egfr/Kras (mouse) |
|---|---|---|
| `hg38_22mer_matches.pickle` | present | absent |
| `mm10_22mer_matches.pickle` | absent | present |

The sets are disjoint, so this is unambiguous. **Filtering was unaffected** — a guide
is dropped when either table matches, and swapping two tables does not change their
union — but every species label the server reported, in the UI and in the
`Human.22mer.match.genes` / `Mouse.22mer.match.genes` CSV columns, was wrong.

`tools/build_offtarget.py` keys its output off the assembly name, so v2 labels these
correctly. Anyone comparing v2 output against archived 2017 CSVs should expect the two
species columns to be exchanged.

## What is not yet carried over

* **Entrez gene IDs** resolve to a single RefSeq transcript. The original designed
  against the *intersection of a gene's isoforms*, using a 112 MB precomputed table.
  The UI states this difference rather than hiding it.
* **Gene names** for off-target hits. The index carries counts; the name lists would
  roughly triple its size and are only needed for the handful of displayed guides.

See `docs/ROADMAP.md`.
