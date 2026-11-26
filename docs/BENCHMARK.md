# Benchmark: the rewrite against the original

Every number here was measured, and the commands to reproduce each are given. Nothing
is estimated or carried over from notes.

Measurements were taken on the machine that hosted the original service
(RHEL 8, 11 GB RAM), with the 2017 server running locally so both sides could be timed
on the same hardware. That server is the only oracle for this comparison and will not
exist forever, which is why the captured responses are committed as fixtures.

## 1. Does it give the same answers?

This is the only question that really matters. A faster tool that designs different
shRNAs is not a rewrite, it is a different tool.

**900 captured responses, 7,975 guides, 100% agreement.**

| Option variant | Cases | Agree |
|---|---|---|
| polyA cut on (the original's default) | 220 | 220 |
| polyA cut off | 220 | 220 |
| CDS only | 220 | 220 |
| Remove restriction sites | 60 | 60 |
| Keep in-gene duplicates | 60 | 60 |
| Keep off-target duplicates | 60 | 60 |
| 15 predictions (bypasses the server's cache) | 60 | 60 |

Agreement means the same guides, in the same order, with the same scores. The largest
score difference across all 7,975 guides was **1.9 × 10⁻⁷**.

That residual is not an approximation in this implementation. Shogun passed kernel
values through a float32 cache (`KERNELCACHE_ELEM`), so the *original* scores carry
about four float32 ulps of noise; this implementation computes in float64 throughout
and is marginally more precise than the thing it reproduces.

```bash
SPLASH_FULL_PARITY=1 npx vitest run tests/parity-sweep.test.ts
```

### How the sample was chosen

Hand-picked genes prove little, so the 220 genes were drawn as a stratified sample
across every axis that changes the code path:

| Stratum | Genes sampled |
|---|---|
| polyA site resolves → transcript is trimmed | 84 |
| No polyA site → falls back to the CDS intersection | 84 |
| No polyA site and no CDS → untrimmed transcript | 52 |
| Minus strand | 112 |
| Short (<1 kb) / mid / long (>5 kb) | all three bands |

The third row matters most: only about 63 genes in the entire table have neither a
polyA annotation nor a CDS, and the sample contains essentially all of them. It is the
rarest branch of `polyaCut()` and the one most likely to hide a mistake.

126 of the 900 cases produced no guides at all — a gene with no CDS asked for CDS-only,
or every candidate removed by a filter. Those count as agreement only when **both**
implementations return nothing, so they test the error paths rather than skipping them.

### A correction worth recording

The first run of this sweep reported mismatches, and the fault was in the measurement.
The extractor pulled any 22-mer out of the server's HTML and deduplicated the result,
which collapsed genuine repeated guides (exactly what "keep in-gene duplicates"
produces) and occasionally picked up a 22-nt run from the displayed transcript.

The page embeds the full result tuples, so the extractor now parses those and recovers
scores as well as guides. The earlier hand-picked fixtures were re-derived with the
precise parser and were unchanged — those genes happened to have no duplicates — but
that was luck. The broad sample is what exposed it.

## 2. Speed

Both sides designing shRNAs for the same gene, on the same machine.

| Gene | Tiles | 2017 server | This implementation |
|---|---|---|---|
| 7157 (TP53) | 1,846 | 0.36 s | **0.07 s** |
| 672 (BRCA2) | 3,130 | 0.52 s | **0.11 s** |
| 1 (A1BG) | — | 0.32 s | — |
| 10 (NAT2) | — | 0.26 s | — |

The 2017 figures are median wall time over three localhost requests at `n=15`, which
forces it to compute rather than serve from its precomputed cache. Timing the cached
path (`n=10`) gave the same numbers to within noise — 0.35 s and 0.52 s — so that
424 MB cache was not buying much.

This is not a like-for-like measurement and should not be read as one: the server
figure includes HTTP handling and HTML rendering, while the v2 figure is the design
pipeline in process. The honest summary is that both are fast enough, and the
interesting difference is not speed but what each needs in order to be fast.

### Engine internals

| | |
|---|---|
| One guide, both classifiers | 23 µs |
| Throughput | ~43,700 guides/s |
| 2,000 distinct guides (a 2 kb transcript) | 0.04 s |
| Model bundle parse and decode | 202 ms, once |
| Gene index load (1.27 MB) | 8 ms |

Guides are distinct pseudo-random sequences. An earlier version of this benchmark
rotated a single sequence and reported roughly three times the throughput; branch
prediction and cache locality made it meaningless.

### The optimisation that makes browser scoring viable

Scoring needs only the weighted total over support vectors, never the individual kernel
values, so the sum collapses into position-indexed lookup tables
(`src/engine/collapsed.ts`).

| | |
|---|---|
| Naive loop over 11,585 support vectors | 3,234 µs |
| Collapsed lookup tables | **5 µs** |
| Speedup | **660×** |
| One-off table construction at load | 247 ms |

The two forms are exactly equal, not approximately: `tests/collapsed.test.ts` asserts
they agree to 10 decimal places.

```bash
SPLASH_BENCH=1 npx vitest run tests/benchmark.test.ts
```

## 3. What each needs to run

This is where the rewrite actually pays off.

| | 2017 server | This implementation |
|---|---|---|
| Resident memory | **2.0 GB** | none — runs in the visitor's browser |
| Startup | ~30 s loading pickles | 202 ms model decode, in the page |
| Server processes | 1 (single-threaded CherryPy) | 0 |
| Runtime dependencies | Shogun 2.1.0 (C++), Python 2.7, CherryPy, BioPython, numpy, scipy | **0** |
| Application source | 2,133 lines Python | 1,961 lines TypeScript |
| Tests | none found | 95, plus a 900-case parity sweep |

The 2 GB is the load-bearing number. The original was fast because it held every
lookup table in RAM on a machine someone had to keep alive — and when that machine
rebooted, the service was down for nine days before anyone restarted it.

## 4. What a visitor downloads

The deployment is 97 MB, but range requests mean nobody fetches it:

| | Raw | Gzipped |
|---|---|---|
| Page, CSS, JavaScript | 33 KB | ~10 KB |
| Scoring model | 409 KB | **112 KB** |
| Off-target block directories (only if filtering) | 11 KB | 10 KB |
| Gene index (only if using a gene ID) | 1.24 MB | 656 KB |
| Each gene, each off-target lookup | a few KB | — |

About a megabyte for a typical session, against a 100 GB/month Pages allowance — room
for roughly 100,000 sessions.

The 58 MB transcript blob and 32 MB CDS blob are never downloaded whole. That depends
on the host honouring HTTP `Range`; GitHub Pages does, verified with ranged GETs at the
start, middle and end of the blob, with the returned bytes checked against the local
file 30 MB deep.

## Reproducing all of it

```bash
npm ci
npm test                                            # 95 tests
SPLASH_FULL_PARITY=1 npx vitest run tests/parity-sweep.test.ts
SPLASH_BENCH=1       npx vitest run tests/benchmark.test.ts
npm run browser-check                               # real Chromium against the build
```

The parity fixtures are committed, so the first three need nothing but this repository.
Re-capturing fresh ground truth would need the 2017 server, which may no longer exist.
