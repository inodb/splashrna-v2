# SplashRNA v2

Prediction of potent miR-E shRNAs, as a **static website that runs entirely in your
browser**. No server, no database, no backend to maintain.

This is a ground-up rewrite of the original [splashrna.mskcc.org](http://splashrna.mskcc.org/)
server (Python 2.7, CherryPy, Shogun 2.1.0), which could no longer be deployed on
current systems. The published classifiers are preserved exactly: the rewrite
reproduces the original scores to within 10⁻⁶ and produces **identical shRNA
rankings**.

> Pelossof R, Fairchild L, Huang CH, et al. Prediction of potent shRNAs with a
> sequential classification algorithm. *Nature Biotechnology* 35, 350–353 (2017).
> [doi:10.1038/nbt.3807](https://doi.org/10.1038/nbt.3807)

## Why a rewrite

The original stack is unreproducible. Shogun 2.1.0 (2013) is unmaintained and will not
build on a modern toolchain; Python 2.7 is end-of-life; the classifiers existed only as
pickles that the original binary could load. The service had already died once.

A static site removes the failure mode entirely. There is no runtime to rot — the whole
application is 440 KB of HTML, CSS, JavaScript and model data that any web host will
still serve in ten years.

## Quick start

```bash
npm install
npm run dev        # development server
npm test           # 30 tests, including ground truth from the original server
npm run build      # static site in dist/
```

`dist/` can be dropped on any static host: GitHub Pages, S3, Netlify, or an
`nginx` document root. Nothing else is required.

## How it works

```
 FASTA input ──▶ resolve ──▶ tile into 22-mers ──▶ score ──▶ filter ──▶ rank
                   │                                 │
                   │ RefSeq / Entrez                 │ miR-30 classifier (gate)
                   ▼                                 ▼
             NCBI E-utilities                  miR-E classifier
             (direct from browser)             combined per the paper
```

Scoring is a two-stage sequential classifier. The miR-30 SVM gates, and the miR-E SVM
contributes only for candidates above threshold:

```
splash = mir30 > 1.1  ?  0.6·mir30 + 0.4·mirE  :  0.6·mir30
```

Each SVM uses a combined kernel over the 22-mer guide: a weighted-degree-position
string kernel with shifts, plus two weighted spectrum kernels over the 1–15 and 16–22
halves. See **[docs/MODEL-ARCHAEOLOGY.md](docs/MODEL-ARCHAEOLOGY.md)** for how these
were recovered from the original binaries and how every constant was verified.

## Layout

```
src/engine/      scoring engine — no dependencies, no DOM, runs in Node or a browser
  kernels.ts       the two string kernels
  collapsed.ts     the optimisation that makes browser scoring fast
  model.ts         model bundle loading
  scorer.ts        the sequential classifier
  design.ts        tiling, filtering and ranking pipeline
  sequence.ts      encoding, tiling, 97-mer constructs, restriction sites
  fasta.ts         input parsing
  export.ts        CSV export
src/data/        NCBI E-utilities client
src/ui/          stylesheet
src/main.ts      browser entry point
tools/           Python: model extraction, reference implementation, bundle build
models/raw/      extracted SVM parameters (committed; the legacy server is not needed)
tests/fixtures/  ground truth captured from the original Shogun models
public/          the browser-ready model bundle
```

The engine has **zero runtime dependencies**. Everything in `package.json` is a dev
tool for building and testing.

## Performance

A 2 kb transcript, ~2,000 candidate guides, scored against 12,601 support vectors:

| | Time |
|---|---|
| Naive per-support-vector loop | 6.2 s |
| With collapsed lookup tables | **0.35 s** |

The speedup is exact, not approximate — `tests/collapsed.test.ts` asserts the two forms
agree to 10 decimal places.

## Fidelity to the original

Verified against the running 2017 server:

| Check | Result |
|---|---|
| 44 guides, both classifiers | max error 3.3 × 10⁻⁷ |
| Every tile of a real transcript (1,979 guides) | max error 3.1 × 10⁻⁷ |
| Top 5 / 10 / 25 / 100 ranking | identical, including order |

The residual is roughly four float32 ulps and comes from Shogun's float32 kernel cache,
so this implementation is marginally *more* precise than the original.

## Not yet carried over

Two features depend on bulk datasets (~800 MB) that have not yet been sharded into
static assets. Both are honest gaps, surfaced in the UI rather than silently skipped:

* **Entrez gene IDs** currently resolve to a single RefSeq transcript. The original
  designed against the *intersection of a gene's isoforms*. Results for multi-isoform
  genes may differ, and the UI says so.
* **Off-target filtering** against the human and mouse transcriptomes. The interface
  exists and is tested against a stub; the data is not bundled.

See [docs/ROADMAP.md](docs/ROADMAP.md).

## Privacy

Sequences you type are scored in your browser and never transmitted. The one exception
is fetching a transcript by RefSeq accession or Entrez gene ID, which queries NCBI
E-utilities directly; only the accession or gene ID is sent. No email address or other
identifying information is included in those requests.

## License

GPL-3.0-or-later, following the original server, which carried GPL-2.0-or-later and
GPL-3.0-or-later headers from Christian Widmer / Max-Planck-Society (2009–2012).

The classifiers are the published SplashRNA models. If you use this tool, please cite
the paper and add to your methods:

> The shRNAs were designed using the SplashRNA algorithm (Pelossof and Fairchild
> et al., 2017).
