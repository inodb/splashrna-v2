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
npm test           # 53 tests, including ground truth from the original server
npm run build      # static site in dist/
npm run browser-check   # drive the built site in a real browser
```

`browser-check` builds the site, serves it on localhost and drives it with Chromium:
it verifies the page renders, the model and indexes load over HTTP, a design run
populates the results table, nothing throws or 404s, dark mode applies, and there is no
horizontal overflow at phone width. Screenshots land in `$SHOT_DIR`.

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
| **Full pipeline, 5 genes by Entrez ID** | **identical ten shRNAs, in order** |

The residual is roughly four float32 ulps and comes from Shogun's float32 kernel cache,
so this implementation is marginally *more* precise than the original.

The last row is the one that matters most. While the 2017 server was still running, it
was queried for genes 1, 10, 672 (BRCA2), 7157 (TP53) and 71458, and the complete
pipeline — isoform intersection, tiling across block boundaries, scoring, in-gene
de-duplication, off-target filtering, ranking — returns exactly the same ten shRNAs in
exactly the same order for all five. Those responses are committed as fixtures
(`tests/fixtures/legacy-genes.json`), so the parity test survives the old server.

## Off-target filtering

Guides that exactly match a second gene are filtered out, against both the human
(hg38) and mouse (mm10) transcriptomes. This is an **optional deployment**: without it
the site works and says so in the results; with it, filtering activates automatically.

```bash
python3 tools/build_offtarget.py /path/to/legacy/dat public/offtarget
```

That reads the legacy 22-mer match tables and produces a 9.7 MB block-structured index
(8.2 MB gzipped). Only the 21 KB block directories are fetched up front; individual
lookups pull a single ~4 KB block with an HTTP Range request, and a query that would
touch most of the file downloads it in one go instead.

The index is derived data and is not committed. What *is* committed is a 750-entry
sample from the original tables, which the test suite checks the index against.

> **A bug in the original, corrected here.** The 2017 server loaded the human table
> into its `mouseMatch` variable and the mouse table into `humanMatch`
> (`server.py:662-663`), so the Human and Mouse columns it reported were swapped.
> Filtering was unaffected — a guide is dropped if either matches — but every species
> label was wrong. This version keys off the assembly name.

## Gene IDs and the isoform intersection

The original designed against the *intersection* of a gene's RefSeq isoforms — the
stretch common to all of them — not any single transcript. That table is also available
as a static index:

```bash
python3 tools/build_genes.py /path/to/legacy/dat public/genes
```

39,760 genes packed 4 bits per symbol into a 58 MB blob, with a 636 KB index. Only the
index is downloaded; a gene costs one Range request of a few kilobytes, so the blob's
size costs the user nothing.

Without it, Entrez IDs fall back to a single RefSeq transcript from NCBI and the UI
says so. With it, gene IDs behave exactly as they did in 2017.

`*` in an intersection sequence marks a break between non-contiguous blocks. Windows
spanning one are skipped and positions are reported in `*`-free coordinates, matching
`entrezID_prepare()` in the original.

## Not yet carried over

Gene *names* for off-target hits. The index reports how many genes a guide matches, not
which; carrying the name lists would roughly triple its size. See
[docs/ROADMAP.md](docs/ROADMAP.md).

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
