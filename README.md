# SplashRNA

Pick potent miR-E shRNAs for a target. Paste a sequence, or look up a human gene by ID
or RefSeq accession, and get back the best-scoring guides with their 97-mer constructs.

**→ [mskcc.github.io/splashrna-v2](https://mskcc.github.io/splashrna-v2/)**

It is a static site. Scoring, off-target filtering and gene lookup all happen in the
browser — there is no server, no database and no backend to keep alive.

> Pelossof R, Fairchild L, Huang CH, et al. Prediction of potent shRNAs with a
> sequential classification algorithm. *Nature Biotechnology* 35, 350–353 (2017).
> [doi:10.1038/nbt.3807](https://doi.org/10.1038/nbt.3807)

This replaces the original [splashrna.mskcc.org](http://splashrna.mskcc.org/) server
(Python 2.7, CherryPy, Shogun 2.1.0), which could no longer be deployed. The published
classifiers are unchanged: scores match the originals to within 10⁻⁶ and the shRNAs
come back in the same order.

## Quick start

```bash
npm install
npm run dev              # development server
npm test                 # 95 tests, including ground truth from the 2017 server
npm run build            # static site in dist/
npm run browser-check    # build it, then drive it in a real Chromium
```

`browser-check` is the one that catches what unit tests cannot: that the page renders,
the model and indexes load over HTTP, a design run fills the results table, nothing
throws or 404s, dark mode applies, and nothing overflows at phone width. Screenshots
land in `$SHOT_DIR`.

## Deploying

`dist/` is the whole application. Drop it on GitHub Pages, S3, Netlify, or an nginx
document root — anything that serves files and honours HTTP `Range`.

Pushing to `main` deploys to GitHub Pages automatically
(`.github/workflows/pages.yml`). The build needs no secrets: the data indexes are
committed, so a deploy is `npm ci && npm run build`.

**`Range` support is not optional.** Gene lookup reads a few kilobytes out of a 58 MB
blob with a ranged request. On a host that ignores `Range` the client still gets the
right answer, but it downloads the whole blob to do it.

### What a visitor actually downloads

The deployment is ~97 MB, but nobody downloads 97 MB:

| | |
|---|---|
| Page, CSS, JavaScript | 33 KB |
| Scoring model | 409 KB |
| Off-target block directories | 21 KB, only if filtering is on |
| Gene index | 1.3 MB, only if you use a gene ID |
| Each gene, each off-target lookup | a few KB, by ranged request |

So roughly a megabyte for a typical session.

## How it works

```
 target ──▶ resolve ──▶ tile into 22-mers ──▶ score ──▶ filter ──▶ rank
              │                                 │
              │ gene ID → isoform intersection  │ miR-30 classifier (gate)
              │ RefSeq  → NCBI E-utilities      ▼
              ▼                           miR-E classifier
        local static index                combined per the paper
```

Scoring is a two-stage sequential classifier. The miR-30 SVM gates, and the miR-E SVM
only contributes for candidates above threshold:

```
splash = mir30 > 1.1  ?  0.6·mir30 + 0.4·mirE  :  0.6·mir30
```

Each SVM uses a combined kernel over the 22-mer guide: a weighted-degree-position
string kernel with shifts, plus two weighted spectrum kernels over the 1–15 and 16–22
halves.

Recovering those kernels from the original binaries was the hard part of this rewrite,
and two constants had to be determined experimentally because the serialized model did
not settle them. **[docs/MODEL-ARCHAEOLOGY.md](docs/MODEL-ARCHAEOLOGY.md)** records the
evidence for each. Read it before touching the scoring path.

## Does it match the original?

Yes, and the tests enforce it on every push rather than taking my word for it.

| Check | Result |
|---|---|
| 44 guides, both classifiers | max error 3.3 × 10⁻⁷ |
| Every tile of a real transcript (1,979 guides) | max error 3.1 × 10⁻⁷ |
| Top 5 / 10 / 25 / 100 ranking | identical, including order |
| **Whole pipeline, 8 genes, 29 cases** | **identical guides, in order** |

The residual is about four float32 ulps and comes from Shogun's float32 kernel cache,
so this implementation is marginally *more* precise than the one it replaces.

The last row matters most. While the 2017 server was still running it was queried for
genes 1, 10, 207, 672 (BRCA2), 4609 (MYC), 5728 (PTEN), 7157 (TP53) and 71458 — with
the polyA cut on and off, and with CDS-only — and the complete pipeline returns the
same guides in the same order every time. Those responses are committed as fixtures, so
the parity tests outlive the old server.

## Performance

A 2 kb transcript, ~2,000 candidate guides, scored against 12,601 support vectors:

| | Time |
|---|---|
| Naive loop over support vectors | 6.2 s |
| Collapsed lookup tables | **0.35 s** |

The sum over support vectors collapses into position-indexed tables, so a guide costs a
few hundred lookups instead of 11,585 kernel evaluations. The speedup is exact, not
approximate — `tests/collapsed.test.ts` asserts the two forms agree to 10 decimal
places.

## The data indexes

Three features need bulk data from the 2017 server: off-target filtering, gene lookup,
and the polyA cut. All three are static indexes, built from the archive in
[mskcc/splashrna-v2-data](https://github.com/mskcc/splashrna-v2-data) (private):

```bash
python3 tools/build_offtarget.py /path/to/dat public/offtarget
python3 tools/build_genes.py     /path/to/dat public/genes
```

The outputs are committed, so you only need to rerun these if the source tables change.

**Off-target filtering** drops guides whose 22-mer exactly matches a second gene, in
human (hg38) or mouse (mm10). 2.1 M keys in 9.7 MB, as blocks of 1024; a lookup is an
in-memory binary search plus one ~4 KB ranged request.

**Gene lookup** resolves an Entrez ID to the intersection of that gene's RefSeq
isoforms — the stretch common to all of them — which is what the original designed
against, not any single transcript. 39,760 genes, 4 bits per symbol, in a 58 MB
transcript blob and a 32 MB CDS blob. `*` marks a break between non-contiguous blocks;
windows spanning one are skipped and positions reported in `*`-free coordinates.

**The polyA cut** trims each transcript at its annotated cleavage site, so no guide
targets 3′ sequence absent from the mature mRNA. The original had this **on by
default**, and so does this. It reduces to one integer per gene, precomputed at build
time. 22,513 of 39,760 genes resolve a cut; the rest fall back to the CDS intersection,
or to the untrimmed transcript when there is no CDS — the same fallbacks, with the same
wording.

That cut is not cosmetic: for gene 10 it changes 4 of the top 10 guides.

> **A bug in the original, corrected here.** The 2017 server loaded the human table
> into its `mouseMatch` variable and the mouse table into `humanMatch`
> (`server.py:662-663`), so the Human and Mouse columns it reported were swapped.
> Filtering was unaffected — a guide is dropped if either matches — but every species
> label was wrong. Archived CSVs from that server have those two columns exchanged.

## Layout

```
src/engine/      scoring engine — no dependencies, no DOM, runs in Node or a browser
  kernels.ts       the two string kernels
  collapsed.ts     the optimisation that makes browser scoring fast
  model.ts         model bundle loading
  scorer.ts        the sequential classifier
  design.ts        tiling, filtering and ranking
  sequence.ts      encoding, tiling, 97-mer constructs, restriction sites
  fasta.ts/export.ts   input parsing, CSV output
src/data/        genes.ts, offtarget.ts (static indexes), ncbi.ts (E-utilities)
src/main.ts      browser entry point; src/ui/ the stylesheet
tools/           Python builders, the Python reference implementation, browser-check
models/raw/      extracted SVM parameters — the legacy server is not needed to build
tests/fixtures/  ground truth captured from the 2017 server
public/          model bundle and the built indexes, served as-is
```

The engine has **zero runtime dependencies**. Everything in `package.json` is a dev
tool.

## Not yet carried over

* Gene *names* for off-target hits. The index says how many genes a guide matches, not
  which; the name lists would roughly triple its size.
* The transcript diagram, and the `/help`, `/contact` and `/cite` pages.

See [docs/ROADMAP.md](docs/ROADMAP.md).

## Privacy

Sequences you paste are scored in your browser and never leave it. Gene IDs resolve
against the local index, so they are not sent anywhere either. The one exception is a
RefSeq accession, which is fetched from NCBI E-utilities — only the accession goes out,
with no email address or other identifying information attached.

## License

GPL-3.0-or-later, following the original server, which carried GPL-2.0-or-later and
GPL-3.0-or-later headers from Christian Widmer / Max-Planck-Society (2009–2012).

If you use these predictions, please cite the paper and add to your methods:

> The shRNAs were designed using the SplashRNA algorithm (Pelossof and Fairchild
> et al., 2017).
