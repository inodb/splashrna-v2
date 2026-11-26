# Roadmap

Phase 1 (done) is a working static site: sequence and RefSeq input, exact scoring,
ranking, CSV export.

Phase 2 (done) covered the bulk datasets the original server held in RAM: off-target
filtering and the isoform intersection are both static indexes now. The goal is to keep
the site fully static, so these become sharded assets fetched on demand rather than a
database behind an API.

## The data

| Dataset | Size | Needed for |
|---|---|---|
| `hg38_22mer_matches.pickle` | 115 MB | Off-target filtering (human) — **done** |
| `mm10_22mer_matches.pickle` | 102 MB | Off-target filtering (mouse) — **done** |
| `intersection_*_seq_updated.txt` | 112 MB | Entrez gene IDs (isoform intersection) — **done** |
| `CDSintersection_seq_*.txt` | 61 MB | "CDS only" option |
| `reference_20151217.pickle` | 22 MB | Isoform coordinates for the gene diagram |
| `predictions/` | 424 MB | Precomputed top-10 per gene (a cache, not required) |

These originally lived only in the `dat/` directory of the legacy server's install,
on a single internal host. They are now archived, gzipped and checksummed, in
[mskcc/splashrna-v2-data](https://github.com/mskcc/splashrna-v2-data) (private, 127 MB)
— the eight files that matter, minus the regenerable `predictions/` cache.

They are deliberately *not* in this repository: the build never reads them at runtime,
and carrying them here would slow every clone and CI checkout.

## Off-target filtering — done

Implemented in `tools/build_offtarget.py` and `src/data/offtarget.ts`.

The pickles are far smaller than their file size suggests: they hold only 22-mers
occurring in two or more genes, which is exactly the set a guide must avoid.

| | Keys | Index |
|---|---|---|
| hg38 (human) | 1,138,230 | 5.16 MB |
| mm10 (mouse) | 990,865 | 4.54 MB |

Blocks of 1024 keys, delta-varint encoded, each block self-contained so it can be
fetched alone. Block directories total 21 KB and are loaded up front; a lookup is an
in-memory binary search plus one Range request. Past a quarter of the file, the client
downloads the whole thing instead, which is what exhaustive mode triggers.

Prefix sharding was rejected: a 2 kb transcript produces ~2,000 tiles whose prefixes
spread across nearly as many shards, degenerating into thousands of requests.

This is also why `offTargetMode` exists. The original filtered every tile before
truncating to N. `"top-candidates"` ranks first and checks only the best few hundred,
which returns the same list in almost every case; `"exhaustive"` reproduces the original
ordering when it matters. The default is the cheap path, and it is a user-visible
setting rather than a silent change.

### Remaining: gene names

The index carries a count, not the matching gene IDs. Adding them means a second file
of uint16 gene indices (~13 MB for both genomes) plus a per-key offset, Range-requested
only for displayed guides. Worth doing if users ask which gene a guide hits; not worth
tripling the download otherwise.

## Entrez gene IDs, polyA cut and CDS-only — done

Implemented in `tools/build_genes.py` and `src/data/genes.ts`. 39,760 genes packed at
4 bits per symbol (A C G T * N) into a 58 MB blob with a 636 KB index; a gene costs one
Range request. Verified against the 2017 server on five genes, which return identical
results through the whole pipeline.

The index also carries each gene's polyA cut as a single prefix length, computed at
build time by replaying `polyaCut()`, so no genomic coordinate lists need shipping.
The CDS intersection rides along in a second blob for the CDS-only option and for the
fallback when a polyA site cannot be resolved.

Genomic coordinates themselves were not carried over — only the intersection sequence,
the CDS, the strand and the cut. They are needed only for the transcript diagram
below.

The `predictions/` directory is a cache of the top 10 per gene under default options. It
is not needed for correctness — the engine recomputes in well under a second — and is
probably not worth carrying over.

## Smaller items

* **Web Worker.** Scoring is fast enough to stay on the main thread, but moving it off
  would keep the UI perfectly responsive on very long transcripts.
* **CDS-only option.** Needs the CDS intersection table; same mechanism as above.
* **Gene diagram — done.** Inline SVG, built from the 2015 isoform annotations; see
  tools/build_annotations.py and src/ui/diagram.ts.
* **Isoform picker.** When a gene has several transcripts, let the user choose rather
  than defaulting to the first.

## Deployment

`npm run build` produces `dist/`, which is the entire application. For
`splashrna.mskcc.org`, point the document root at it and serve over HTTPS; no
application server, no Python, no ports to open. Keep the old host available until
phase 2 lands, since it is the only source of the bulk data and of regenerated ground
truth.
