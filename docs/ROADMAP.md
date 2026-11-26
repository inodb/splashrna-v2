# Roadmap

Phase 1 (done) is a working static site: sequence and RefSeq input, exact scoring,
filtering, ranking, CSV export.

Phase 2 is about the bulk datasets the original server held in RAM. The goal is to keep
the site fully static, so these become sharded assets fetched on demand rather than a
database behind an API.

## The data

| Dataset | Size | Needed for |
|---|---|---|
| `hg38_22mer_matches.pickle` | 115 MB | Off-target filtering (human) |
| `mm10_22mer_matches.pickle` | 102 MB | Off-target filtering (mouse) |
| `intersection_*_seq_updated.txt` | 112 MB | Entrez gene IDs (isoform intersection) |
| `CDSintersection_seq_*.txt` | 61 MB | "CDS only" option |
| `reference_20151217.pickle` | 22 MB | Isoform coordinates for the gene diagram |
| `predictions/` | 424 MB | Precomputed top-10 per gene (a cache, not required) |

All of these live on the old host under
`plvsplashrna1-copy/thankamv/softwares/splashrna/server/server/dat/`. They are *not* in
this repository and must not be committed.

## Off-target filtering

Each table maps a 22-mer to the genes it matches. The design pipeline already takes an
`OffTargetIndex` interface, which is implemented and tested against a stub, so only the
data layer is missing.

Proposed approach — a sorted static index queried with HTTP range requests:

1. Sort all 22-mers, pack them 2 bits per base (6 bytes each) with a gene-list offset.
2. Emit one sorted key file plus a payload file per genome.
3. Ship a sparse in-memory prefix table (say every 4096th key, ~100 KB) so a lookup
   needs one or two range requests rather than a binary search from scratch.

Roughly 25 range requests for a default query. Both `Range` support and CORS are
standard on static hosts.

Why not shard by prefix: a 2 kb transcript produces ~2,000 tiles whose prefixes spread
across nearly as many shards, so prefix sharding degenerates into thousands of requests.
Range requests over one sorted file avoid that.

This is also why `offTargetMode` exists. The original filtered every tile before
truncating to N. `"top-candidates"` ranks first and checks only the best few hundred,
which returns the same list in almost every case; `"exhaustive"` reproduces the original
ordering when it matters. The default is the cheap path, and it is a user-visible
setting rather than a silent change.

## Entrez gene IDs

The original designed against the intersection of a gene's isoforms, not a single
transcript. Phase 1 resolves a gene ID to one RefSeq transcript and says so in the UI.

To restore the original behaviour, export the intersection table as one JSON file per
gene (~20 k small files, or grouped into range-addressable bundles) holding the
intersection sequence, genomic positions and strand. A gene ID then costs one fetch.

The `predictions/` directory is a cache of the top 10 per gene under default options. It
is not needed for correctness — the engine recomputes in well under a second — and is
probably not worth carrying over.

## Smaller items

* **Web Worker.** Scoring is fast enough to stay on the main thread, but moving it off
  would keep the UI perfectly responsive on very long transcripts.
* **CDS-only option.** Needs the CDS intersection table; same mechanism as above.
* **Gene diagram.** The original drew shRNA positions along the transcript. Worth
  rebuilding as inline SVG once isoform coordinates are available.
* **Isoform picker.** When a gene has several transcripts, let the user choose rather
  than defaulting to the first.

## Deployment

`npm run build` produces `dist/`, which is the entire application. For
`splashrna.mskcc.org`, point the document root at it and serve over HTTPS; no
application server, no Python, no ports to open. Keep the old host available until
phase 2 lands, since it is the only source of the bulk data and of regenerated ground
truth.
