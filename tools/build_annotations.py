"""
Build the static annotation index that backs the transcript diagram.

    python3 tools/build_annotations.py /path/to/legacy/dat public/annotations

Sources on the legacy server:
    reference_20151217.pickle          per-gene isoform exon structures
    intersection_20151217_updated.txt  genomic ranges of the isoform intersection

Why ship our own annotations rather than embed a genome browser
---------------------------------------------------------------
igv.js and JBrowse 2 both run without a backend and can load hg38 RefSeq tracks from
hosted files, which would save building this. They were rejected for three reasons:

  1. Correctness. Those tracks carry *current* RefSeq annotations. The predictions were
     computed against the 2015 annotation frozen in these files, so a diagram drawn
     from today's RefSeq could show isoforms that are not the ones the intersection was
     derived from. A picture that disagrees with the numbers beside it is worse than no
     picture.
  2. Durability. Hosted annotation files go down -- the IGV forum has reports of
     exactly that taking out the RefSeq track. This rewrite exists because the previous
     version depended on things that stopped working.
  3. Privacy. The site scores sequences locally; sending the genes a user is looking at
     to a third-party host would quietly undo that.

Format
------
`annotations.bin`  per gene, two newline-terminated records, verbatim in the legacy
  encoding so the browser parses exactly what the server did:

      chrom,strand:REFSEQ,annot,start,end;annot,start,end:REFSEQ2,...
      start,end*;start,end;...            (intersection ranges)

`annotations.idx`  one 20-byte record per gene, sorted by gene ID:

      uint32  Entrez gene ID
      uint32  byte offset of the isoform record
      uint32  length of the isoform record
      uint32  byte offset of the intersection record
      uint32  length of the intersection record
"""

from __future__ import annotations

import json
import pathlib
import pickle
import struct
import sys

RECORD = struct.Struct("<IIIII")
REFERENCE = "reference_20151217.pickle"
INTERSECTION = "intersection_20151217_updated.txt"


def main() -> int:
    if len(sys.argv) < 2:
        print(__doc__)
        return 2

    dat = pathlib.Path(sys.argv[1])
    out_dir = pathlib.Path(sys.argv[2] if len(sys.argv) > 2 else "public/annotations")
    for name in (REFERENCE, INTERSECTION):
        if not (dat / name).exists():
            print(f"missing: {dat / name}", file=sys.stderr)
            return 1

    with (dat / REFERENCE).open("rb") as handle:
        isoforms: dict[str, str] = pickle.load(handle, encoding="latin1")

    ranges: dict[str, str] = {}
    for line in (dat / INTERSECTION).open():
        parts = line.strip().split()
        if len(parts) == 4 and parts[0].isdigit():
            ranges[parts[0]] = parts[3]

    gene_ids = sorted(
        {g for g in isoforms if g.isdigit()} | {g for g in ranges if g.isdigit()},
        key=int,
    )

    blob = bytearray()
    index = bytearray()
    both = 0

    for gene_id in gene_ids:
        iso = isoforms.get(gene_id, "")
        rng = ranges.get(gene_id, "")
        if iso and rng:
            both += 1

        iso_bytes = (iso + "\n").encode()
        iso_offset = len(blob)
        blob += iso_bytes

        rng_bytes = (rng + "\n").encode()
        rng_offset = len(blob)
        blob += rng_bytes

        index += RECORD.pack(
            int(gene_id), iso_offset, len(iso_bytes) - 1, rng_offset, len(rng_bytes) - 1
        )

    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "annotations.bin").write_bytes(bytes(blob))
    (out_dir / "annotations.idx").write_bytes(bytes(index))
    (out_dir / "manifest.json").write_text(
        json.dumps(
            {
                "formatVersion": 1,
                "note": (
                    "Isoform exon structures and intersection ranges, in the legacy "
                    "encoding, for the transcript diagram. These are the 2015 "
                    "annotations the predictions were computed against, not current "
                    "RefSeq."
                ),
                "source": [REFERENCE, INTERSECTION],
                "count": len(gene_ids),
                "withBoth": both,
                "blobBytes": len(blob),
                "indexBytes": len(index),
                "recordBytes": RECORD.size,
                "blob": "annotations.bin",
                "index": "annotations.idx",
            },
            indent=2,
        )
    )

    print(f"genes            : {len(gene_ids):,}  ({both:,} with both records)")
    print(f"annotations.bin  : {len(blob)/1e6:.1f} MB")
    print(f"annotations.idx  : {len(index)/1e6:.2f} MB  (fetched only when a diagram is opened)")
    print(f"mean per gene    : {len(blob)/max(1,len(gene_ids)):.0f} bytes -> one range request")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
