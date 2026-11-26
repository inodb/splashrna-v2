"""
Build the static gene index: Entrez gene ID -> isoform-intersection transcript.

    python3 tools/build_genes.py /path/to/legacy/dat public/genes

The source is `intersection_20151217_seq_updated.txt` on the legacy server, a 112 MB
table of `geneId <TAB> strand <TAB> sequence` for 39,760 genes. The sequence is the
intersection of a gene's RefSeq isoforms -- the region common to all of them -- which
is what the original server designed against. `*` marks a boundary between
non-contiguous blocks; a 22-mer spanning one is not a real sequence and is skipped,
exactly as `entrezID_prepare()` did.

Output
------
`genes.bin`  Sequences packed at 4 bits per symbol (A C G T * N), two per byte, each
  gene starting on a byte boundary so it can be fetched with a single Range request.

`genes.idx`  One 16-byte record per gene, sorted by gene ID for binary search:

      uint32  Entrez gene ID
      uint32  byte offset into genes.bin
      uint32  length in symbols
      uint8   strand, 0 for + and 1 for -
      3 bytes padding

`manifest.json`  Counts and file sizes.

Clients download only `genes.idx` (~620 KB) and Range-request the one gene they need,
so the size of `genes.bin` costs the user nothing.
"""

from __future__ import annotations

import json
import pathlib
import struct
import sys

SYMBOLS = "ACGT*N"  # N occurs once in the whole table, in gene 73047
CODE = {c: i for i, c in enumerate(SYMBOLS)}
SOURCE_NAME = "intersection_20151217_seq_updated.txt"


def pack(sequence: str) -> bytes:
    """Two symbols per byte, high nibble first."""
    out = bytearray((len(sequence) + 1) // 2)
    for i, char in enumerate(sequence):
        value = CODE[char]
        if i % 2 == 0:
            out[i // 2] = value << 4
        else:
            out[i // 2] |= value
    return bytes(out)


def main() -> int:
    if len(sys.argv) < 2:
        print(__doc__)
        return 2

    source = pathlib.Path(sys.argv[1]) / SOURCE_NAME
    out_dir = pathlib.Path(sys.argv[2] if len(sys.argv) > 2 else "public/genes")
    if not source.exists():
        print(f"missing: {source}", file=sys.stderr)
        return 1

    records = []
    skipped = 0
    alphabet: set[str] = set()

    with source.open() as handle:
        for line in handle:
            parts = line.rstrip("\n").split()
            if len(parts) != 3:
                skipped += 1
                continue
            gene_id, strand, sequence = parts
            sequence = sequence.upper()
            alphabet.update(sequence)
            if not gene_id.isdigit() or sequence in ("", "NONE"):
                skipped += 1
                continue
            if not set(sequence) <= set(SYMBOLS):
                skipped += 1  # unexpected symbols: leave it out rather than guess
                continue
            records.append((int(gene_id), strand, sequence))

    records.sort(key=lambda r: r[0])

    blob = bytearray()
    index = bytearray()
    for gene_id, strand, sequence in records:
        index += struct.pack("<IIIBxxx", gene_id, len(blob), len(sequence), 1 if strand == "-" else 0)
        blob += pack(sequence)

    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "genes.bin").write_bytes(bytes(blob))
    (out_dir / "genes.idx").write_bytes(bytes(index))
    (out_dir / "manifest.json").write_text(
        json.dumps(
            {
                "formatVersion": 1,
                "source": SOURCE_NAME,
                "note": (
                    "Entrez gene ID to isoform-intersection transcript. '*' marks a "
                    "boundary between non-contiguous blocks; 22-mers spanning one are "
                    "not valid guides."
                ),
                "symbols": SYMBOLS,
                "count": len(records),
                "blobBytes": len(blob),
                "indexBytes": len(index),
                "blob": "genes.bin",
                "index": "genes.idx",
            },
            indent=2,
        )
    )

    lengths = [len(r[2]) for r in records]
    print(f"genes      : {len(records):,} ({skipped} skipped)")
    print(f"alphabet   : {''.join(sorted(alphabet))}")
    print(f"genes.bin  : {len(blob)/1e6:.1f} MB")
    print(f"genes.idx  : {len(index)/1e3:.0f} KB  (downloaded up front)")
    print(f"sequence   : min {min(lengths):,}  median {sorted(lengths)[len(lengths)//2]:,}  max {max(lengths):,}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
