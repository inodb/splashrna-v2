"""
Build the static gene index: Entrez gene ID -> isoform-intersection transcript,
CDS intersection, and the polyadenylation cut point.

    python3 tools/build_genes.py /path/to/legacy/dat public/genes

Sources on the legacy server:
    intersection_20151217_seq_updated.txt       39,760 genes: strand + intersection
    CDSintersection_seq_20151217_updated.txt    CDS-only intersection
    intersection_20151217_updated.txt           genomic coordinates of the intersection
    pA_20160114.bed                             polyadenylation sites

The polyA cut
-------------
The original server trimmed each transcript at its annotated cleavage site before
designing, and had that switched ON by default, so it never proposed guides against
3' sequence absent from the mature mRNA. Reproducing `polyaCut()`:

    position = genomic_positions.index(polyA_coordinate)
    kept     = sequence[:position]                     on the + strand
    kept     = sequence[:len(sequence) - position]     on the - strand

Both branches take a *prefix* of the transcript, so the whole thing reduces to one
integer per gene — computed here rather than shipping genomic coordinate lists.

When the cut cannot be resolved (no polyA annotation, or the coordinate is not inside
the intersection) the original fell back to the CDS intersection, or to the untrimmed
transcript when there was no CDS either. `apaCut = -1` marks that case; the client
then uses the CDS sequence if one exists.

Coverage: 22,513 of 39,760 genes resolve a cut; the remaining 17,247 fall back.

Output
------
`genes.bin`  intersection sequences, 4 bits per symbol (A C G T * N), two per byte
`cds.bin`    CDS intersection sequences, same encoding
`genes.idx`  one 32-byte record per gene, sorted by gene ID:

      uint32  Entrez gene ID
      uint32  byte offset into genes.bin
      uint32  intersection length in symbols
      int32   polyA cut index, or -1
      uint32  byte offset into cds.bin
      uint32  CDS length in symbols, or 0
      uint8   strand, 0 for + and 1 for -
      7 bytes padding
"""

from __future__ import annotations

import json
import pathlib
import struct
import sys

SYMBOLS = "ACGT*N"  # N occurs once in the whole table, in gene 73047
CODE = {c: i for i, c in enumerate(SYMBOLS)}
RECORD = struct.Struct("<IIIiIIB7x")

SEQ_FILE = "intersection_20151217_seq_updated.txt"
CDS_FILE = "CDSintersection_seq_20151217_updated.txt"
POS_FILE = "intersection_20151217_updated.txt"
PA_FILE = "pA_20160114.bed"


def pack(sequence: str) -> bytes:
    out = bytearray((len(sequence) + 1) // 2)
    for i, char in enumerate(sequence):
        if i % 2 == 0:
            out[i // 2] = CODE[char] << 4
        else:
            out[i // 2] |= CODE[char]
    return bytes(out)


def read_sequences(path: pathlib.Path) -> dict[str, tuple[str, str]]:
    """geneId -> (strand, sequence), skipping rows with unexpected symbols."""
    table: dict[str, tuple[str, str]] = {}
    for line in path.open():
        parts = line.rstrip("\n").split()
        if len(parts) != 3 or not parts[0].isdigit():
            continue
        sequence = parts[2].upper()
        if sequence and sequence != "NONE" and set(sequence) <= set(SYMBOLS):
            table[parts[0]] = (parts[1], sequence)
    return table


def expand_positions(ranges: str) -> list[int] | None:
    """"start,end*;start,end" -> every genomic coordinate it covers."""
    out: list[int] = []
    for item in ranges.strip(";").split(";"):
        parts = item.replace(":", ",").split(",")
        if len(parts) == 1 or parts[1] == ".":
            return None
        out.extend(range(int(parts[0]), int(parts[1].rstrip("*")) + 1))
    return out


def polya_cuts(dat: pathlib.Path, sequences: dict[str, tuple[str, str]]) -> dict[str, int]:
    """geneId -> prefix length to keep. Absent when the cut cannot be resolved."""
    sites: dict[str, str] = {}
    for line in (dat / PA_FILE).open():
        parts = line.split()
        if len(parts) >= 4:
            sites[parts[3]] = parts[1]

    ranges: dict[str, str] = {}
    for line in (dat / POS_FILE).open():
        parts = line.strip().split()
        if len(parts) == 4:
            ranges[parts[0]] = parts[3]

    cuts: dict[str, int] = {}
    for gene_id, (strand, sequence) in sequences.items():
        if gene_id not in sites or gene_id not in ranges:
            continue
        positions = expand_positions(ranges[gene_id])
        if positions is None:
            continue
        try:
            index = positions.index(int(sites[gene_id]))
        except ValueError:
            continue
        cuts[gene_id] = index if strand == "+" else len(sequence) - index
    return cuts


def main() -> int:
    if len(sys.argv) < 2:
        print(__doc__)
        return 2

    dat = pathlib.Path(sys.argv[1])
    out_dir = pathlib.Path(sys.argv[2] if len(sys.argv) > 2 else "public/genes")
    for name in (SEQ_FILE, CDS_FILE, POS_FILE, PA_FILE):
        if not (dat / name).exists():
            print(f"missing: {dat / name}", file=sys.stderr)
            return 1

    sequences = read_sequences(dat / SEQ_FILE)
    cds = read_sequences(dat / CDS_FILE)
    cuts = polya_cuts(dat, sequences)

    blob = bytearray()
    cds_blob = bytearray()
    index = bytearray()

    for gene_id in sorted(sequences, key=int):
        strand, sequence = sequences[gene_id]
        cds_sequence = cds.get(gene_id, ("", ""))[1]

        seq_offset = len(blob)
        blob += pack(sequence)
        cds_offset = len(cds_blob)
        if cds_sequence:
            cds_blob += pack(cds_sequence)

        index += RECORD.pack(
            int(gene_id),
            seq_offset,
            len(sequence),
            cuts.get(gene_id, -1),
            cds_offset,
            len(cds_sequence),
            1 if strand == "-" else 0,
        )

    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "genes.bin").write_bytes(bytes(blob))
    (out_dir / "cds.bin").write_bytes(bytes(cds_blob))
    (out_dir / "genes.idx").write_bytes(bytes(index))
    (out_dir / "manifest.json").write_text(
        json.dumps(
            {
                "formatVersion": 2,
                "source": SEQ_FILE,
                "note": (
                    "Entrez gene ID to isoform-intersection transcript, CDS intersection "
                    "and polyA cut index. '*' marks a boundary between non-contiguous "
                    "blocks; 22-mers spanning one are not valid guides. apaCut is the "
                    "prefix length to keep, or -1 when the cut could not be resolved."
                ),
                "symbols": SYMBOLS,
                "count": len(sequences),
                "withPolyACut": len(cuts),
                "withCds": sum(1 for g in sequences if cds.get(g, ("", ""))[1]),
                "blobBytes": len(blob),
                "cdsBytes": len(cds_blob),
                "indexBytes": len(index),
                "recordBytes": RECORD.size,
                "blob": "genes.bin",
                "cds": "cds.bin",
                "index": "genes.idx",
            },
            indent=2,
        )
    )

    print(f"genes        : {len(sequences):,}")
    print(f"  polyA cut  : {len(cuts):,}  ({100*len(cuts)/len(sequences):.0f}%; rest fall back)")
    print(f"  with CDS   : {sum(1 for g in sequences if cds.get(g, ('', ''))[1]):,}")
    print(f"genes.bin    : {len(blob)/1e6:.1f} MB")
    print(f"cds.bin      : {len(cds_blob)/1e6:.1f} MB")
    print(f"genes.idx    : {len(index)/1e6:.2f} MB  (downloaded up front)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
