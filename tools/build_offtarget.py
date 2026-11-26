"""
Build the static off-target index from the legacy 22-mer match tables.

    python3 tools/build_offtarget.py /path/to/dat public/offtarget

The source tables live only on the legacy server, at
`.../splashrna/server/server/dat/{hg38,mm10}_22mer_matches.pickle`. They are Python 2
pickles of `{22-mer: [entrez gene id, ...]}` and hold *only* 22-mers that occur in two
or more genes, which is exactly the set a guide must avoid.

    hg38: 1,138,230 keys, 12,937 distinct genes
    mm10:   990,865 keys, 13,588 distinct genes

NOTE ON SPECIES: the original server loaded these into swapped variables --
`hg38` (human) into `self.mouseMatch` and `mm10` (mouse) into `self.humanMatch`
(server.py:662-663), so every species label it reported was wrong. Filtering was
unaffected, since a guide is dropped when either table matches. This builder keys the
output off the assembly name, which is correct.

Output format, per genome
-------------------------
`<genome>-keys.bin`  Blocks of BLOCK_SIZE entries. Each block is self-contained:

      6 bytes   first 22-mer, big-endian, 2 bits per base (A=0 C=1 G=2 T=3)
      varint    gene count for that first key
      then, for each remaining entry in the block:
        varint  delta from the previous key
        varint  gene count

  Self-contained blocks mean a client can fetch one block with a Range request and
  decode it without touching the rest of the file.

`<genome>-keys.dir`  One 14-byte record per block, so a lookup is a binary search in
  memory followed by at most one range request:

      uint32    byte offset of the block within the keys file
      uint32    index of the block's first gene reference in the genes file
      6 bytes   first 22-mer of the block

`<genome>-genes.bin`  The matching genes themselves: uint16 indices into `geneIds`,
  concatenated in key order. A key's slice starts at its block's gene index plus the
  counts of the keys before it in that block, and runs for its own count — so no
  per-key offset table is needed. Fetched only for guides actually shown.

`<genome>-genes.json`  The Entrez gene IDs those uint16 values index.

`manifest.json`  Sizes, counts, and the block size.
"""

from __future__ import annotations

import json
import pathlib
import pickle
import struct
import sys

BASE = {"A": 0, "C": 1, "G": 2, "T": 3}
BLOCK_SIZE = 1024
GENOMES = {"human": "hg38", "mouse": "mm10"}
FORMAT_VERSION = 2


def pack_kmer(seq: str) -> int:
    value = 0
    for char in seq:
        value = value * 4 + BASE[char]
    return value


def varint(value: int) -> bytes:
    out = bytearray()
    while value >= 0x80:
        out.append((value & 0x7F) | 0x80)
        value >>= 7
    out.append(value)
    return bytes(out)


def build_genome(source: pathlib.Path, out_dir: pathlib.Path, genome: str) -> dict:
    with source.open("rb") as handle:
        table = pickle.load(handle, encoding="latin1")

    # Map every Entrez ID to a uint16 index, so the payload is two bytes per reference.
    gene_ids = sorted({int(g) for genes in table.values() for g in genes})
    gene_index = {g: i for i, g in enumerate(gene_ids)}
    if len(gene_ids) > 0xFFFF:
        raise SystemExit(f"{genome}: {len(gene_ids)} genes exceeds uint16")

    entries = []
    skipped = 0
    for kmer, genes in table.items():
        try:
            entries.append((pack_kmer(kmer), [gene_index[int(g)] for g in genes]))
        except KeyError:
            skipped += 1  # ambiguous bases; the scorer cannot handle these anyway
    entries.sort()
    del table

    keys = bytearray()
    directory = bytearray()
    genes_payload = bytearray()

    for start in range(0, len(entries), BLOCK_SIZE):
        block = entries[start : start + BLOCK_SIZE]
        directory += struct.pack(">I", len(keys))
        directory += struct.pack(">I", len(genes_payload) // 2)
        directory += block[0][0].to_bytes(6, "big")

        keys += block[0][0].to_bytes(6, "big")
        keys += varint(len(block[0][1]))
        genes_payload += struct.pack(f">{len(block[0][1])}H", *block[0][1])
        previous = block[0][0]
        for value, genes in block[1:]:
            keys += varint(value - previous)
            keys += varint(len(genes))
            genes_payload += struct.pack(f">{len(genes)}H", *genes)
            previous = value

    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / f"{genome}-keys.bin").write_bytes(bytes(keys))
    (out_dir / f"{genome}-keys.dir").write_bytes(bytes(directory))
    (out_dir / f"{genome}-genes.bin").write_bytes(bytes(genes_payload))
    (out_dir / f"{genome}-genes.json").write_text(json.dumps(gene_ids, separators=(",", ":")))

    blocks = (len(entries) + BLOCK_SIZE - 1) // BLOCK_SIZE
    print(
        f"{genome:6s} {source.name:28s} {len(entries):>9,} keys  "
        f"{len(keys)/1e6:5.2f} MB keys  {len(genes_payload)/1e6:5.2f} MB genes  "
        f"{len(gene_ids):,} ids  {blocks:,} blocks"
        + (f"  ({skipped} skipped)" if skipped else "")
    )

    return {
        "assembly": GENOMES[genome],
        "count": len(entries),
        "blocks": blocks,
        "keysBytes": len(keys),
        "dirBytes": len(directory),
        "keys": f"{genome}-keys.bin",
        "directory": f"{genome}-keys.dir",
        "genes": f"{genome}-genes.bin",
        "geneIds": f"{genome}-genes.json",
        "geneCount": len(gene_ids),
        "geneRefs": len(genes_payload) // 2,
        "dirRecordBytes": 14,
    }


def main() -> int:
    if len(sys.argv) < 2:
        print(__doc__)
        return 2

    dat = pathlib.Path(sys.argv[1])
    out_dir = pathlib.Path(sys.argv[2] if len(sys.argv) > 2 else "public/offtarget")

    manifest = {
        "formatVersion": FORMAT_VERSION,
        "blockSize": BLOCK_SIZE,
        "note": (
            "22-mers occurring in two or more genes. Built from the legacy SplashRNA "
            "match tables; keyed by assembly, correcting the swapped species labels in "
            "the original server."
        ),
        "genomes": {},
    }

    for genome, assembly in GENOMES.items():
        source = dat / f"{assembly}_22mer_matches.pickle"
        if not source.exists():
            print(f"missing: {source}", file=sys.stderr)
            return 1
        manifest["genomes"][genome] = build_genome(source, out_dir, genome)

    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))
    total = sum(g["keysBytes"] + g["dirBytes"] for g in manifest["genomes"].values())
    print(f"wrote {out_dir}/manifest.json — {total/1e6:.1f} MB total")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
