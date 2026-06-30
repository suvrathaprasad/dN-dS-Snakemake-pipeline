#!/usr/bin/env python3
"""
extract_gene_pair.py — Extract one gene pair's FAA and FNA sequences.

Uses indexed FASTA lookup (single pass to build index, then seek to offset)
so large proteomes are not loaded entirely into memory for every gene.

Called by the Snakemake rule extract_gene_sequences.
"""

# =============================================================================
# Author:  Suvratha Jayaprasad
# Contact: https://suvrathaprasad.github.io/index.html
# License: CC BY-NC-ND 4.0 — you may use and share with attribution,
#          but you may not modify or use for commercial purposes.
#          © 2025 Suvratha Jayaprasad. All rights reserved.
#          Full terms: https://creativecommons.org/licenses/by-nc-nd/4.0/
# =============================================================================

import argparse
import sys
from pathlib import Path


def build_index(path: str) -> dict:
    """
    Single-pass index: {sequence_id: (byte_offset_of_header, byte_offset_after_header)}
    Allows random access into large FASTA files without loading them into memory.
    """
    index = {}
    with open(path, "rb") as fh:
        while True:
            offset = fh.tell()
            line = fh.readline()
            if not line:
                break
            if line.startswith(b">"):
                seq_id = line[1:].split()[0].decode()
                index[seq_id] = offset
    return index


def fetch_sequence(path: str, offset: int) -> tuple:
    """
    Seek to offset in FASTA file and read header + sequence.
    Returns (header_without_gt, sequence_string).
    """
    with open(path, "rb") as fh:
        fh.seek(offset)
        header_line = fh.readline().decode().rstrip()
        header = header_line[1:]  # strip leading >
        seq_parts = []
        while True:
            line = fh.readline()
            if not line or line.startswith(b">"):
                break
            seq_parts.append(line.decode().rstrip())
    return header, "".join(seq_parts)


def write_fasta(path: str, records: list) -> None:
    with open(path, "w") as fh:
        for header, seq in records:
            fh.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i+60] + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Extract FAA and FNA sequences for one RBH gene pair."
    )
    parser.add_argument("--pairs",   required=True, help="RBH TSV (query TAB target ...)")
    parser.add_argument("--idx",     required=True, type=int, help="1-based row index")
    parser.add_argument("--qfna",    required=True, help="Query CDS nucleotide FASTA")
    parser.add_argument("--tfna",    required=True, help="Target CDS nucleotide FASTA")
    parser.add_argument("--qfaa",    required=True, help="Query protein FASTA")
    parser.add_argument("--tfaa",    required=True, help="Target protein FASTA")
    parser.add_argument("--out_fna", required=True, help="Output nucleotide FASTA")
    parser.add_argument("--out_faa", required=True, help="Output protein FASTA")
    args = parser.parse_args()

    # ── Read the gene pair at 1-based index (skip header) ────────────────────
    with open(args.pairs) as fh:
        lines = [l.strip() for l in fh
                 if l.strip() and not l.startswith("Query")]
    if args.idx < 1 or args.idx > len(lines):
        sys.exit(f"Index {args.idx} out of range (1–{len(lines)})")
    query_id, target_id = lines[args.idx - 1].split("\t")[:2]

    # ── Build indexes (fast single pass per file) ─────────────────────────────
    idx_qfaa = build_index(args.qfaa)
    idx_tfaa = build_index(args.tfaa)
    idx_qfna = build_index(args.qfna)
    idx_tfna = build_index(args.tfna)

    # ── Validate IDs exist ────────────────────────────────────────────────────
    for label, store, key in [
        ("query FAA",   idx_qfaa, query_id),
        ("target FAA",  idx_tfaa, target_id),
        ("query FNA",   idx_qfna, query_id),
        ("target FNA",  idx_tfna, target_id),
    ]:
        if key not in store:
            sys.exit(f"ID '{key}' not found in {label}")

    # ── Fetch and write ───────────────────────────────────────────────────────
    Path(args.out_faa).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out_fna).parent.mkdir(parents=True, exist_ok=True)

    write_fasta(args.out_faa, [
        fetch_sequence(args.qfaa, idx_qfaa[query_id]),
        fetch_sequence(args.tfaa, idx_tfaa[target_id]),
    ])
    write_fasta(args.out_fna, [
        fetch_sequence(args.qfna, idx_qfna[query_id]),
        fetch_sequence(args.tfna, idx_tfna[target_id]),
    ])


if __name__ == "__main__":
    main()
