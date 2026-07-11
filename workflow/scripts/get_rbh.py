#!/usr/bin/env python3
"""
get_rbh.py — Reciprocal Best Hits from two protein FASTA files.

Supports two search methods (set via --method):
  blast   — NCBI BLAST+ blastp (default, gold standard)
  diamond — DIAMOND blastp (100x faster for large proteomes, same output format)

Steps:
  1. Build search database for both FAA files
  2. Search query → target  (forward)
  3. Search target → query  (reverse)
  4. Parse both outputs, extract best hit per query by bitscore
  5. Find reciprocal best hits (RBH)
  6. Filter by e-value and minimum coverage of the shorter sequence
  7. Write TSV: Query TAB Target TAB eValue TAB bitScore TAB pident

Usage:
  python3 get_rbh.py \
      --query      species1.faa \
      --target     species2.faa \
      --out        rbh_pairs.tsv \
      --blast_fwd  blast_fwd.tsv \
      --blast_rev  blast_rev.tsv \
      --db_dir     blast_dbs/ \
      --method     blast \
      --evalue     1e-5 \
      --min_cov    60 \
      --threads    8
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
import os
import subprocess
import sys
from pathlib import Path


# ── Shared output format (identical for BLAST+ and DIAMOND) ──────────────────

OUTFMT = "6 qseqid sseqid evalue bitscore pident qstart qend qlen sstart send slen"


# ── BLAST+ ────────────────────────────────────────────────────────────────────

def make_blastdb(faa: str, db_path: str) -> None:
    cmd = ["makeblastdb", "-in", faa, "-dbtype", "prot", "-out", db_path]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        sys.exit(f"makeblastdb failed for {faa}:\n{result.stderr}")


def run_blastp(query: str, db: str, out: str, evalue: float, threads: int) -> None:
    cmd = [
        "blastp",
        "-query", query,
        "-db", db,
        "-out", out,
        "-evalue", str(evalue),
        "-num_threads", str(threads),
        "-seg", "yes",
        "-soft_masking", "true",
        "-comp_based_stats", "0",
        "-outfmt", OUTFMT,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        sys.exit(f"blastp failed:\n{result.stderr}")


# ── DIAMOND ───────────────────────────────────────────────────────────────────

def make_diamond_db(faa: str, db_path: str, threads: int) -> None:
    cmd = [
        "diamond", "makedb",
        "--in", faa,
        "--db", db_path,
        "--threads", str(threads),
        "--quiet",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        sys.exit(f"diamond makedb failed for {faa}:\n{result.stderr}")


def run_diamond(query: str, db: str, out: str, evalue: float, threads: int) -> None:
    cmd = [
        "diamond", "blastp",
        "--query", query,
        "--db", db,
        "--out", out,
        "--evalue", str(evalue),
        "--threads", str(threads),
        "--outfmt", *OUTFMT.split(),   # diamond takes outfmt tokens as separate args
        "--more-sensitive",            # comparable sensitivity to BLAST+
        "--quiet",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        sys.exit(f"diamond blastp failed:\n{result.stderr}")


# ── Parse search output (identical format for both tools) ────────────────────

def parse_hits(path: str, min_cov: float) -> dict:
    """
    Return best_hits[query] = (target, evalue, bitscore, pident, qcov, scov)
    keeping only the hit with the highest bitscore per query, subject to
    minimum coverage of the shorter sequence.
    """
    best: dict = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 11:
                continue
            (qid, sid, evalue, bitscore, pident,
             qstart, qend, qlen, sstart, send, slen) = parts[:11]

            evalue   = float(evalue)
            bitscore = float(bitscore)
            pident   = float(pident)
            qstart, qend, qlen = int(qstart), int(qend), int(qlen)
            sstart, send, slen = int(sstart), int(send), int(slen)

            qcov = (qend - qstart + 1) / qlen * 100 if qlen > 0 else 0
            scov = (send - sstart + 1) / slen * 100 if slen > 0 else 0

            # Coverage requirement applies to whichever sequence is shorter —
            # that's the one a partial/domain-only alignment would most
            # inflate. Using OR here (pass if *either* side clears min_cov)
            # would let a small alignable domain on a much longer sequence
            # through even when the short sequence itself is barely covered,
            # which is exactly the spurious-hit case this filter exists to
            # catch.
            shorter_cov = scov if slen <= qlen else qcov
            if shorter_cov < min_cov:
                continue

            if qid not in best or bitscore > best[qid][2]:
                best[qid] = (sid, evalue, bitscore, pident, qcov, scov)

    return best


# ── RBH logic ─────────────────────────────────────────────────────────────────

def find_rbh(fwd: dict, rev: dict, max_evalue: float) -> list:
    """
    Return RBH pairs sorted by bitscore descending.
    """
    rbh = []
    for query, (target, evalue, bitscore, pident, qcov, scov) in fwd.items():
        if target in rev and rev[target][0] == query:
            if evalue <= max_evalue:
                rbh.append((query, target, evalue, bitscore, pident))
    rbh.sort(key=lambda x: -x[3])
    return rbh


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Find Reciprocal Best Hits between two protein FASTA files."
    )
    parser.add_argument("--query",     required=True, help="Query protein FASTA")
    parser.add_argument("--target",    required=True, help="Target protein FASTA")
    parser.add_argument("--out",       required=True, help="Output RBH TSV file")
    parser.add_argument("--blast_fwd", required=True, help="Forward search output TSV (query→target)")
    parser.add_argument("--blast_rev", required=True, help="Reverse search output TSV (target→query)")
    parser.add_argument("--db_dir",    required=True, help="Directory to store search databases")
    parser.add_argument("--method",    default="blast", choices=["blast", "diamond"],
                        help="Search method: blast (default) or diamond")
    parser.add_argument("--evalue",    type=float, default=1e-5)
    parser.add_argument("--min_cov",   type=float, default=60.0)
    parser.add_argument("--threads",   type=int,   default=4)
    args = parser.parse_args()

    for f in [args.query, args.target]:
        if not os.path.isfile(f):
            sys.exit(f"Input file not found: {f}")

    for path in [args.out, args.blast_fwd, args.blast_rev]:
        Path(path).parent.mkdir(parents=True, exist_ok=True)

    db_dir = Path(args.db_dir)
    db_dir.mkdir(parents=True, exist_ok=True)

    print(f"[get_rbh] Method: {args.method}")

    if args.method == "blast":
        db_query  = str(db_dir / "query_db")
        db_target = str(db_dir / "target_db")
        print(f"[get_rbh] Building BLAST databases in {db_dir}")
        make_blastdb(args.query,  db_query)
        make_blastdb(args.target, db_target)
        print("[get_rbh] Running blastp: query → target")
        run_blastp(args.query,  db_target, args.blast_fwd, args.evalue, args.threads)
        print("[get_rbh] Running blastp: target → query")
        run_blastp(args.target, db_query,  args.blast_rev, args.evalue, args.threads)

    else:  # diamond
        db_query  = str(db_dir / "query_db")
        db_target = str(db_dir / "target_db")
        print(f"[get_rbh] Building DIAMOND databases in {db_dir}")
        make_diamond_db(args.query,  db_query,  args.threads)
        make_diamond_db(args.target, db_target, args.threads)
        print("[get_rbh] Running diamond blastp: query → target")
        run_diamond(args.query,  db_target, args.blast_fwd, args.evalue, args.threads)
        print("[get_rbh] Running diamond blastp: target → query")
        run_diamond(args.target, db_query,  args.blast_rev, args.evalue, args.threads)

    print("[get_rbh] Parsing search results...")
    fwd = parse_hits(args.blast_fwd, args.min_cov)
    rev = parse_hits(args.blast_rev, args.min_cov)
    print(f"[get_rbh]   Forward best hits: {len(fwd)}")
    print(f"[get_rbh]   Reverse best hits: {len(rev)}")

    rbh = find_rbh(fwd, rev, args.evalue)
    print(f"[get_rbh] Reciprocal best hits found: {len(rbh)}")

    if len(rbh) == 0:
        sys.exit(
            "[get_rbh] ERROR: No reciprocal best hits found. "
            "Check your input files and try relaxing --evalue or --min_cov."
        )

    with open(args.out, "w") as fh:
        fh.write("Query\tTarget\teValue\tbitScore\tpident\n")
        for query, target, evalue, bitscore, pident in rbh:
            fh.write(f"{query}\t{target}\t{evalue}\t{bitscore}\t{pident:.2f}\n")

    print(f"[get_rbh] Written: {args.out}")


if __name__ == "__main__":
    main()
