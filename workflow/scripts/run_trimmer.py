#!/usr/bin/env python3
"""
run_trimmer.py — Snakemake script directive wrapper for alignment trimming.

Supports two trimming tools (set via config tools.trimmer):
  gblocks — Gblocks 0.91b (default, conservative codon-aware trimming)
  trimal  — trimAl 1.4    (faster, less aggressive, explicit output path)

Both tools take a codon alignment FASTA and write a trimmed FASTA.
The output filename is always the Snakemake-declared output path regardless
of which tool is used, so downstream rules are unaffected by the choice.

Gblocks quirks handled here:
  - Appends suffix to input filename rather than writing to a declared path
  - Exits non-zero when all positions are trimmed (writes nothing)
  - Writes unwanted .ps and .html sidecar files

trimAl writes directly to --out so none of the above apply.

In both cases an empty sentinel file is created if no positions pass
trimming, so codeml can skip the gene gracefully.
"""

# =============================================================================
# Author:  Suvratha Jayaprasad
# Contact: https://suvrathaprasad.github.io/index.html
# License: CC BY-NC-ND 4.0 — you may use and share with attribution,
#          but you may not modify or use for commercial purposes.
#          © 2025 Suvratha Jayaprasad. All rights reserved.
#          Full terms: https://creativecommons.org/licenses/by-nc-nd/4.0/
# =============================================================================

import subprocess
from pathlib import Path

# ── Snakemake-injected objects ────────────────────────────────────────────────
input_aln = Path(snakemake.input[0])
output_fa = Path(snakemake.output[0])
log_path  = Path(snakemake.log[0])
tool      = snakemake.params.get("trimmer", "gblocks")

log_path.parent.mkdir(parents=True, exist_ok=True)


def write_log(msg: str, mode: str = "a") -> None:
    with open(log_path, mode) as f:
        f.write(msg + "\n")


def sentinel_if_empty(label: str) -> None:
    """Create empty sentinel file if output is missing or empty."""
    if not output_fa.exists() or output_fa.stat().st_size == 0:
        output_fa.touch()
        write_log(f"{label} trimmed all positions — empty sentinel created.")


# ── Gblocks ───────────────────────────────────────────────────────────────────
if tool == "gblocks":
    write_log(f"Trimmer: Gblocks", mode="w")
    cmd = [
        "Gblocks", str(input_aln),
        "-t=c",
        "-e=-gb1",
        "-b4=10",
        "-d=y",
        "-b1=2",
        "-b2=2",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    write_log(result.stdout)
    if result.stderr:
        write_log(result.stderr)

    # Gblocks writes {input}-gb1.fa — move it to the declared output path
    gblocks_out = Path(str(input_aln) + "-gb1.fa")
    if gblocks_out.exists() and gblocks_out.stat().st_size > 0:
        gblocks_out.rename(output_fa)
    else:
        # Also clean up empty sidecar files Gblocks may leave behind
        for suffix in ("-gb1.fa", "-gb1.ps", "-gb1.html"):
            sidecar = Path(str(input_aln) + suffix)
            if sidecar.exists():
                sidecar.unlink()
        sentinel_if_empty("Gblocks")

# ── trimAl ────────────────────────────────────────────────────────────────────
elif tool == "trimal":
    write_log(f"Trimmer: trimAl", mode="w")
    cmd = [
        "trimal",
        "-in",       str(input_aln),
        "-out",      str(output_fa),
        "-automated1",   # heuristic selection of best trimming method
        "-fasta",        # ensure FASTA output
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    write_log(result.stdout)
    if result.stderr:
        write_log(result.stderr)
    sentinel_if_empty("trimAl")

else:
    raise ValueError(
        f"Unknown trimmer '{tool}'. "
        f"Valid options: gblocks, trimal. "
        f"Set tools.trimmer in config/config.yaml."
    )
