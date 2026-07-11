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
trimming, so codeml can skip the gene gracefully downstream. That "all
positions trimmed" outcome is a normal, expected result for some genes and
is not treated as an error. To keep it distinguishable from a genuine tool
failure (missing binary, crash, bad install, etc.) that happens to also
leave no output behind, the tool's exit code and full stdout/stderr are
always written to the log next to the sentinel message — this is a pure
logging addition and does not change what happens to any gene's output
either way.

Unwanted Gblocks sidecar files (.ps, .html, and any leftover empty -gb1.fa)
are moved — not deleted — into a per-gene folder next to the alignment
(e.g. aligns/gene1/) so nothing is silently discarded, while the wanted
trimmed-alignment file stays at the Snakemake-declared output path used by
downstream rules.
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

# Per-gene folder for unwanted sidecar files, e.g. aligns/gene1/
# (input_aln.stem strips just the last suffix: "gene1.pal2nal" -> "gene1")
gene_sidecar_dir = input_aln.parent / input_aln.stem


def write_log(msg: str, mode: str = "a") -> None:
    with open(log_path, mode) as f:
        f.write(msg + "\n")


def stash_sidecars(paths: list) -> None:
    """
    Move any of the given files (if they exist) into gene_sidecar_dir instead
    of deleting them, so nothing produced by the trimmer is silently lost.
    """
    existing = [p for p in paths if p.exists()]
    if not existing:
        return
    gene_sidecar_dir.mkdir(parents=True, exist_ok=True)
    for p in existing:
        dest = gene_sidecar_dir / p.name
        p.rename(dest)
        write_log(f"Moved unwanted sidecar file {p.name} → {dest}")


def sentinel_if_empty(tool_name: str, returncode: int) -> bool:
    """
    Create empty sentinel file if output is missing or empty.
    Returns True if sentinel was created (i.e. all positions were trimmed).

    Always logs the tool's exit code alongside the sentinel message. This
    is expected to be non-zero for Gblocks specifically when it trims
    everything (see module docstring) — logging it either way just makes
    a genuinely unexpected failure (wrong exit code for the tool in use,
    or a crash) traceable later without changing what happens to the gene.
    """
    if not output_fa.exists() or output_fa.stat().st_size == 0:
        output_fa.touch()
        write_log(
            f"{tool_name} trimmed all positions for "
            f"{input_aln.name} (exit code {returncode}) — empty sentinel "
            f"created. This gene will be skipped by codeml."
        )
        return True
    return False


# ── Gblocks ───────────────────────────────────────────────────────────────────
if tool == "gblocks":
    write_log("Trimmer: Gblocks", mode="w")
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
    write_log(f"Gblocks exit code: {result.returncode}")
    write_log(result.stdout)
    if result.stderr:
        write_log(result.stderr)

    # Gblocks writes {input}-gb1.fa — move it to the declared output path
    gblocks_out = Path(str(input_aln) + "-gb1.fa")
    sidecars = [Path(str(input_aln) + suffix) for suffix in ("-gb1.ps", "-gb1.html")]

    if gblocks_out.exists() and gblocks_out.stat().st_size > 0:
        gblocks_out.rename(output_fa)
        write_log(f"Gblocks trimmed successfully → {output_fa.name}")
        stash_sidecars(sidecars)
    else:
        # Nothing usable was produced — stash whatever Gblocks did leave
        # behind (including a leftover empty -gb1.fa, if any) rather than
        # deleting it, then fall through to the normal sentinel handling.
        stash_sidecars(sidecars + [gblocks_out])
        sentinel_if_empty("Gblocks", result.returncode)

# ── trimAl ────────────────────────────────────────────────────────────────────
elif tool == "trimal":
    write_log("Trimmer: trimAl", mode="w")
    cmd = [
        "trimal",
        "-in",       str(input_aln),
        "-out",      str(output_fa),
        "-automated1",   # heuristic selection of best trimming method
        "-fasta",        # ensure FASTA output
    ]
    # capture_output=True prevents trimAl's stdout from appearing on the
    # Snakemake terminal — all output goes to the log file only
    result = subprocess.run(cmd, capture_output=True, text=True)
    write_log(f"trimAl exit code: {result.returncode}")
    write_log(result.stdout)
    if result.stderr:
        write_log(result.stderr)

    trimmed_all = sentinel_if_empty("trimAl", result.returncode)
    if not trimmed_all:
        write_log(f"trimAl trimmed successfully → {output_fa.name}")

else:
    raise ValueError(
        f"Unknown trimmer '{tool}'. "
        f"Valid options: gblocks, trimal. "
        f"Set tools.trimmer in config/config.yaml."
    )
