"""
run_codeml.py — Snakemake script directive wrapper for PAML codeml.

PAML has a hardcoded filename buffer of ~96 characters. Absolute paths
on deep directory trees exceed this limit and produce the cryptic error
"err: option file. add space around the equal sign?"

The fix: run codeml from inside its working directory and use short
relative symlinks for the input alignment and output file so codeml
only ever sees filenames like "input.fa" and "output.txt".
"""

# =============================================================================
# Author:  Suvratha Jayaprasad
# Contact: https://suvrathaprasad.github.io/index.html
# License: CC BY-NC-ND 4.0 — you may use and share with attribution,
#          but you may not modify or use for commercial purposes.
#          © 2025 Suvratha Jayaprasad. All rights reserved.
#          Full terms: https://creativecommons.org/licenses/by-nc-nd/4.0/
# =============================================================================


import os
import shutil
import subprocess
import sys
from pathlib import Path

# ── Resolve all paths to absolute ─────────────────────────────────────────────
aln_path    = Path(snakemake.input.aln).resolve()
ctl_template = Path(snakemake.input.ctl).resolve()
out_txt     = Path(snakemake.output.txt).resolve()
log_path    = Path(snakemake.log[0]).resolve()
workdir     = Path(snakemake.params.workdir).resolve()

# ── Ensure log directory exists ───────────────────────────────────────────────
log_path.parent.mkdir(parents=True, exist_ok=True)

def log(msg):
    with open(log_path, "a") as f:
        f.write(msg + "\n")

# ── Skip empty alignments (Gblocks trimmed everything) ───────────────────────
if aln_path.stat().st_size == 0:
    log("Skipped: Gblocks produced empty alignment")
    out_txt.parent.mkdir(parents=True, exist_ok=True)
    out_txt.write_text("Skipped: Gblocks produced empty alignment\n")
    sys.exit(0)

# ── Clean workdir to remove stale codeml auxiliary files from prior runs ──────
if workdir.exists():
    shutil.rmtree(workdir)
workdir.mkdir(parents=True)

# ── Use short local names inside the workdir to stay under PAML's 96-char limit
# PAML reads seqfile and outfile from the ctl as-is — keep them short.
local_aln = "input.fa"
local_out = "output.txt"

# Symlink the alignment in; copy avoids cross-device issues
shutil.copy(aln_path, workdir / local_aln)

# ── Write a clean codeml.ctl using short local paths, no inline comments ──────
ctl_out = workdir / "codeml.ctl"
with open(ctl_template) as fh:
    lines = fh.readlines()

with open(ctl_out, "w") as fh:
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("seqfile"):
            fh.write(f"seqfile = {local_aln}\n")
        elif stripped.startswith("outfile"):
            fh.write(f"outfile = {local_out}\n")
        else:
            # All other lines (including PAML's own inline '*' comments)
            # are written through unchanged — codeml ignores anything
            # after '*' itself, so no extra stripping is needed here.
            fh.write(line)

log(f"Running codeml in {workdir}")
log(f"  seqfile = {local_aln}  (-> {aln_path})")
log(f"  outfile = {local_out}  (-> {out_txt})")

# ── Run codeml from inside its workdir ───────────────────────────────────────
result = subprocess.run(
    ["codeml", "codeml.ctl"],
    cwd=workdir,
    stdout=subprocess.PIPE,
    stderr=subprocess.STDOUT,
    text=True
)

log(result.stdout)

if result.returncode != 0:
    log(f"ERROR: codeml exited with code {result.returncode}")
    sys.exit(result.returncode)

# ── Copy local output to the declared Snakemake output path ──────────────────
local_out_path = workdir / local_out
if local_out_path.exists():
    shutil.copy(local_out_path, out_txt)
    log(f"Output written to {out_txt}")
else:
    log("ERROR: codeml did not produce output file")
    sys.exit(1)

log("codeml completed successfully")
