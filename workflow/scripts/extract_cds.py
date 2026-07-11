#!/usr/bin/env python3
"""
extract_cds.py — Snakemake script directive for CDS extraction.

Supports two tools (set via config tools.cds_extraction):
  anchorwave — Anchorwave gff2seq (default, whole-genome alignment aware)
  gffread    — gffread -x (faster, simpler, widely used for CDS extraction)

Both take a genome FASTA + GFF and write a CDS nucleotide FASTA.
Output is then passed to clean_fna.py for header sanitisation.
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
import sys
from pathlib import Path

# ── Snakemake-injected objects ────────────────────────────────────────────────
fasta_in  = snakemake.input.fasta
gff_in    = snakemake.input.gff
output    = snakemake.output[0]
log_path  = snakemake.log[0]
tool      = snakemake.params.get("cds_tool", "anchorwave")

Path(output).parent.mkdir(parents=True, exist_ok=True)
Path(log_path).parent.mkdir(parents=True, exist_ok=True)

# ── Anchorwave ────────────────────────────────────────────────────────────────
if tool == "anchorwave":
    cmd = [
        "anchorwave", "gff2seq",
        "-r", fasta_in,
        "-i", gff_in,
        "-o", output,
    ]

# ── gffread ───────────────────────────────────────────────────────────────────
elif tool == "gffread":
    cmd = [
        "gffread",
        gff_in,
        "-g", fasta_in,
        "-x", output,   # -x writes CDS sequences
    ]

else:
    sys.exit(
        f"Unknown cds_extraction tool '{tool}'. "
        f"Valid options: anchorwave, gffread. "
        f"Set tools.cds_extraction in config/config.yaml."
    )

result = subprocess.run(cmd, capture_output=True, text=True)

with open(log_path, "w") as log:
    log.write(f"Tool: {tool}\n")
    log.write(f"Command: {' '.join(cmd)}\n\n")
    log.write(result.stdout)
    if result.stderr:
        log.write(result.stderr)

if result.returncode != 0:
    sys.exit(f"{tool} failed. Check {log_path} for details.")

if not Path(output).exists() or Path(output).stat().st_size == 0:
    sys.exit(
        f"{tool} exited successfully (code 0) but produced no output "
        f"(or an empty file) at {output}. Check {log_path} for details — "
        f"this usually means the GFF and FASTA don't match (wrong "
        f"chromosome/contig names, mismatched assembly versions, etc.)."
    )

print(f"[extract_cds] {tool} completed → {output}")
