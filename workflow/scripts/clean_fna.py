#!/usr/bin/env python3
"""
clean_fna.py — Sanitise CDS FASTA headers after extraction.

Anchorwave appends trailing annotations like " gene_123" to sequence
headers. This script removes them and prepends the species prefix.

gffread produces clean headers already, so only the prefix is added
when gffread is the CDS extraction tool.

Called by the Snakemake rules clean_fna_query and clean_fna_target.
"""

# =============================================================================
# Author:  Suvratha Jayaprasad
# Contact: https://suvrathaprasad.github.io/index.html
# License: CC BY-NC-ND 4.0 — you may use and share with attribution,
#          but you may not modify or use for commercial purposes.
#          © 2025 Suvratha Jayaprasad. All rights reserved.
#          Full terms: https://creativecommons.org/licenses/by-nc-nd/4.0/
# =============================================================================

import re
from pathlib import Path

# ── Snakemake-injected objects ────────────────────────────────────────────────
input_fna  = snakemake.input[0]
output_fna = snakemake.output[0]
prefix     = snakemake.params.prefix
cds_tool   = snakemake.params.get("cds_tool", "anchorwave")

Path(output_fna).parent.mkdir(parents=True, exist_ok=True)

with open(input_fna) as fin, open(output_fna, "w") as fout:
    for line in fin:
        if line.startswith(">"):
            header = line.rstrip()[1:]  # strip leading >
            if cds_tool == "anchorwave":
                # Anchorwave appends " gene_42" — strip it
                header = re.sub(r"\s+gene_\d+\s*$", "", header)
            # gffread headers are already clean — just add the prefix
            fout.write(f">{prefix}_{header}\n")
        else:
            fout.write(line)
