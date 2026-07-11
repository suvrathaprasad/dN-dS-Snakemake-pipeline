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

ANCHORWAVE_SUFFIX = re.compile(r"\s+gene_\d+\s*$")
n_headers = 0
n_unmatched = 0
unmatched_examples = []

with open(input_fna) as fin, open(output_fna, "w") as fout:
    for line in fin:
        if line.startswith(">"):
            n_headers += 1
            header = line.rstrip()[1:]  # strip leading >
            if cds_tool == "anchorwave":
                # Anchorwave appends " gene_42" — strip it
                new_header = ANCHORWAVE_SUFFIX.sub("", header)
                if new_header == header:
                    n_unmatched += 1
                    if len(unmatched_examples) < 5:
                        unmatched_examples.append(header)
                header = new_header
            # gffread headers are already clean — just add the prefix
            fout.write(f">{prefix}_{header}\n")
        else:
            fout.write(line)

# Every header should match the expected Anchorwave suffix pattern; if a
# chunk of them don't, that's most likely a future Anchorwave version
# changing its annotation format rather than genuinely clean headers —
# worth a visible trace back to here rather than a silent no-op per header.
if cds_tool == "anchorwave" and n_unmatched > 0:
    print(
        f"[clean_fna] WARNING: {n_unmatched}/{n_headers} header(s) did not "
        f"match the expected Anchorwave ' gene_<digits>' suffix pattern and "
        f"were left unchanged (only the '{prefix}_' prefix was added). "
        f"This may mean Anchorwave's output format has changed. "
        f"Example unmatched header(s): {unmatched_examples}"
    )
