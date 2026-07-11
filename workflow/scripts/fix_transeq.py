#!/usr/bin/env python3
"""
fix_transeq.py — Run EMBOSS transeq and post-process its protein output.

transeq writes sequences across multiple lines with non-standard wrapping
and appends "_1" to every header (the frame number). This script runs
transeq itself (capturing its stdout directly, rather than receiving it
via a pipe) and then:
  1. Joins multi-line sequences into single-line records
  2. Removes the trailing "_1" frame suffix from headers
  3. Drops any empty records (logging how many, and their headers, so a
     downstream "ID not found" error in extract_gene_pair.py can be traced
     back to a translation that silently produced nothing here)

Called by the Snakemake rules translate_query and translate_target.
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
import subprocess
import sys
from pathlib import Path

output_faa = snakemake.output[0]
log_path   = snakemake.log[0]
input_fna  = snakemake.input[0]

Path(output_faa).parent.mkdir(parents=True, exist_ok=True)
Path(log_path).parent.mkdir(parents=True, exist_ok=True)

# Run transeq, capturing stdout
result = subprocess.run(
    ["transeq", input_fna, "stdout", "-trim", "Y", "-clean", "Y"],
    capture_output=True, text=True
)

with open(log_path, "w") as log:
    log.write(result.stderr)

if result.returncode != 0:
    sys.exit(f"transeq failed:\n{result.stderr}")

# Parse and clean the output
records = []
dropped_empty_headers = []
current_header = None
current_seq = []


def _flush():
    """Append the current record if non-empty, else note it as dropped."""
    if current_header is None:
        return
    if current_seq:
        records.append((current_header, "".join(current_seq)))
    else:
        dropped_empty_headers.append(current_header)


for line in result.stdout.splitlines():
    line = line.strip()
    if not line:
        continue
    if line.startswith(">"):
        _flush()
        # Remove trailing "_1" frame suffix added by transeq
        current_header = re.sub(r"_1$", "", line[1:])
        current_seq = []
    else:
        current_seq.append(line)
_flush()

with open(output_faa, "w") as fout:
    for header, seq in records:
        fout.write(f">{header}\n")
        for i in range(0, len(seq), 60):
            fout.write(seq[i:i+60] + "\n")

print(f"[fix_transeq] Written {len(records)} sequences to {output_faa}")
if dropped_empty_headers:
    print(
        f"[fix_transeq] WARNING: dropped {len(dropped_empty_headers)} empty "
        f"record(s) (transeq produced no residues for these — usually a "
        f"CDS that's all-N, all-gap, or too short to translate). "
        f"IDs will NOT appear in {output_faa} and will fail downstream "
        f"'ID not found' checks if they're referenced later: "
        f"{dropped_empty_headers[:10]}"
        f"{' ...' if len(dropped_empty_headers) > 10 else ''}"
    )
if len(records) == 0:
    print(
        f"[fix_transeq] WARNING: {output_faa} contains zero sequences. "
        f"Check {input_fna} and {log_path}."
    )
