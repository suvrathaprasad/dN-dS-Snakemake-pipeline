#!/usr/bin/env python3
"""
fix_transeq.py — Post-process EMBOSS transeq protein output.

transeq writes sequences across multiple lines with non-standard wrapping
and appends "_1" to every header (the frame number). This script:
  1. Joins multi-line sequences into single-line records
  2. Removes the trailing "_1" frame suffix from headers
  3. Drops any empty records

Called by the Snakemake rules translate_query and translate_target,
which pipe transeq stdout directly into this script via stdin.
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
import sys

output_faa = snakemake.output[0]
log_path   = snakemake.log[0]
input_fna  = snakemake.input[0]

import subprocess
from pathlib import Path

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
current_header = None
current_seq = []

for line in result.stdout.splitlines():
    line = line.strip()
    if not line:
        continue
    if line.startswith(">"):
        if current_header and current_seq:
            records.append((current_header, "".join(current_seq)))
        # Remove trailing "_1" frame suffix added by transeq
        header = re.sub(r"_1$", "", line[1:])
        current_header = header
        current_seq = []
    else:
        current_seq.append(line)

if current_header and current_seq:
    records.append((current_header, "".join(current_seq)))

with open(output_faa, "w") as fout:
    for header, seq in records:
        fout.write(f">{header}\n")
        for i in range(0, len(seq), 60):
            fout.write(seq[i:i+60] + "\n")

print(f"[fix_transeq] Written {len(records)} sequences to {output_faa}")
