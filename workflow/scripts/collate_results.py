"""
collate_results.py — Snakemake script directive for collating codeml outputs.

Parses all per-gene codeml output files and writes a single TSV.
Skips genes where Gblocks trimmed the entire alignment (sentinel files)
without failing the pipeline.
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
input_files = list(snakemake.input)
output_file = Path(snakemake.output[0])

output_file.parent.mkdir(parents=True, exist_ok=True)


def parse_codeml(path: str):
    """
    Parse a single codeml output file.
    Returns list of (gene1, gene2, t, dN, dS, dNdS) tuples,
    or empty list if the file is a skip sentinel or unparseable.
    """
    try:
        with open(path) as fh:
            text = fh.read()
    except OSError:
        return []

    # Skip sentinel files written when Gblocks trimmed everything
    if text.strip().startswith("Skipped:"):
        return []

    results = []
    status = False
    gene1 = gene2 = None

    for line in text.splitlines():
        if line.startswith("pairwise comparison, codon frequencies"):
            status = True

        if status:
            if line and line[0].isdigit():
                # Gene pair line: "1 (gene1...) vs 2 (gene2...)"
                clean = re.sub(r"[()]", "", line)
                parts = clean.split()
                try:
                    gene1 = parts[1].split("...")[0]
                    gene2 = parts[4].split("...")[0]
                except IndexError:
                    continue

            if line.startswith("t=") and gene1 and gene2:
                clean = re.sub(r"=", "", line)
                tokens = re.split(r"\s+", clean.strip())
                try:
                    t    = tokens[1]
                    dnds = tokens[7]
                    dn   = tokens[9]
                    ds   = tokens[11]
                    results.append((gene1, gene2, t, dn, ds, dnds))
                except IndexError:
                    continue

    return results


# ── Collate ───────────────────────────────────────────────────────────────────
rows = []
skipped = 0

for f in input_files:
    parsed = parse_codeml(f)
    if not parsed:
        skipped += 1
    else:
        rows.extend(parsed)

print(f"[collate_results] Parsed {len(rows)} gene pairs, skipped {skipped} (empty after Gblocks)")

with open(output_file, "w") as fh:
    fh.write("Gene_query\tGene_target\tt\tdN\tdS\tdNdS\n")
    for gene1, gene2, t, dn, ds, dnds in rows:
        fh.write(f"{gene1}\t{gene2}\t{t}\t{dn}\t{ds}\t{dnds}\n")

print(f"[collate_results] Written: {output_file}")
