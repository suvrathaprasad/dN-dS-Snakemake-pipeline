#!/usr/bin/env python3
"""
check_pseudogenes.py — Sequence-level pseudogene evidence for dN/dS candidates.

For each gene pair in genes_degenerate.tsv (ω ≥ 1), checks the CDS nucleotide
sequences for two classic hallmarks of pseudogenisation:

  1. Premature stop codons — a stop codon (TAA, TAG, TGA) appearing before
     the expected final codon position, indicating a disrupted reading frame.

  2. Frameshifts — CDS length not divisible by 3, indicating an insertion
     or deletion that shifts the reading frame.

Sequences are read from the per-gene-pair FNA files in intermediate/input/
(e.g. gene1.fna, gene2.fna ...) where each file contains exactly two
sequences — one per species — with headers matching the gene IDs in
genes_degenerate.tsv. This is the correct source because:

  - The headers are guaranteed to match genes_degenerate.tsv exactly
  - The sequences are already trimmed to the gene of interest
  - The files exist regardless of which input mode (A, B, C) was used

Each gene is checked independently in both species. Results are written to
genes_degenerate_annotated.tsv and pseudogene_evidence.pdf.

Note: this is a first-pass sequence screen, not a replacement for dedicated
pseudogene detection pipelines. Treat strong evidence genes as priority
candidates for manual inspection or follow-up with dedicated tools.
"""

# =============================================================================
# Author:  Suvratha Jayaprasad
# Contact: https://suvrathaprasad.github.io/index.html
# License: CC BY-NC-ND 4.0 — you may use and share with attribution,
#          but you may not modify or use for commercial purposes.
#          © 2025 Suvratha Jayaprasad. All rights reserved.
#          Full terms: https://creativecommons.org/licenses/by-nc-nd/4.0/
# =============================================================================

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

# ── Snakemake-injected objects ────────────────────────────────────────────────
degenerate_tsv = snakemake.input.degenerate
input_dir      = Path(snakemake.params.input_dir)  # intermediate/input/
out_tsv        = snakemake.output.annotated_tsv
out_pdf        = snakemake.output.evidence_pdf
log_path       = Path(snakemake.log[0])

log_path.parent.mkdir(parents=True, exist_ok=True)
Path(out_tsv).parent.mkdir(parents=True, exist_ok=True)
Path(out_pdf).parent.mkdir(parents=True, exist_ok=True)

_log_fh = open(log_path, "w")

def log(msg: str) -> None:
    _log_fh.write(msg + "\n")
    _log_fh.flush()
    print(msg)

# ── Codon table ───────────────────────────────────────────────────────────────
STOP_CODONS = {"TAA", "TAG", "TGA"}


def parse_fasta_file(path: Path) -> dict:
    """
    Parse a two-sequence FASTA file (one gene pair).
    Returns {header: sequence} for both sequences in the file.
    """
    seqs = {}
    curr, parts = None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if curr is not None:
                    seqs[curr] = "".join(parts)
                curr = line[1:].split()[0]
                parts = []
            else:
                parts.append(line.upper())
        if curr is not None:
            seqs[curr] = "".join(parts)
    return seqs


def check_frameshift(seq: str) -> bool:
    """True if sequence length is not divisible by 3."""
    return len(seq) % 3 != 0


def check_premature_stop(seq: str) -> bool:
    """
    True if a stop codon appears before the final codon.
    The last codon is excluded (legitimate terminal stop).
    Only complete codons are considered.
    """
    n_complete = (len(seq) // 3) * 3
    codons = [seq[i:i+3] for i in range(0, n_complete, 3)]
    if not codons:
        return False
    # All codons except the last
    for codon in codons[:-1]:
        if codon in STOP_CODONS:
            return True
    return False


def evidence_level(ps_q: bool, fs_q: bool, ps_t: bool, fs_t: bool) -> str:
    """
    strong — at least one hallmark detected in either species
    weak   — ω ≥ 1 only, no sequence-level evidence
    """
    if ps_q or fs_q or ps_t or fs_t:
        return "strong"
    return "weak"


# ── Load degenerate candidates ────────────────────────────────────────────────
log(f"Loading degenerate candidates: {degenerate_tsv}")
df = pd.read_csv(degenerate_tsv, sep="\t")
log(f"  {len(df)} gene pairs to check")

if len(df) == 0:
    log("No degenerate candidates found. Writing empty outputs.")
    df_out = df.copy()
    for col in ["premature_stop_query", "premature_stop_target",
                "frameshift_query", "frameshift_target", "pseudogene_evidence"]:
        df_out[col] = []
    df_out.to_csv(out_tsv, sep="\t", index=False)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.text(0.5, 0.5, "No degenerate gene pairs (ω ≥ 1) found.",
            ha="center", va="center", fontsize=12, color="#888888",
            transform=ax.transAxes)
    ax.set_axis_off()
    fig.savefig(out_pdf, format="pdf", bbox_inches="tight")
    plt.close()
    _log_fh.close()
    sys.exit(0)

# ── Build a map from gene_name → fna file path ────────────────────────────────
# intermediate/input/ contains gene1.fna, gene2.fna, ... Each file has
# exactly two sequences: one query and one target for that gene pair.
log(f"\nScanning input directory: {input_dir}")
fna_files = sorted(input_dir.glob("*.fna"))
log(f"  Found {len(fna_files)} FNA files")

# Map gene ID → (gene_file_path) by loading headers from all fna files
# This is the correct approach: match on the actual sequence headers
# which are guaranteed to match the gene IDs in genes_degenerate.tsv
header_to_file = {}
for fna_path in fna_files:
    seqs = parse_fasta_file(fna_path)
    for header in seqs:
        header_to_file[header] = (fna_path, seqs)

log(f"  Indexed {len(header_to_file)} sequence headers")

# ── Check each candidate ──────────────────────────────────────────────────────
log("\nChecking for premature stop codons and frameshifts...")

results = []
not_found = 0

for _, row in df.iterrows():
    qid = row["Gene_query"]
    tid = row["Gene_target"]

    q_entry = header_to_file.get(qid)
    t_entry = header_to_file.get(tid)

    if q_entry is None or t_entry is None:
        missing = []
        if q_entry is None:
            missing.append(qid)
        if t_entry is None:
            missing.append(tid)
        log(f"  WARNING: sequence(s) not found: {', '.join(missing)}")
        not_found += 1
        results.append({
            "premature_stop_query":  "unknown",
            "premature_stop_target": "unknown",
            "frameshift_query":      "unknown",
            "frameshift_target":     "unknown",
            "pseudogene_evidence":   "unknown",
        })
        continue

    _, q_seqs = q_entry
    _, t_seqs = t_entry
    q_seq = q_seqs[qid]
    t_seq = t_seqs[tid]

    ps_q = check_premature_stop(q_seq)
    fs_q = check_frameshift(q_seq)
    ps_t = check_premature_stop(t_seq)
    fs_t = check_frameshift(t_seq)
    ev   = evidence_level(ps_q, fs_q, ps_t, fs_t)

    results.append({
        "premature_stop_query":  "yes" if ps_q else "no",
        "premature_stop_target": "yes" if ps_t else "no",
        "frameshift_query":      "yes" if fs_q else "no",
        "frameshift_target":     "yes" if fs_t else "no",
        "pseudogene_evidence":   ev,
    })

# ── Write annotated TSV ───────────────────────────────────────────────────────
annot  = pd.DataFrame(results)
df_out = pd.concat([df.reset_index(drop=True), annot], axis=1)
df_out.to_csv(out_tsv, sep="\t", index=False)
log(f"\nWritten: {out_tsv}")

# ── Summary ───────────────────────────────────────────────────────────────────
n_strong  = (annot["pseudogene_evidence"] == "strong").sum()
n_weak    = (annot["pseudogene_evidence"] == "weak").sum()
n_unknown = (annot["pseudogene_evidence"] == "unknown").sum()

log(f"\nPseudogene evidence summary:")
log(f"  Strong (stop codon or frameshift detected): {n_strong}")
log(f"  Weak   (ω ≥ 1 only, no sequence evidence): {n_weak}")
if n_unknown > 0:
    log(f"  Unknown (sequence not found):              {n_unknown}")

n_ps_q = (annot["premature_stop_query"]  == "yes").sum()
n_ps_t = (annot["premature_stop_target"] == "yes").sum()
n_fs_q = (annot["frameshift_query"]      == "yes").sum()
n_fs_t = (annot["frameshift_target"]     == "yes").sum()
log(f"\n  Premature stop — query:  {n_ps_q}")
log(f"  Premature stop — target: {n_ps_t}")
log(f"  Frameshift     — query:  {n_fs_q}")
log(f"  Frameshift     — target: {n_fs_t}")

# ── Plot ──────────────────────────────────────────────────────────────────────
log("\nGenerating pseudogene evidence plot...")

COL_STRONG  = "#F44336"
COL_WEAK    = "#FF9800"
COL_UNKNOWN = "#9E9E9E"

fig, axes = plt.subplots(1, 2, figsize=(13, 6),
                          gridspec_kw={"wspace": 0.4})
fig.suptitle("Pseudogene Evidence for Degenerate Candidates (ω ≥ 1)",
             fontsize=14, fontweight="bold", y=1.02)

# ── Left: evidence level bar chart ───────────────────────────────────────────
ax = axes[0]
categories = ["Strong\n(stop or frameshift)", "Weak\n(ω ≥ 1 only)"]
counts     = [int(n_strong), int(n_weak)]
colours    = [COL_STRONG, COL_WEAK]
if n_unknown > 0:
    categories.append("Unknown\n(seq not found)")
    counts.append(int(n_unknown))
    colours.append(COL_UNKNOWN)

total = sum(counts)
bars = ax.bar(categories, counts, color=colours,
              edgecolor="white", linewidth=0.8, width=0.5, zorder=3)
for bar, count in zip(bars, counts):
    pct = count / total * 100 if total > 0 else 0
    ax.text(bar.get_x() + bar.get_width() / 2,
            bar.get_height() + max(counts) * 0.01,
            f"{count}\n({pct:.1f}%)",
            ha="center", va="bottom", fontsize=10, fontweight="bold")

ax.set_ylabel("Number of gene pairs", fontsize=12)
ax.set_title("Evidence level", fontsize=12, fontweight="bold")
ax.set_ylim(top=max(counts) * 1.25 if max(counts) > 0 else 1)
ax.spines[["top", "right"]].set_visible(False)
ax.yaxis.grid(True, linestyle="--", alpha=0.4, zorder=0)
ax.set_axisbelow(True)
ax.set_facecolor("#fafafa")
ax.text(0.98, 0.97, f"Total degenerate candidates: {total}",
        transform=ax.transAxes, fontsize=9,
        ha="right", va="top", color="#555555")

# ── Right: hallmark type breakdown ───────────────────────────────────────────
ax = axes[1]
hallmarks = [
    ("Premature stop\n(query)",  int(n_ps_q), "#EF5350"),
    ("Premature stop\n(target)", int(n_ps_t), "#E57373"),
    ("Frameshift\n(query)",      int(n_fs_q), "#7B1FA2"),
    ("Frameshift\n(target)",     int(n_fs_t), "#BA68C8"),
]
h_labels  = [h[0] for h in hallmarks]
h_counts  = [h[1] for h in hallmarks]
h_colours = [h[2] for h in hallmarks]
ymax = max(h_counts + [1]) * 1.25

bars2 = ax.bar(h_labels, h_counts, color=h_colours,
               edgecolor="white", linewidth=0.8, width=0.5, zorder=3)
for bar, count in zip(bars2, h_counts):
    ax.text(bar.get_x() + bar.get_width() / 2,
            bar.get_height() + ymax * 0.01,
            str(count),
            ha="center", va="bottom", fontsize=10, fontweight="bold")

ax.set_ylabel("Number of gene pairs", fontsize=12)
ax.set_title("Hallmark type breakdown\n(strong evidence genes only)",
             fontsize=12, fontweight="bold")
ax.set_ylim(top=ymax)
ax.spines[["top", "right"]].set_visible(False)
ax.yaxis.grid(True, linestyle="--", alpha=0.4, zorder=0)
ax.set_axisbelow(True)
ax.set_facecolor("#fafafa")
ax.text(0.98, 0.97,
        "Note: one gene pair may show\nmultiple hallmarks",
        transform=ax.transAxes, fontsize=8,
        ha="right", va="top", color="#777777", style="italic")

fig.text(
    0.5, -0.04,
    "Strong evidence = premature stop codon or frameshift detected in the CDS of either species.\n"
    "This is a first-pass screen. Confirm candidates with a dedicated pseudogene tool "
    "before drawing biological conclusions.",
    ha="center", fontsize=8, color="#666666", style="italic"
)

plt.tight_layout()
fig.savefig(out_pdf, format="pdf", bbox_inches="tight", dpi=300)
plt.close()
log(f"Written: {out_pdf}")

log("\nPseudogene evidence check complete.")
_log_fh.close()
