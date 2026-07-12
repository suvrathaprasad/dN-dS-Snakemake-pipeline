#!/usr/bin/env python3
"""
write_summary.py — One-page PDF run summary for the dN/dS pipeline.

Produces results/run_summary.pdf: a single-page report covering run metadata,
gene pair counts, dN/dS summary statistics, functional classification,
pseudogene evidence, the dS saturation check, and a checklist of every output
file the pipeline is expected to produce.

Layout strategy (guaranteed single-page fit):
  Every line of content is tagged with a "kind" (title / header / body /
  small, etc.), and each kind has a relative weight. Before drawing anything,
  the script sums the weights of all lines it is about to render and divides
  the available vertical space by that total, so the per-line height always
  shrinks or grows to fit exactly — there is no fixed per-line increment that
  can silently overflow the page as the number of gene pairs, output files,
  or evidence categories changes from run to run.

  The file checklist is the one section whose length is genuinely
  unbounded (it grows with the number of declared pipeline outputs), so it
  additionally switches from one column to two once it has more than
  MAX_SINGLE_COLUMN_FILES entries, halving its vertical footprint, and it
  always renders in a smaller font than the rest of the report.

  The footer is drawn separately, pinned near the bottom margin reserved by
  fig.subplots_adjust(bottom=...), and the figure is saved with
  bbox_inches=None so matplotlib respects that margin instead of cropping to
  content (bbox_inches="tight" was what turned the footer into a watermark).
"""

# =============================================================================
# Author:  Suvratha Jayaprasad
# Contact: https://suvrathaprasad.github.io/index.html
# License: CC BY-NC-ND 4.0 — you may use and share with attribution,
#          but you may not modify or use for commercial purposes.
#          © 2025 Suvratha Jayaprasad. All rights reserved.
#          Full terms: https://creativecommons.org/licenses/by-nc-nd/4.0/
# =============================================================================

import json
import math
import shutil
import subprocess
import sys
import textwrap
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ── Snakemake-injected objects ────────────────────────────────────────────────
dnds_tsv      = snakemake.input.dnds_tsv
saturated_tsv = getattr(snakemake.input, "saturated_tsv", None)
annotated_tsv = getattr(snakemake.input, "annotated_tsv", None)
rbh_tsv       = getattr(snakemake.input, "rbh_tsv", None)

out_pdf     = Path(snakemake.output[0])
results_dir = out_pdf.parent
log_path    = Path(snakemake.log[0])

p = snakemake.params
query_prefix     = p.query_prefix
target_prefix    = p.target_prefix
query_mode       = p.query_mode
target_mode      = p.target_mode
cds_tool         = p.cds_tool
search_method    = p.search_method
trimmer          = p.trimmer
evalue           = p.evalue
ds_sat_threshold = float(p.ds_sat_threshold)

out_pdf.parent.mkdir(parents=True, exist_ok=True)
log_path.parent.mkdir(parents=True, exist_ok=True)

# config.yaml's query/target dicts (fasta/gff/faa/fna, per input mode) aren't
# among this rule's declared params, but snakemake.config is always available
# in a script directive regardless of what the rule declares, so we read the
# resolved input file path(s) straight from there rather than needing a
# Snakefile change.
cfg = snakemake.config
query_spec  = cfg.get("query", {})
target_spec = cfg.get("target", {})


def get_git_commit(start_dir) -> str:
    """
    Return the short git commit hash of the pipeline repo, flagged if the
    working tree has uncommitted changes, or a short explanatory string
    if this isn't a git checkout (e.g. a downloaded tarball) or git isn't
    installed. Never raises.
    """
    if shutil.which("git") is None:
        return "unavailable (git not installed)"
    try:
        result = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=str(start_dir), capture_output=True, text=True, timeout=10,
        )
        if result.returncode != 0:
            return "unavailable (not a git checkout)"
        commit = result.stdout.strip()
        dirty = subprocess.run(
            ["git", "status", "--porcelain"],
            cwd=str(start_dir), capture_output=True, text=True, timeout=10,
        )
        suffix = " (+ uncommitted changes)" if dirty.stdout.strip() else ""
        return f"{commit}{suffix}"
    except Exception as exc:
        return f"unavailable ({exc.__class__.__name__})"


def species_path_lines(spec: dict, mode: str) -> list:
    """
    Return [(label, path), ...] describing the actual input file(s) used for
    a species, based on its resolved mode. Falls back to auto-detecting from
    whichever keys are present (mirroring the Mode B > Mode C > Mode A
    priority documented in config.yaml) if mode wasn't resolved.
    """
    mode = (mode or "").upper()
    if mode == "B" or ("faa" in spec and "fna" in spec and mode not in ("A", "C")):
        return [("faa", spec.get("faa", "(not set)")),
                ("fna", spec.get("fna", "(not set)"))]
    if mode == "C" or ("fna" in spec and "faa" not in spec and mode != "A"):
        return [("fna", spec.get("fna", "(not set)"))]
    if mode == "A" or ("fasta" in spec and "gff" in spec):
        return [("fasta", spec.get("fasta", "(not set)")),
                ("gff",   spec.get("gff", "(not set)"))]
    return [("path", "(unresolved — check config.yaml)")]


def wrap_path_line(label: str, path: str, width: int = 100) -> list:
    """Wrap a long 'label: path' line so it never overflows the page width."""
    lead = f"  {label}: "
    full = f"{lead}{path}"
    if len(full) <= width:
        return [full]
    wrapped = textwrap.wrap(str(path), width=max(width - len(lead), 20))
    indent = " " * len(lead)
    return [f"{lead}{wrapped[0]}"] + [f"{indent}{w}" for w in wrapped[1:]]

_log_fh = open(log_path, "w")


def log(msg: str = "") -> None:
    _log_fh.write(msg + "\n")
    _log_fh.flush()
    print(msg)


# ── Load dN/dS results ────────────────────────────────────────────────────────
log(f"Loading dN/dS results: {dnds_tsv}")
df = pd.read_csv(dnds_tsv, sep="\t", comment="#")
for col in ("dN", "dS", "dNdS"):
    if col in df.columns:
        df[col] = pd.to_numeric(df[col], errors="coerce")

total_pairs = len(df)
valid = df[(df["dS"] > 0) & df["dNdS"].notna() & np.isfinite(df["dNdS"])].copy()
undefined = total_pairs - len(valid)
log(f"  {total_pairs} total gene pairs, {len(valid)} valid (dS > 0), {undefined} undefined (dS = 0)")

# ── RBH pairs (upstream filter stage, before codeml/collate) ──────────────────
n_rbh = None
if rbh_tsv and Path(rbh_tsv).exists():
    try:
        with open(rbh_tsv) as fh:
            n_rbh = sum(1 for line in fh if line.strip() and not line.startswith("Query"))
    except OSError as exc:
        log(f"  WARNING: could not read {rbh_tsv}: {exc}")
n_skipped_alignment = (n_rbh - total_pairs) if n_rbh is not None else None

# ── Functional classification (matches plot_results.py thresholds) ───────────
conserved  = valid[valid["dNdS"] < 0.5]
relaxed    = valid[(valid["dNdS"] >= 0.5) & (valid["dNdS"] < 1.0)]
degenerate = valid[valid["dNdS"] >= 1.0]

n_valid      = len(valid)
n_conserved  = len(conserved)
n_relaxed    = len(relaxed)
n_degenerate = len(degenerate)


def pct(n: int, total: int) -> str:
    return f"{(n / total * 100):.1f}%" if total > 0 else "n/a"


# ── dN/dS summary statistics ──────────────────────────────────────────────────
def stat_row(col: str) -> dict:
    s = valid[col]
    if len(s) == 0:
        return {"mean": float("nan"), "median": float("nan"),
                "min": float("nan"), "max": float("nan")}
    return {"mean": s.mean(), "median": s.median(), "min": s.min(), "max": s.max()}


stats = {col: stat_row(col) for col in ("dN", "dS", "dNdS")}

# ── Pseudogene evidence summary ───────────────────────────────────────────────
# genes_degenerate.tsv isn't a direct Snakemake input of this rule (it's
# transitively guaranteed via annotated_tsv -> check_pseudogenes -> plot_results),
# so it's looked up on disk relative to the output directory instead of
# declaring a second dependency the rule doesn't otherwise need.
degenerate_tsv = results_dir / "tables" / "genes_degenerate.tsv"
n_degenerate_declared = None
if degenerate_tsv.exists():
    try:
        n_degenerate_declared = len(pd.read_csv(degenerate_tsv, sep="\t"))
    except Exception:
        n_degenerate_declared = None

pseudo_available = False
n_strong = n_weak = n_unknown = 0
if annotated_tsv and Path(annotated_tsv).exists():
    try:
        annot = pd.read_csv(annotated_tsv, sep="\t")
        if "pseudogene_evidence" in annot.columns:
            pseudo_available = True
            n_strong  = int((annot["pseudogene_evidence"] == "strong").sum())
            n_weak    = int((annot["pseudogene_evidence"] == "weak").sum())
            n_unknown = int((annot["pseudogene_evidence"] == "unknown").sum())
    except Exception as exc:
        log(f"  WARNING: could not read {annotated_tsv}: {exc}")

n_pseudo_total = n_strong + n_weak + n_unknown

# ── dS saturation ──────────────────────────────────────────────────────────────
# ds_sat_threshold already comes from snakemake.params (see top of script).
if saturated_tsv and Path(saturated_tsv).exists():
    try:
        n_saturated = len(pd.read_csv(saturated_tsv, sep="\t"))
    except Exception:
        n_saturated = int((valid["dS"] > ds_sat_threshold).sum())
else:
    n_saturated = int((valid["dS"] > ds_sat_threshold).sum())

# ── Run metadata ───────────────────────────────────────────────────────────────
generated_at = datetime.now().strftime("%Y-%m-%d %H:%M")

# ── Provenance ─────────────────────────────────────────────────────────────────
# Which exact version of the pipeline and its tools produced this report —
# the config-declared tool *choice* (e.g. "blast") is already shown in Run
# Metadata; this is the actually-installed version behind that choice,
# useful for a methods section or reproducing a run later. Repo root is
# taken from Snakemake's own workflow.basedir (the directory containing
# the Snakefile) when available, falling back to the current directory —
# git works from any subdirectory of a checkout regardless.
_repo_dir = getattr(getattr(snakemake, "workflow", None), "basedir", ".")
git_commit = get_git_commit(_repo_dir)

python_version = sys.version.split()[0]

# Snakemake's own version is captured once in the Snakefile itself (the
# master process actually running Snakemake, where it's guaranteed to be
# installed) and passed in as a param — NOT looked up here via subprocess.
# Every other tool version below is read from a shared JSON file that each
# tool's own rule wrote to at the point it actually ran. Neither of these
# can be discovered by asking this rule's own environment: under
# --use-conda, this script runs in its own isolated (typically
# plotting-only) conda env, which was never going to have blastp, mafft,
# Gblocks, codeml, or even Snakemake itself on PATH — those tools were
# only ever reachable from inside their own rule's environment, which no
# longer exists by the time this report is written. See
# record_tool_version.py for the write side of this.
snakemake_version = getattr(snakemake.params, "snakemake_version", "unknown")

tool_versions = {}
_version_file = getattr(snakemake.params, "version_file", None)
if _version_file and Path(_version_file).exists():
    try:
        tool_versions = json.loads(Path(_version_file).read_text())
    except (json.JSONDecodeError, OSError):
        tool_versions = {}


def _tool_version(key: str) -> str:
    return tool_versions.get(key, "unavailable (not recorded — rerun this pipeline "
                                    "with the current Snakefile to populate it)")


search_tool_version = _tool_version(search_method)   # key: "blast" or "diamond"
trimmer_version     = _tool_version(trimmer)          # key: "gblocks" or "trimal"
mafft_version        = _tool_version("mafft")
codeml_version       = _tool_version("codeml")

# ── Output file checklist ─────────────────────────────────────────────────────
# (relative to RESULTS) — mirrors the tree documented in the Snakefile.
CHECKLIST = [
    ("dnds_output.tsv",                      "Final dN/dS table"),
    ("plots/dnds_boxplot.pdf",                "Boxplot of dN, dS, dN/dS"),
    ("plots/dnds_violin.pdf",                 "Violin plot of dN, dS, dN/dS"),
    ("plots/dnds_scatter.pdf",                "dN vs dS scatter"),
    ("plots/functional_summary.pdf",          "Gene category bar chart"),
    ("plots/pseudogene_evidence.pdf",         "Pseudogene evidence chart"),
    ("run_summary.pdf",                       "This report"),
    ("tables/genes_conserved.tsv",            "ω < 0.5"),
    ("tables/genes_relaxed.tsv",              "0.5 ≤ ω < 1"),
    ("tables/genes_degenerate.tsv",           "ω ≥ 1"),
    ("tables/genes_degenerate_annotated.tsv", "ω ≥ 1 + sequence evidence"),
    ("tables/genes_ds_saturated.tsv",         f"dS > {ds_sat_threshold:g}"),
    ("tables/genes_undefined_ds.tsv",         "dS = 0"),
    ("tables/genes_high_dn.tsv",              "dN above median"),
    ("tables/genes_high_ds.tsv",              "dS above median"),
    ("tables/genes_high_dnds.tsv",            "dN/dS above median"),
    ("tables/genes_above_diagonal.tsv",       "dN > dS"),
    ("tables/genes_on_diagonal.tsv",          "dN ≈ dS (±10%)"),
    ("tables/genes_below_diagonal.tsv",       "dN < dS"),
]
checklist_status = []
for rel, desc in CHECKLIST:
    full_path = results_dir / rel
    if full_path == out_pdf:
        # This is the report being written right now — it doesn't exist on
        # disk yet at the moment this check runs (matplotlib hasn't saved
        # the figure), but it deterministically will by the time anyone
        # opens the PDF to read this checklist. Reporting it as missing
        # would be misleading, not accurate.
        ok = True
    else:
        ok = full_path.exists()
    checklist_status.append((rel, desc, ok))
n_present = sum(1 for _, _, ok in checklist_status if ok)
log(f"\nOutput checklist: {n_present}/{len(checklist_status)} files present under {results_dir}")

# =============================================================================
# Layout — build a flat list of (kind, content) lines, then render with a
# dynamically computed line height so everything fits on one page no matter
# how many gene pairs or output files this particular run produced.
# =============================================================================

KIND_WEIGHT = {
    "title":   2.4,
    "meta":    1.05,
    "header":  1.75,
    "sub":     1.2,
    "body":    1.0,
    "small":   0.62,   # file checklist — deliberately compact
    "gap":     0.55,
}
KIND_FONTSIZE = {
    "title":  16,
    "meta":   8.5,
    "header": 11,
    "sub":    9.5,
    "body":   8.5,
    "small":  6.6,
}
MAX_SINGLE_COLUMN_FILES = 10  # switch the checklist to two columns beyond this

lines = []  # each item: (kind, content) where content is str, or for
            # "small" checklist rows a tuple (left_text, right_text_or_None)

lines.append(("title", "dN/dS Pipeline — Run Summary"))
lines.append(("meta", f"Generated {generated_at}    |    Results directory: {results_dir}"))
lines.append(("gap", ""))

# ── Run metadata ──────────────────────────────────────────────────────────
lines.append(("header", "Run Metadata"))
lines.append(("body", f"Query species:  {query_prefix}  (Mode {query_mode})"))
for label, path in species_path_lines(query_spec, query_mode):
    for wrapped in wrap_path_line(label, path):
        lines.append(("body", wrapped))
lines.append(("body", f"Target species: {target_prefix}  (Mode {target_mode})"))
for label, path in species_path_lines(target_spec, target_mode):
    for wrapped in wrap_path_line(label, path):
        lines.append(("body", wrapped))
lines.append(("body", f"CDS extraction: {cds_tool}    "
                       f"Search method: {search_method}    "
                       f"Trimmer: {trimmer}"))
lines.append(("body", f"BLAST e-value cutoff: {evalue}"))
lines.append(("body", f"dS saturation threshold: {ds_sat_threshold:g}"))
lines.append(("gap", ""))

# ── Provenance ──────────────────────────────────────────────────────────────
lines.append(("header", "Provenance"))
lines.append(("body", f"Pipeline commit: {git_commit}"))
lines.append(("body", f"Python: {python_version}    Snakemake: {snakemake_version}"))
lines.append(("body", f"{search_method}: {search_tool_version}"))
lines.append(("body", f"{trimmer}: {trimmer_version}"))
lines.append(("body", f"mafft: {mafft_version}"))
lines.append(("body", f"codeml: {codeml_version}"))
lines.append(("gap", ""))

# ── Gene pair counts ────────────────────────────────────────────────────────
lines.append(("header", "Gene Pair Counts"))
if n_rbh is not None:
    lines.append(("body", f"RBH gene pairs (pre-alignment):  {n_rbh}"))
    if n_skipped_alignment:
        lines.append(("body", f"  Skipped at alignment/trimming stage: {n_skipped_alignment}"))
lines.append(("body", f"Gene pairs with dN/dS results: {total_pairs}"))
lines.append(("body", f"Valid pairs (dS > 0):          {n_valid}   ({pct(n_valid, total_pairs)})"))
lines.append(("body", f"Undefined (dS = 0):            {undefined}   ({pct(undefined, total_pairs)})"))
lines.append(("gap", ""))

# ── dN/dS summary statistics ────────────────────────────────────────────────
lines.append(("header", "dN/dS Summary Statistics"))
lines.append(("sub", f"{'':<8}{'mean':>10}{'median':>10}{'min':>10}{'max':>10}"))
for col in ("dN", "dS", "dNdS"):
    s = stats[col]
    lines.append(("body",
        f"{col:<8}{s['mean']:>10.4f}{s['median']:>10.4f}{s['min']:>10.4f}{s['max']:>10.4f}"))
lines.append(("gap", ""))

# ── Functional classification ───────────────────────────────────────────────
lines.append(("header", "Functional Classification"))
lines.append(("body", f"Conserved  (ω < 0.5):        {n_conserved:>6}   ({pct(n_conserved, n_valid)})"))
lines.append(("body", f"Relaxed    (0.5 ≤ ω < 1):    {n_relaxed:>6}   ({pct(n_relaxed, n_valid)})"))
lines.append(("body", f"Degenerate (ω ≥ 1):          {n_degenerate:>6}   ({pct(n_degenerate, n_valid)})"))
if n_degenerate_declared is not None and n_degenerate_declared != n_degenerate:
    lines.append(("body", f"  (genes_degenerate.tsv reports {n_degenerate_declared} rows)"))
lines.append(("gap", ""))

# ── Pseudogene evidence ──────────────────────────────────────────────────────
lines.append(("header", "Pseudogene Evidence (degenerate candidates, ω ≥ 1)"))
if pseudo_available:
    lines.append(("body", f"Strong evidence (stop codon or frameshift): {n_strong:>4}   ({pct(n_strong, n_pseudo_total)})"))
    lines.append(("body", f"Weak evidence (ω ≥ 1 only):                 {n_weak:>4}   ({pct(n_weak, n_pseudo_total)})"))
    if n_unknown:
        lines.append(("body", f"Unknown (sequence not found):               {n_unknown:>4}   ({pct(n_unknown, n_pseudo_total)})"))
else:
    lines.append(("body", "Not available for this run (genes_degenerate_annotated.tsv not found)."))
lines.append(("gap", ""))

# ── dS saturation ────────────────────────────────────────────────────────────
lines.append(("header", "dS Saturation"))
lines.append(("body", f"Gene pairs with dS > {ds_sat_threshold:g}: {n_saturated}   ({pct(n_saturated, n_valid)})"))
if n_saturated > 0:
    lines.append(("body", "Note: high dS values may reflect substitution saturation rather than"))
    lines.append(("body", "true divergence — interpret dN/dS for these pairs with caution."))
lines.append(("gap", ""))

# ── Output file checklist ────────────────────────────────────────────────────
lines.append(("header", f"Output File Checklist ({n_present}/{len(checklist_status)} present)"))

two_col = len(checklist_status) > MAX_SINGLE_COLUMN_FILES
if two_col:
    half = math.ceil(len(checklist_status) / 2)
    left_col, right_col = checklist_status[:half], checklist_status[half:]
    for i in range(half):
        rel_l, desc_l, ok_l = left_col[i]
        mark_l = "✓" if ok_l else "✗"
        left_txt = f"{mark_l} {rel_l} — {desc_l}"
        if i < len(right_col):
            rel_r, desc_r, ok_r = right_col[i]
            mark_r = "✓" if ok_r else "✗"
            right_txt = f"{mark_r} {rel_r} — {desc_r}"
        else:
            right_txt = None
        lines.append(("small", (left_txt, right_txt)))
else:
    for rel, desc, ok in checklist_status:
        mark = "✓" if ok else "✗"
        lines.append(("small", (f"{mark} {rel} — {desc}", None)))

# ── Render ────────────────────────────────────────────────────────────────────
log("\nRendering run_summary.pdf ...")

fig, ax = plt.subplots(figsize=(8.5, 11))
ax.set_axis_off()
fig.subplots_adjust(left=0.08, right=0.94, top=0.97, bottom=0.06)

TOP_Y, BOTTOM_Y = 0.985, 0.03
available = TOP_Y - BOTTOM_Y
total_weight = sum(KIND_WEIGHT[kind] for kind, _ in lines)
unit = available / total_weight

LEFT_X, RIGHT_X = 0.0, 0.5
cursor = TOP_Y

for kind, content in lines:
    h = KIND_WEIGHT[kind] * unit
    fontsize = KIND_FONTSIZE.get(kind)

    if kind == "title":
        ax.text(0.5, cursor, content, transform=ax.transAxes, ha="center", va="top",
                fontsize=fontsize, fontweight="bold", color="#212121")
    elif kind == "meta":
        ax.text(0.5, cursor, content, transform=ax.transAxes, ha="center", va="top",
                fontsize=fontsize, color="#666666", style="italic")
    elif kind == "header":
        ax.text(LEFT_X, cursor, content, transform=ax.transAxes, ha="left", va="top",
                fontsize=fontsize, fontweight="bold", color="#1B5E20")
        ax.plot([LEFT_X, 1.0], [cursor - h * 0.92, cursor - h * 0.92],
                transform=ax.transAxes, color="#C8E6C9", linewidth=1.0, zorder=1)
    elif kind == "sub":
        ax.text(LEFT_X, cursor, content, transform=ax.transAxes, ha="left", va="top",
                fontsize=fontsize, fontweight="bold", family="monospace", color="#424242")
    elif kind == "body":
        ax.text(LEFT_X, cursor, content, transform=ax.transAxes, ha="left", va="top",
                fontsize=fontsize, family="monospace", color="#424242")
    elif kind == "small":
        left_txt, right_txt = content
        ax.text(LEFT_X, cursor, left_txt, transform=ax.transAxes, ha="left", va="top",
                fontsize=fontsize, family="monospace", color="#555555")
        if right_txt is not None:
            ax.text(RIGHT_X, cursor, right_txt, transform=ax.transAxes, ha="left", va="top",
                    fontsize=fontsize, family="monospace", color="#555555")
    # "gap" — blank spacer, nothing to draw

    cursor -= h

# ── Footer (fixed position, drawn after the flowed content) ─────────────────
footer_text = (
    "dN/dS Pipeline · Suvratha Jayaprasad · "
    "https://github.com/suvrathaprasad/dN-dS-Snakemake-pipeline · "
    "© 2025 Suvratha Jayaprasad · CC BY-NC-ND 4.0"
)
fig.text(0.5, 0.015, footer_text, ha="center", va="bottom",
          fontsize=7, color="#888888", style="italic")

fig.savefig(out_pdf, format="pdf", bbox_inches=None, dpi=300)
plt.close(fig)

log(f"Written: {out_pdf}")
log("\nRun summary complete.")
_log_fh.close()
