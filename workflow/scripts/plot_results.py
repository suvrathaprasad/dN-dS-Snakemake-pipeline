#!/usr/bin/env python3
"""
plot_results.py — Visualisation and functional classification of dN/dS results.

Produces four publication-quality PDF plots and ten accompanying TSV tables.

Outlier handling:
  All data is kept in tables. For display purposes, axes are clipped at the
  99th percentile (minimum cap: ω=2 for dN/dS, adaptive for dN and dS).
  Points above the display limit are shown as upward triangles (▲) at the
  top of the axis with a note on the plot. This keeps the bulk of the
  distribution visible without removing any data from the analysis.

Plots:
  1. dnds_boxplot.pdf       — Boxplot of dN, dS, and dN/dS (independent axes)
  2. dnds_violin.pdf        — Violin plot of dN, dS, and dN/dS (independent axes)
  3. dnds_scatter.pdf       — dN vs dS scatter, colour-coded by diagonal position
  4. functional_summary.pdf — Bar chart of gene counts per functional category

Tables (output/results/tables/):
  genes_conserved.tsv        — ω < 0.5
  genes_relaxed.tsv          — 0.5 ≤ ω < 1
  genes_degenerate.tsv       — ω ≥ 1
  genes_undefined_ds.tsv     — dS = 0 (excluded from plots)
  genes_high_dn.tsv          — dN above median
  genes_high_ds.tsv          — dS above median
  genes_high_dnds.tsv        — dN/dS above median
  genes_above_diagonal.tsv   — dN > dS
  genes_on_diagonal.tsv      — dN ≈ dS (within ±10%)
  genes_below_diagonal.tsv   — dN < dS
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
matplotlib.use("Agg")   # non-interactive backend for HPC
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages   # used for placeholder plots
import numpy as np
import pandas as pd

# ── Snakemake-injected objects ────────────────────────────────────────────────
input_tsv  = snakemake.input.tsv
plots_dir  = Path(snakemake.params.plots_dir)
tables_dir = Path(snakemake.params.tables_dir)
diag_tol   = float(snakemake.params.get("diagonal_tolerance", 0.10))
log_path   = Path(snakemake.log[0])

out_box        = snakemake.output.boxplot
out_violin     = snakemake.output.violin
out_scatter    = snakemake.output.scatter
out_summary    = snakemake.output.summary
out_degenerate = snakemake.output.degenerate

plots_dir.mkdir(parents=True, exist_ok=True)
tables_dir.mkdir(parents=True, exist_ok=True)
log_path.parent.mkdir(parents=True, exist_ok=True)
Path(out_degenerate).parent.mkdir(parents=True, exist_ok=True)

# ── Colour palette ────────────────────────────────────────────────────────────
# Note: COL["conserved"] and COL["dn"] share the same blue intentionally —
# conserved genes (low dN/dS) map naturally to the dN colour.
COL = {
    "dn":         "#2196F3",   # blue   — dN metric and conserved genes
    "ds":         "#4CAF50",   # green  — dS metric and neutral diagonal
    "dnds":       "#7B1FA2",   # purple — dN/dS metric
    "conserved":  "#2196F3",   # blue   — ω < 0.5
    "relaxed":    "#FF9800",   # amber  — 0.5 ≤ ω < 1
    "degenerate": "#F44336",   # red    — ω ≥ 1
    "undefined":  "#9E9E9E",   # grey   — dS = 0
    "neutral":    "#4CAF50",   # green  — on diagonal
}

# ── Logging — single open handle per run ─────────────────────────────────────
_log_fh = open(log_path, "w")

def log(msg: str) -> None:
    _log_fh.write(msg + "\n")
    _log_fh.flush()
    print(msg)

# ── Load and validate ─────────────────────────────────────────────────────────
log(f"Loading: {input_tsv}")
df = pd.read_csv(input_tsv, sep="\t")
log(f"Total gene pairs loaded: {len(df)}")

required = {"Gene_query", "Gene_target", "dN", "dS", "dNdS"}
missing  = required - set(df.columns)
if missing:
    _log_fh.close()
    sys.exit(f"Missing columns in {input_tsv}: {missing}")

for col in ["dN", "dS", "dNdS"]:
    df[col] = pd.to_numeric(df[col], errors="coerce")

# ── Separate undefined dS ─────────────────────────────────────────────────────
undef = df[df["dS"] == 0].copy()
valid = df[df["dS"] >  0].copy()

log(f"  Valid (dS > 0):   {len(valid)}")
log(f"  Undefined (dS=0): {len(undef)}")

# ── Placeholder PDFs when no valid data (e.g. back-translated test data) ──────
if len(valid) == 0:
    log("WARNING: No valid gene pairs with dS > 0. Generating placeholder plots.")
    log("This is expected with back-translated test data.")
    for out in [out_box, out_violin, out_scatter, out_summary]:
        with PdfPages(out) as pdf:
            fig, ax = plt.subplots(figsize=(6, 4))
            ax.text(0.5, 0.5,
                    "No data with dS > 0\n(expected with test data)",
                    ha="center", va="center", fontsize=12, color="#888888",
                    transform=ax.transAxes)
            ax.set_axis_off()
            pdf.savefig(fig, bbox_inches="tight")
            plt.close()
    for fname in ["genes_conserved.tsv", "genes_relaxed.tsv", "genes_degenerate.tsv",
                  "genes_high_dn.tsv", "genes_high_ds.tsv", "genes_high_dnds.tsv",
                  "genes_above_diagonal.tsv", "genes_on_diagonal.tsv",
                  "genes_below_diagonal.tsv"]:
        (tables_dir / fname).write_text("\t".join(df.columns) + "\n")
    undef.to_csv(tables_dir / "genes_undefined_ds.tsv", sep="\t", index=False)
    _log_fh.close()
    sys.exit(0)

# ── Functional classification ─────────────────────────────────────────────────
conserved  = valid[valid["dNdS"] <  0.5].copy()
relaxed    = valid[(valid["dNdS"] >= 0.5) & (valid["dNdS"] < 1.0)].copy()
degenerate = valid[valid["dNdS"] >= 1.0].copy()

log(f"\nFunctional classification:")
log(f"  Conserved  (ω < 0.5):       {len(conserved)}")
log(f"  Relaxed    (0.5 ≤ ω < 1.0): {len(relaxed)}")
log(f"  Degenerate (ω ≥ 1.0):       {len(degenerate)}")
log(f"  Undefined  (dS = 0):        {len(undef)}")

# ── Distribution tables ───────────────────────────────────────────────────────
high_dn   = valid[valid["dN"]   > valid["dN"].median()].copy()
high_ds   = valid[valid["dS"]   > valid["dS"].median()].copy()
high_dnds = valid[valid["dNdS"] > valid["dNdS"].median()].copy()

# ── Diagonal classification ───────────────────────────────────────────────────
ratio      = valid["dN"] / valid["dS"]
above_diag = valid[ratio >  (1 + diag_tol)].copy()
on_diag    = valid[(ratio >= (1 - diag_tol)) & (ratio <= (1 + diag_tol))].copy()
below_diag = valid[ratio <  (1 - diag_tol)].copy()

log(f"\nDiagonal classification (tolerance ±{int(diag_tol*100)}%):")
log(f"  Above diagonal (dN > dS): {len(above_diag)}")
log(f"  On diagonal   (dN ≈ dS): {len(on_diag)}")
log(f"  Below diagonal (dN < dS): {len(below_diag)}")

# ── Write all tables ──────────────────────────────────────────────────────────
table_map = {
    "genes_conserved.tsv":      conserved,
    "genes_relaxed.tsv":        relaxed,
    "genes_degenerate.tsv":     degenerate,
    "genes_undefined_ds.tsv":   undef,
    "genes_high_dn.tsv":        high_dn,
    "genes_high_ds.tsv":        high_ds,
    "genes_high_dnds.tsv":      high_dnds,
    "genes_above_diagonal.tsv": above_diag,
    "genes_on_diagonal.tsv":    on_diag,
    "genes_below_diagonal.tsv": below_diag,
}
for fname, subset in table_map.items():
    path = tables_dir / fname
    subset.to_csv(path, sep="\t", index=False)
    log(f"Written: {path}  ({len(subset)} genes)")


# ─────────────────────────────────────────────────────────────────────────────
# Shared helpers
# ─────────────────────────────────────────────────────────────────────────────

def clip_limit(series, col_name, min_cap=2.0):
    """
    Display upper limit: 99th percentile × 1.1.
    dN/dS has a hard minimum cap of min_cap (default 2.0) so the
    neutral line at ω=1 always has breathing room above it.
    dN and dS use a floor of 3× median so the axis is never too tight.
    """
    p99 = series.quantile(0.99)
    if col_name == "dNdS":
        return max(p99 * 1.1, min_cap)
    return max(p99 * 1.1, series.median() * 3, 0.01)


def jitter(data, pos=1, width=0.06, seed=42):
    rng = np.random.default_rng(seed)
    return pos + rng.uniform(-width, width, size=len(data))


def add_datapoints(ax, data, pos=1, alpha=0.35, size=18):
    """Jittered black strip overlay for individual data points."""
    ax.scatter(jitter(data, pos=pos), data,
               color="black", alpha=alpha, s=size, zorder=5, linewidths=0)


def draw_outlier_markers(ax, data, ymax, pos=1):
    """Pin out-of-range points as upward black triangles at top of axis."""
    outside = data[data > ymax]
    if len(outside) > 0:
        ax.scatter(jitter(outside.values, pos=pos),
                   np.full(len(outside), ymax * 0.975),
                   marker="^", color="black", alpha=0.55, s=24, zorder=6)
    return len(outside)


def outlier_note(ax, n_clipped, ymax, metric_label):
    """Add italic note below axis about clipped points."""
    if n_clipped == 0:
        return
    ax.text(0.5, -0.09,
            f"{n_clipped} gene(s) above display limit "
            f"({metric_label} > {ymax:.3g}) shown as ▲ — see tables for full data",
            transform=ax.transAxes, fontsize=7.5, ha="center",
            color="#777777", style="italic", va="top")


def add_stats_box(ax, values):
    """Median, mean, n summary in top-right corner."""
    ax.text(0.97, 0.97,
            f"median = {values.median():.4f}\nmean = {values.mean():.4f}\nn = {len(values)}",
            transform=ax.transAxes, fontsize=8, va="top", ha="right",
            color="#444444",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                      edgecolor="#cccccc", alpha=0.8))


def neutral_line(ax, ymax):
    """Dashed red line at ω=1 with shaded regions."""
    ax.axhline(y=1.0, color="#F44336", linewidth=1.3,
               linestyle="--", zorder=4)
    ax.axhspan(0,   1.0,  alpha=0.03, color=COL["conserved"],  zorder=0)
    ax.axhspan(1.0, ymax, alpha=0.03, color=COL["degenerate"], zorder=0)


def make_boxplot_panel(ax, col, vals, colour, label):
    """Draw one boxplot panel with data points and outlier handling."""
    data = vals.dropna()
    ymax = clip_limits[col]
    bp = ax.boxplot(
        data, patch_artist=True, widths=0.45,
        medianprops=dict(color="white", linewidth=2.0),
        flierprops=dict(marker="", markersize=0),  # hide default fliers
    )
    bp["boxes"][0].set_facecolor(colour)
    bp["boxes"][0].set_alpha(0.75)
    add_datapoints(ax, data[data <= ymax].values)
    n_clip = draw_outlier_markers(ax, data, ymax)
    if col == "dNdS":
        neutral_line(ax, ymax)
    ax.set_ylim(bottom=0, top=ymax)
    ax.set_xticks([])
    ax.set_xlabel(label, fontsize=11, color=colour, fontweight="bold")
    ax.set_ylabel(label, fontsize=11)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_facecolor("#fafafa")
    add_stats_box(ax, data)
    outlier_note(ax, n_clip, ymax, label)


def make_violin_panel(ax, col, vals, colour, label, pos=1):
    """Draw one violin panel with data points and outlier handling."""
    data = vals.dropna()
    ymax = clip_limits[col]
    # Violin fitted to all data — shape reflects true distribution including
    # the outlier tail, but y-axis is clipped for readability
    parts = ax.violinplot(data, positions=[pos],
                          showmedians=True, showextrema=True)
    for pc in parts["bodies"]:
        pc.set_facecolor(colour)
        pc.set_alpha(0.72)
    parts["cmedians"].set_color("white")
    parts["cmedians"].set_linewidth(2.0)
    for part in ["cmins", "cmaxes", "cbars"]:
        parts[part].set_color(colour)
        parts[part].set_alpha(0.6)
    add_datapoints(ax, data[data <= ymax].values, pos=pos)
    n_clip = draw_outlier_markers(ax, data, ymax, pos=pos)
    return n_clip


# Pre-compute per-metric display clip limits
metrics = [
    ("dN",   valid["dN"],   COL["dn"],   "dN"),
    ("dS",   valid["dS"],   COL["ds"],   "dS"),
    ("dNdS", valid["dNdS"], COL["dnds"], "dN/dS (ω)"),
]
clip_limits = {col: clip_limit(vals.dropna(), col)
               for col, vals, _, _ in metrics}

log(f"\nDisplay clip limits:")
for col, lim in clip_limits.items():
    n_above = int((valid[col].dropna() > lim).sum())
    log(f"  {col}: display ≤ {lim:.4g}  ({n_above} gene(s) above limit shown as ▲)")

# Shared legend handles for distribution plots
legend_handles = [
    mpatches.Patch(color=COL["conserved"],  alpha=0.4,
                   label="Purifying selection (ω < 1)"),
    mpatches.Patch(color=COL["degenerate"], alpha=0.4,
                   label="Relaxed / positive selection (ω ≥ 1)"),
    plt.Line2D([0], [0], color="#F44336", linewidth=1.3,
               linestyle="--", label="Neutral expectation (ω = 1)"),
]

# ─────────────────────────────────────────────────────────────────────────────
# PLOT 1 — Boxplot: dN, dS, dN/dS
# ─────────────────────────────────────────────────────────────────────────────
log("\nGenerating Plot 1: Boxplot")

fig, axes = plt.subplots(1, 3, figsize=(13, 6),
                         gridspec_kw={"wspace": 0.35})
fig.suptitle("Distribution of dN, dS and dN/dS",
             fontsize=14, fontweight="bold", y=1.02)
for ax, (col, vals, colour, label) in zip(axes, metrics):
    make_boxplot_panel(ax, col, vals, colour, label)
fig.legend(handles=legend_handles, loc="lower center", ncol=3,
           fontsize=9, bbox_to_anchor=(0.5, -0.06))
plt.tight_layout()
fig.savefig(out_box, format="pdf", bbox_inches="tight", dpi=300)
plt.close()
log(f"Written: {out_box}")

# ─────────────────────────────────────────────────────────────────────────────
# PLOT 2 — Violin plot: dN, dS, dN/dS
# ─────────────────────────────────────────────────────────────────────────────
log("\nGenerating Plot 2: Violin plot")

fig, axes = plt.subplots(1, 3, figsize=(13, 6),
                         gridspec_kw={"wspace": 0.35})
fig.suptitle("Distribution of dN, dS and dN/dS",
             fontsize=14, fontweight="bold", y=1.02)
for ax, (col, vals, colour, label) in zip(axes, metrics):
    n_clip = make_violin_panel(ax, col, vals, colour, label)
    ymax   = clip_limits[col]
    if col == "dNdS":
        neutral_line(ax, ymax)
    ax.set_ylim(bottom=0, top=ymax)
    ax.set_xticks([])
    ax.set_xlabel(label, fontsize=11, color=colour, fontweight="bold")
    ax.set_ylabel(label, fontsize=11)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_facecolor("#fafafa")
    add_stats_box(ax, vals.dropna())
    outlier_note(ax, n_clip, ymax, label)
fig.legend(handles=legend_handles, loc="lower center", ncol=3,
           fontsize=9, bbox_to_anchor=(0.5, -0.06))
plt.tight_layout()
fig.savefig(out_violin, format="pdf", bbox_inches="tight", dpi=300)
plt.close()
log(f"Written: {out_violin}")

# ─────────────────────────────────────────────────────────────────────────────
# PLOT 3 — dN vs dS scatter
#
# Both axes clipped at min(mean + 2SD, 5×median), floor of 3×median.
# Points outside the display window are pinned as same-colour triangles
# on the nearest axis edge.
# ─────────────────────────────────────────────────────────────────────────────
log("\nGenerating Plot 3: dN vs dS scatter")

def scatter_clip(series):
    """min(mean + 2SD, 5×median) with a floor of 3×median."""
    mean, sd, med = series.mean(), series.std(), series.median()
    return max(min(mean + 2 * sd, 5 * med), med * 3, 0.01)

ds_clip = scatter_clip(valid["dS"])
dn_clip = scatter_clip(valid["dN"])

n_clipped_scatter = int(((valid["dS"] > ds_clip) | (valid["dN"] > dn_clip)).sum())
log(f"  Scatter clip: dS ≤ {ds_clip:.4g}, dN ≤ {dn_clip:.4g}  "
    f"({n_clipped_scatter} gene(s) outside display window)")

def make_scatter_subsets(df_cat):
    """Split a diagonal category into in-window and out-of-window subsets."""
    in_win = (df_cat["dS"] <= ds_clip) & (df_cat["dN"] <= dn_clip)
    return df_cat[in_win], df_cat[~in_win]

fig, ax = plt.subplots(figsize=(9, 8))

scatter_cats = [
    (below_diag, COL["conserved"],
     f"Below diagonal — dN < dS (purifying, n={len(below_diag)})"),
    (on_diag,    COL["neutral"],
     f"On diagonal — dN ≈ dS (±{int(diag_tol*100)}%, neutral, n={len(on_diag)})"),
    (above_diag, COL["degenerate"],
     f"Above diagonal — dN > dS (degenerate / positive, n={len(above_diag)})"),
]

for subset, colour, label in scatter_cats:
    in_win, out_win = make_scatter_subsets(subset)
    if len(in_win) > 0:
        ax.scatter(in_win["dS"], in_win["dN"],
                   c=colour, alpha=0.75, s=45,
                   edgecolors="white", linewidths=0.3,
                   label=label, zorder=3)
    if len(out_win) > 0:
        ax.scatter(np.clip(out_win["dS"].values, 0, ds_clip * 0.975),
                   np.clip(out_win["dN"].values, 0, dn_clip * 0.975),
                   marker="^", c=colour, alpha=0.55, s=35,
                   edgecolors="white", linewidths=0.3, zorder=4)

ax.plot([0, ds_clip], [0, ds_clip], color="#555555",
        linewidth=1.2, linestyle="--", label="dN = dS (neutral)", zorder=2)
ax.fill_between(
    [0, ds_clip],
    [0, ds_clip * (1 - diag_tol)],
    [0, min(ds_clip * (1 + diag_tol), dn_clip)],
    alpha=0.07, color=COL["neutral"], zorder=1,
    label=f"±{int(diag_tol*100)}% tolerance band",
)
ax.set_xlim(left=0, right=ds_clip)
ax.set_ylim(bottom=0, top=dn_clip)
ax.set_xlabel("dS (synonymous substitution rate)", fontsize=12)
ax.set_ylabel("dN (non-synonymous substitution rate)", fontsize=12)
ax.set_title("dN vs dS — Scatter Plot", fontsize=14, fontweight="bold")
ax.legend(fontsize=9, loc="upper left", framealpha=0.9)
ax.spines[["top", "right"]].set_visible(False)

notes = []
if n_clipped_scatter > 0:
    notes.append(f"{n_clipped_scatter} gene(s) outside display window "
                 f"(dS>{ds_clip:.3g} or dN>{dn_clip:.3g}) shown as ▲")
if len(undef) > 0:
    notes.append(f"{len(undef)} gene(s) with dS=0 excluded (see genes_undefined_ds.tsv)")
if notes:
    ax.text(0.98, 0.02, "\n".join(notes),
            transform=ax.transAxes, fontsize=7.5, color="#777777",
            ha="right", va="bottom", style="italic")

plt.tight_layout()
fig.savefig(out_scatter, format="pdf", bbox_inches="tight", dpi=300)
plt.close()
log(f"Written: {out_scatter}")

# ─────────────────────────────────────────────────────────────────────────────
# PLOT 4 — Functional classification bar chart
# ─────────────────────────────────────────────────────────────────────────────
log("\nGenerating Plot 4: Functional classification summary")

categories = ["Conserved\n(ω < 0.5)",
              "Relaxed\n(0.5 ≤ ω < 1)",
              "Degenerate\n(ω ≥ 1)",
              "Undefined\n(dS = 0)"]
counts  = [len(conserved), len(relaxed), len(degenerate), len(undef)]
colours = [COL["conserved"], COL["relaxed"], COL["degenerate"], COL["undefined"]]
total   = sum(counts)

fig, ax = plt.subplots(figsize=(9, 6))
bars = ax.bar(categories, counts, color=colours,
              edgecolor="white", linewidth=0.8, width=0.55, zorder=3)
for bar, count in zip(bars, counts):
    pct = count / total * 100 if total > 0 else 0
    ax.text(bar.get_x() + bar.get_width() / 2,
            bar.get_height() + max(counts) * 0.01,
            f"{count}\n({pct:.1f}%)",
            ha="center", va="bottom", fontsize=10, fontweight="bold")
ax.set_ylabel("Number of gene pairs", fontsize=12)
ax.set_title("Functional Classification of Gene Pairs",
             fontsize=14, fontweight="bold")
ax.set_ylim(top=max(counts) * 1.18)
ax.spines[["top", "right"]].set_visible(False)
ax.yaxis.grid(True, linestyle="--", alpha=0.4, zorder=0)
ax.set_axisbelow(True)
ax.text(0.98, 0.97, f"Total gene pairs: {total}",
        transform=ax.transAxes, fontsize=10,
        ha="right", va="top", color="#555555")

plt.tight_layout()
fig.savefig(out_summary, format="pdf", bbox_inches="tight", dpi=300)
plt.close()
log(f"Written: {out_summary}")

log("\nAll plots and tables written successfully.")
_log_fh.close()
