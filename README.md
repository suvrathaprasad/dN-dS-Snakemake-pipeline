# dN/dS Pipeline

A reproducible Snakemake pipeline for estimating **pairwise synonymous (dS) and non-synonymous (dN) substitution rates** between two species, starting from genome assemblies, CDS nucleotide sequences, or pre-translated protein files.

**Author:** Suvratha Jayaprasad  
**Contact:** [suvrathaprasad.github.io](https://suvrathaprasad.github.io/index.html) — for questions, bug reports, or collaboration enquiries, please use the contact form on the website.  
**License:** [CC BY-NC-ND 4.0](https://creativecommons.org/licenses/by-nc-nd/4.0/) — you may use and share this pipeline with attribution, but you may not modify it or use it for commercial purposes. See `LICENSE` for full terms.  
**Citation:** If you use this pipeline, please cite this repository and the tools it wraps (see [Dependencies](#dependencies)):
> Jayaprasad, S. (2025). *dN/dS Pipeline*. GitHub. https://github.com/suvrathaprasad/dnds_pipeline

---

## Overview

```
 ┌─ Mode A: Genome + GFF ──────┐   ┌─ Mode B: FAA + FNA ──────┐   ┌─ Mode C: FNA only ───────┐
 │                             │   │                          │   │                          │
 │  Anchorwave or gffread      │   │  (skip extraction and    │   │  EMBOSS transeq          │
 │    Extract CDS per species  │   │     translation)         │   │    Translate → protein   │
 │      │                      │   │                          │   │                          │
 │      ▼                      │   │                          │   │                          │
 │  EMBOSS transeq → Protein   │   │                          │   │                          │
 └──────────────┬──────────────┘   └────────────┬─────────────┘   └────────────┬─────────────┘
                │                               │                              │
                └──────────────────────────┬────┴──────────────────────────────┘
                                           │
                                           ▼
          BLAST+ or DIAMOND   Bidirectional search + Reciprocal Best Hits
                                           │
                                           ▼
          extract_gene_pair.py   Per-gene-pair FAA + FNA extraction
                                           │
                                           ▼
          MAFFT                  Protein alignment per gene pair
                                           │
                                           ▼
          pal2nal                Back-translate → codon alignment
                                           │
                                           ▼
          Gblocks or trimAl      Trim poorly aligned regions
                                           │
                                           ▼
          PAML codeml            Pairwise dN/dS (runmode = -2)
                                           │
                                           ▼
          collate_results.py     Parse codeml output + dS saturation check
                                           │
                                           ▼
          output/results/dnds_output.tsv + plots/ + tables/
                                           │
                                           ▼
          check_pseudogenes.py   Premature stop & frameshift screen
                                           │
                                           ▼
          genes_degenerate_annotated.tsv + pseudogene_evidence.pdf
                                           │
                                           ▼
          write_summary.py       Collate everything into one report
                                           │
                                           ▼
                             output/results/run_summary.pdf
```

**Why RBH and not OrthoFinder?** For strict pairwise dN/dS between two species, RBH identifies one-to-one orthologs more directly and with less overhead. OrthoFinder is designed for multi-species ortholog inference and produces many-to-many groups that require additional filtering for pairwise analyses.

---

## Output structure

```
output/
├── intermediate/
│   ├── cds/           CDS and protein sequences (Modes A and C)
│   ├── blast/         Search databases, forward and reverse results
│   ├── rbh_pairs.tsv  RBH gene pair list (query TAB target + stats)
│   ├── input/         Per-gene-pair FAA and FNA files, plus cached
│   │                  `*.offset_idx.json` index files (auto-generated
│   │                  next to each source FASTA to speed up per-gene
│   │                  extraction; safe to delete, rebuilt automatically)
│   ├── aligns/        MAFFT, pal2nal, and trimmed alignments, plus a
│   │                  per-gene subfolder (e.g. `gene1/`) holding any
│   │                  Gblocks `.ps`/`.html` sidecar files for that gene
│   └── codeml/        Per-gene codeml output directories
├── results/
│   ├── dnds_output.tsv              ← FINAL dN/dS TABLE  (all gene pairs)
│   ├── plots/
│   │   ├── dnds_boxplot.pdf              Boxplot of dN, dS, and dN/dS
│   │   ├── dnds_violin.pdf               Violin plot of dN, dS, and dN/dS
│   │   ├── dnds_scatter.pdf              dN vs dS scatter
│   │   ├── functional_summary.pdf        Gene category bar chart
│   │   └── pseudogene_evidence.pdf       Pseudogene evidence chart
│   ├── tables/
│   │   ├── genes_conserved.tsv               ω < 0.5
│   │   ├── genes_relaxed.tsv                 0.5 ≤ ω < 1
│   │   ├── genes_degenerate.tsv              ω ≥ 1
│   │   ├── genes_degenerate_annotated.tsv    ω ≥ 1 + sequence evidence
│   │   ├── genes_ds_saturated.tsv            dS > saturation threshold
│   │   ├── genes_undefined_ds.tsv            dS = 0 (excluded from plots)
│   │   ├── genes_high_dn.tsv                 dN above median
│   │   ├── genes_high_ds.tsv                 dS above median
│   │   ├── genes_high_dnds.tsv               dN/dS above median
│   │   ├── genes_above_diagonal.tsv          dN > dS
│   │   ├── genes_on_diagonal.tsv             dN ≈ dS (±10%)
│   │   └── genes_below_diagonal.tsv          dN < dS
│   └── run_summary.pdf              ← FULL RUN SUMMARY REPORT
└── logs/              Log files for every step
```

### `dnds_output.tsv` format

Tab-separated, one row per orthologous gene pair:

| Column      | Description                          |
|-------------|--------------------------------------|
| Gene_query  | Gene ID from the query species       |
| Gene_target | Gene ID from the target species      |
| t           | Divergence (branch length)           |
| dN          | Non-synonymous substitution rate     |
| dS          | Synonymous substitution rate         |
| dNdS        | ω = dN/dS                            |

---

## Requirements

- [Snakemake](https://snakemake.readthedocs.io) ≥ 7.0
- [Conda](https://docs.conda.io) or [Mamba](https://mamba.readthedocs.io) (recommended)

All other dependencies are installed automatically per-rule via the `envs/` conda environments.

---

## Installation

```bash
# 1. Clone the repository
git clone https://github.com/suvrathaprasad/dnds_pipeline.git
cd dnds_pipeline

# 2. Install Snakemake (if not already installed)
conda create -n snakemake -c conda-forge -c bioconda snakemake mamba
conda activate snakemake
```

---

## Quick Start

### 1. Configure inputs

Edit `config/config.yaml`. Three input modes are supported and can be mixed between species:

**Mode A — genome + GFF** (CDS extracted automatically, EMBOSS translates):
```yaml
query:
  fasta: "path/to/species1.fasta"
  gff:   "path/to/species1.gff"
  prefix: "species1"
```

**Mode B — pre-made FAA + FNA** (skip extraction and translation entirely):
```yaml
query:
  faa: "path/to/species1.faa"   # protein sequences
  fna: "path/to/species1.fna"   # CDS nucleotide sequences
  prefix: "species1"
```

**Mode C — CDS nucleotide FNA only** (skip extraction, transeq runs automatically):
```yaml
query:
  fna: "path/to/species1.fna"   # CDS nucleotide sequences
  prefix: "species1"
```

The pipeline detects which mode to use automatically. Priority when multiple keys are present: Mode B > Mode C > Mode A.

By default all output is written under `output/` (paths throughout this
README assume that default). This is configurable via `outdir` in
`config/config.yaml`:

```yaml
outdir: "output"   # change to write intermediate/ and results/ elsewhere
```

### 2. Configure tool switches (optional)

Three tool choices are available in `config/config.yaml`. All default to the original tools used in the published analysis:

```yaml
tools:
  cds_extraction: "anchorwave"  # anchorwave (default) | gffread
  search_method:  "blast"       # blast (default)      | diamond
  trimmer:        "gblocks"     # gblocks (default)    | trimal
```

| Switch | Default | Alternative | When to switch |
|--------|---------|-------------|----------------|
| `cds_extraction` | Anchorwave | gffread | Faster CDS extraction when whole-genome alignment is not needed |
| `search_method` | BLAST+ | DIAMOND | Large proteomes (tens of thousands of proteins) where speed matters |
| `trimmer` | Gblocks | trimAl | Less aggressive trimming; preferable for divergent sequences |

### 3. Run

```bash
# Dry run first (recommended)
snakemake --snakefile workflow/Snakefile \
          --use-conda \
          --cores 16 \
          -n

# Full run
snakemake --snakefile workflow/Snakefile \
          --use-conda \
          --cores 16
```

### 4. HPC / SLURM

```bash
snakemake --snakefile workflow/Snakefile \
          --use-conda \
          --executor slurm \
          --jobs 100 \
          --default-resources slurm_partition=standard mem_mb=8000 runtime=120
```

---

## Test dataset

A minimal test dataset (20 gene pairs, Mode B) is provided in `test/`.

```bash
snakemake --snakefile workflow/Snakefile \
          --configfile test/config_test.yaml \
          --use-conda \
          --cores 4
```

Expected output: `test/output/results/dnds_output.tsv`

> **Note:** The test nucleotide sequences are back-translations of the protein sequences, so dS will be near zero. The test is for toolchain verification, not biological interpretation.

---

## Visualisation and functional classification

Four PDF plots and ten TSV gene-list tables are generated automatically
after `collate_results` and written to `output/results/plots/` and
`output/results/tables/`.

### Plot 1 — Boxplot (`dnds_boxplot.pdf`)
Side-by-side boxplots of dN, dS, and dN/dS with independent y-axes. The
dN/dS panel includes a dashed line at ω = 1 and background shading
separating purifying selection (ω < 1) from relaxed or positive selection
(ω ≥ 1). Individual data points are overlaid as a jittered strip.
Outliers above the display limit are shown as triangles (▲) with a note.

### Plot 2 — Violin plot (`dnds_violin.pdf`)
The same three metrics as Plot 1, shown as violin plots. The violin shape
is fitted to all data including outliers, but the y-axis is clipped for
readability. Saved as a separate PDF so boxplot and violin can be used
independently in publications.

### Plot 3 — dN vs dS scatter (`dnds_scatter.pdf`)
Classical molecular evolution scatter plot with genes colour-coded by
their position relative to the neutral diagonal (dN = dS):

| Colour | Category | Interpretation |
|--------|----------|----------------|
| Blue   | Below diagonal (dN < dS) | Under purifying selection |
| Green  | On diagonal (dN ≈ dS, ±10%) | Evolving neutrally |
| Red    | Above diagonal (dN > dS) | Degenerate / positive selection candidates |

The dS saturation threshold (see [below](#ds-saturation-warning)) is also
drawn on this plot: a vertical dotted orange line at `dS = dS_saturation_threshold`
if that value falls within the displayed axis range, or a text annotation
noting the threshold if it falls outside the range shown.

### Plot 4 — Functional summary (`functional_summary.pdf`)
Bar chart showing gene counts per functional category with percentage
labels. Categories follow standard dN/dS interpretation thresholds:

| Category | ω range | Interpretation |
|----------|---------|----------------|
| Conserved | ω < 0.5 | Strong purifying selection, gene likely functional |
| Relaxed | 0.5 ≤ ω < 1 | Weakened purifying selection |
| Degenerate | ω ≥ 1 | Relaxed or positive selection; pseudogene candidates |
| Undefined | dS = 0 | No synonymous divergence; excluded from plots |

> **Important note on the "Degenerate" category:** genes classified here
> have ω ≥ 1, which is consistent with relaxed purifying selection or
> pseudogenisation, but **this classification alone is not sufficient to
> confirm pseudogene status**. Elevated dN/dS can also result from
> positive selection, short alignments, or noisy estimates in genes with
> low dS. We strongly recommend treating `genes_degenerate.tsv` as a
> candidate list and following up with dedicated pseudogene detection tools
> (e.g. [PGAP](https://github.com/ncbi/pgap),
> [PseudoPipe](http://pseudogene.org/pseudopipe/), or manual inspection of
> the alignments) before drawing biological conclusions about gene
> functionality.

### TSV gene lists
Ten TSV files are written to `output/results/tables/`, each containing
the full gene pair records from `dnds_output.tsv` for that category:

| File | Contents |
|------|----------|
| `genes_conserved.tsv` | ω < 0.5 — strong purifying selection |
| `genes_relaxed.tsv` | 0.5 ≤ ω < 1 — weakened purifying selection |
| `genes_degenerate.tsv` | ω ≥ 1 — pseudogene candidates (see note above) |
| `genes_undefined_ds.tsv` | dS = 0 — excluded from plots |
| `genes_high_dn.tsv` | dN above median |
| `genes_high_ds.tsv` | dS above median |
| `genes_high_dnds.tsv` | dN/dS above median |
| `genes_above_diagonal.tsv` | dN > dS (scatter plot above diagonal) |
| `genes_on_diagonal.tsv` | dN ≈ dS within ±10% |
| `genes_below_diagonal.tsv` | dN < dS (scatter plot below diagonal) |

> `genes_ds_saturated.tsv` also lives in `output/results/tables/`, but is
> produced by `collate_results.py` rather than this step — see
> [dS saturation warning](#ds-saturation-warning) below.

---

## Pseudogene sequence evidence

For each gene pair in `genes_degenerate.tsv` (ω ≥ 1), the pipeline
automatically checks the CDS nucleotide sequences for two classic
hallmarks of pseudogenisation:

**Premature stop codons** — a stop codon (TAA, TAG, TGA) appearing
before the final codon position, indicating a disrupted reading frame.

**Frameshifts** — CDS length not divisible by 3, indicating an insertion
or deletion that shifts the reading frame.

Each gene is checked independently in both species. Results are written to
`output/results/tables/genes_degenerate_annotated.tsv`, which extends
`genes_degenerate.tsv` with five additional columns:

| Column | Values | Description |
|--------|--------|-------------|
| `premature_stop_query` | yes / no / unknown | Premature stop in query species CDS |
| `premature_stop_target` | yes / no / unknown | Premature stop in target species CDS |
| `frameshift_query` | yes / no / unknown | Frameshift in query species CDS |
| `frameshift_target` | yes / no / unknown | Frameshift in target species CDS |
| `pseudogene_evidence` | strong / weak / unknown | Strong = at least one hallmark detected; Weak = ω ≥ 1 only; Unknown = the CDS sequence for this gene ID couldn't be located under `intermediate/input/` |

A two-panel PDF (`pseudogene_evidence.pdf`) summarises the counts:
- Left panel: gene pairs by evidence level (strong / weak, plus an
  Unknown bar if any gene IDs couldn't be matched to a sequence)
- Right panel: breakdown by hallmark type (premature stop vs frameshift, per species)

> **Important:** this is a first-pass sequence screen, not a replacement
> for dedicated pseudogene detection tools. False positives can arise from
> sequencing errors, incomplete annotations, or genuine positive selection.
> Treat `strong` evidence genes as priority candidates for manual
> inspection or follow-up with tools such as
> [PGAP](https://github.com/ncbi/pgap) or
> [PseudoPipe](http://pseudogene.org/pseudopipe/).

**Why not a full pseudogene detection pipeline?**

Dedicated pseudogene tools such as PseudoPipe and PGAP operate on a
single genome assembly and GFF — they scan for disrupted open reading
frames, truncated proteins, and homology to known functional genes
across the whole genome. This pipeline is fundamentally different: it
is a *comparative* analysis working on *pairs* of orthologous genes
between two species, and by the time the pipeline reaches this step the
sequences have already been through Gblocks or trimAl, which can remove
or alter the very frameshifts and stop codons that pseudogene tools rely
on. A full pseudogene pipeline would also require Mode A input (genome +
GFF) for both species, immediately excluding users running Modes B or C.

The `check_pseudogenes.py` step is therefore intentionally scoped to
what is reliably detectable from the per-gene CDS sequences that the
pipeline already has: premature stop codons and frameshifts in the
pre-trimmed FNA files. This provides a rapid, dependency-free first
screen that works in all input modes and flags the candidates most worth
following up, without making claims the data cannot support.

---

## dS saturation warning

Synonymous sites approach mutational saturation when dS is high — multiple
substitutions occur at the same site, making dS an unreliable measure of
divergence time and inflating or deflating dN/dS estimates unpredictably.
The standard threshold in the field is dS > 2.0.

Gene pairs exceeding the threshold are written to
`output/results/tables/genes_ds_saturated.tsv` and flagged in the run
summary. They are **not** removed from `dnds_output.tsv` — the decision
of whether to exclude them is left to the user — but they are clearly
identified so downstream interpretation can account for them.

The threshold is configurable in `config/config.yaml`:

```yaml
dS_saturation_threshold: 2.0   # standard field value
```

---

## Run summary report

After all steps complete, the pipeline writes
`output/results/run_summary.pdf` — a single-page PDF summarising the
entire run. It includes:

- Run metadata: species names, resolved input mode and file path(s) for
  each species (the actual `fasta`/`gff`/`faa`/`fna` path(s) used, not
  just the config prefix), tool switches, thresholds
- Gene pair counts at each filter stage
- dN, dS, and dN/dS summary statistics (median, mean, min, max)
- Functional classification breakdown with percentages
- dS saturation count with warning if any genes are flagged
- Pseudogene evidence summary
- A checklist of all output files with existence verification

This report is designed to be attached to a manuscript as a supplementary
methods summary or shared with collaborators for a quick overview of
the analysis.

---

## Codeml settings

The PAML codeml control file is at `config/codeml.ctl`. Default settings:

| Parameter  | Value | Meaning                              |
|------------|-------|--------------------------------------|
| runmode    | -2    | Pairwise comparison                  |
| seqtype    | 1     | Codons                               |
| CodonFreq  | 0     | Equal codon frequencies              |
| model      | 0     | One ω for all branches               |
| NSsites    | 0     | No site-specific variation           |
| fix_kappa  | 1     | κ (ts/tv) fixed at 1                 |
| fix_omega  | 0     | ω estimated freely                   |

Modify `config/codeml.ctl` to use different codon frequency models (F1X4: `CodonFreq = 1`, F3X4: `2`, F61: `3`) or to free κ (`fix_kappa = 0`).

---

## Input requirements

- **Genome FASTA (Mode A):** standard multi-FASTA assembly; chromosome/scaffold names must match the GFF `seqname` field.
- **GFF (Mode A):** GFF3 format with `gene` and `CDS` features. Tested with outputs from [MAKER](https://www.yandell-lab.org/software/maker.html) and [BRAKER](https://github.com/Gaius-Augustus/BRAKER).
- **FAA (Mode B):** standard protein FASTA. Headers must be unique and consistent with the FNA file.
- **FNA (Modes B and C):** CDS nucleotide FASTA. Each sequence must be in-frame. Headers must match those in the FAA file (Mode B) or be unique (Mode C).

---

## Dependencies

| Tool       | Version tested | Used for                        | Reference                                                                  |
|------------|----------------|---------------------------------|----------------------------------------------------------------------------|
| Anchorwave | 1.2.3          | CDS extraction (default)        | Song & Zhu (2022) *PNAS* doi:10.1073/pnas.2113075119                      |
| gffread    | 0.12.7         | CDS extraction (alternative)    | Pertea & Pertea (2020) *F1000Research* doi:10.12688/f1000research.23297.2  |
| EMBOSS     | 6.6.0          | CDS translation                 | Rice *et al.* (2000) *Trends Genet* doi:10.1016/S0168-9525(00)02024-2     |
| BLAST+     | 2.15.0         | RBH search (default)            | Camacho *et al.* (2009) *BMC Bioinformatics* doi:10.1186/1471-2105-10-421 |
| DIAMOND    | 2.1.0          | RBH search (alternative)        | Buchfink *et al.* (2021) *Nature Methods* doi:10.1038/s41592-021-01101-x  |
| MAFFT      | 7.520          | Protein alignment               | Katoh & Standley (2013) *MBE* doi:10.1093/molbev/mst010                   |
| pal2nal    | 14             | Codon alignment                 | Suyama *et al.* (2006) *NAR* doi:10.1093/nar/gkl315                       |
| Gblocks    | 0.91b          | Alignment trimming (default)    | Castresana (2000) *MBE* doi:10.1093/oxfordjournals.molbev.a026334          |
| trimAl     | 1.4.1          | Alignment trimming (alternative)| Capella-Gutiérrez *et al.* (2009) *Bioinformatics* doi:10.1093/bioinformatics/btp348 |
| PAML       | 4.9j           | dN/dS estimation                | Yang (2007) *MBE* doi:10.1093/molbev/msm088                               |
| Snakemake  | ≥7.0           | Workflow management             | Mölder *et al.* (2021) *F1000Research* doi:10.12688/f1000research.29032.2  |
| matplotlib | ≥3.7           | Result visualisation            | Hunter (2007) *CSE* doi:10.1109/MCSE.2007.55                               |
| pandas     | ≥1.5           | Data handling                   | McKinney (2010) *Proc. SciPy* doi:10.25080/Majora-92bf1922-00a             |

---

## Troubleshooting

**Conda environment installation fails on first run**  
Snakemake builds conda environments automatically before each rule runs.
If your cluster's conda solver is set to `classic`, this can fail with
solver conflicts, particularly for environments mixing bioconda and
conda-forge packages. The workaround is to switch to `libmamba` as the
solver and pre-install all environments manually before launching the
pipeline:

```bash
# Install libmamba solver if not already available
conda install -n base conda-libmamba-solver
conda config --set solver libmamba

# Pre-install all pipeline environments manually
mamba env create -f envs/anchorwave.yaml
mamba env create -f envs/emboss.yaml
mamba env create -f envs/rbh.yaml
mamba env create -f envs/mafft.yaml
mamba env create -f envs/pal2nal.yaml
mamba env create -f envs/trimmer.yaml
mamba env create -f envs/paml.yaml
mamba env create -f envs/plotting.yaml
```

Once all environments are installed, run the pipeline as normal with
`--use-conda`. Snakemake will detect the existing environments and skip
reinstallation.

**No RBH pairs found**  
Try relaxing `blast.evalue` or `blast.min_cov` in `config/config.yaml`.
Check `output/logs/get_rbh.log` for details on how many hits were found
in each direction.

**Many genes lost at the trimming step**  
Switch to `trimmer: trimal` in `config/config.yaml` for less aggressive
trimming. Genes trimmed to zero length are skipped gracefully in both
cases rather than failing the pipeline.

**codeml output is empty**  
Check that `output/intermediate/aligns/{gene}.pal2nal-gb1.fa` is
non-empty. Very short or highly divergent gene pairs can produce
alignments that are fully trimmed by Gblocks or trimAl.

**Plots not produced / visualisation errors**  
Check `output/logs/plot_results.log`. If all dS values are zero (common
with back-translated test data), placeholder PDFs are written with a
message explaining this. Run with real biological data to produce
meaningful plots.

**Pseudogene check shows all genes as weak evidence**  
This is common when the FNA sequences are complete CDS with intact
reading frames — the ω ≥ 1 signal may reflect positive selection
rather than degeneration, or the sequences may have been repaired
during annotation. Check the raw alignments and consider running
a dedicated pseudogene tool for confirmation.

**Re-running from an intermediate step**  
All intermediate files are kept under `output/intermediate/`. Snakemake's
dependency tracking means if you delete a file and rerun, only steps
downstream of that file are recomputed.

**Questions or issues?**  
Please use the contact form at [suvrathaprasad.github.io](https://suvrathaprasad.github.io/index.html).

---

## License

This pipeline is released under the [Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License (CC BY-NC-ND 4.0)](https://creativecommons.org/licenses/by-nc-nd/4.0/).

© 2025 Suvratha Jayaprasad. You may use and share this pipeline with attribution, but you **may not modify it** and you **may not use it for commercial purposes**. See `LICENSE` for full terms.
