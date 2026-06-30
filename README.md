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
 ┌─ Mode A: Genome + GFF ───────┐   ┌─ Mode B: FAA + FNA ──────┐   ┌─ Mode C: FNA only ───────┐
 │                              │   │                          │   │                          │
 │  Anchorwave or gffread       │   │  (skip extraction and    │   │  EMBOSS transeq          │
 │    Extract CDS per species   │   │     translation)         │   │    Translate → protein   │
 │      │                       │   │                          │   │                          │
 │      ▼                       │   │                          │   │                          │
 │  EMBOSS transeq → Protein    │   │                          │   │                          │
 └──────────────┬───────────────┘   └────────────┬─────────────┘   └────────────┬─────────────┘
                │                                │                              │
                └──────────────────────────┬─────┴──────────────────────────────┘
                                           │
                                           ▼
          BLAST+ or DIAMOND   Bidirectional search + Reciprocal Best Hits
                                           │
                                           ▼
          MAFFT                Protein alignment per gene pair
                                           │
                                           ▼
          pal2nal              Back-translate → codon alignment
                                           │
                                           ▼
          Gblocks or trimAl    Trim poorly aligned regions
                                           │
                                           ▼
          PAML codeml          Pairwise dN/dS (runmode = -2)
                                           │
                                           ▼
                             output/results/dnds_output.tsv
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
│   ├── input/         Per-gene-pair FAA and FNA files
│   ├── aligns/        MAFFT, pal2nal, and trimmed alignments
│   └── codeml/        Per-gene codeml output directories
├── results/
│   └── dnds_output.tsv   ← FINAL RESULT
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
| Anchorwave | 1.2.3          | CDS extraction (default)        | Song & Zhu (2022) *PNAS* doi:10.1073/pnas.2113075119                       |
| gffread    | 0.12.7         | CDS extraction (alternative)    | Pertea & Pertea (2020) *F1000Research* doi:10.12688/f1000research.23297.2  |
| EMBOSS     | 6.6.0          | CDS translation                 | Rice *et al.* (2000) *Trends Genet* doi:10.1016/S0168-9525(00)02024-2      |
| BLAST+     | 2.15.0         | RBH search (default)            | Camacho *et al.* (2009) *BMC Bioinformatics* doi:10.1186/1471-2105-10-421  |
| DIAMOND    | 2.1.0          | RBH search (alternative)        | Buchfink *et al.* (2021) *Nature Methods* doi:10.1038/s41592-021-01101-x   |
| MAFFT      | 7.520          | Protein alignment               | Katoh & Standley (2013) *MBE* doi:10.1093/molbev/mst010                    |
| pal2nal    | 14             | Codon alignment                 | Suyama *et al.* (2006) *NAR* doi:10.1093/nar/gkl315                        |
| Gblocks    | 0.91b          | Alignment trimming (default)    | Castresana (2000) *MBE* doi:10.1093/oxfordjournals.molbev.a026334          |
| trimAl     | 1.4.1          | Alignment trimming (alternative)| Capella-Gutiérrez *et al.* (2009) *Bioinformatics* doi:10.1093/bioinformatics/btp348 |
| PAML       | 4.9j           | dN/dS estimation                | Yang (2007) *MBE* doi:10.1093/molbev/msm088                                |
| Snakemake  | ≥7.0           | Workflow management             | Mölder *et al.* (2021) *F1000Research* doi:10.12688/f1000research.29032.2  |

---

## Troubleshooting

**No RBH pairs found**  
Try relaxing `blast.evalue` or `blast.min_cov` in `config/config.yaml`. Check `output/logs/get_rbh.log` for details.

**Many genes lost at the trimming step**  
Switch to `trimmer: trimal` in `config/config.yaml` for less aggressive trimming. Genes trimmed to zero length are skipped gracefully in both cases.

**codeml output is empty**  
Check that `intermediate/aligns/{gene}.pal2nal-gb1.fa` is non-empty. Very short or highly divergent gene pairs can produce alignments that are fully trimmed.

**Re-running from an intermediate step**  
All intermediate files are kept under `output/intermediate/`. Snakemake's dependency tracking means if you delete a file and rerun, only steps downstream of that file are recomputed.

**Questions or issues?**  
Please use the contact form at [suvrathaprasad.github.io](https://suvrathaprasad.github.io/index.html).

---

## License

This pipeline is released under the [Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License (CC BY-NC-ND 4.0)](https://creativecommons.org/licenses/by-nc-nd/4.0/).

© 2025 Suvratha Jayaprasad. You may use and share this pipeline with attribution, but you **may not modify it** and you **may not use it for commercial purposes**. See `LICENSE` for full terms.
