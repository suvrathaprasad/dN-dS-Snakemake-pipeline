# Test Dataset

A minimal test dataset for verifying the pipeline before running on your own data.

## Contents

```
test/
├── data/
│   ├── species1.faa   # 20 protein sequences (P24XY grasshopper, scaffold 6)
│   ├── species1.fna   # Corresponding CDS nucleotide sequences
│   ├── species2.faa   # 20 diverged protein sequences (~5% amino acid divergence)
│   └── species2.fna   # Corresponding back-translated CDS nucleotide sequences
└── config_test.yaml   # Config in Mode B (pre-made FAA + FNA)
```

Species2 sequences are derived from species1 by introducing ~5% random amino
acid substitutions, with independent back-translation to nucleotide. This
guarantees reciprocal BLAST hits exist while simulating real between-species
divergence. dS values will be near zero because there is no synonymous
divergence in back-translated sequences — the test is for toolchain
verification only, not biological interpretation.

## Running the test

From the **repository root**:

```bash
snakemake --snakefile workflow/Snakefile \
          --configfile test/config_test.yaml \
          --use-conda \
          --cores 4
```

Expected outputs:
- `test/output/results/dnds_output.tsv` — dN/dS table
- `test/output/results/plots/` — four PDF figures (boxplot, violin, scatter, summary)
- `test/output/results/tables/` — ten TSV gene-list files

> **Note on plots with test data:** the test sequences are back-translations
> of proteins so dS values will be near zero. All genes will fall into the
> "Undefined (dS=0)" category and the pipeline will generate placeholder
> PDFs with a message explaining this. This is expected — the test validates
> the toolchain end to end, not the biological output. Run with real data
> to produce meaningful plots and tables.

The test uses **Mode B** (pre-made FAA + FNA), so CDS extraction and
translation are skipped automatically. No genome assembly or GFF file is needed.

## Testing alternative tools

`config_test.yaml` includes a `tools:` section with the same switches available
in the main pipeline config. Since this test runs in Mode B, `cds_extraction`
is not used, but `search_method` and `trimmer` can both be tested on the
small dataset before committing to a choice for your real data:

```yaml
tools:
  search_method: "diamond"   # try DIAMOND instead of BLAST+
  trimmer:       "trimal"    # try trimAl instead of Gblocks
```

Edit `config_test.yaml` directly, or pass an override on the command line:

```bash
snakemake --snakefile workflow/Snakefile \
          --configfile test/config_test.yaml \
          --config tools='{"search_method": "diamond", "trimmer": "trimal"}' \
          --use-conda \
          --cores 4
```

Both runs should produce the same set of RBH gene pairs and comparable
dN/dS estimates, since the alternative tools are designed as drop-in
replacements rather than different analyses.

---

© 2025 Suvratha Jayaprasad — CC BY-NC-ND 4.0
Contact: https://suvrathaprasad.github.io/index.html
