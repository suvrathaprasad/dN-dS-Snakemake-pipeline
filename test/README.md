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

Running the pipeline will also create `*.offset_idx.json` cache files
alongside the files in `data/` (cached FASTA indexes, rebuilt automatically
if the source file changes). These are gitignored — safe to ignore or
delete.

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
- `test/output/results/plots/` — five PDF figures (boxplot, violin, scatter, functional summary, pseudogene evidence)
- `test/output/results/tables/` — eleven TSV gene-list files (ten from
  functional classification, plus `genes_ds_saturated.tsv`)
- `test/output/results/tables/genes_degenerate_annotated.tsv` — pseudogene sequence evidence
- `test/output/results/run_summary.pdf` — one-page summary of the whole run;
  the fastest way to confirm the toolchain ran end to end without opening
  every individual file

> **Note on plots with test data:** the test sequences are back-translations
> of proteins so dS values will be near zero. All genes will fall into the
> "Undefined (dS=0)" category and the pipeline will generate placeholder
> PDFs with a message explaining this. This is expected — the test validates
> the toolchain end to end, not the biological output. Run with real data
> to produce meaningful plots and tables.
>
> The same applies downstream: with no genes classified as "Degenerate"
> (ω ≥ 1), `genes_degenerate.tsv` and `genes_degenerate_annotated.tsv` will
> be empty and `pseudogene_evidence.pdf` will show a placeholder message
> rather than real evidence counts. `run_summary.pdf` will still generate
> normally and correctly report zero degenerate candidates — it's designed
> to handle this case gracefully, not treat it as a failure.

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
