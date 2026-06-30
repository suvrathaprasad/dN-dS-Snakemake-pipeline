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

Expected output: `test/output/results/dnds_output.txt`

The test uses **Mode B** (pre-made FAA + FNA), so Anchorwave and transeq
are skipped automatically. No genome assembly or GFF file is needed.

---

© 2025 Suvratha Jayaprasad — CC BY-NC-ND 4.0
Contact: https://suvrathaprasad.github.io/index.html
