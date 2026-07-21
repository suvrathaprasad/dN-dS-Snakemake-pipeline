"""
collate_results.py — Snakemake script directive for collating codeml outputs.

Parses all per-gene codeml output files and writes a single TSV.
Skips genes where Gblocks trimmed the entire alignment (sentinel files)
without failing the pipeline.

Gene pairs are excluded from dnds_output.tsv for three distinct reasons,
which are counted and reported separately rather than lumped into one
"skipped" total:

  - sentinel   — the trimmer legitimately removed the whole alignment
                 (expected, benign, see run_trimmer.py)
  - unparsed   — the codeml output file didn't have the expected structure
                 (no pairwise-comparison section, or no gene-pair/t= line
                 found) — usually means codeml itself failed or produced
                 something unexpected
  - bad_labels — a "t=" line was found, but the tokens immediately before
                 the dN/dS, dN and dS values weren't the labels we expect
                 at those fixed positions. The value extraction below
                 relies on codeml's pairwise output always tokenising to
                 the same fixed positions (a known PAML fragility point);
                 this check makes sure that assumption still holds before
                 trusting the numbers, rather than silently writing
                 misaligned values into the results table.

Only "sentinel" is an expected, silent outcome. "unparsed" and
"bad_labels" are logged with the offending file/line so a real problem
(bad install, PAML version drift, etc.) doesn't get mistaken for normal
trimming attrition. All three categories are also written out to
genes_dnds_skipped.tsv (see below), rather than only being visible as a
bare count in the log.

dnds_output.tsv also carries three alignment-quality columns (pident,
query_coverage, target_coverage) joined in from rbh_pairs.tsv by gene ID.
These were already computed during the RBH search and coverage filter in
get_rbh.py but previously discarded afterwards — keeping them lets a
reviewer sanity-check whether an outlier ω is riding on a well-covered,
high-identity alignment or a marginal one, without having to separately
dig through the RBH intermediate file by hand.

genes_dnds_skipped.tsv lists every gene pair that never made it into
dnds_output.tsv, across all three skip categories above, with its real
(Gene_query, Gene_target) IDs. A skipped codeml output file has no gene
names in it (there was no alignment left to name, for a sentinel skip),
so those IDs are recovered by cross-referencing the file's position in
the gene list back to the matching row of rbh_pairs.tsv — the same
1-based numbering the Snakefile itself uses when creating each gene's
{gene}/{gene}_codeml.txt directory (see get_genes() and
extract_gene_index() below).
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
input_files  = list(snakemake.input.codeml_outputs)
rbh_tsv      = snakemake.input.rbh_tsv
output_file  = Path(snakemake.output[0])
sat_tsv      = Path(snakemake.output[1])
skip_tsv     = Path(snakemake.output[2])
trimmer      = snakemake.config.get("tools", {}).get("trimmer", "gblocks")
ds_sat_threshold = float(
    snakemake.config.get("dS_saturation_threshold", 2.0)
)

output_file.parent.mkdir(parents=True, exist_ok=True)


def load_rbh_table(path: str):
    """
    Read rbh_pairs.tsv once, returning two lookups built from the same
    pass over the file:

      by_pair:  {(query, target): (pident, query_coverage, target_coverage)}
                — used to join alignment-quality columns onto genes that
                DID make it into dnds_output.tsv.

      by_index: {1-based row number: (query, target)}
                — used to recover the real gene IDs for genes that were
                SKIPPED before ever reaching a parseable codeml output,
                by matching a skipped file's position in the gene list
                back to this same row number.

    The numbering in by_index deliberately mirrors the Snakefile's own
    get_genes() exactly — every non-header line counts, including any
    blank line, since get_genes() doesn't special-case blank lines
    either. Building both lookups from a single read (rather than two
    separate passes, or worse, two different line-selection rules that
    could quietly drift apart) is what keeps this numbering guaranteed
    consistent with the Snakefile's.
    """
    with open(path) as fh:
        all_lines = fh.readlines()

    header = None
    data_lines = []
    for line in all_lines:
        if line.startswith("Query"):
            header = line.rstrip("\n").split("\t")
            continue
        data_lines.append(line)

    col_idx = {name: i for i, name in enumerate(header)} if header else {}
    has_cov = "query_coverage" in col_idx and "target_coverage" in col_idx

    by_pair = {}
    by_index = {}
    for i, line in enumerate(data_lines):
        gene_index = i + 1  # matches get_genes()'s f"gene{i+1}" numbering exactly
        parts = line.rstrip("\n").split("\t") if line.strip() else []
        query = target = None
        if header:
            try:
                query = parts[col_idx["Query"]]
                target = parts[col_idx["Target"]]
            except (IndexError, KeyError):
                pass
        by_index[gene_index] = (query, target)
        if query is not None and target is not None:
            try:
                pident = parts[col_idx["pident"]]
                qcov = parts[col_idx["query_coverage"]] if has_cov else "NA"
                tcov = parts[col_idx["target_coverage"]] if has_cov else "NA"
            except (IndexError, KeyError):
                pident = qcov = tcov = "NA"
            by_pair[(query, target)] = (pident, qcov, tcov)

    return by_pair, by_index


def extract_gene_index(path: str):
    """
    Recover the 1-based gene index from a codeml output file's path
    (".../codeml/gene42/gene42_codeml.txt" -> 42) — the parent directory
    name is exactly the Snakefile's {gene} wildcard value. Returns None
    if the path doesn't match the expected "geneN" pattern, rather than
    raising, so a genuinely unexpected path just becomes an "unknown"
    gene in the output table instead of crashing the whole collation.
    """
    match = re.match(r"^gene(\d+)$", Path(path).parent.name)
    return int(match.group(1)) if match else None


def split_skip_reason(reason: str):
    """
    Split parse_codeml()'s reason string into (category, detail) for a
    cleaner TSV — e.g. "unparsed (no 'pairwise comparison' section
    found)" -> ("unparsed", "no 'pairwise comparison' section found").
    "sentinel" has no detail.
    """
    if reason == "sentinel":
        return "sentinel", ""
    category, _, rest = reason.partition(" ")
    detail = rest.strip()
    if detail.startswith("(") and detail.endswith(")"):
        detail = detail[1:-1]
    return category, detail


rbh_lookup, rbh_by_index = load_rbh_table(rbh_tsv)

# Expected labels immediately preceding the values we pull out of codeml's
# "t= ... dN/dS= ... dN = ... dS = ..." line, once all "=" signs are
# stripped and the line is split on whitespace. See parse_codeml().
EXPECTED_LABELS = {6: "dN/dS", 8: "dN", 10: "dS"}


def parse_codeml(path: str):
    """
    Parse a single codeml output file.

    Returns (rows, reason) where rows is a list of
    (gene1, gene2, t, dN, dS, dNdS) tuples (empty if nothing was parsed),
    and reason is None on success or one of "sentinel", "unparsed",
    "bad_labels" describing why nothing was parsed.
    """
    try:
        with open(path) as fh:
            text = fh.read()
    except OSError as exc:
        return [], f"unparsed (could not read file: {exc})"

    # Skip sentinel files written when the trimmer removed all positions —
    # this is the normal, expected "nothing to parse" case.
    if text.strip().startswith("Skipped:"):
        return [], "sentinel"

    results = []
    status = False
    gene1 = gene2 = None
    saw_pairwise_section = False
    saw_t_line = False
    bad_label_line = None

    for line in text.splitlines():
        if line.startswith("pairwise comparison, codon frequencies"):
            status = True
            saw_pairwise_section = True

        if status:
            if line and line[0].isdigit():
                # Gene pair header line. For a 2-sequence pairwise alignment,
                # PAML/codeml always prints the HIGHER-indexed sequence
                # first: "2 (seq2) ... 1 (seq1)", not "1 (seq1) ... 2 (seq2)"
                # as the token positions might suggest at a glance. Since
                # extract_gene_pair.py always writes the query as sequence 1
                # and the target as sequence 2 into the alignment codeml
                # sees, that means the token right after "1" (parts[4]) is
                # always the query, and the token right after "2" (parts[1])
                # is always the target — the reverse of their positions in
                # the line. Confirmed against real codeml output; this is a
                # deterministic PAML convention for pairwise mode, not a
                # per-file quirk needing a conditional check.
                clean = re.sub(r"[()]", "", line)
                parts = clean.split()
                try:
                    gene1 = parts[4].split("...")[0]  # query  (sequence 1)
                    gene2 = parts[1].split("...")[0]  # target (sequence 2)
                except IndexError:
                    continue

            if line.startswith("t=") and gene1 and gene2:
                saw_t_line = True
                clean = re.sub(r"=", "", line)
                tokens = re.split(r"\s+", clean.strip())

                # Validate the fixed-position labels before trusting the
                # values next to them — protects against silently writing
                # misaligned numbers if codeml's output format ever shifts.
                mismatches = [
                    f"tokens[{idx}]='{tokens[idx] if idx < len(tokens) else '<missing>'}' "
                    f"(expected '{label}')"
                    for idx, label in EXPECTED_LABELS.items()
                    if idx >= len(tokens) or tokens[idx] != label
                ]
                if mismatches:
                    bad_label_line = f"{line.strip()}  [{'; '.join(mismatches)}]"
                    continue

                try:
                    t    = tokens[1]
                    dnds = tokens[7]
                    dn   = tokens[9]
                    ds   = tokens[11]
                    results.append((gene1, gene2, t, dn, ds, dnds))
                except IndexError:
                    continue

    if results:
        return results, None
    if bad_label_line:
        return [], f"bad_labels ({bad_label_line})"
    if not saw_pairwise_section:
        return [], "unparsed (no 'pairwise comparison' section found)"
    if not saw_t_line:
        return [], "unparsed (pairwise section found, but no 't=' line)"
    return [], "unparsed (unknown reason)"


# ── Collate ───────────────────────────────────────────────────────────────────
rows = []
skipped_rows = []  # (gene_index, Gene_query, Gene_target, category, detail, file)
n_sentinel = 0
n_unparsed = 0
n_no_rbh_match = 0
n_skip_no_gene_index = 0
n_skip_no_rbh_row = 0
unparsed_details = []

for f in input_files:
    parsed, reason = parse_codeml(f)
    if parsed:
        rows.extend(parsed)
        continue

    if reason == "sentinel":
        n_sentinel += 1
    else:
        n_unparsed += 1
        unparsed_details.append((f, reason))

    gene_index = extract_gene_index(f)
    if gene_index is None:
        n_skip_no_gene_index += 1
        gq, gt = "unknown", "unknown"
    else:
        gq, gt = rbh_by_index.get(gene_index, (None, None))
        if gq is None or gt is None:
            n_skip_no_rbh_row += 1
            gq, gt = "unknown", "unknown"

    category, detail = split_skip_reason(reason)
    skipped_rows.append((
        str(gene_index) if gene_index is not None else "unknown",
        gq, gt, category, detail, f,
    ))

# Join in alignment-quality columns from rbh_pairs.tsv by (gene1, gene2).
# This should always match — every gene pair reaching codeml originated
# from an RBH row — so a miss here would indicate an ID mismatch
# somewhere upstream (e.g. header mangling) rather than a normal outcome,
# and is worth surfacing rather than silently writing blank columns.
enriched_rows = []
for gene1, gene2, t, dn, ds, dnds in rows:
    quality = rbh_lookup.get((gene1, gene2))
    if quality is None:
        n_no_rbh_match += 1
        pident, qcov, tcov = "NA", "NA", "NA"
    else:
        pident, qcov, tcov = quality
    enriched_rows.append((gene1, gene2, t, dn, ds, dnds, pident, qcov, tcov))

tool_label = "Gblocks" if trimmer == "gblocks" else "trimAl"
print(f"[collate_results] Parsed {len(rows)} gene pairs")
print(f"[collate_results]   {n_sentinel} skipped — empty alignment after {tool_label} (expected)")
if n_unparsed:
    print(f"[collate_results]   {n_unparsed} skipped — codeml output did NOT match the expected "
          f"format (unexpected — see details below):")
    for f, reason in unparsed_details:
        print(f"[collate_results]     {f}: {reason}")
if n_no_rbh_match:
    print(f"[collate_results]   WARNING: {n_no_rbh_match} gene pair(s) parsed from codeml output "
          f"had no matching row in {rbh_tsv} — alignment quality columns are 'NA' for these. "
          f"This is unexpected; check for gene ID mismatches upstream.")
if n_skip_no_gene_index:
    print(f"[collate_results]   WARNING: {n_skip_no_gene_index} skipped file(s) had a path that "
          f"didn't match the expected 'geneN' pattern — their IDs are 'unknown' in "
          f"{skip_tsv.name}. This is unexpected; check for a Snakefile path convention change.")
if n_skip_no_rbh_row:
    print(f"[collate_results]   WARNING: {n_skip_no_rbh_row} skipped file(s) had a gene index with "
          f"no matching row in {rbh_tsv} — their IDs are 'unknown' in {skip_tsv.name}. "
          f"This is unexpected; check that rbh_pairs.tsv wasn't modified after the gene list "
          f"was fixed for this run.")

with open(output_file, "w") as fh:
    fh.write("Gene_query\tGene_target\tt\tdN\tdS\tdNdS\tpident\tquery_coverage\ttarget_coverage\n")
    for gene1, gene2, t, dn, ds, dnds, pident, qcov, tcov in enriched_rows:
        fh.write(f"{gene1}\t{gene2}\t{t}\t{dn}\t{ds}\t{dnds}\t{pident}\t{qcov}\t{tcov}\n")

print(f"[collate_results] Written: {output_file}")

# ── dS saturation check ───────────────────────────────────────────────────────
# Synonymous sites approach saturation at dS > threshold (default 2.0).
# At high dS, multiple substitutions at the same site become likely,
# making dS an unreliable proxy for divergence time and dN/dS estimates
# potentially unreliable. These genes are flagged but kept in dnds_output.tsv.
saturated = []
for row in enriched_rows:
    ds = row[4]
    try:
        if float(ds) > ds_sat_threshold:
            saturated.append(row)
    except ValueError:
        pass

sat_tsv.parent.mkdir(parents=True, exist_ok=True)
with open(sat_tsv, "w") as fh:
    fh.write("Gene_query\tGene_target\tt\tdN\tdS\tdNdS\tpident\tquery_coverage\ttarget_coverage\n")
    for row in saturated:
        fh.write("\t".join(row) + "\n")

if saturated:
    print(
        f"[collate_results] WARNING: {len(saturated)} gene pair(s) have "
        f"dS > {ds_sat_threshold} (synonymous saturation threshold). "
        f"dN/dS estimates for these genes may be unreliable. "
        f"See: {sat_tsv}"
    )
else:
    print(f"[collate_results] No genes exceed dS saturation threshold ({ds_sat_threshold})")

print(f"[collate_results] Saturation check written: {sat_tsv}")

# ── Skipped genes table ───────────────────────────────────────────────────────
# Every gene pair that never made it into dnds_output.tsv, across all
# three skip categories from parse_codeml() above — sentinel (expected,
# benign trimming attrition), unparsed, and bad_labels (both genuinely
# unexpected). Sorted by gene_index purely for readability; "unknown"
# sorts after all numeric indices via the key function below.
skip_tsv.parent.mkdir(parents=True, exist_ok=True)


def _skip_sort_key(row):
    idx = row[0]
    return (0, int(idx)) if idx.isdigit() else (1, idx)


with open(skip_tsv, "w") as fh:
    fh.write("gene_index\tGene_query\tGene_target\tskip_category\tskip_detail\tcodeml_output_file\n")
    for gene_index, gq, gt, category, detail, f in sorted(skipped_rows, key=_skip_sort_key):
        fh.write(f"{gene_index}\t{gq}\t{gt}\t{category}\t{detail}\t{f}\n")

print(f"[collate_results] {len(skipped_rows)} skipped gene pair(s) written: {skip_tsv}")
