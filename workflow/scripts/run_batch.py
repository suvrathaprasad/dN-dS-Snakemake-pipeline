#!/usr/bin/env python3
"""
run_batch.py — Run the dN/dS pipeline for one reference species against
multiple target species.

This is a thin wrapper, not a Snakemake rule or a second Snakefile: for
each target, it generates an ordinary single-pair config.yaml (reference
as "query", that target as "target") and invokes the real, unmodified
Snakefile once, sequentially, exactly as if you'd run that pairwise
comparison by hand. The Snakefile itself has no batch-mode awareness at
all — every existing rule, and the mode-resolution logic (Mode B > C > A,
independently per species) that already lives there, works completely
unchanged, because from the Snakefile's point of view each invocation is
just a normal single-pair run.

Why a subprocess-per-pair wrapper rather than a native Snakemake batch
rule: the comparisons are genuinely independent (nothing to share in one
DAG beyond the reference's CDS extraction, a small and bounded cost),
failures stay isolated and debuggable (if pair 2 fails, `cd` into its own
output directory and rerun that exact single-pair snakemake command by
hand — it's an ordinary run, unrelated to the other pairs), and it adds
zero risk to the already-tested core Snakefile, since this script never
modifies or imports it — only shells out to the real `snakemake`
executable, the same way a person would.

Usage:
  python3 workflow/scripts/run_batch.py \
      --configfile batch_config.yaml \
      --cores 16

Run from the repository root (same requirement as running the Snakefile
directly) — paths in batch_config.yaml, and the "workflow/Snakefile"
path this script invokes, are both relative to the current directory.
"""

# =============================================================================
# Author:  Suvratha Jayaprasad
# Contact: https://suvrathaprasad.github.io/index.html
# License: CC BY-NC-ND 4.0 — you may use and share with attribution,
#          but you may not modify or use for commercial purposes.
#          © 2025 Suvratha Jayaprasad. All rights reserved.
#          Full terms: https://creativecommons.org/licenses/by-nc-nd/4.0/
# =============================================================================

import argparse
import subprocess
import sys
from datetime import datetime
from pathlib import Path

try:
    import yaml
except ImportError:
    sys.exit(
        "run_batch.py requires PyYAML, which isn't installed in this "
        "environment. Install it with:\n"
        "  pip install pyyaml\n"
        "or add it to whichever conda environment you're running "
        "run_batch.py from."
    )


# =============================================================================
# Config loading and validation
# =============================================================================

def resolve_mode(spec: dict) -> str:
    """
    Same Mode B > Mode C > Mode A priority as the main Snakefile, applied
    to one species' config block (reference, or one target).
    """
    if "faa" in spec and "fna" in spec:
        return "B"
    if "fna" in spec:
        return "C"
    if "fasta" in spec and "gff" in spec:
        return "A"
    return None


def mode_file_checks(spec: dict, mode: str) -> list:
    """Return [(label, path), ...] of the raw, user-supplied files to check
    for this species' resolved mode."""
    if mode == "A":
        return [("fasta", spec.get("fasta")), ("gff", spec.get("gff"))]
    if mode == "B":
        return [("faa", spec.get("faa")), ("fna", spec.get("fna"))]
    if mode == "C":
        return [("fna", spec.get("fna"))]
    return []


def validate_prefix(label: str, prefix, errors: list) -> None:
    if not prefix or any(c in prefix for c in " /\\\t\n"):
        errors.append(
            f"  {label} = '{prefix}' — must be non-empty and contain "
            f"no spaces or path separators (used directly in output "
            f"directory/file names)"
        )


def validate_species_files(label: str, spec: dict, errors: list) -> None:
    mode = resolve_mode(spec)
    if mode is None:
        errors.append(
            f"  {label}: no valid input mode found — needs either "
            f"fasta+gff (Mode A), faa+fna (Mode B), or fna (Mode C)"
        )
        return

    for file_label, path in mode_file_checks(spec, mode):
        if not path:
            errors.append(f"  {label}.{file_label} is not set (required for Mode {mode})")
        elif not Path(path).is_file():
            errors.append(f"  {label}.{file_label} does not exist: '{path}'")


def validate_species(label: str, spec: dict, errors: list) -> None:
    """Validate one species block (a target, or the reference when no
    target overrides it), appending any problems found to errors rather
    than raising immediately — so every problem across the whole batch
    can be reported together."""
    validate_prefix(f"{label}.prefix", spec.get("prefix"), errors)
    validate_species_files(label, spec, errors)


def validate_batch_config(cfg: dict) -> list:
    """
    Validate the whole batch config, collecting every problem found
    rather than stopping at the first one — the same principle as the
    main Snakefile's own config validation, applied here across the
    whole target list before any pair is launched. Returns a list of
    error strings; empty means the config is valid.
    """
    errors = []

    reference = cfg.get("reference")
    if not reference:
        errors.append("  'reference' block is missing")
        reference = {}
    else:
        # Prefix is always required — it names every comparison's output
        # directory and is used as query.prefix in every generated
        # per-pair config, even for a target whose reference_override
        # points at a completely different file (e.g. a scaffold-
        # restricted extract of the same reference species).
        validate_prefix("reference.prefix", reference.get("prefix"), errors)

    targets = cfg.get("targets")
    if not targets:
        errors.append("  'targets' list is missing or empty")
        targets = []

    # The top-level reference's own files are only actually required if
    # at least one target doesn't supply its own reference_override —
    # if every target overrides it, the top-level files (if any) are
    # never used by anything and shouldn't be required to exist.
    reference_files_needed = False

    for i, target in enumerate(targets):
        validate_species(f"targets[{i}]", target, errors)

        override = target.get("reference_override")
        if override:
            validate_species_files(f"targets[{i}].reference_override", override, errors)
        else:
            reference_files_needed = True

    if reference and reference_files_needed:
        validate_species_files("reference", reference, errors)

    # Prefix collisions — checked even if individual species validation
    # above already failed for one of them, since this is a cheap,
    # independent check worth surfacing regardless. reference_override
    # doesn't introduce a separate prefix to check here: it's still the
    # same reference species for naming purposes, just a different file
    # for that one comparison.
    if reference and targets:
        all_prefixes = [("reference", reference.get("prefix"))]
        all_prefixes += [(f"targets[{i}]", t.get("prefix")) for i, t in enumerate(targets)]
        seen = {}
        for label, prefix in all_prefixes:
            if not prefix:
                continue  # already reported above
            if prefix in seen:
                errors.append(
                    f"  Prefix collision: '{prefix}' is used by both "
                    f"{seen[prefix]} and {label} — every target's prefix "
                    f"must be unique, and different from the reference's, "
                    f"since it's used directly in the output directory "
                    f"name for that comparison"
                )
            else:
                seen[prefix] = label

    if not cfg.get("outdir"):
        errors.append("  'outdir' is missing")

    blast_cfg = cfg.get("blast", {})
    min_cov = blast_cfg.get("min_cov", 60)
    try:
        if not (0 <= float(min_cov) <= 100):
            errors.append(f"  blast.min_cov = {min_cov} — must be between 0 and 100")
    except (TypeError, ValueError):
        errors.append(f"  blast.min_cov = '{min_cov}' — must be numeric")

    evalue = blast_cfg.get("evalue")
    try:
        if evalue is None or float(evalue) <= 0:
            errors.append(f"  blast.evalue = {evalue} — must be a positive number")
    except (TypeError, ValueError):
        errors.append(f"  blast.evalue = '{evalue}' — must be numeric")

    ds_threshold = cfg.get("dS_saturation_threshold", 2.0)
    try:
        if float(ds_threshold) <= 0:
            errors.append(
                f"  dS_saturation_threshold = {ds_threshold} — must be a positive number"
            )
    except (TypeError, ValueError):
        errors.append(f"  dS_saturation_threshold = '{ds_threshold}' — must be numeric")

    for label, val in (
        ("blast.threads", blast_cfg.get("threads", 1)),
        ("mafft.threads", cfg.get("mafft", {}).get("threads", 1)),
    ):
        if not isinstance(val, int) or val < 1:
            errors.append(f"  {label} = {val} — must be a positive integer")

    return errors


# =============================================================================
# Per-pair config generation
# =============================================================================

def build_pair_config(batch_cfg: dict, reference: dict, target: dict, pair_outdir: str) -> dict:
    """
    Build one ordinary single-pair config (same shape as config.yaml)
    for this (reference, target) comparison. Everything under blast/
    mafft/tools/dS_saturation_threshold is copied straight from the
    batch config — shared across every comparison in the batch.

    If this target has a reference_override block, its files are used
    for query instead of the top-level reference's — e.g. when each
    comparison needs a different, pre-restricted extract of the
    reference (a single scaffold/chromosome relevant to that particular
    target), rather than the same reference file for every comparison.
    The reference's prefix is always used regardless, since it's still
    logically the same species in every comparison — only which file(s)
    represent it for this one comparison changes.
    """
    override = target.get("reference_override")
    query_spec = dict(override) if override else dict(reference)
    query_spec["prefix"] = reference["prefix"]

    target_spec = {k: v for k, v in target.items() if k != "reference_override"}

    return {
        "query":  query_spec,
        "target": target_spec,
        "outdir": pair_outdir,
        "blast":  batch_cfg.get("blast", {}),
        "mafft":  batch_cfg.get("mafft", {}),
        "dS_saturation_threshold": batch_cfg.get("dS_saturation_threshold", 2.0),
        "tools":  batch_cfg.get("tools", {}),
    }


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Run the dN/dS pipeline for one reference against multiple targets."
    )
    parser.add_argument("--configfile", required=True, help="Path to batch_config.yaml")
    parser.add_argument("--cores", required=True, type=int, help="Cores passed to each snakemake run")
    parser.add_argument("--no-use-conda", action="store_true",
                         help="Don't pass --use-conda to snakemake (default: pass it)")
    parser.add_argument("-n", "--dry-run", action="store_true",
                         help="Pass -n/--dry-run through to snakemake for every pair")
    parser.add_argument("--snakemake-args", default="",
                         help="Extra arguments appended verbatim to every snakemake "
                              "invocation, e.g. --snakemake-args=\"--rerun-incomplete\"")
    args = parser.parse_args()

    configfile_path = Path(args.configfile)
    if not configfile_path.is_file():
        sys.exit(f"Batch config file not found: {args.configfile}")

    with open(configfile_path) as fh:
        batch_cfg = yaml.safe_load(fh)

    errors = validate_batch_config(batch_cfg)
    if errors:
        sys.exit(
            "Batch config validation failed — fix the following in "
            f"{args.configfile} before running:\n" + "\n".join(errors)
        )

    reference = batch_cfg["reference"]
    targets = batch_cfg["targets"]
    batch_outdir = Path(batch_cfg["outdir"])
    generated_dir = batch_outdir / "_generated_configs"
    generated_dir.mkdir(parents=True, exist_ok=True)

    batch_log_path = batch_outdir / "batch_run.log"
    batch_outdir.mkdir(parents=True, exist_ok=True)
    batch_log = open(batch_log_path, "a")

    def log(msg: str) -> None:
        print(msg)
        batch_log.write(msg + "\n")
        batch_log.flush()

    log(f"\n{'=' * 80}")
    log(f"[run_batch] Batch started {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    log(f"[run_batch] Reference: {reference['prefix']}")
    log(f"[run_batch] Targets ({len(targets)}): {', '.join(t['prefix'] for t in targets)}")
    log(f"[run_batch] Output root: {batch_outdir}")
    log(f"{'=' * 80}\n")

    results = []  # (pair_name, pair_outdir, success: bool, generated_config: Path)

    for i, target in enumerate(targets, start=1):
        pair_name = f"{reference['prefix']}_vs_{target['prefix']}"
        pair_outdir = str(batch_outdir / pair_name)
        generated_config_path = generated_dir / f"{pair_name}.yaml"

        pair_config = build_pair_config(batch_cfg, reference, target, pair_outdir)
        with open(generated_config_path, "w") as fh:
            yaml.safe_dump(pair_config, fh, sort_keys=False)

        log(f"{'-' * 80}")
        log(f"[run_batch] [{i}/{len(targets)}] Starting: {pair_name}")
        log(f"[run_batch]   Config:    {generated_config_path}")
        log(f"[run_batch]   Output:    {pair_outdir}")
        log(f"[run_batch]   Started:   {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        log(f"{'-' * 80}")

        cmd = [
            "snakemake",
            "--snakefile", "workflow/Snakefile",
            "--configfile", str(generated_config_path),
            "--cores", str(args.cores),
        ]
        if not args.no_use_conda:
            cmd.append("--use-conda")
        if args.dry_run:
            cmd.append("-n")
        if args.snakemake_args:
            cmd.extend(args.snakemake_args.split())

        batch_log.flush()
        # Deliberately NOT capturing stdout/stderr — a real batch run can
        # take hours per pair, and streaming snakemake's own live progress
        # straight through (rather than silently buffering it until the
        # whole pair finishes) is the only usable experience at that
        # timescale. Only the orchestration-level messages above/below go
        # into batch_run.log; each pair's own rule-level logs already land
        # under its own output directory exactly as in a normal run.
        result = subprocess.run(cmd)

        success = result.returncode == 0
        results.append((pair_name, pair_outdir, success, generated_config_path))

        status = "SUCCEEDED" if success else "FAILED"
        log(f"[run_batch]   Finished:  {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} — {status}")
        if not success:
            log(f"[run_batch]   To debug this pair on its own, rerun:")
            log(f"[run_batch]     snakemake --snakefile workflow/Snakefile "
                f"--configfile {generated_config_path} --use-conda --cores {args.cores}")
        log("")

    # ── Summary ───────────────────────────────────────────────────────────────
    n_ok = sum(1 for _, _, ok, _ in results if ok)
    log(f"{'=' * 80}")
    log(f"[run_batch] Batch complete: {n_ok}/{len(results)} comparisons succeeded")
    log(f"{'=' * 80}")
    for pair_name, pair_outdir, ok, _ in results:
        mark = "✓" if ok else "✗"
        log(f"  {mark} {pair_name}  ->  {pair_outdir}")
    log("")

    batch_log.close()

    if n_ok < len(results):
        sys.exit(1)


if __name__ == "__main__":
    main()
