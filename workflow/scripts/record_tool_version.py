#!/usr/bin/env python3
"""
record_tool_version.py — Idempotently record one tool's version string
into a shared JSON file, for later use in run_summary.pdf's Provenance
section.

Why this exists: with Snakemake's --use-conda, each rule runs in its own
isolated conda environment. By the time write_summary.py runs (in its
own, unrelated plotting environment), it has no way to see whether
blastp/mafft/Gblocks/codeml are even on PATH, let alone which version
ran — those tools were only ever reachable from inside their own rule's
environment, which no longer exists by the time the report is written.

The fix is to capture each tool's version at the point its own rule
actually runs — where it genuinely is on PATH — rather than trying to
rediscover it afterwards from an unrelated environment. This script is
that capture point: call it from within (or right after) the rule that
actually invokes the tool, and write_summary.py just reads the result
back from disk afterward instead of querying PATH itself.

Safe to call redundantly from thousands of parallel per-gene jobs (e.g.
every single mafft invocation calls this): skips re-running the version
command if already recorded, and uses an atomic write (temp file +
rename) so concurrent first-time writers can't corrupt the shared file —
worst case a few of them redundantly query the version before one wins
the write race, never data loss or corruption.

Usage:
  # Run a command and record its (first line of) output:
  record_tool_version.py --file tool_versions.json --tool mafft \
      --cmd "mafft --version"

  # Record a fixed string directly (e.g. a tool with no --version flag,
  # or a version string already extracted from the tool's own output
  # elsewhere, such as PAML's banner in a codeml output file):
  record_tool_version.py --file tool_versions.json --tool codeml \
      --literal "CODONML (in paml version 4.9j, February 2020)"
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
import json
import os
import shutil
import subprocess
from pathlib import Path


def get_version_from_cmd(cmd: list, stdin_input: str = None) -> str:
    """
    Run a tool's version-check command and return its first line of
    output, or a short explanatory string. Never raises.
    """
    if shutil.which(cmd[0]) is None:
        return "not found on PATH"
    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=10,
            input=stdin_input if stdin_input is not None else None,
        )
        output = ((result.stdout or "") + (result.stderr or "")).strip()
        return output.splitlines()[0] if output else "unavailable (no version output)"
    except Exception as exc:
        return f"unavailable ({exc.__class__.__name__})"


def record(path: Path, tool: str, version: str) -> None:
    """
    Add {tool: version} to the shared JSON file at path, preserving
    whatever other tools are already recorded there, via an atomic
    write (temp file + rename) so concurrent first-time writers from
    parallel Snakemake jobs can't corrupt each other's output.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    existing = {}
    if path.exists():
        try:
            existing = json.loads(path.read_text())
        except (json.JSONDecodeError, OSError):
            existing = {}
    if tool in existing:
        return  # already recorded (by us or a concurrent job) — nothing to do
    existing[tool] = version
    tmp = path.with_name(path.name + f".tmp{os.getpid()}")
    try:
        tmp.write_text(json.dumps(existing, indent=2, sort_keys=True))
        os.replace(tmp, path)
    except OSError:
        try:
            tmp.unlink(missing_ok=True)
        except OSError:
            pass


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--file", required=True, help="Shared tool_versions.json path")
    parser.add_argument("--tool", required=True, help="Key to store this version under")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--cmd",
        help="Version-check command as one quoted string, e.g. "
             "--cmd \"mafft --version\" (kept as a single string rather than "
             "nargs='+' because argparse otherwise mistakes a flag like "
             "--version for a new top-level option instead of part of --cmd)"
    )
    group.add_argument("--literal", help="Store this exact string instead of running a command")
    parser.add_argument(
        "--stdin-input", default=None,
        help="Sent as stdin to --cmd (needed for tools that would otherwise "
             "wait on interactive input, e.g. PAML's codeml with no args)"
    )
    args = parser.parse_args()

    path = Path(args.file)

    # Fast path: already recorded (by an earlier job), nothing to do — this
    # is what keeps thousands of redundant per-gene calls (e.g. every mafft
    # invocation) cheap after the very first one.
    if path.exists():
        try:
            existing = json.loads(path.read_text())
            if args.tool in existing:
                return
        except (json.JSONDecodeError, OSError):
            pass  # fall through and try to (re)record

    if args.literal is not None:
        version = args.literal
    else:
        version = get_version_from_cmd(args.cmd.split(), args.stdin_input)

    record(path, args.tool, version)


if __name__ == "__main__":
    main()
