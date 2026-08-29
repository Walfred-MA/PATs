#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from tooling import resolve_executable, run


def cleanup_selector_scratch(output: Path, target_fasta: Path) -> None:
    """Remove kmer_selector sidecars without touching its declared output."""

    if output.parent.is_dir():
        for candidate in output.parent.iterdir():
            if candidate == output or not candidate.name.startswith(output.name):
                continue
            if candidate.is_file() or candidate.is_symlink():
                candidate.unlink(missing_ok=True)

    Path(f"{target_fasta}_allkmer.cache").unlink(missing_ok=True)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Select k-mers exclusive to one graph/locus FASTA among all query assemblies."
    )
    parser.add_argument("-i", "--input", required=True, help="target locus FASTA")
    parser.add_argument("-q", "--queries", required=True, help="two-column C++ query table")
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--scripts-dir", required=True)
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args()

    executable = resolve_executable(
        args.scripts_dir,
        ("kmer_selector", "kmerselector"),
        label="kmer_selector",
    )
    output = Path(args.output)
    target_fasta = Path(args.input)
    output.parent.mkdir(parents=True, exist_ok=True)
    succeeded = False
    try:
        run(
            [
                executable,
                "-I",
                args.queries,
                "-t",
                target_fasta,
                "-o",
                output,
                "-n",
                str(args.threads),
            ],
            label="select exclusive target k-mers across all assemblies",
        )
        if not output.is_file() or output.stat().st_size == 0:
            raise RuntimeError(f"kmer_selector produced no exclusive-kmer file: {output}")
        succeeded = True
    finally:
        cleanup_selector_scratch(output, target_fasta)
        if not succeeded:
            output.unlink(missing_ok=True)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
