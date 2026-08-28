#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from tooling import resolve_executable, run


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
    output.parent.mkdir(parents=True, exist_ok=True)
    run(
        [
            executable,
            "-I",
            args.queries,
            "-t",
            args.input,
            "-o",
            output,
            "-n",
            str(args.threads),
        ],
        label="select exclusive target k-mers across all assemblies",
    )
    if not output.is_file() or output.stat().st_size == 0:
        raise RuntimeError(f"kmer_selector produced no exclusive-kmer file: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

