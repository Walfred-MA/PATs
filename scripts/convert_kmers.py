#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from tooling import convert_fasta_to_kmers


def fasta_has_records(path: str) -> bool:
    with Path(path).open() as handle:
        return any(line.startswith(">") for line in handle)


def main() -> int:
    parser = argparse.ArgumentParser(description="Convert unmasked FASTA sequence to PATs k-mers.")
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("-k", "--kmer-size", type=int, default=31)
    parser.add_argument("--scripts-dir", required=True)
    parser.add_argument("--allow-empty", action="store_true")
    args = parser.parse_args()
    if args.allow_empty and not fasta_has_records(args.input):
        output = Path(args.output)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text("")
        return 0
    convert_fasta_to_kmers(
        args.input,
        args.output,
        scripts_dir=args.scripts_dir,
        kmer_size=args.kmer_size,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
