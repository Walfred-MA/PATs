#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from patslib.core import ensure_fasta_index, parse_query_table  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description="Ensure every PATs query FASTA has a samtools index.")
    parser.add_argument("-q", "--queries", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--samtools", default="samtools")
    args = parser.parse_args()
    records = parse_query_table(args.queries)
    for record in records:
        ensure_fasta_index(record.fasta, samtools=args.samtools)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(f"{record.sample}\t{record.fasta}.fai" for record in records) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

