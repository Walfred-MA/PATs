#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
from pathlib import Path
import tempfile


def copy_fasta(source: Path, destination) -> None:
    if not source.is_file() or source.stat().st_size == 0:
        raise FileNotFoundError(f"FASTA shard is absent or empty: {source}")

    with source.open("rb") as handle:
        first = handle.read(1)
        if first != b">":
            raise ValueError(f"FASTA shard does not start with a header: {source}")
        destination.write(first)

        last = first
        while True:
            chunk = handle.read(1024 * 1024)
            if not chunk:
                break
            destination.write(chunk)
            last = chunk[-1:]

        if last != b"\n":
            destination.write(b"\n")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Concatenate final PATs allele FASTA shards in target-group order."
    )
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("inputs", nargs="+")
    args = parser.parse_args()

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{output.name}.", suffix=".tmp", dir=output.parent
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as destination:
            for value in args.inputs:
                copy_fasta(Path(value), destination)

        if temporary.stat().st_size == 0:
            raise RuntimeError("combined FASTA is empty")
        os.replace(temporary, output)
    finally:
        temporary.unlink(missing_ok=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
