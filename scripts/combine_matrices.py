#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description="Concatenate PATs matrix shards in group order.")
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("inputs", nargs="+")
    args = parser.parse_args()
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("wb") as destination:
        for value in args.inputs:
            source = Path(value)
            if not source.is_file() or source.stat().st_size == 0:
                raise FileNotFoundError(f"matrix shard is absent or empty: {source}")
            with source.open("rb") as handle:
                shutil_copyfileobj(handle, destination)
    if output.stat().st_size == 0:
        raise RuntimeError("combined matrix is empty")
    return 0


def shutil_copyfileobj(source, destination, length: int = 1024 * 1024) -> None:
    while True:
        chunk = source.read(length)
        if not chunk:
            return
        destination.write(chunk)


if __name__ == "__main__":
    raise SystemExit(main())

