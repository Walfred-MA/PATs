#!/usr/bin/env python3
"""Convert primary CHM13 GFF3 seqids from chr* names to RefSeq accessions."""

from __future__ import annotations

import argparse
from pathlib import Path


CHM13_PRIMARY_CONTIGS = {
    "chr1": "NC_060925.1",
    "chr2": "NC_060926.1",
    "chr3": "NC_060927.1",
    "chr4": "NC_060928.1",
    "chr5": "NC_060929.1",
    "chr6": "NC_060930.1",
    "chr7": "NC_060931.1",
    "chr8": "NC_060932.1",
    "chr9": "NC_060933.1",
    "chr10": "NC_060934.1",
    "chr11": "NC_060935.1",
    "chr12": "NC_060936.1",
    "chr13": "NC_060937.1",
    "chr14": "NC_060938.1",
    "chr15": "NC_060939.1",
    "chr16": "NC_060940.1",
    "chr17": "NC_060941.1",
    "chr18": "NC_060942.1",
    "chr19": "NC_060943.1",
    "chr20": "NC_060944.1",
    "chr21": "NC_060945.1",
    "chr22": "NC_060946.1",
    "chrX": "NC_060947.1",
    "chrY": "NC_060948.1",
}


def convert_gff3(source: Path, destination: Path) -> tuple[int, int]:
    if source.resolve() == destination.resolve():
        raise ValueError("input and output must differ; refusing to overwrite the GFF3 input")
    if not source.is_file():
        raise FileNotFoundError(f"input GFF3 does not exist: {source}")

    destination.parent.mkdir(parents=True, exist_ok=True)
    features = 0
    converted = 0
    with source.open() as reader, destination.open("w") as writer:
        for line_number, raw in enumerate(reader, 1):
            if raw.startswith("#") or not raw.strip():
                writer.write(raw)
                continue
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                raise ValueError(
                    f"{source}:{line_number}: expected nine tab-separated GFF3 columns"
                )
            features += 1
            replacement = CHM13_PRIMARY_CONTIGS.get(fields[0])
            if replacement is not None:
                fields[0] = replacement
                converted += 1
            writer.write("\t".join(fields) + "\n")
    return features, converted


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Convert chr1..chr22/chrX/chrY seqids in a CHM13 CAT/LiftOff GFF3 "
            "to the NC_060925.1..NC_060948.1 names used by the RefSeq CHM13v2 FASTA."
        )
    )
    parser.add_argument("-i", "--input", required=True, type=Path, help="input GFF3")
    parser.add_argument("-o", "--output", required=True, type=Path, help="new output GFF3")
    args = parser.parse_args()

    try:
        features, converted = convert_gff3(args.input.expanduser(), args.output.expanduser())
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    print(f"[PATs] wrote {args.output}: converted {converted} of {features} feature seqids")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
