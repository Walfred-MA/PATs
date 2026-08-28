#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from patslib.core import FastaRecord, read_fasta, sanitize_name, write_fasta  # noqa: E402
from tooling import convert_fasta_to_kmers, resolve_executable, run  # noqa: E402


REGION_RE = re.compile(r"^([^\s:]+):(-?\d+)-(\d+)(?:[+-])?$")


@dataclass
class TargetGroup:
    raw_name: str
    records: list[FastaRecord]


def numeric_partition_key(path: Path) -> tuple[int, str]:
    match = re.search(r"p(\d+)\.fa$", path.name)
    return (int(match.group(1)) if match else 10**9, path.name)


def parse_region(description: str) -> tuple[str, int, int] | None:
    if not description:
        return None
    token = description.split()[0]
    match = REGION_RE.match(token)
    if not match:
        return None
    return match.group(1), int(match.group(2)), int(match.group(3))


def invoke_partitioner(args: argparse.Namespace, raw_dir: Path) -> list[Path]:
    executable = resolve_executable(
        args.scripts_dir,
        ("kmerpartition", "kmer_partition"),
        label="kmerpartition",
    )
    raw_dir.mkdir(parents=True, exist_ok=True)
    output_prefix = str(raw_dir) + "/"
    run(
        [
            executable,
            "-i",
            args.input,
            "-k",
            args.kmers,
            "-o",
            output_prefix,
            "-n",
            str(args.threads),
            "-c",
            str(args.min_kmers),
            "-s",
            str(args.similarity),
            "-m",
            "1",
        ],
        label="partition target records by k-mer similarity",
    )
    return sorted(raw_dir.glob("p*.fa"), key=numeric_partition_key)


def build_groups(args: argparse.Namespace, original: list[FastaRecord], raw_dir: Path) -> list[TargetGroup]:
    by_name = {record.name: record for record in original}
    assigned: set[str] = set()
    groups: list[TargetGroup] = []

    if len(original) > 1:
        raw_fastas = invoke_partitioner(args, raw_dir)
        for raw_fasta in raw_fastas:
            records: list[FastaRecord] = []
            for raw_record in read_fasta(raw_fasta):
                record = by_name.get(raw_record.name, raw_record)
                if record.name in assigned:
                    continue
                records.append(record)
                assigned.add(record.name)
            if records:
                groups.append(TargetGroup(raw_fasta.stem, records))

    for record in original:
        if record.name not in assigned:
            groups.append(TargetGroup(f"singleton_{record.name}", [record]))
            assigned.add(record.name)

    if not groups:
        raise RuntimeError("target partitioning produced no groups")
    return groups


def main() -> int:
    parser = argparse.ArgumentParser(description="Partition and deterministically name target groups.")
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-k", "--kmers", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--mode", choices=["bed", "gene", "fasta"], required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--scripts-dir", required=True)
    parser.add_argument("--kmer-size", type=int, default=31)
    parser.add_argument("--min-kmers", type=int, default=500)
    parser.add_argument("--similarity", type=float, default=0.20)
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args()

    original = list(read_fasta(args.input))
    if not original:
        raise ValueError(f"target FASTA has no records: {args.input}")
    if len({record.name for record in original}) != len(original):
        raise ValueError("target FASTA record names must be unique")

    output_dir = Path(args.output)
    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True)
    groups = build_groups(args, original, output_dir / "raw")
    name_prefix = sanitize_name(args.prefix, fallback="PATs")

    manifest = Path(args.manifest)
    manifest.parent.mkdir(parents=True, exist_ok=True)
    with manifest.open("w") as manifest_handle:
        manifest_handle.write("group\tfasta\tkmers\tinitial_bed\n")
        for group_number, group in enumerate(groups, 1):
            group_id = f"group{group_number}"
            group_dir = output_dir / group_id
            group_dir.mkdir(parents=True)
            group_fasta = group_dir / "targets.fa"
            group_kmers = group_dir / "targets.kmers.txt"
            initial_bed = group_dir / "initial.bed"

            renamed = [
                FastaRecord(
                    f"{name_prefix}group{group_number}_{record.name}",
                    record.sequence,
                    record.description,
                )
                for record in group.records
            ]
            write_fasta(renamed, group_fasta)
            convert_fasta_to_kmers(
                group_fasta,
                group_kmers,
                scripts_dir=args.scripts_dir,
                kmer_size=args.kmer_size,
            )
            with initial_bed.open("w") as bed_handle:
                if args.mode in {"bed", "gene"}:
                    for record in renamed:
                        region = parse_region(record.description)
                        if region is None:
                            raise ValueError(
                                f"target {record.name!r} lacks a genomic region in its FASTA description"
                            )
                        contig, start, end = region
                        bed_handle.write(f"{contig}\t{start}\t{end}\t{record.name}\n")

            manifest_handle.write(
                "\t".join(
                    [
                        group_id,
                        str(group_fasta.resolve()),
                        str(group_kmers.resolve()),
                        str(initial_bed.resolve()),
                    ]
                )
                + "\n"
            )

    print(f"[PATs] target partitioning produced {len(groups)} group(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

