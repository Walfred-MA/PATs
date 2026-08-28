#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import shutil
from pathlib import Path

from tooling import resolve_executable, run


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"target-group manifest is empty: {path}")
    required = {"group", "fasta", "kmers", "initial_bed"}
    if not required.issubset(rows[0]):
        raise ValueError(f"target-group manifest lacks columns: {', '.join(sorted(required))}")
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description="Search target or pooled exon k-mers in all assemblies.")
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--queries", required=True, help="two-column C++ search input table")
    parser.add_argument("--exon-kmers", default="")
    parser.add_argument("--output", required=True)
    parser.add_argument("--scripts-dir", required=True)
    parser.add_argument("--kmer-size", type=int, default=31)
    parser.add_argument("--hit-kmers", type=int, default=1_000)
    parser.add_argument("--exon-hit-kmers", type=int, default=300)
    parser.add_argument("--window", type=int, default=2_000)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--sample-threads", type=int, default=1)
    args = parser.parse_args()

    rows = read_manifest(Path(args.manifest))
    output_dir = Path(args.output)
    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True)
    exon_kmers = Path(args.exon_kmers) if args.exon_kmers else None
    use_exons = bool(exon_kmers and exon_kmers.is_file() and exon_kmers.stat().st_size > 0)

    target_list = output_dir / "targets.list"
    output_list = output_dir / "outputs.list"
    with target_list.open("w") as targets, output_list.open("w") as outputs:
        for row in rows:
            group_output = output_dir / row["group"]
            group_output.mkdir()
            targets.write(f"{exon_kmers if use_exons else row['kmers']}\n")
            outputs.write(str(group_output) + "/\n")

    executable = resolve_executable(
        args.scripts_dir,
        ("kmer_searcher", "KmerSearcher", "kmersearcher"),
        label="KmerSearcher",
    )
    cutoff = args.exon_hit_kmers if use_exons else args.hit_kmers
    run(
        [
            executable,
            "-T",
            target_list,
            "-O",
            output_list,
            "-I",
            args.queries,
            "-k",
            str(args.kmer_size),
            "-n",
            str(args.threads),
            "-N",
            str(args.sample_threads),
            "-c",
            str(cutoff),
            "-w",
            str(args.window),
            "--cache",
            "0",
        ],
        label=(
            f"KmerSearcher ({cutoff} pooled exon k-mers/{args.window} bp)"
            if use_exons
            else f"KmerSearcher ({cutoff} target k-mers/{args.window} bp)"
        ),
    )
    (output_dir / ".done").write_text(
        f"mode={'exon' if use_exons else 'target'}\tcutoff={cutoff}\twindow={args.window}\n"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

