#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Sequence

from patslib.core import (
    QueryRecord,
    parse_query_table,
    paths_refer_to_same_file,
    sanitize_name,
    write_query_tables,
    write_text_if_changed,
)


ROOT = Path(__file__).resolve().parent
DEFAULT_SNAKEFILE = ROOT / "Snakefile"


def existing_file(value: str) -> str:
    path = Path(value).expanduser().resolve()
    if not path.is_file():
        raise argparse.ArgumentTypeError(f"file does not exist: {path}")
    return str(path)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="pats.py",
        description="Build a PATs allele FASTA and indexed matrix from BED, gene, or FASTA targets.",
    )
    targets = parser.add_mutually_exclusive_group(required=True)
    targets.add_argument("-b", "--bed", type=existing_file, help="0-based half-open target BED")
    targets.add_argument(
        "-g", "--gene", action="append", help="gene-name prefix, comma list, or a file of prefixes"
    )
    targets.add_argument("-f", "--fasta", type=existing_file, help="target multi-FASTA")

    parser.add_argument("--gff3", type=existing_file, help="GFF3 required with --gene")
    parser.add_argument(
        "--exon",
        type=existing_file,
        help="pooled exon BED for BED mode or pooled exon FASTA for FASTA mode; gene mode derives exons",
    )
    parser.add_argument("-r", "--reference", required=True, type=existing_file, help="reference FASTA")
    parser.add_argument("-q", "--queries", required=True, type=existing_file, help="query TSV")
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="output prefix; writes PREFIX.fa, PREFIX.matrix.txt, PREFIX.matrix.txt.index, and PREFIX.work/",
    )

    parser.add_argument("--gene-extension", "--geneExtension", type=int, default=5_000)
    parser.add_argument("--gene-merge-distance", type=int, default=10_000)
    parser.add_argument("--anchor-size", type=int, default=5_000)
    parser.add_argument("--edge-distance", type=int, default=20_000)
    parser.add_argument("--kmer-size", type=int, default=31)
    parser.add_argument("--partition-min-kmers", type=int, default=500)
    parser.add_argument("--partition-similarity", type=float, default=0.20)
    parser.add_argument("--hit-kmers", type=int, default=1_000)
    parser.add_argument("--exon-hit-kmers", type=int, default=300)
    parser.add_argument("--hit-window", type=int, default=2_000)
    parser.add_argument("--reference-prefix", default="", help="override inferred reference-contig prefix")
    parser.add_argument("--threads", "--cores", type=int, default=max(1, os.cpu_count() or 1))
    parser.add_argument(
        "--sample-threads",
        type=int,
        default=1,
        help="threads used within one indexed FASTA by KmerSearcher",
    )
    parser.add_argument(
        "--scripts-dir",
        default=str(ROOT / "scripts"),
        help="directory containing compiled PATs binaries (also searched via PATH)",
    )
    parser.add_argument("--snakefile", default=str(DEFAULT_SNAKEFILE))
    parser.add_argument("--snakemake", default="snakemake")
    parser.add_argument("--profile", default="")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument(
        "--prepare-only",
        action="store_true",
        help="validate inputs and write normalized config without launching Snakemake",
    )
    parser.add_argument(
        "--snakemake-arg",
        action="append",
        default=[],
        help="additional single argument passed to Snakemake (repeatable)",
    )
    return parser


def flatten_gene_arguments(values: Sequence[str] | None) -> list[str]:
    prefixes: list[str] = []
    for value in values or []:
        candidate = Path(value).expanduser()
        tokens: list[str]
        if candidate.is_file():
            tokens = [
                token
                for raw in candidate.read_text().splitlines()
                if raw.strip() and not raw.lstrip().startswith("#")
                for token in raw.replace(",", " ").split()
            ]
        else:
            tokens = [token.strip() for token in value.split(",") if token.strip()]
        prefixes.extend(tokens)
    deduplicated = list(dict.fromkeys(prefixes))
    if not deduplicated:
        raise ValueError("--gene did not provide any non-empty prefixes")
    return deduplicated


def validate_numeric_args(args: argparse.Namespace) -> None:
    positive = {
        "anchor size": args.anchor_size,
        "k-mer size": args.kmer_size,
        "partition minimum": args.partition_min_kmers,
        "hit k-mers": args.hit_kmers,
        "exon hit k-mers": args.exon_hit_kmers,
        "hit window": args.hit_window,
        "threads": args.threads,
        "sample threads": args.sample_threads,
    }
    for label, value in positive.items():
        if value < 1:
            raise ValueError(f"{label} must be greater than zero")
    nonnegative = {
        "gene extension": args.gene_extension,
        "gene merge distance": args.gene_merge_distance,
        "edge distance": args.edge_distance,
    }
    for label, value in nonnegative.items():
        if value < 0:
            raise ValueError(f"{label} must be zero or greater")
    if not 0.0 <= args.partition_similarity <= 1.0:
        raise ValueError("--partition-similarity must be between 0 and 1")


def infer_mode(args: argparse.Namespace) -> tuple[str, str, list[str]]:
    if args.bed:
        if args.gff3:
            raise ValueError("--gff3 is only valid with --gene")
        return "bed", args.bed, []
    if args.fasta:
        if args.gff3:
            raise ValueError("--gff3 is only valid with --gene")
        return "fasta", args.fasta, []
    if not args.gff3:
        raise ValueError("--gene requires --gff3")
    if args.exon:
        raise ValueError("gene mode derives exons from --gff3; do not also pass --exon")
    return "gene", "", flatten_gene_arguments(args.gene)


def normalize_queries(args: argparse.Namespace, work_dir: Path) -> tuple[list[QueryRecord], Path, Path]:
    # The user-supplied query table is read-only.  Downstream tools consume
    # generated manifests in the output work directory instead.
    records = parse_query_table(args.queries)
    reference = Path(args.reference).resolve()
    first = records[0]
    if not paths_refer_to_same_file(reference, first.fasta):
        raise ValueError(
            "the reference FASTA must be the first query-table row: "
            f"--reference={reference}, first-row FASTA={first.fasta}"
        )
    if not first.whitelist:
        print(
            f"[PATs] reference row {first.sample!r} is forced to whitelist=YES",
            file=sys.stderr,
        )
        records[0] = QueryRecord(first.sample, first.fasta, True)

    full = work_dir / "inputs" / "queries.tsv"
    search = work_dir / "inputs" / "queries.search.tsv"
    for generated in (full, search):
        if paths_refer_to_same_file(args.queries, generated):
            raise ValueError(
                "refusing to overwrite the input query table with a generated manifest: "
                f"{generated}"
            )
    write_query_tables(records, full, search)
    return records, full, search


def build_config(args: argparse.Namespace) -> tuple[dict, Path]:
    validate_numeric_args(args)
    mode, target_input, gene_prefixes = infer_mode(args)
    output_prefix = Path(args.output).expanduser().resolve()
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    work_dir = Path(f"{output_prefix}.work")
    work_dir.mkdir(parents=True, exist_ok=True)
    records, query_table, search_table = normalize_queries(args, work_dir)

    fasta_output = Path(f"{output_prefix}.fa")
    protected_inputs = [Path(args.reference), Path(args.queries)]
    protected_inputs.extend(record.fasta for record in records)
    for value in (args.bed, args.fasta, args.gff3, args.exon):
        if value:
            protected_inputs.append(Path(value))
    for protected in protected_inputs:
        if paths_refer_to_same_file(fasta_output, protected):
            raise ValueError(
                f"refusing to overwrite input file {protected} with final FASTA {fasta_output}; "
                "choose a different --output prefix"
            )

    name_prefix = sanitize_name(output_prefix.name, fallback="PATs")
    config = {
        "version": 1,
        "mode": mode,
        "target_input": target_input,
        "gene_prefixes": gene_prefixes,
        "gff3": args.gff3 or "",
        "exon": args.exon or "",
        "reference": str(Path(args.reference).resolve()),
        "reference_sample": records[0].sample,
        "query_table": str(query_table),
        "search_table": str(search_table),
        "output_prefix": str(output_prefix),
        "name_prefix": name_prefix,
        "work_dir": str(work_dir),
        "tools_dir": str(Path(args.scripts_dir).expanduser().resolve()),
        "gene_extension": args.gene_extension,
        "gene_merge_distance": args.gene_merge_distance,
        "anchor_size": args.anchor_size,
        "edge_distance": args.edge_distance,
        "kmer_size": args.kmer_size,
        "partition_min_kmers": args.partition_min_kmers,
        "partition_similarity": args.partition_similarity,
        "hit_kmers": args.hit_kmers,
        "exon_hit_kmers": args.exon_hit_kmers,
        "hit_window": args.hit_window,
        "reference_prefix": args.reference_prefix,
        "threads": args.threads,
        "sample_threads": args.sample_threads,
        "fasta": str(fasta_output),
        "matrix": f"{output_prefix}.matrix.txt",
        "matrix_index": f"{output_prefix}.matrix.txt.index",
    }
    config_path = work_dir / "config.json"
    config_text = json.dumps(config, indent=2, sort_keys=True) + "\n"
    write_text_if_changed(config_path, config_text)
    return config, config_path


def snakemake_command(args: argparse.Namespace, config_path: Path) -> list[str]:
    snakefile = Path(args.snakefile).expanduser().resolve()
    if not snakefile.is_file():
        raise FileNotFoundError(f"Snakefile does not exist: {snakefile}")
    executable = shutil.which(args.snakemake) or args.snakemake
    command = [
        executable,
        "--snakefile",
        str(snakefile),
        "--configfile",
        str(config_path),
        "--cores",
        str(args.threads),
        "--rerun-incomplete",
        "--printshellcmds",
    ]
    if args.profile:
        command.extend(["--profile", args.profile])
    if args.dry_run:
        command.append("--dry-run")
    command.extend(args.snakemake_arg)
    return command


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        config, config_path = build_config(args)
        command = snakemake_command(args, config_path)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))

    print(f"[PATs] normalized configuration: {config_path}")
    print(f"[PATs] output FASTA: {config['fasta']}")
    print(f"[PATs] output matrix: {config['matrix']}")
    print(f"[PATs] command: {shlex.join(command)}")
    if args.prepare_only:
        return 0
    try:
        return subprocess.run(command).returncode
    except FileNotFoundError:
        print(
            f"[PATs] ERROR: Snakemake executable not found: {args.snakemake!r}",
            file=sys.stderr,
        )
        return 127


if __name__ == "__main__":
    raise SystemExit(main())
