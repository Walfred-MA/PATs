#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

from tooling import resolve_executable, run


def run_python(script: Path, arguments: list[object], *, stdout: Path | None = None) -> None:
    command = [sys.executable, str(script), *(str(value) for value in arguments)]
    print(f"[PATs] packedrun: {' '.join(command)}", flush=True)
    if stdout is None:
        subprocess.run(command, check=True)
    else:
        stdout.parent.mkdir(parents=True, exist_ok=True)
        with stdout.open("w") as handle:
            subprocess.run(command, check=True, stdout=handle)


def link_graph(source: Path, destination: Path) -> None:
    if destination.is_symlink() or destination.exists():
        destination.unlink()
    try:
        destination.symlink_to(source.resolve())
    except OSError:
        shutil.copyfile(source, destination)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Normalize a final locus FASTA and compile one PATs matrix shard."
    )
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-k", "--kmers", required=True)
    parser.add_argument("-r", "--reference", required=True, help="target-group reference FASTA")
    parser.add_argument("-g", "--graph", default="", help="gfixbreaks graph FASTA base")
    parser.add_argument("-a", "--alignment", default="", help="gfixbreaks full graph alignment")
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--scripts-dir", required=True, help="redesigned PATs Python scripts")
    parser.add_argument("--tools-dir", required=True, help="compiled PATs binaries")
    parser.add_argument("--reference-sample", required=True)
    parser.add_argument("--kmer-cutoff", type=int, default=1_000)
    parser.add_argument("--kmer-ratio", type=float, default=0.10)
    args = parser.parse_args()

    scripts_dir = Path(args.scripts_dir)
    tools_dir = Path(args.tools_dir)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    prefix = output.parent / output.stem
    filtered_fasta = Path(f"{prefix}.filtered.fa")
    filtered_kmers = Path(f"{filtered_fasta}.kmer.list")
    norm = Path(f"{filtered_fasta}.norm.gz")
    tree = Path(f"{prefix}.tree.fa")
    annotated_kmers = Path(f"{tree}.kmer.list")

    query_filter = scripts_dir / "querylist_filterbykmer.py"
    run_python(
        query_filter,
        [
            "-c",
            args.kmer_cutoff,
            "-r",
            args.kmer_ratio,
            "-i",
            args.input,
            "-k",
            args.kmers,
            "-o",
            filtered_fasta,
            "-l",
            filtered_kmers,
            "-p",
            "1",
        ],
        stdout=Path(f"{filtered_fasta}.info.txt"),
    )
    if not filtered_fasta.is_file() or filtered_fasta.stat().st_size == 0:
        raise RuntimeError("packedrun filtering removed every locus sequence")

    kmernorm = resolve_executable(tools_dir, ("kmernorm",), label="kmernorm")
    run(
        [
            kmernorm,
            "-i",
            filtered_fasta,
            "-k",
            filtered_kmers,
            "-o",
            norm,
            "-w",
            "1",
            "-h",
            "1",
        ],
        label="packedrun k-mer normalization",
    )
    kmertree = resolve_executable(tools_dir, ("kmertree",), label="kmertree")
    run([kmertree, "-i", filtered_fasta, "-n", norm, "-o", tree], label="packedrun tree")

    kmerstrd = resolve_executable(
        tools_dir,
        ("kmerstrd", "KmerStrd"),
        label="kmerstrd",
    )
    oriented_tree = Path(f"{tree}.oriented")
    run(
        [kmerstrd, "-i", tree, "-r", args.reference, "-o", oriented_tree],
        label="packedrun orient tree",
    )
    os.replace(oriented_tree, tree)

    graph = Path(args.graph) if args.graph else None
    alignment = Path(args.alignment) if args.alignment else None
    annotate_args: list[object] = [
        "-k",
        filtered_kmers,
        "-o",
        annotated_kmers,
        "-p",
        args.reference_sample,
    ]
    if graph and alignment and graph.is_file() and alignment.is_file() and alignment.stat().st_size > 0:
        local_alignment = Path(f"{tree}.allgraphalign.out")
        run_python(
            scripts_dir / "getLocalCigar.py",
            ["-i", tree, "-s", graph, "-a", alignment, "-o", local_alignment],
        )
        graph_fasta = Path(f"{graph}_graph.FA")
        if graph_fasta.is_file() and local_alignment.is_file() and local_alignment.stat().st_size > 0:
            tree_graph = Path(f"{tree}_graph.FA")
            link_graph(graph_fasta, tree_graph)
            annotate_args = [
                "-k",
                filtered_kmers,
                "-s",
                tree,
                "-g",
                tree_graph,
                "-a",
                local_alignment,
                "-o",
                annotated_kmers,
                "-p",
                args.reference_sample,
            ]

    run_python(scripts_dir / "kmerannotate.py", annotate_args)
    run_python(
        scripts_dir / "matrixcompile.py",
        ["-s", tree, "-k", annotated_kmers, "-o", output],
    )
    if not output.is_file() or output.stat().st_size == 0:
        raise RuntimeError(f"matrixcompile.py produced no matrix shard: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
