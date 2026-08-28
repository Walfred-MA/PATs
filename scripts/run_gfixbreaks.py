#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path


def infer_prefix(reference: Path) -> str:
    fai = Path(f"{reference}.fai")
    contigs: list[str] = []
    with fai.open() as handle:
        for raw in handle:
            if raw.strip():
                contigs.append(raw.split()[0])
            if len(contigs) >= 24:
                break
    if not contigs:
        raise ValueError(f"reference index contains no contigs: {fai}")
    prefix = os.path.commonprefix(contigs)
    if prefix:
        return prefix
    first = contigs[0]
    for index, character in enumerate(first):
        if character.isdigit():
            return first[:index] or first
    return first


def main() -> int:
    parser = argparse.ArgumentParser(description="Run gfixbreaks.py without graph caching.")
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("-q", "--queries", required=True)
    parser.add_argument("-r", "--reference", required=True)
    parser.add_argument("--target-fasta", required=True)
    parser.add_argument("--reference-prefix", default="")
    parser.add_argument("--group-prefix", default="")
    parser.add_argument("--anchor-size", type=int, default=5_000)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--scripts-dir", required=True, help="redesigned PATs Python scripts")
    parser.add_argument("--tools-dir", required=True, help="compiled PATs binaries")
    args = parser.parse_args()
    prefix = args.reference_prefix or infer_prefix(Path(args.reference))
    command = [
        sys.executable,
        str(Path(args.scripts_dir) / "gfixbreaks.py"),
        "-i",
        args.input,
        "-o",
        args.output,
        "-q",
        args.queries,
        "-t",
        str(args.threads),
        "-r",
        prefix,
        "--target-fasta",
        args.target_fasta,
        "--anchor-size",
        str(args.anchor_size),
        "--tools-dir",
        args.tools_dir,
    ]
    if args.group_prefix:
        command.extend(["--group-prefix", args.group_prefix])
    print(f"[PATs] gfixbreaks reference prefix: {prefix}", flush=True)
    print(f"[PATs] gfixbreaks (cache disabled): {' '.join(command)}", flush=True)
    subprocess.run(command, check=True)
    output = Path(args.output)
    expected = [
        output,
        Path(f"{output}_loci.txt"),
        Path(f"{output}_loci.txt.fasta"),
        Path(f"{output}_loci.txt.fasta_graph.FA"),
        Path(f"{output}_loci.txt.fasta_allgraphalign.out"),
    ]
    missing = [path for path in expected if not path.is_file() or path.stat().st_size == 0]
    if missing:
        raise RuntimeError(
            "gfixbreaks.py did not produce required graph output(s): "
            + ", ".join(str(path) for path in missing)
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
