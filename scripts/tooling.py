from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path
from typing import Iterable, Sequence


def resolve_executable(
    scripts_dir: str | Path,
    candidates: Sequence[str],
    *,
    label: str,
) -> Path:
    scripts = Path(scripts_dir)
    checked: list[str] = []
    for name in candidates:
        candidate = Path(name)
        if not candidate.is_absolute():
            candidate = scripts / name
        checked.append(str(candidate))
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return candidate.resolve()
    for name in candidates:
        found = shutil.which(Path(name).name)
        if found:
            return Path(found).resolve()
    raise FileNotFoundError(
        f"{label} executable was not found; checked {', '.join(checked)} and PATH"
    )


def run(command: Iterable[str], *, label: str = "command") -> None:
    command_list = [str(value) for value in command]
    print(f"[PATs] {label}: {' '.join(command_list)}", flush=True)
    subprocess.run(command_list, check=True)


def convert_fasta_to_kmers(
    fasta: str | Path,
    output: str | Path,
    *,
    scripts_dir: str | Path,
    kmer_size: int,
) -> Path:
    executable = resolve_executable(
        scripts_dir,
        ("kmer_convertor8", "kmer_converter", "kmer_convertor", "kmercounter8"),
        label="FASTA-to-kmer converter",
    )
    output_path = Path(output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    run(
        [
            executable,
            "-i",
            Path(fasta),
            "-o",
            output_path,
            "-k",
            str(kmer_size),
            "-s",
            "0",
            "-m",
            "1",
        ],
        label="convert unmasked FASTA sequence to canonical k-mers",
    )
    if not output_path.is_file() or output_path.stat().st_size == 0:
        raise RuntimeError(f"k-mer converter did not produce a non-empty file: {output_path}")
    return output_path
