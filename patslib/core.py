from __future__ import annotations

import os
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, Mapping, Sequence


TRUE_VALUES = {"yes", "y", "true", "1"}
FALSE_VALUES = {"no", "n", "false", "0", ""}


@dataclass(frozen=True)
class QueryRecord:
    sample: str
    fasta: Path
    whitelist: bool = False


@dataclass(frozen=True)
class FastaRecord:
    name: str
    sequence: str
    description: str = ""

    @property
    def header(self) -> str:
        return self.name if not self.description else f"{self.name}\t{self.description}"


def _clean_line(raw: str) -> str:
    return raw.strip()


def normalize_yes_no(value: str, *, line_number: int | None = None) -> bool:
    normalized = value.strip().lower()
    if normalized in TRUE_VALUES:
        return True
    if normalized in FALSE_VALUES:
        return False
    where = f" on line {line_number}" if line_number is not None else ""
    raise ValueError(
        f"invalid whitelist value{where}: {value!r}; expected YES/NO (legacy 1/0 is accepted)"
    )


def parse_query_table(path: str | Path, *, require_files: bool = True) -> list[QueryRecord]:
    table_path = Path(path).expanduser().resolve()
    records: list[QueryRecord] = []
    seen: set[str] = set()

    with table_path.open() as handle:
        for line_number, raw in enumerate(handle, 1):
            stripped = _clean_line(raw)
            if not stripped or stripped.startswith("#"):
                continue
            fields = stripped.split()
            if len(fields) < 2:
                raise ValueError(
                    f"{table_path}:{line_number}: expected SAMPLE FASTA [YES|NO]"
                )
            if len(fields) > 3:
                raise ValueError(
                    f"{table_path}:{line_number}: expected at most three columns"
                )
            sample = fields[0]
            if "/" in sample or "\\" in sample:
                raise ValueError(
                    f"{table_path}:{line_number}: sample names cannot contain path separators"
                )
            if sample in seen:
                raise ValueError(f"{table_path}:{line_number}: duplicate sample {sample!r}")
            seen.add(sample)

            fasta = Path(fields[1]).expanduser()
            if not fasta.is_absolute():
                fasta = (table_path.parent / fasta).resolve()
            else:
                fasta = fasta.resolve()
            if require_files and not fasta.is_file():
                raise FileNotFoundError(f"query FASTA does not exist: {fasta}")

            whitelist = normalize_yes_no(fields[2], line_number=line_number) if len(fields) >= 3 else False
            records.append(QueryRecord(sample=sample, fasta=fasta, whitelist=whitelist))

    if not records:
        raise ValueError(f"query table has no data rows: {table_path}")
    return records


def write_query_tables(
    records: Sequence[QueryRecord],
    full_path: str | Path,
    search_path: str | Path,
) -> None:
    full = Path(full_path)
    search = Path(search_path)
    full.parent.mkdir(parents=True, exist_ok=True)
    search.parent.mkdir(parents=True, exist_ok=True)
    full_text = "".join(
        f"{record.sample}\t{record.fasta}\t{'YES' if record.whitelist else 'NO'}\n"
        for record in records
    )
    search_text = "".join(
        f"{record.sample}\t{record.fasta}\n" for record in records
    )
    write_text_if_changed(full, full_text)
    write_text_if_changed(search, search_text)


def write_text_if_changed(path: str | Path, text: str) -> bool:
    """Write text only when its content changed, preserving useful mtimes."""

    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.is_file() and output.read_text() == text:
        return False
    output.write_text(text)
    return True


def read_fasta(path: str | Path) -> Iterator[FastaRecord]:
    fasta_path = Path(path)
    name: str | None = None
    description = ""
    sequence: list[str] = []

    with fasta_path.open() as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.rstrip("\r\n")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield FastaRecord(name, "".join(sequence), description)
                header = line[1:].strip()
                if not header:
                    raise ValueError(f"{fasta_path}:{line_number}: empty FASTA header")
                parts = header.split(maxsplit=1)
                name = parts[0]
                description = parts[1] if len(parts) == 2 else ""
                sequence = []
            else:
                if name is None:
                    raise ValueError(f"{fasta_path}:{line_number}: sequence before first header")
                sequence.append("".join(line.split()))
    if name is not None:
        yield FastaRecord(name, "".join(sequence), description)


def write_fasta(
    records: Iterable[FastaRecord],
    path: str | Path,
    *,
    width: int = 80,
) -> int:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with output.open("w") as handle:
        for record in records:
            count += 1
            handle.write(f">{record.header}\n")
            sequence = record.sequence
            if width <= 0:
                handle.write(sequence + "\n")
            else:
                for start in range(0, len(sequence), width):
                    handle.write(sequence[start : start + width] + "\n")
    return count


def fasta_record_count(path: str | Path) -> int:
    return sum(1 for _ in read_fasta(path))


def read_fai(path: str | Path) -> dict[str, int]:
    sizes: dict[str, int] = {}
    with Path(path).open() as handle:
        for line_number, raw in enumerate(handle, 1):
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise ValueError(f"{path}:{line_number}: malformed FASTA index row")
            sizes[fields[0]] = int(fields[1])
    return sizes


def fasta_index_path(fasta: str | Path) -> Path:
    return Path(f"{Path(fasta)}.fai")


def ensure_fasta_index(fasta: str | Path, *, samtools: str = "samtools") -> Path:
    fasta_path = Path(fasta)
    index = fasta_index_path(fasta_path)
    if (
        index.is_file()
        and index.stat().st_size > 0
        and index.stat().st_mtime_ns >= fasta_path.stat().st_mtime_ns
    ):
        return index
    subprocess.run([samtools, "faidx", str(fasta_path)], check=True)
    if not index.is_file():
        raise RuntimeError(f"samtools faidx did not create {index}")
    return index


def extract_fasta_region(
    fasta: str | Path,
    contig: str,
    start: int,
    end: int,
    *,
    samtools: str = "samtools",
) -> str:
    if start < 0 or end <= start:
        raise ValueError(f"invalid zero-based half-open interval: {contig}:{start}-{end}")
    region = f"{contig}:{start + 1}-{end}"
    try:
        result = subprocess.run(
            [samtools, "faidx", str(fasta), region],
            check=True,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    except subprocess.CalledProcessError as exc:
        detail = (exc.stderr or "").strip()
        message = f"samtools faidx failed for reference region {region!r} in {fasta}"
        if detail:
            message += f": {detail}"
        raise RuntimeError(message) from exc
    return "".join(line.strip() for line in result.stdout.splitlines() if not line.startswith(">"))


def merge_intervals(
    intervals: Iterable[tuple[int, int]],
    *,
    max_gap: int = 0,
) -> list[tuple[int, int]]:
    ordered = sorted((int(start), int(end)) for start, end in intervals if int(end) > int(start))
    if not ordered:
        return []
    merged: list[list[int]] = [[ordered[0][0], ordered[0][1]]]
    for start, end in ordered[1:]:
        previous = merged[-1]
        if start - previous[1] <= max_gap:
            previous[1] = max(previous[1], end)
        else:
            merged.append([start, end])
    return [(start, end) for start, end in merged]


def sanitize_name(text: str, *, fallback: str = "target") -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", text.strip()).strip("_.-")
    return cleaned or fallback


def paths_refer_to_same_file(left: str | Path, right: str | Path) -> bool:
    left_path = Path(left)
    right_path = Path(right)
    try:
        return os.path.samefile(left_path, right_path)
    except OSError:
        return left_path.resolve() == right_path.resolve()


def find_sample_for_contig(
    contig: str,
    contig_to_samples: Mapping[str, Sequence[str]],
    *,
    preferred: str | None = None,
) -> str:
    candidates = list(contig_to_samples.get(contig, ()))
    if preferred and preferred in candidates:
        return preferred
    if len(candidates) == 1:
        return candidates[0]
    if not candidates:
        raise KeyError(f"contig {contig!r} is absent from all query FASTA indexes")
    raise ValueError(
        f"contig {contig!r} occurs in multiple query FASTAs ({', '.join(candidates)}); "
        "the FASTA header must encode the sample or the caller must provide it explicitly"
    )
