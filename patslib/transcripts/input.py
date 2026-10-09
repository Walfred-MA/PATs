"""Read native ParaTyper calls and extract their target assembly intervals."""
from __future__ import annotations

import csv
import gzip
import math
import re
from dataclasses import dataclass
from pathlib import Path


REQUIRED = {
    "GENE_index", "model_type", "transcript_id", "gene_id", "gene_name",
    "query_contig", "query_start", "query_end", "strand",
    "exon_query_coordinates", "found_exon_numbers", "exon_coverages",
    "call_status", "assignment_status", "tie_group_id", "pipeline_version",
}
ARRAYS = {
    "exon_query_coordinates", "exon_reference_coordinates", "exon_ids",
    "exon_coverages", "exon_identities", "alignment_AS", "found_exon_numbers",
}


@dataclass(frozen=True)
class Sample:
    name: str
    assembly: Path
    calls: Path
    fragments: Path | None


def read_manifest(path: Path) -> list[Sample]:
    """Paths are relative to the manifest. ParaTyper tables are passed unchanged."""
    samples = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"sample", "assembly", "calls"} <= set(reader.fieldnames or []):
            raise ValueError("manifest requires tab-separated sample, assembly, calls columns")
        for row in reader:
            if None in row or any(row.get(key) is None for key in ("sample", "assembly", "calls")):
                raise ValueError(f"{path}:{reader.line_num}: malformed manifest row")
            name = row["sample"].strip()
            if not name or any(s.name == name for s in samples):
                raise ValueError(f"{path}:{reader.line_num}: empty or duplicate sample")
            paths = []
            for key in ("assembly", "calls", "fragments"):
                value = (row.get(key) or "").strip()
                if not value and key == "fragments":
                    paths.append(None)
                    continue
                p = Path(value).expanduser()
                p = (path.parent / p).resolve() if not p.is_absolute() else p.resolve()
                if not p.is_file():
                    raise ValueError(f"{path}:{reader.line_num}: missing {key}: {p}")
                paths.append(p)
            samples.append(Sample(name, *paths))
    if not samples:
        raise ValueError("manifest has no samples")
    return samples


def alternatives(row: dict[str, str]) -> list[dict[str, str]]:
    """Expand aligned alternatives; shared scalar fields may be broadcast."""
    n = len(row["gene_id"].split(";"))
    result = [dict() for _ in range(n)]
    for key, value in row.items():
        parts = value.split(";")
        if len(parts) not in (1, n):
            raise ValueError(f"{key}: alternatives do not match gene_id")
        if key in ARRAYS and n > 1 and value and len(parts) != n:
            raise ValueError(f"{key}: tied exon arrays must retain one entry per alternative")
        for i in range(n):
            result[i][key] = parts[i] if len(parts) == n else value
    return result


def intervals(text: str) -> list[tuple[int, int]]:
    result = []
    for value in text.split(","):
        match = re.fullmatch(r"(\d+)-(\d+)", value)
        if not match:
            raise ValueError(f"invalid exon interval: {value!r}")
        start, end = map(int, match.groups())
        if end <= start:
            raise ValueError(f"empty or reversed exon interval: {value}")
        result.append((start, end))
    if result != sorted(result):
        raise ValueError("exons must be in increasing target coordinate order")
    return result


def validate_call(row: dict[str, str]) -> None:
    if row["model_type"] not in {"transcript", "full_gene"}:
        raise ValueError(f"unknown model_type: {row['model_type']}")
    if not row["gene_id"] or not row["transcript_id"] or not row["query_contig"]:
        raise ValueError("gene_id, transcript_id and query_contig must be nonempty")
    if row["strand"] not in {"+", "-"}:
        raise ValueError("strand must be + or -")
    if row["call_status"] not in {"complete", "partial"}:
        raise ValueError("call_status must be complete or partial")
    if row["assignment_status"] not in {"unique", "tied"}:
        raise ValueError("assignment_status must be unique or tied")
    start, end = int(row["query_start"]), int(row["query_end"])
    exons = intervals(row["exon_query_coordinates"])
    if start < 0 or end <= start or any(s < start or e > end for s, e in exons):
        raise ValueError("exons fall outside the call bounds")
    for key in ARRAYS & row.keys():
        if row[key] and len(row[key].split(",")) != len(exons):
            raise ValueError(f"{key}: length does not match exon coordinates")
    if len(row["found_exon_numbers"].split(",")) != len(exons):
        raise ValueError("missing exon numbers")
    if any(int(x) < 1 for x in row["found_exon_numbers"].split(",")):
        raise ValueError("exon numbers must be positive")
    cov = [float(x) for x in row["exon_coverages"].split(",")]
    if len(cov) != len(exons) or any(not math.isfinite(x) or not 0 <= x <= 100 for x in cov):
        raise ValueError("invalid exon coverages")


def read_calls(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        missing = REQUIRED - set(fields)
        if missing or len(set(fields)) != len(fields):
            raise ValueError(f"{path}: missing/duplicate ParaTyper columns: {sorted(missing)}")
        for raw in reader:
            try:
                if None in raw or any(value is None for value in raw.values()):
                    raise ValueError("malformed TSV row")
                expanded = alternatives(raw)
                for row in expanded:
                    validate_call(row)
                yield reader.line_num, raw, expanded
            except ValueError as exc:
                raise ValueError(f"{path}:{reader.line_num}: {exc}") from exc


class Fasta:
    """Seek a plain FASTA using its FAI, or scan to build an in-memory contig index.

    Assemblies are never loaded into RAM or modified. Regular wrapped FASTAs
    with LF/CRLF and an optional shorter last line are supported.
    """

    def __init__(self, path: Path):
        self.path = path
        self.handle = path.open("rb")
        self.index = {}
        try:
            if self.handle.read(2) == b"\x1f\x8b":
                raise ValueError(f"use an uncompressed assembly FASTA: {path}")
            self.handle.seek(0)
            fai = Path(str(path) + ".fai")
            if fai.is_file() and fai.stat().st_mtime_ns >= path.stat().st_mtime_ns:
                for line in fai.read_text().splitlines():
                    fields = line.split("\t")
                    if len(fields) < 5 or fields[0] in self.index:
                        raise ValueError(f"invalid FASTA index: {fai}")
                    values = tuple(map(int, fields[1:5]))
                    if values[0] < 0 or values[1] < 0 or values[2] < 1 or values[3] < values[2]:
                        raise ValueError(f"invalid FASTA index: {fai}")
                    self.index[fields[0]] = values
            else:
                self._index()
        except Exception:
            self.handle.close()
            raise

    def _index(self):
        name = None
        length = offset = bases = width = 0
        previous = None
        while True:
            position = self.handle.tell()
            raw = self.handle.readline()
            if not raw or raw.startswith(b">"):
                if name is not None:
                    if not length:
                        raise ValueError(f"empty FASTA contig: {name}")
                    self.index[name] = (length, offset, bases, width)
                if not raw:
                    break
                tokens = raw[1:].split()
                if not tokens:
                    raise ValueError("empty FASTA header")
                name = tokens[0].decode("utf-8")
                if name in self.index:
                    raise ValueError(f"duplicate FASTA contig: {name}")
                length = bases = width = 0
                previous = None
                offset = self.handle.tell()
                continue
            sequence = raw.rstrip(b"\r\n")
            if name is None or not sequence or re.search(rb"[^A-Za-z]", sequence):
                raise ValueError(f"invalid FASTA sequence at byte {position}: {self.path}")
            if previous is not None and previous != (bases, width):
                raise ValueError(f"irregular FASTA wrapping in {name}")
            if not bases:
                bases, width = len(sequence), len(raw)
            elif len(sequence) > bases:
                raise ValueError(f"irregular FASTA wrapping in {name}")
            previous = (len(sequence), len(raw))
            length += len(sequence)
        if not self.index:
            raise ValueError(f"empty FASTA: {self.path}")

    def size(self, contig: str) -> int:
        if contig not in self.index:
            raise ValueError(f"contig {contig!r} is absent from {self.path}")
        return self.index[contig][0]

    def fetch(self, contig: str, start: int, end: int) -> str:
        size = self.size(contig)
        if not 0 <= start < end <= size:
            raise ValueError(f"interval outside assembly: {contig}:{start}-{end} (length {size})")
        _, offset, bases, width = self.index[contig]
        first = offset + start // bases * width + start % bases
        last = offset + (end - 1) // bases * width + (end - 1) % bases + 1
        self.handle.seek(first)
        sequence = self.handle.read(last - first).replace(b"\r", b"").replace(b"\n", b"")
        if len(sequence) != end - start or re.search(rb"[^A-Za-z]", sequence):
            raise ValueError(f"FASTA/index mismatch: {self.path}, {contig}:{start}-{end}")
        return sequence.decode("ascii").upper()

    def close(self):
        self.handle.close()
