#!/usr/bin/env python3
from __future__ import annotations

import argparse
import collections
import re
import sys
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from patslib.core import (  # noqa: E402
    FastaRecord,
    extract_fasta_region,
    parse_query_table,
    read_fai,
    read_fasta,
    sanitize_name,
    write_fasta,
)


REGION_RE = re.compile(r"^([^\s:]+):(-?\d+)-(\d+)([+-])?$")


@dataclass(frozen=True)
class LocatedRecord:
    record: FastaRecord
    sample: str
    contig: str
    start: int
    end: int


def description_region(description: str) -> tuple[str, int, int] | None:
    if not description:
        return None
    match = REGION_RE.match(description.split()[0])
    if not match:
        return None
    return match.group(1), int(match.group(2)), int(match.group(3))


def explicit_sample(description: str) -> str:
    for token in description.split():
        if token.startswith("Sample="):
            return token.split("=", 1)[1]
    return ""


def resolve_contig_alias(
    contig: str, contig_to_samples: dict[str, list[str]]
) -> str:
    """Resolve a gfixbreaks ``contig_N`` locus key using query FASTA indexes."""

    if contig in contig_to_samples:
        return contig
    match = re.fullmatch(r"(.+)_([0-9]+)", contig)
    if match and match.group(1) in contig_to_samples:
        return match.group(1)
    return contig


def identify_sample(
    record: FastaRecord,
    contig: str,
    sample_names: list[str],
    contig_to_samples: dict[str, list[str]],
) -> str:
    sample = explicit_sample(record.description)
    if sample:
        return sample
    by_name = [sample for sample in sample_names if sample in record.name]
    if len(by_name) == 1:
        return by_name[0]
    candidates = contig_to_samples.get(contig, [])
    if len(candidates) == 1:
        return candidates[0]
    raise ValueError(
        f"cannot determine assembly for FASTA record {record.name!r} on contig {contig!r}"
    )


def merge_located(records: list[LocatedRecord]) -> list[tuple[str, str, int, int]]:
    by_key: dict[tuple[str, str], list[tuple[int, int]]] = collections.defaultdict(list)
    for located in records:
        by_key[(located.sample, located.contig)].append((located.start, located.end))
    merged: list[tuple[str, str, int, int]] = []
    for (sample, contig), intervals in sorted(by_key.items()):
        intervals.sort()
        start, end = intervals[0]
        for next_start, next_end in intervals[1:]:
            if next_start <= end:
                end = max(end, next_end)
            else:
                merged.append((sample, contig, start, end))
                start, end = next_start, next_end
        merged.append((sample, contig, start, end))
    return merged


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Add QC-filtered locus regions back to gfixbreaks output before k-mer selection."
    )
    parser.add_argument("-i", "--input", required=True, help="gfixbreaks FASTA")
    parser.add_argument("-a", "--add", required=True, help="QC-filtered FASTA")
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("-q", "--queries", required=True)
    parser.add_argument("--prefix", default="PATs", help="output FASTA record prefix")
    parser.add_argument("--samtools", default="samtools")
    args = parser.parse_args()

    queries = parse_query_table(args.queries)
    query_by_sample = {query.sample: query for query in queries}
    contig_to_samples: dict[str, list[str]] = collections.defaultdict(list)
    for query in queries:
        for contig in read_fai(Path(f"{query.fasta}.fai")):
            contig_to_samples[contig].append(query.sample)

    located: list[LocatedRecord] = []
    opaque: list[FastaRecord] = []
    reported_aliases: set[tuple[str, str]] = set()
    for fasta_path in (Path(args.input), Path(args.add)):
        if not fasta_path.is_file() or fasta_path.stat().st_size == 0:
            continue
        for record in read_fasta(fasta_path):
            region = description_region(record.description)
            if region is None:
                opaque.append(record)
                continue
            contig, start, end = region
            resolved_contig = resolve_contig_alias(contig, contig_to_samples)
            if resolved_contig != contig and (contig, resolved_contig) not in reported_aliases:
                print(
                    f"[PATs] internal locus contig {contig!r} resolved to FASTA contig "
                    f"{resolved_contig!r}",
                    file=sys.stderr,
                )
                reported_aliases.add((contig, resolved_contig))
            contig = resolved_contig
            sample = identify_sample(
                record,
                contig,
                list(query_by_sample),
                contig_to_samples,
            )
            located.append(LocatedRecord(record, sample, contig, start, end))

    output_records: list[FastaRecord] = list(opaque)
    for index, (sample, contig, start, end) in enumerate(merge_located(located), 1):
        sequence = extract_fasta_region(
            query_by_sample[sample].fasta,
            contig,
            start,
            end,
            samtools=args.samtools,
        )
        output_records.append(
            FastaRecord(
                f"{sanitize_name(args.prefix)}_{sample}_{index}",
                sequence,
                f"{contig}:{start}-{end}\tSample={sample}",
            )
        )
    if not output_records:
        raise RuntimeError("addregion.py received no FASTA records")
    write_fasta(output_records, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
