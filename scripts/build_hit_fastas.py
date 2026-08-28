#!/usr/bin/env python3
from __future__ import annotations

import argparse
import collections
import json
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from patslib.core import (  # noqa: E402
    FastaRecord,
    extract_fasta_region,
    merge_intervals,
    parse_query_table,
    read_fai,
    read_fasta,
    sanitize_name,
    write_fasta,
)


@dataclass(frozen=True)
class Region:
    sample: str
    contig: str
    start: int
    end: int
    initial: bool = False


def read_hotspots(path: Path) -> dict[str, list[tuple[int, int]]]:
    regions: dict[str, list[tuple[int, int]]] = collections.defaultdict(list)
    if not path.is_file():
        return regions
    with path.open() as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.split()
            if len(fields) < 4:
                raise ValueError(f"{path}:{line_number}: malformed KmerSearcher hotspot row")
            contig = fields[0]
            start, end = sorted((int(fields[-2]), int(fields[-1])))
            if end > start:
                regions[contig].append((start, end))
    return regions


def merged_anchored_hits(
    sample: str,
    hotspots: dict[str, list[tuple[int, int]]],
    sizes: dict[str, int],
    anchor_size: int,
) -> list[Region]:
    result: list[Region] = []
    for contig, intervals in hotspots.items():
        if contig not in sizes:
            print(
                f"[PATs] warning: hotspot contig {contig!r} is absent from {sample}'s FASTA index",
                file=sys.stderr,
            )
            continue
        for start, end in merge_intervals(intervals, max_gap=2 * anchor_size):
            result.append(
                Region(
                    sample,
                    contig,
                    max(0, start - anchor_size),
                    min(sizes[contig], end + anchor_size),
                )
            )
    return result


def read_initial_bed(path: Path, reference_sample: str, sizes: dict[str, int]) -> list[Region]:
    result: list[Region] = []
    if not path.is_file():
        return result
    with path.open() as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.split()
            if len(fields) < 3:
                raise ValueError(f"{path}:{line_number}: malformed initial BED row")
            contig, start, end = fields[0], int(fields[1]), int(fields[2])
            if contig not in sizes:
                raise ValueError(f"initial contig is absent from reference index: {contig}")
            start, end = max(0, start), min(sizes[contig], end)
            if end > start:
                result.append(Region(reference_sample, contig, start, end, True))
    return result


def merge_regions(regions: list[Region]) -> list[Region]:
    by_key: dict[tuple[str, str], list[Region]] = collections.defaultdict(list)
    for region in regions:
        by_key[(region.sample, region.contig)].append(region)
    merged: list[Region] = []
    for (sample, contig), values in sorted(by_key.items()):
        values.sort(key=lambda item: (item.start, item.end))
        current = values[0]
        for region in values[1:]:
            if region.start <= current.end:
                current = Region(
                    sample,
                    contig,
                    current.start,
                    max(current.end, region.end),
                    current.initial or region.initial,
                )
            else:
                merged.append(current)
                current = region
        merged.append(current)
    return merged


def consistent_units_for_reference(
    target_fasta: Path,
    reference_regions: list[Region],
    reference_fasta: Path,
    consistent_script: Path,
    threads: int,
    samtools: str,
) -> list[Region]:
    if not reference_regions:
        raise RuntimeError(
            "FASTA target mode found no KmerSearcher hits in the first-row reference; "
            "find_consistent_units.py has no reference loci to initialize"
        )
    with tempfile.TemporaryDirectory(prefix="pats.consistent.") as temp_name:
        temp = Path(temp_name)
        hit_fasta = temp / "reference_hits.fa"
        output_bed = temp / "consistent.bed"
        mapping: dict[str, Region] = {}
        records: list[FastaRecord] = []
        for index, region in enumerate(reference_regions, 1):
            record_name = f"refhit_{index}"
            mapping[record_name] = region
            sequence = extract_fasta_region(
                reference_fasta,
                region.contig,
                region.start,
                region.end,
                samtools=samtools,
            )
            records.append(
                FastaRecord(record_name, sequence, f"{region.contig}:{region.start}-{region.end}")
            )
        write_fasta(records, hit_fasta)
        subprocess.run(
            [
                sys.executable,
                str(consistent_script),
                "-a",
                str(target_fasta),
                "-b",
                str(hit_fasta),
                "-o",
                str(output_bed),
                "-t",
                str(threads),
            ],
            check=True,
        )
        result: list[Region] = []
        if not output_bed.is_file():
            return result
        with output_bed.open() as handle:
            for raw in handle:
                if not raw.strip() or raw.startswith("#"):
                    continue
                fields = raw.split()
                if len(fields) < 3 or fields[0] not in mapping:
                    continue
                parent = mapping[fields[0]]
                local_start, local_end = int(fields[1]), int(fields[2])
                start = max(parent.start, parent.start + local_start)
                end = min(parent.end, parent.start + local_end)
                if end > start:
                    result.append(
                        Region(parent.sample, parent.contig, start, end, initial=True)
                    )
        if not result:
            raise RuntimeError(
                f"find_consistent_units.py produced no usable reference initial units for {target_fasta}"
            )
        return result


def write_outputs(
    regions: list[Region],
    query_by_sample,
    sizes_by_sample,
    name_prefix: str,
    edge_distance: int,
    passing_path: Path,
    filtered_path: Path,
    bed_path: Path,
    report_path: Path,
    samtools: str,
) -> None:
    passing: list[FastaRecord] = []
    filtered: list[FastaRecord] = []
    reasons: collections.Counter[str] = collections.Counter()
    bed_path.parent.mkdir(parents=True, exist_ok=True)

    with bed_path.open("w") as bed_handle:
        for index, region in enumerate(regions, 1):
            query = query_by_sample[region.sample]
            sizes = sizes_by_sample[region.sample]
            sequence = extract_fasta_region(
                query.fasta,
                region.contig,
                region.start,
                region.end,
                samtools=samtools,
            )
            reason = ""
            if "N" in sequence.upper():
                reason = "N_gap"
            elif not query.whitelist:
                edge = min(region.start, sizes[region.contig] - region.end)
                if edge < edge_distance:
                    reason = "non_whitelist_edge"

            record_name = f"{name_prefix}_{region.sample}_{index}"
            description = (
                f"{region.contig}:{region.start}-{region.end}\tSample={region.sample}"
                f"\tInitial={'YES' if region.initial else 'NO'}"
            )
            record = FastaRecord(record_name, sequence, description)
            status = "FILTERED" if reason else "PASS"
            bed_handle.write(
                "\t".join(
                    [
                        region.contig,
                        str(region.start),
                        str(region.end),
                        record_name,
                        region.sample,
                        status,
                        reason or ".",
                        "YES" if region.initial else "NO",
                    ]
                )
                + "\n"
            )
            if reason:
                filtered.append(record)
                reasons[reason] += 1
            else:
                passing.append(record)

    write_fasta(passing, passing_path)
    write_fasta(filtered, filtered_path)
    if not passing:
        raise RuntimeError("all merged hit regions failed QC; no loci remain for gfixbreaks.py")
    report = {
        "total": len(regions),
        "passing": len(passing),
        "filtered": len(filtered),
        "filtered_reasons": dict(reasons),
    }
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description="Merge KmerSearcher hits, anchor, QC, and emit locus FASTAs.")
    parser.add_argument("--group", required=True)
    parser.add_argument("--name-prefix", required=True)
    parser.add_argument("--mode", choices=["bed", "gene", "fasta"], required=True)
    parser.add_argument("--search-dir", required=True)
    parser.add_argument("--queries", required=True)
    parser.add_argument("--reference-sample", required=True)
    parser.add_argument("--target-fasta", required=True)
    parser.add_argument("--initial-bed", required=True)
    parser.add_argument("--consistent-script", required=True)
    parser.add_argument("--anchor-size", type=int, default=5_000)
    parser.add_argument("--edge-distance", type=int, default=20_000)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--passing", required=True)
    parser.add_argument("--filtered", required=True)
    parser.add_argument("--regions-bed", required=True)
    parser.add_argument("--report", required=True)
    parser.add_argument("--samtools", default="samtools")
    args = parser.parse_args()

    queries = parse_query_table(args.queries)
    query_by_sample = {record.sample: record for record in queries}
    if args.reference_sample not in query_by_sample:
        raise ValueError(f"reference sample is absent from query table: {args.reference_sample}")
    sizes_by_sample = {
        record.sample: read_fai(Path(f"{record.fasta}.fai")) for record in queries
    }

    regions: list[Region] = []
    group_search = Path(args.search_dir) / args.group
    for query in queries:
        hotspots = read_hotspots(group_search / f"{query.sample}_hotspot.txt")
        regions.extend(
            merged_anchored_hits(
                query.sample,
                hotspots,
                sizes_by_sample[query.sample],
                args.anchor_size,
            )
        )

    if args.mode in {"bed", "gene"}:
        regions.extend(
            read_initial_bed(
                Path(args.initial_bed),
                args.reference_sample,
                sizes_by_sample[args.reference_sample],
            )
        )
    else:
        reference_hits = [
            region for region in regions if region.sample == args.reference_sample
        ]
        regions.extend(
            consistent_units_for_reference(
                Path(args.target_fasta),
                reference_hits,
                query_by_sample[args.reference_sample].fasta,
                Path(args.consistent_script),
                args.threads,
                args.samtools,
            )
        )

    regions = merge_regions(regions)
    prefix = sanitize_name(f"{args.name_prefix}{args.group}")
    write_outputs(
        regions,
        query_by_sample,
        sizes_by_sample,
        prefix,
        args.edge_distance,
        Path(args.passing),
        Path(args.filtered),
        Path(args.regions_bed),
        Path(args.report),
        args.samtools,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

