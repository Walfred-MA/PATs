#!/usr/bin/env python3
from __future__ import annotations

import argparse
import collections
import json
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import unquote

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from patslib.core import (  # noqa: E402
    FastaRecord,
    ensure_fasta_index,
    extract_fasta_region,
    merge_intervals,
    read_fai,
    read_fasta,
    sanitize_name,
    write_fasta,
)


@dataclass(frozen=True)
class BedTarget:
    contig: str
    start: int
    end: int
    name: str
    genes: tuple[str, ...] = ()


@dataclass(frozen=True)
class GffFeature:
    contig: str
    feature: str
    start: int
    end: int
    strand: str
    attributes: dict[str, str]

    @property
    def gene_name(self) -> str:
        return self.attributes.get("gene_name", "")


def parse_attributes(text: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for item in text.strip().strip(";").split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
        elif " " in item:
            key, value = item.split(None, 1)
            value = value.strip('"')
        else:
            continue
        attributes[key.strip()] = unquote(value.strip())
    return attributes


def read_gff(path: Path) -> list[GffFeature]:
    features: list[GffFeature] = []
    with path.open() as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"{path}:{line_number}: expected nine GFF3 columns")
            start = int(fields[3]) - 1
            end = int(fields[4])
            if start < 0 or end <= start:
                raise ValueError(f"{path}:{line_number}: invalid GFF3 coordinates")
            features.append(
                GffFeature(
                    contig=fields[0],
                    feature=fields[2].lower(),
                    start=start,
                    end=end,
                    strand=fields[6],
                    attributes=parse_attributes(fields[8]),
                )
            )
    return features


def normalize_gff_contigs(
    features: list[GffFeature], sizes: dict[str, int]
) -> tuple[list[GffFeature], dict[str, str]]:
    """Resolve numbered LiftOff aliases such as NC_060946.1_1.

    A RefSeq accession already has its assembly-version suffix after the dot.
    Some lifted annotations append another ``_N`` copy suffix.  When the
    underlying accession exists in the reference FASTA, use that accession
    consistently for both gene anchors and exon extraction.
    """

    remapping: dict[str, str] = {}
    normalized: list[GffFeature] = []
    for feature in features:
        contig = feature.contig
        match = re.fullmatch(r"(NC_\d+\.\d+)_\d+", contig)
        if match and match.group(1) in sizes:
            resolved = match.group(1)
            remapping[contig] = resolved
            feature = GffFeature(
                contig=resolved,
                feature=feature.feature,
                start=feature.start,
                end=feature.end,
                strand=feature.strand,
                attributes=feature.attributes,
            )
        normalized.append(feature)
    return normalized, remapping


def parse_bed(path: Path, *, require_name: bool = False) -> list[BedTarget]:
    targets: list[BedTarget] = []
    used_names: collections.Counter[str] = collections.Counter()
    with path.open() as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.split()
            if len(fields) < 3:
                raise ValueError(f"{path}:{line_number}: expected at least three BED columns")
            contig, start, end = fields[0], int(fields[1]), int(fields[2])
            if start < 0 or end <= start:
                raise ValueError(f"{path}:{line_number}: invalid BED interval")
            if require_name and len(fields) < 4:
                raise ValueError(f"{path}:{line_number}: target BED requires column 4 names")
            base_name = sanitize_name(fields[3] if len(fields) >= 4 else f"target{len(targets) + 1}")
            used_names[base_name] += 1
            name = base_name if used_names[base_name] == 1 else f"{base_name}_{used_names[base_name]}"
            targets.append(BedTarget(contig, start, end, name))
    if not targets:
        raise ValueError(f"BED contains no intervals: {path}")
    return targets


def clip_targets(targets: list[BedTarget], sizes: dict[str, int]) -> list[BedTarget]:
    clipped: list[BedTarget] = []
    for target in targets:
        if target.contig not in sizes:
            raise ValueError(f"target contig is absent from the reference index: {target.contig}")
        start = max(0, target.start)
        end = min(sizes[target.contig], target.end)
        if end <= start:
            raise ValueError(f"target becomes empty after clipping: {target.contig}:{target.start}-{target.end}")
        clipped.append(BedTarget(target.contig, start, end, target.name, target.genes))
    return clipped


def gene_targets(
    features: list[GffFeature],
    prefixes: list[str],
    extension: int,
    merge_distance: int,
    sizes: dict[str, int],
) -> tuple[list[BedTarget], dict[str, list[str]], set[str]]:
    genes = [feature for feature in features if feature.feature == "gene" and feature.gene_name]
    found: dict[str, list[str]] = {prefix: [] for prefix in prefixes}
    selected: list[GffFeature] = []
    selected_keys: set[tuple[str, int, int, str]] = set()

    for gene in genes:
        matched = [prefix for prefix in prefixes if gene.gene_name.startswith(prefix)]
        if not matched:
            continue
        for prefix in matched:
            if gene.gene_name not in found[prefix]:
                found[prefix].append(gene.gene_name)
        key = (gene.contig, gene.start, gene.end, gene.gene_name)
        if key not in selected_keys:
            selected.append(gene)
            selected_keys.add(key)

    missing = [prefix for prefix, names in found.items() if not names]
    if missing:
        raise ValueError(f"gene prefix(es) not found in GFF3: {', '.join(missing)}")

    by_contig: dict[str, list[tuple[int, int, str]]] = collections.defaultdict(list)
    for gene in selected:
        if gene.contig not in sizes:
            raise ValueError(f"GFF3 gene contig is absent from reference: {gene.contig}")
        by_contig[gene.contig].append((gene.start, gene.end, gene.gene_name))

    targets: list[BedTarget] = []
    for contig in sorted(by_contig):
        intervals = sorted(by_contig[contig], key=lambda item: (item[0], item[1], item[2]))
        current_start, current_end, first_gene = intervals[0]
        current_genes = [first_gene]
        for start, end, gene_name in intervals[1:]:
            if start - current_end <= merge_distance:
                current_end = max(current_end, end)
                if gene_name not in current_genes:
                    current_genes.append(gene_name)
            else:
                name = sanitize_name("__".join(current_genes))
                targets.append(
                    BedTarget(
                        contig,
                        max(0, current_start - extension),
                        min(sizes[contig], current_end + extension),
                        name,
                        tuple(current_genes),
                    )
                )
                current_start, current_end, current_genes = start, end, [gene_name]
        name = sanitize_name("__".join(current_genes))
        targets.append(
            BedTarget(
                contig,
                max(0, current_start - extension),
                min(sizes[contig], current_end + extension),
                name,
                tuple(current_genes),
            )
        )

    # The same gene label can occur on alternate/reference contigs.  FASTA
    # record identifiers must still be unique for partitioning and manifests.
    name_counts: collections.Counter[str] = collections.Counter()
    unique_targets: list[BedTarget] = []
    for target in targets:
        name_counts[target.name] += 1
        name = (
            target.name
            if name_counts[target.name] == 1
            else f"{target.name}_{name_counts[target.name]}"
        )
        unique_targets.append(
            BedTarget(target.contig, target.start, target.end, name, target.genes)
        )

    return unique_targets, found, {gene.gene_name for gene in selected}


def write_bed(targets: list[BedTarget], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for target in targets:
            fields = [target.contig, str(target.start), str(target.end), target.name]
            if target.genes:
                fields.append("Genes=" + ",".join(target.genes))
            handle.write("\t".join(fields) + "\n")


def target_fasta_records(targets: list[BedTarget], reference: Path, samtools: str):
    for target in targets:
        sequence = extract_fasta_region(
            reference, target.contig, target.start, target.end, samtools=samtools
        )
        description = f"{target.contig}:{target.start}-{target.end}"
        if target.genes:
            description += "\tGenes=" + ",".join(target.genes)
        yield FastaRecord(target.name, sequence, description)


def extract_exon_bed(exon_bed: Path, reference: Path, sizes: dict[str, int], samtools: str):
    targets = clip_targets(parse_bed(exon_bed), sizes)
    for index, target in enumerate(targets, 1):
        sequence = extract_fasta_region(
            reference, target.contig, target.start, target.end, samtools=samtools
        )
        yield FastaRecord(
            f"Exon_{index}", sequence, f"{target.contig}:{target.start}-{target.end}"
        )


def extract_gene_exons(
    features: list[GffFeature],
    selected_gene_names: set[str],
    reference: Path,
    sizes: dict[str, int],
    samtools: str,
):
    by_contig: dict[str, list[tuple[int, int]]] = collections.defaultdict(list)
    interval_genes: dict[tuple[str, int, int], set[str]] = collections.defaultdict(set)
    for feature in features:
        if feature.feature != "exon" or feature.gene_name not in selected_gene_names:
            continue
        if feature.contig not in sizes:
            continue
        start = max(0, feature.start)
        end = min(sizes[feature.contig], feature.end)
        if end > start:
            by_contig[feature.contig].append((start, end))
            interval_genes[(feature.contig, start, end)].add(feature.gene_name)

    index = 0
    for contig in sorted(by_contig):
        for start, end in merge_intervals(by_contig[contig]):
            index += 1
            genes = sorted(
                {
                    gene
                    for (one_contig, one_start, one_end), names in interval_genes.items()
                    if one_contig == contig and one_start < end and one_end > start
                    for gene in names
                }
            )
            sequence = extract_fasta_region(reference, contig, start, end, samtools=samtools)
            yield FastaRecord(
                f"Exon_{index}", sequence, f"{contig}:{start}-{end}\tGenes={','.join(genes)}"
            )


def copy_normalized_fasta(source: Path, destination: Path) -> int:
    records = list(read_fasta(source))
    if not records:
        raise ValueError(f"FASTA contains no records: {source}")
    names: set[str] = set()
    for record in records:
        if record.name in names:
            raise ValueError(f"duplicate FASTA record name: {record.name}")
        names.add(record.name)
    return write_fasta(records, destination)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Normalize PATs targets and pooled exon sequences.")
    parser.add_argument("--mode", choices=["bed", "gene", "fasta"], required=True)
    parser.add_argument("--input", default="")
    parser.add_argument("--genes-json", default="[]")
    parser.add_argument("--gff3", default="")
    parser.add_argument("--exon", default="")
    parser.add_argument("--reference", required=True)
    parser.add_argument("--gene-extension", type=int, default=5_000)
    parser.add_argument("--gene-merge-distance", type=int, default=10_000)
    parser.add_argument("--targets-bed", required=True)
    parser.add_argument("--targets-fasta", required=True)
    parser.add_argument("--exons-fasta", required=True)
    parser.add_argument("--report", required=True)
    parser.add_argument("--samtools", default="samtools")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    reference = Path(args.reference)
    index = ensure_fasta_index(reference, samtools=args.samtools)
    sizes = read_fai(index)
    targets_bed = Path(args.targets_bed)
    targets_fasta = Path(args.targets_fasta)
    exons_fasta = Path(args.exons_fasta)
    report_path = Path(args.report)
    report_path.parent.mkdir(parents=True, exist_ok=True)

    report: dict = {"mode": args.mode, "reference": str(reference), "matched_genes": {}}
    if args.mode == "bed":
        targets = clip_targets(parse_bed(Path(args.input)), sizes)
        write_bed(targets, targets_bed)
        write_fasta(target_fasta_records(targets, reference, args.samtools), targets_fasta)
        if args.exon:
            write_fasta(
                extract_exon_bed(Path(args.exon), reference, sizes, args.samtools), exons_fasta
            )
        else:
            write_fasta([], exons_fasta)
    elif args.mode == "gene":
        prefixes = json.loads(args.genes_json)
        features = read_gff(Path(args.gff3))
        features, contig_remapping = normalize_gff_contigs(features, sizes)
        report["gff3_contig_remapping"] = contig_remapping
        for source, resolved in sorted(contig_remapping.items()):
            print(
                f"[PATs] GFF3 contig {source!r} resolved to reference contig {resolved!r}",
                file=sys.stderr,
            )
        targets, found, selected_gene_names = gene_targets(
            features,
            prefixes,
            args.gene_extension,
            args.gene_merge_distance,
            sizes,
        )
        report["matched_genes"] = found
        write_bed(targets, targets_bed)
        write_fasta(target_fasta_records(targets, reference, args.samtools), targets_fasta)
        write_fasta(
            extract_gene_exons(features, selected_gene_names, reference, sizes, args.samtools),
            exons_fasta,
        )
        for prefix in prefixes:
            print(f"[PATs] gene prefix {prefix!r}: {', '.join(found[prefix])}", file=sys.stderr)
    else:
        copy_normalized_fasta(Path(args.input), targets_fasta)
        targets_bed.parent.mkdir(parents=True, exist_ok=True)
        targets_bed.write_text("")
        if args.exon:
            first_nonempty = next(
                (line for line in Path(args.exon).read_text().splitlines() if line.strip()), ""
            )
            if not first_nonempty.startswith(">"):
                raise ValueError("FASTA target mode requires --exon to be FASTA")
            copy_normalized_fasta(Path(args.exon), exons_fasta)
        else:
            write_fasta([], exons_fasta)

    report["target_records"] = sum(1 for _ in read_fasta(targets_fasta))
    report["exon_records"] = sum(1 for _ in read_fasta(exons_fasta))
    if report["target_records"] == 0:
        raise ValueError("target preparation produced no FASTA records")
    with report_path.open("w") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
