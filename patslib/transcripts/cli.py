"""Opt-in `pats.py transcripts run` command. The existing workflow is untouched."""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import json
import os
from pathlib import Path
import tempfile

from patslib.core import sanitize_name
from .catalog import add_sample, create_catalog, partition_genes
from .input import read_manifest
from .matrix import INDEX_VERSION_LINE, write_partition


SUFFIXES = (".fa", ".matrix.txt", ".matrix.txt.index", ".annotations.jsonl",
            ".partitions.tsv", ".merges.tsv", ".trees.tsv", ".report.json")


def build_parser():
    parser = argparse.ArgumentParser(prog="pats.py transcripts", description="Experimental ParaTyper matrix builder")
    sub = parser.add_subparsers(dest="command", required=True)
    run = sub.add_parser("run", help="build a native Ctyper/PATs WES or RNA matrix")
    run.add_argument("--mode", choices=("wes", "rna"), required=True)
    run.add_argument("--manifest", type=Path, required=True, help="TSV: sample, assembly, calls, optional fragments")
    run.add_argument("-o", "--output", type=Path, required=True, help="output prefix")
    run.add_argument("--reference-sample", default="", help="optional sample used for native index reference intervals")
    run.add_argument("--max-partition-candidates", type=int, default=5000,
                     help="limit unique profiles before allocating the quadratic tree distances (default 5000)")
    run.add_argument("--scratch-dir", type=Path, help="temporary SQLite catalog location (requires sufficient disk)")
    run.add_argument("--force", action="store_true", help="replace existing outputs after a successful build")
    return parser


def execute(args):
    if not 1 <= args.max_partition_candidates <= 32768:
        raise ValueError("--max-partition-candidates must be between 1 and 32768")
    manifest = args.manifest.expanduser().resolve()
    samples = read_manifest(manifest)
    if args.reference_sample and args.reference_sample not in {s.name for s in samples}:
        raise ValueError("--reference-sample is absent from the manifest")
    output = args.output.expanduser().resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    destinations = [Path(str(output) + suffix) for suffix in SUFFIXES]
    inputs = {manifest} | {p for s in samples for p in (s.assembly, s.calls, s.fragments) if p is not None}
    if any(p.resolve() in inputs for p in destinations):
        raise ValueError("an output path would replace an input")
    if not args.force and any(p.exists() for p in destinations):
        raise ValueError("output already exists; choose a new prefix or pass --force")
    stats = Counter()
    with tempfile.TemporaryDirectory(prefix="pats-transcripts.", dir=args.scratch_dir) as scratch, \
            tempfile.TemporaryDirectory(prefix=".pats-transcripts.", dir=output.parent) as staging:
        stage = Path(staging)
        db = create_catalog(Path(scratch) / "catalog.sqlite")
        try:
            for sample in samples:
                print(f"[PATs transcripts] importing {sample.name}", flush=True)
                add_sample(db, sample, args.mode, stats)
            if not stats["candidate_models"]:
                raise ValueError("no candidate models in the supplied annotations")
            print("[PATs transcripts] merging genes with >100 shared canonical 31-mers", flush=True)
            groups = partition_genes(db)
            name_prefix = sanitize_name(output.name, fallback="transcripts")
            results = []
            with (stage / "matrix.txt").open("wb") as matrix, \
                    (stage / "matrix.txt.index").open("w") as index, \
                    (stage / "trees.tsv").open("w") as trees:
                index.write(INDEX_VERSION_LINE + "\n")
                trees.write("partition\tnewick\n")
                db.execute("CREATE TABLE assignments (candidate INTEGER PRIMARY KEY, data TEXT)")
                for number in range(1, len(groups) + 1):
                    name = f"{name_prefix}_group{number}"
                    print(f"[PATs transcripts] matrix {number}/{len(groups)}: {name}", flush=True)
                    assignments, report = write_partition(db, number, name, matrix, index, trees,
                                                          args.reference_sample, args.max_partition_candidates)
                    db.executemany("INSERT INTO assignments VALUES (?,?)", ((c, json.dumps(a)) for c, a in assignments.items()))
                    results.append(report)
            if not any(r["matrix_candidates"] for r in results):
                raise ValueError("no informative 31-mers remain after excluding k-mers shared across partitions")
            with (stage / "fa").open("w") as fasta, (stage / "annotations.jsonl").open("w") as annotations:
                for candidate, metadata, parts, assignment in db.execute("""SELECT c.id,c.metadata,c.parts,a.data
                    FROM candidates c LEFT JOIN assignments a ON a.candidate=c.id ORDER BY c.id"""):
                    metadata, parts = json.loads(metadata), json.loads(parts)
                    metadata["candidate_id"] = f"tx_{candidate:09d}"
                    metadata["mode"] = args.mode
                    metadata["matrix"] = json.loads(assignment) if assignment else None
                    metadata["matrix_exclusion"] = None if assignment else "no_partition_exclusive_31mers"
                    fasta_ids = []
                    for i, sequence in enumerate(parts, 1):
                        name = metadata["candidate_id"] + (f".segment{i}" if len(parts) > 1 or args.mode == "wes" else "")
                        fasta_ids.append(name)
                        fasta.write(f">{name}\n")
                        for start in range(0, len(sequence), 80):
                            fasta.write(sequence[start:start+80] + "\n")
                    metadata["fasta_records"] = fasta_ids
                    annotations.write(json.dumps(metadata, sort_keys=True) + "\n")
            with (stage / "partitions.tsv").open("w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(("gene_unit", "partition"))
                writer.writerows((gene, f"{name_prefix}_group{number}") for gene, number in db.execute("SELECT name,partition FROM genes ORDER BY name"))
            with (stage / "merges.tsv").open("w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(("gene_a", "gene_b", "shared_31mers"))
                writer.writerows(db.execute("SELECT a.name,b.name,e.shared FROM edges e JOIN genes a ON a.id=e.a JOIN genes b ON b.id=e.b ORDER BY a.name,b.name"))
            report = {
                "schema_version": 1, "mode": args.mode, "kmer_size": 31,
                "wes_anchor": 30 if args.mode == "wes" else 0,
                "matrix_format": INDEX_VERSION_LINE, "graph_paths": False,
                "partition_rule": "connected components of original gene units sharing >100 distinct canonical 31-mers",
                "specificity_scope": "supplied ParaTyper candidate catalog; no whole-genome off-target screen",
                "shared_partition_kmers_excluded": db.execute("SELECT count(*) FROM shared_kmers").fetchone()[0],
                "tree": "average linkage / weighted Jaccard of retained k-mer count profiles",
                "manifest": str(manifest), "samples": [s.name for s in samples],
                "reference_sample": args.reference_sample or None,
                "stats": dict(stats), "partitions": results,
            }
            (stage / "report.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
            # Publish only after all extraction, encoding and range checks succeed.
            for suffix, destination in zip(SUFFIXES, destinations):
                os.replace(stage / suffix[1:], destination)
        finally:
            db.close()
    print(f"[PATs transcripts] wrote {output}.matrix.txt and index", flush=True)


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        execute(args)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    return 0
