"""Build assay candidates and partition genes using disk-backed k-mer postings."""
from __future__ import annotations

from collections import Counter, defaultdict
import json
import sqlite3

from patslib.core import merge_intervals
from .input import Fasta, Sample, intervals, read_calls


K = 31
COMPLEMENT = str.maketrans("ACGTRYMKSWBDHVN", "TGCAYRKMSWVHDBN")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def count_kmers(parts: list[str]) -> Counter:
    """PATs uses the numerically larger strand as its canonical 31-mer."""
    counts = Counter()
    mask = (1 << (2 * K)) - 1
    bases = {"A": 0, "C": 1, "G": 2, "T": 3}
    for sequence in parts:
        forward = reverse = valid = 0
        for base in sequence.upper():
            value = bases.get(base)
            if value is None:
                forward = reverse = valid = 0
                continue
            forward = ((forward << 2) | value) & mask
            reverse = (reverse >> 2) | ((3 - value) << (2 * (K - 1)))
            valid += 1
            if valid >= K:
                counts[max(forward, reverse)] += 1
    return counts


def locus_groups(path, used_parents: set[str]):
    parent = None
    children = []
    for line, raw, rows in read_calls(path):
        call = {"line": line, "source": str(path), "raw": raw, "alternatives": rows}
        kinds = {row["model_type"] for row in rows}
        if len(kinds) != 1:
            raise ValueError(f"{path}:{line}: mixed parent/transcript alternatives")
        if rows[0]["model_type"] == "full_gene":
            if parent is not None or children:
                yield parent, children
            index = raw["GENE_index"]
            if not index or index in used_parents:
                raise ValueError(f"{path}:{line}: empty/reused parent GENE_index")
            used_parents.add(index)
            parent, children = call, []
        elif raw["GENE_index"]:
            if parent is None or raw["GENE_index"] != parent["raw"]["GENE_index"]:
                raise ValueError(f"{path}:{line}: transcript has no preceding matching parent")
            children.append(call)
        else:
            if parent is not None or children:
                yield parent, children
            parent, children = None, []
            yield None, [call]
    if parent is not None or children:
        yield parent, children


def source_info(call):
    return {key: call[key] for key in ("source", "line", "raw")}


def belongs(child, parent):
    return set(child["gene_id"].split("&")) <= set(parent["gene_id"].split("&"))


def candidates(sample: Sample, mode: str, stats: Counter):
    fasta = Fasta(sample.assembly)
    used = set()
    try:
        for path, fragment in ((sample.calls, False), (sample.fragments, True)):
            if path is None:
                continue
            for parent, children in locus_groups(path, used):
                stats["parent_rows"] += parent is not None
                stats["transcript_rows"] += len(children)
                if mode == "wes" and parent is None:
                    raise ValueError(f"{path}: WES requires full_gene parents and GENE_index")
                for child in children:
                    for row in child["alternatives"]:
                        if parent and not any(belongs(row, p) for p in parent["alternatives"]):
                            raise ValueError(f"{path}:{child['line']}: child gene is absent from its parent")
                if mode == "wes":
                    models = [(p, parent) for p in parent["alternatives"]]
                else:
                    models = [(row, child) for child in children for row in child["alternatives"]]
                for row, call in models:
                    contig, strand = row["query_contig"], row["strand"]
                    size = fasta.size(contig)
                    exon_intervals = intervals(row["exon_query_coordinates"])
                    gene = row["gene_id"]
                    gene_name = row["gene_name"]
                    sources = [source_info(call)]
                    parent_index = row["GENE_index"] or f"row:{path.name}:{call['line']}"
                    if parent and mode == "rna":
                        matching = [p for p in parent["alternatives"] if belongs(row, p)]
                        if len(matching) != 1:
                            raise ValueError(f"{path}:{call['line']}: ambiguous child-to-parent unit mapping")
                        if matching[0]["query_contig"] != contig or matching[0]["strand"] != strand:
                            raise ValueError(f"{path}:{call['line']}: child contig/strand differs from parent")
                        gene, gene_name = matching[0]["gene_id"], matching[0]["gene_name"]
                        sources.append(source_info(parent))
                    if mode == "wes":
                        # Parent representatives can be shorter than the owned union of
                        # aligned exons. Include those blocks, never the intronic span.
                        for field in ("merged_interval_coordinates", "inserted_exon_query_coordinates"):
                            if row.get(field):
                                exon_intervals.extend(intervals(row[field]))
                        for child in children:
                            matches = [c for c in child["alternatives"] if belongs(c, row)]
                            if matches:
                                sources.append(source_info(child))
                            for c in matches:
                                if c["query_contig"] != contig or c["strand"] != strand:
                                    raise ValueError(f"{path}:{child['line']}: child contig/strand differs from parent")
                                exon_intervals.extend(intervals(c["exon_query_coordinates"]))
                        if any(e > size for _, e in exon_intervals):
                            raise ValueError(f"{path}:{call['line']}: exon outside assembly")
                        spans = merge_intervals((max(0, s - 30), min(size, e + 30)) for s, e in exon_intervals)
                        parts = [fasta.fetch(contig, s, e) for s, e in spans]
                        segments = [[span] for span in spans]
                        reconstruction = "anchored_exon_union"
                    else:
                        numbers = [int(x) for x in row["found_exon_numbers"].split(",")]
                        coverages = [float(x) for x in row["exon_coverages"].split(",")]
                        parts, segments = [], []
                        for i, (s, e) in enumerate(exon_intervals):
                            if i and s < exon_intervals[i - 1][1]:
                                raise ValueError(f"{path}:{call['line']}: overlapping exons within one transcript")
                            sequence = fasta.fetch(contig, s, e)
                            # Unknown exon ends or missing exon numbers must not create splice junctions.
                            adjacent = i and numbers[i] - numbers[i - 1] == (1 if strand == "+" else -1)
                            join = adjacent and coverages[i] == 100 and coverages[i - 1] == 100
                            if join:
                                parts[-1] += sequence
                                segments[-1].append((s, e))
                            else:
                                parts.append(sequence)
                                segments.append([(s, e)])
                        if strand == "-":
                            parts = [reverse_complement(p) for p in reversed(parts)]
                            segments = [list(reversed(s)) for s in reversed(segments)]
                        reconstruction = "spliced" if len(parts) == 1 else "segmented"
                        stats["segmented_rna_models"] += len(parts) > 1
                    stats["candidate_models"] += 1
                    stats["partial_models"] += row["call_status"] == "partial"
                    yield parts, {
                        "sample": sample.name, "assembly": str(sample.assembly),
                        "gene_unit": gene, "gene_name": gene_name,
                        "source_gene_id": row["gene_id"],
                        "transcript_id": row["transcript_id"] if mode == "rna" else None,
                        "locus_id": parent_index, "contig": contig, "strand": strand,
                        "call_status": row["call_status"],
                        "assignment_status": row["assignment_status"],
                        "tie_group_id": row["tie_group_id"],
                        "parent_assignment_status": parent["raw"]["assignment_status"] if parent else None,
                        "candidate_fragment": fragment,
                        "reconstruction": reconstruction, "segments": segments,
                        "sources": sources,
                    }
    finally:
        fasta.close()


def create_catalog(path):
    db = sqlite3.connect(path)
    db.executescript("""
        PRAGMA temp_store=FILE;
        CREATE TABLE genes (id INTEGER PRIMARY KEY, name TEXT UNIQUE, partition INTEGER);
        CREATE TABLE candidates (id INTEGER PRIMARY KEY, gene INTEGER, metadata TEXT, parts TEXT);
        CREATE TABLE counts (candidate INTEGER, kmer INTEGER, n INTEGER,
                             PRIMARY KEY(candidate, kmer)) WITHOUT ROWID;
        CREATE TABLE gene_kmers (kmer INTEGER, gene INTEGER,
                                 PRIMARY KEY(kmer, gene)) WITHOUT ROWID;
    """)
    return db


def add_sample(db, sample, mode, stats):
    for parts, metadata in candidates(sample, mode, stats):
        db.execute("INSERT OR IGNORE INTO genes(name) VALUES (?)", (metadata["gene_unit"],))
        gene = db.execute("SELECT id FROM genes WHERE name=?", (metadata["gene_unit"],)).fetchone()[0]
        cursor = db.execute("INSERT INTO candidates(gene, metadata, parts) VALUES (?,?,?)",
                            (gene, json.dumps(metadata), json.dumps(parts)))
        candidate = cursor.lastrowid
        counts = count_kmers(parts)
        stats["models_without_31mers"] += not counts
        db.executemany("INSERT INTO counts VALUES (?,?,?)", ((candidate, k, n) for k, n in counts.items()))
        db.executemany("INSERT OR IGNORE INTO gene_kmers VALUES (?,?)", ((k, gene) for k in counts))
    db.commit()


def partition_genes(db):
    genes = dict(db.execute("SELECT id,name FROM genes ORDER BY name"))
    parent = {gene: gene for gene in genes}

    def root(gene):
        while parent[gene] != gene:
            parent[gene] = parent[parent[gene]]
            gene = parent[gene]
        return gene

    db.execute("""CREATE TABLE edges AS
        SELECT a.gene AS a, b.gene AS b, count(*) AS shared
        FROM gene_kmers a JOIN gene_kmers b ON a.kmer=b.kmer AND a.gene<b.gene
        GROUP BY a.gene,b.gene HAVING count(*) > 100""")
    for a, b in db.execute("SELECT a,b FROM edges ORDER BY a,b"):
        parent[root(b)] = root(a)
    components = defaultdict(list)
    for gene in genes:
        components[root(gene)].append(gene)
    groups = sorted(components.values(), key=lambda group: tuple(sorted(genes[g] for g in group)))
    for number, group in enumerate(groups, 1):
        db.executemany("UPDATE genes SET partition=? WHERE id=?", ((number, g) for g in group))
    db.executescript("""
        CREATE INDEX candidates_gene ON candidates(gene);
        CREATE TABLE shared_kmers (kmer INTEGER PRIMARY KEY);
        INSERT INTO shared_kmers SELECT kmer FROM gene_kmers k JOIN genes g ON k.gene=g.id
            GROUP BY kmer HAVING count(DISTINCT g.partition)>1;
    """)
    db.commit()
    return groups
