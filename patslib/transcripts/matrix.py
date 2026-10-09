"""Graph-free adapter to the existing PATs/Ctyper v2.0.1 wire format."""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
import json
from urllib.parse import quote

from scripts.matrixformat import INDEX_VERSION_LINE, encode_fixed, encode_metadata


def integer(value):
    return encode_fixed(value, max(1, (int(value).bit_length() + 5) // 6), "integer")


@dataclass
class Node:
    name: str = ""
    height: float = 0.0
    children: list = field(default_factory=list)
    parent: object = None
    leaf: int = -1


def make_tree(profiles, names):
    """Average-linkage tree over weighted Jaccard distances of count profiles."""
    import numpy as np
    from scipy.cluster.hierarchy import linkage

    nodes = [Node(name=name, leaf=i) for i, name in enumerate(names)]
    if len(nodes) > 1:
        distances = np.empty(len(nodes) * (len(nodes) - 1) // 2, dtype=float)
        totals = [sum(p.values()) for p in profiles]
        pos = 0
        for i, left in enumerate(profiles):
            for j in range(i + 1, len(profiles)):
                right = profiles[j]
                small, large = (left, right) if len(left) <= len(right) else (right, left)
                overlap = sum(min(n, large.get(k, 0)) for k, n in small.items())
                distances[pos] = 1 - overlap / (totals[i] + totals[j] - overlap)
                pos += 1
        for a, b, distance, _ in linkage(distances, method="average"):
            node = Node(height=float(distance) / 2, children=[nodes[int(a)], nodes[int(b)]])
            for child in node.children:
                child.parent = node
            nodes.append(node)
    root = nodes[-1]
    preorder, order, stack = [], [], [root]
    while stack:
        node = stack.pop()
        preorder.append(node)
        if node.children:
            stack.extend(reversed(node.children))
        else:
            order.append(node.leaf)
            node.leaf = len(order) - 1
    rendered = {}
    for node in reversed(preorder):
        label = "(" + ",".join(rendered[id(c)] for c in node.children) + ")" if node.children else node.name
        distance = max(0.0, node.parent.height - node.height) if node.parent else 0.0
        rendered[id(node)] = f"{label}:{distance:.7f}"
    return rendered[id(root)] + ";", preorder, order


def annotate_tree(nodes, counts):
    """Native PATs node deltas, using the reader's preorder node numbering."""
    totals, sizes = {}, {}
    for node in reversed(nodes):
        key = id(node)
        if node.children:
            totals[key] = sum(totals[id(c)] for c in node.children)
            sizes[key] = sum(sizes[id(c)] for c in node.children)
        else:
            totals[key], sizes[key] = counts.get(node.leaf, 0), 1
    rounded = {key: (2 * total + sizes[key]) // (2 * sizes[key]) for key, total in totals.items()}
    positive, negative = [], []
    for index, node in enumerate(nodes):
        delta = rounded[id(node)] - (rounded[id(node.parent)] if node.parent else 0)
        (positive if delta >= 0 else negative).extend([index] * abs(delta))
    return positive, negative


def allele_title(candidate, metadata):
    # Native title separators must not be introduced by arbitrary assembly IDs.
    contig = quote(metadata["sample"], safe="") + "#" + quote(metadata["contig"], safe="")
    spans = [span for segment in metadata["segments"] for span in segment]
    start, end = min(s for s, _ in spans), max(e for _, e in spans)
    return f">tx_{candidate:09d}\t{contig}:{start}-{end}{metadata['strand']}\tNA:NA"


def load_partition(db, number):
    profiles, members, signatures = [], [], {}
    omitted = []
    rows = db.execute("""SELECT c.id,c.metadata FROM candidates c JOIN genes g ON c.gene=g.id
                         WHERE g.partition=? ORDER BY c.id""", (number,)).fetchall()
    for candidate, metadata in rows:
        profile = dict(db.execute("""SELECT kmer,n FROM counts WHERE candidate=?
            AND NOT EXISTS (SELECT 1 FROM shared_kmers s WHERE s.kmer=counts.kmer) ORDER BY kmer""", (candidate,)))
        if not profile:
            omitted.append(candidate)
            continue
        signature = tuple(profile.items())
        if signature not in signatures:
            signatures[signature] = len(profiles)
            profiles.append(profile)
            members.append([])
        members[signatures[signature]].append((candidate, json.loads(metadata)))
    return profiles, members, omitted


def write_partition(db, number, name, matrix, index, trees, reference_sample, max_candidates):
    profiles, members, omitted = load_partition(db, number)
    if not profiles:
        return {}, {"partition": name, "matrix_candidates": 0, "omitted_candidates": omitted}
    if len(profiles) > min(max_candidates, 32768):
        raise ValueError(f"{name}: {len(profiles)} candidate profiles exceed the partition limit ({max_candidates}); "
                         "the tree uses quadratic memory and Ctyper has 16-bit node indexes")
    names = [f"tx_{group[0][0]:09d}" for group in members]
    newick, nodes, order = make_tree(profiles, names)
    profiles = [profiles[i] for i in order]
    members = [members[i] for i in order]
    kmers = defaultdict(dict)
    for i, profile in enumerate(profiles):
        if sum(profile.values()) > 2**32 - 1:
            raise ValueError(f"{name}: candidate k-mer count exceeds Ctyper uint32 capacity")
        for kmer, n in profile.items():
            kmers[kmer][i] = n
    counts_total = sum(sum(p.values()) for p in profiles)
    start = matrix.tell()

    def emit(text):
        encoded = (text + "\n").encode("ascii")
        if len(encoded) >= 2_000_000:
            raise ValueError(f"{name}: row exceeds Ctyper's line buffer")
        matrix.write(encoded)

    emit(f"#{name}\t{len(kmers)}\t{len(profiles)}\t{counts_total}\t" + " " * 30)
    emit("%" + " ".join(str(sum(p.values())) for p in profiles))
    emit("!" + newick)
    assignments = {}
    for column, group in enumerate(members):
        emit(";".join(allele_title(candidate, metadata) for candidate, metadata in group))
        for candidate, _ in group:
            assignments[candidate] = {"partition": name, "column": column, "representative": f"tx_{group[0][0]:09d}"}
    # Preserve native group records, without silently clustering distinguishable isoforms.
    for i, profile in enumerate(profiles):
        emit(f"@{len(profile)}\t{i},")
    for i in range(len(profiles)):
        emit(f"?{i},")
    matrix_size = phylo_size = 0
    # Ctyper's post-correction tag: bias=0, variance=100, reciprocal=100.
    # Raw encode_tag(0,1) would decode to variance 26 and be discarded by
    # Ctyper's variance threshold. These neutral values imply no calibration.
    metadata_fields = ("!{S", *encode_metadata(0, 0, 0, 0, 0, 0, 1.0)[1:])
    for kmer, counts in sorted(kmers.items()):
        occurrences = sum(counts.values())
        if occurrences + 6 > 65535:
            raise ValueError(f"{name}: k-mer multiplicity exceeds Ctyper's uint16 row capacity")
        positive, negative = annotate_tree(nodes, counts)
        entries = ",".join(integer(i) + ("*" + integer(n) if n > 1 else "") for i, n in sorted(counts.items())) + ","
        phylo = " ".join(",".join(integer(i) for i in values) + ("," if values else "") for values in (positive, negative))
        emit("\t".join(("&+", *metadata_fields, encode_fixed(kmer, 11, "k-mer"), entries + " " + phylo)))
        matrix_size += occurrences + 15
        phylo_size += len(positive) + len(negative) + 2
    stop = matrix.tell()
    labels, reference_regions = set(), set()
    for group in members:
        for _, metadata in group:
            labels.update(metadata["gene_unit"].split("&"))
            labels.update(metadata["gene_name"].split("&"))
            if metadata["sample"] == reference_sample:
                contig = metadata["contig"]
                if any(c in contig for c in ":;\t\r\n"):
                    raise ValueError("reference contig cannot be represented in native index regions")
                for segment in metadata["segments"]:
                    for s, e in segment:
                        reference_regions.add(f"{contig}:{s+1}-{e}{metadata['strand']};")
    label_text = "".join(quote(label, safe="_.-&") + ";" for label in sorted(labels) if label)
    index.write(f"#{name}\t{start}\t{stop}\t{len(kmers)}\t{matrix_size}\t{phylo_size}\t0\t{len(profiles)}\t"
                + "".join(sorted(reference_regions)) + "\t" + label_text + "\n")
    trees.write(name + "\t" + newick + "\n")
    return assignments, {"partition": name, "matrix_candidates": len(profiles), "kmers": len(kmers),
                         "omitted_candidates": omitted, "equivalent_models": sum(len(g) - 1 for g in members)}
