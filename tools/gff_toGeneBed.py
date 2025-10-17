#!/usr/bin/env python3
import argparse
from pathlib import Path
from collections import defaultdict

def parse_attrs(attr_str):
    out = {}
    for part in attr_str.strip().split(';'):
        if not part:
            continue
        if '=' in part:
            k, v = part.split('=', 1)
            out[k] = v
        else:
            kv = part.strip().split()
            if len(kv) == 2:
                out[kv[0]] = kv[1]
    return out

def merge_intervals(intervals):
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = [list(intervals[0])]
    for s, e in intervals[1:]:
        last = merged[-1]
        if s <= last[1]:
            if e > last[1]:
                last[1] = e
        else:
            merged.append([s, e])
    return [(s, e) for s, e in merged]

def _read_gene_terms_from_file(path: Path):
    terms = []
    with path.open() as f:
        for line in f:
            t = line.strip()
            if t:
                terms.append(t)
    return terms

def _split_exact_and_prefix(terms):
    """
    Separate exact terms and prefix patterns (ending with '*').
    Returns (exact_set, prefix_list) where prefix_list contains strings with the '*' removed.
    """
    exact = set()
    pref = []
    for t in terms:
        if t.endswith('*'):
            stem = t[:-1]
            if stem:
                pref.append(stem)
        else:
            exact.add(t)
    return exact, pref

def build_gene_selector(sel_arg):
    """
    If sel_arg is None -> return (None, None) meaning 'no filter'.
    If file exists -> read one per line; otherwise split by comma.
    Supports prefix patterns ending with '*' (e.g., 'SMN*').
    """
    if sel_arg is None:
        return None, None
    p = Path(sel_arg)
    if p.exists() and p.is_file():
        terms = _read_gene_terms_from_file(p)
    else:
        terms = [t for t in sel_arg.split(',') if t]
    return _split_exact_and_prefix(terms)

def gene_keys_for_match(gene_id, gene_name):
    """
    Candidate identifiers to match selectors against.
    """
    keys = set()
    if gene_id:   keys.add(gene_id)
    if gene_name: keys.add(gene_name)
    return keys

def matches_selection(candidates, exact_set, prefix_list):
    """
    True if any candidate is in exact_set OR startswith any prefix in prefix_list.
    """
    if exact_set is None and prefix_list is None:
        return True
    # exact
    if exact_set and any(c in exact_set for c in candidates):
        return True
    # prefix
    if prefix_list:
        for c in candidates:
            for pre in prefix_list:
                if c.startswith(pre):
                    return True
    return False

def main():
    ap = argparse.ArgumentParser(
        prog="gff_toGeneBed.py",
        description=(
            "Convert a GFF3 to a BED of per-gene intervals.\n\n"
            "• No -g: include ALL genes.\n"
            "• -g X: filter genes. If X is a file, read one per line; else treat X as comma-separated.\n"
            "        Items ending with '*' are prefixes (e.g., 'SMN*'). Matches against gene ID and Name.\n"
            "• -e: add merged exon intervals as 'Target=' in column 5 (absolute genomic, 0-based, half-open).\n"
            "• -a N: expand BED start/end by N bp on both sides (start floored at 0; end not capped).\n\n"
            "Coordinate conventions: GFF3 is 1-based inclusive; BED is 0-based half-open."
        ),
        epilog=(
            "Examples:\n"
            "  # All genes, with exon Targets and 1kb anchors\n"
            "  gff_toGeneBed.py -i genes.gff3 -o genes.bed -e -a 1000\n\n"
            "  # Only selected genes (exact + prefixes)\n"
            "  gff_toGeneBed.py -i genes.gff3 -o sel.bed -g 'SMN*,BRCA1' -e\n\n"
            "  # Selected genes from a file (one per line; supports prefixes like HLA*)\n"
            "  gff_toGeneBed.py -i genes.gff3 -o hla.bed -g my_genes.txt -e\n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("-i", "--input", required=True, help="Input GFF3 file")
    ap.add_argument("-o", "--output", required=True, help="Output BED file")
    ap.add_argument("-g", "--genes", default=None,
                    help="Gene selector: file path (one per line) or comma-separated list. "
                         "Items ending with '*' are treated as prefixes (e.g., 'SMN*').")
    ap.add_argument("-e", "--exon-targets", action="store_true",
                    help="If set, add merged exon intervals as Target=... in 5th column (absolute, 0-based half-open)")
    ap.add_argument("-a", "--anchor", type=int, default=0,
                    help="Expand each gene interval by this many bp upstream and downstream (columns 2 and 3). Default: 0")
    args = ap.parse_args()

    gff_path = Path(args.input)
    out_path = Path(args.output)
    anchor = max(0, int(args.anchor or 0))

    # Build selection (exact, prefixes)
    exact_sel, prefix_sel = build_gene_selector(args.genes)

    gene_info = {}                 # gid -> dict(chrom,start0,end0,name,strand)
    transcript_to_gene = {}        # tid -> gid
    exon_by_gene = defaultdict(list)  # gid -> [(s0,e), ...]

    with gff_path.open() as fin:
        for line in fin:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, source, ftype, start1_str, end1_str, score, strand, phase, attrs_str = parts
            try:
                start1 = int(start1_str)
                end1 = int(end1_str)
            except ValueError:
                continue
            attrs = parse_attrs(attrs_str)
            start0 = start1 - 1  # GFF3: 1-based inclusive -> BED: 0-based half-open

            ftype_l = ftype.lower()
            if ftype_l == "gene":
                gid = attrs.get("ID") or attrs.get("gene_id") or attrs.get("Name")
                if not gid:
                    continue
                gname = attrs.get("gene_name") or attrs.get("Name") or attrs.get("gene") or gid
                gene_info[gid] = {
                    "chrom": chrom,
                    "start0": start0,
                    "end0": end1,   # end1 is 1-based inclusive; BED end is half-open → ok to use end1 directly
                    "name": gname,
                    "strand": strand
                }

            elif ftype_l in ("mrna", "transcript"):
                tid = attrs.get("ID")
                parent = attrs.get("Parent")
                if tid and parent:
                    transcript_to_gene[tid] = parent.split(",")[0]

            elif ftype_l == "exon":
                parents = attrs.get("Parent")
                if not parents:
                    continue
                for pid in parents.split(","):
                    gid = transcript_to_gene.get(pid)
                    if not gid:
                        # sometimes exon parents the gene directly
                        if pid in gene_info:
                            gid = pid
                    if gid:
                        exon_by_gene[gid].append((start0, end1))

    with out_path.open("w") as fout:
        for gid, info in gene_info.items():
            chrom = info["chrom"]
            gstart = info["start0"]
            gend = info["end0"]
            gname = info["name"]

            # Selection filter (exact and/or prefix*)
            if exact_sel is not None or prefix_sel is not None:
                if not matches_selection(
                    gene_keys_for_match(gid, gname),
                    exact_sel, prefix_sel
                ):
                    continue

            # Anchor expansion (symmetric)
            adj_start = gstart - anchor
            if adj_start < 0:
                adj_start = 0
            adj_end = gend + anchor  # not capped (no chrom sizes here)

            if args.exon_targets:
                merged = merge_intervals(exon_by_gene.get(gid, []))
                target_str = "Target=" + ",".join(f"{s}-{e}" for s, e in merged) if merged else "Target="
                fout.write(f"{chrom}\t{adj_start}\t{adj_end}\t{gname}\t{target_str}\n")
            else:
                fout.write(f"{chrom}\t{adj_start}\t{adj_end}\t{gname}\n")

if __name__ == "__main__":
    main()
