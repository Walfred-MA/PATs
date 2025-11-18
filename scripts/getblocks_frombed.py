#!/usr/bin/env python3

import argparse
import subprocess
import collections as cl
from pathlib import Path

nametochr = {
	'NC_060925.1': 'chr1', 'NC_060926.1': 'chr2', 'NC_060927.1': 'chr3', 'NC_060928.1': 'chr4', 'NC_060929.1': 'chr5',
	'NC_060930.1': 'chr6', 'NC_060931.1': 'chr7', 'NC_060932.1': 'chr8', 'NC_060933.1': 'chr9', 'NC_060934.1': 'chr10',
	'NC_060935.1': 'chr11', 'NC_060936.1': 'chr12', 'NC_060937.1': 'chr13', 'NC_060938.1': 'chr14', 'NC_060939.1': 'chr15',
	'NC_060940.1': 'chr16', 'NC_060941.1': 'chr17', 'NC_060942.1': 'chr18', 'NC_060943.1': 'chr19', 'NC_060944.1': 'chr20',
	'NC_060945.1': 'chr21', 'NC_060946.1': 'chr22', 'NC_060947.1': 'chrX', 'NC_060948.1': 'chrY'
	}

def parse_target_absolute_col5(fields):
	"""ONLY parse Target= from column 5 (index 4); absolute 0-based half-open."""
	if len(fields) >= 5:
		col5 = fields[4]
		if '=' in col5:
			k, v = col5.split('=', 1)
			if k.lower() == 'target':
				out = []
				for seg in v.split(','):
					seg = seg.strip()
					if not seg:
						continue
					a, b = seg.split('-')
					a, b = int(a), int(b)
					if b > a:
						out.append((a, b))
				return out
	return []

def select_query_and_chr(chr_raw, queries, singleref):
	"""Choose reference path and normalize chr name (NC_0609xx -> chrN)."""
	if singleref:
		query = singleref
	else:
		if "#" not in chr_raw:
			haplo = "CHM13_h1" if chr_raw.startswith("NC_0609") else "HG38_h1"
		else:
			haplo = "_h".join(chr_raw.split("#")[:2])
		query = queries.get(haplo)
		if not query:
			raise ValueError(f"No reference found for haplotype key: {haplo}")
	chr_out = nametochr.get(chr_raw, chr_raw) if chr_raw.startswith("NC_0609") else chr_raw
	return query, chr_out

def faidx_get_seq(query, chr_name, start0, end0):
	"""Return raw sequence (no header) as string; empty if faidx fails."""
	cmd = f"samtools faidx {query} {chr_name}:{start0+1}-{end0}"
	p = subprocess.run(cmd, shell=True, capture_output=True, text=True)
	if p.returncode != 0:
		return ""
	lines = p.stdout.splitlines()
	return "".join(ln.strip() for ln in lines if not ln.startswith('>'))

def merge_intervals(intervals):
	if not intervals: return []
	intervals.sort()
	merged = [list(intervals[0])]
	for s, e in intervals[1:]:
		last = merged[-1]
		if s <= last[1]:
			if e > last[1]:
				last[1] = e
		else:
			merged.append([s, e])
	return [(s, e) for s, e in merged]

def clip_abs_to_rel(abs_targets, start0, end0):
	"""Convert absolute targets to relative [0, L) and clip; merge overlaps."""
	L = end0 - start0
	rel = []
	for a, b in abs_targets:
		rs, re = a - start0, b - start0
		if re <= 0 or rs >= L:
			continue
		if rs < 0: rs = 0
		if re > L: re = L
		if rs < re:
			rel.append((rs, re))
	return merge_intervals(rel)

def uppercase_ranges(seq_lower, rel_ranges):
	if not rel_ranges: return seq_lower
	arr = bytearray(seq_lower, 'ascii')
	n = len(arr)
	for rs, re in rel_ranges:
		if rs < 0: rs = 0
		if re > n: re = n
		if rs >= re: continue
		for i in range(rs, re):
			c = arr[i]
			if 97 <= c <= 122:  # a..z
				arr[i] = c - 32  # -> A..Z
	return arr.decode('ascii')

def wrap80(seq):
	return "\n".join(seq[i:i+80] for i in range(0, len(seq), 80)) + ("\n" if seq else "")

# --- NEW: .fai loader with caching ---
_fai_cache = {}  # query_path -> {chrom: length}

def load_fai_lengths(query_path):
	"""
	Load .fai for a reference FASTA. Caches per query_path.
	Returns dict {chrom: length}. Raises FileNotFoundError if .fai missing.
	"""
	if query_path in _fai_cache:
		return _fai_cache[query_path]
	q = Path(query_path)
	# Try exact .fai; if FASTA is gz, samtools still creates .fai next to it.
	fai_path = q.with_suffix(q.suffix + ".fai") if q.suffix == ".gz" else q.with_suffix(q.suffix + ".fai")
	if not fai_path.exists():
		# Fallback: append .fai directly
		alt = Path(str(q) + ".fai")
		if alt.exists():
			fai_path = alt
		else:
			raise FileNotFoundError(f"No .fai found for {query_path}")
	lengths = {}
	with open(fai_path) as f:
		for line in f:
			parts = line.strip().split("\t")
			if len(parts) >= 2:
				chrom = parts[0]
				try:
					length = int(parts[1])
				except ValueError:
					continue
				lengths[chrom] = length
	_fai_cache[query_path] = lengths
	return lengths

def clamp_to_fai(query, chr_name, start0, end0):
	"""
	Clamp [start0, end0) to chromosome length from the reference's .fai.
	Returns (ok, start0_new, end0_new). ok=False if chr not in .fai or invalid range.
	"""
	try:
		lengths = load_fai_lengths(query)
	except FileNotFoundError:
		print("Missing: ", query)
		# If no .fai, we cannot validate; treat as not ok.
		return (False, start0, end0)
	if chr_name not in lengths:
		return (False, start0, end0)
	chrom_len = lengths[chr_name]
	if start0 < 0:
		start0 = 0
	if start0 >= chrom_len:
		return (False, start0, end0)
	if end0 > chrom_len:
		end0 = chrom_len
	if end0 <= start0:
		return (False, start0, end0)
	return (True, start0, end0)

def main():
	ap = argparse.ArgumentParser(
		description=(
			"Extract regions with samtools. "
			"If -o is an existing directory, write per-name files (<out>/<name>/<name>.fa). "
			"Otherwise, write ONE multi-FASTA to -o. "
			"If column 5 is Target= (absolute, 0-based, half-open), lowercase all and uppercase only targets. "
			"Validates coordinates against the reference .fai: skips missing chroms and clamps end to chrom length."
		)
	)
	ap.add_argument("-i", "--input", required=True, help="BED-like: chrom start end name [Target=a-b,c-d,...] (Target only in column 5)")
	ap.add_argument("-r", "--ref", required=True, help="Reference FASTA or mapping file (hap_key\\tpath)")
	ap.add_argument("-o", "--output", required=True, help="Output path: existing dir => per-name files; else a single multi-FASTA")
	args = ap.parse_args()
	
	# Load reference(s)
	singleref = ""
	queries = {}
	refp = args.ref
	if refp.endswith((".fa", ".fasta", ".fna", ".fa.gz", ".fasta.gz", ".fna.gz")):
		singleref = refp
	else:
		with open(refp) as rf:
			for line in rf:
				parts = line.strip().split()
				if len(parts) >= 2:
					queries[parts[0]] = parts[1]
					
	output_path = Path(args.output)
	dir_mode = output_path.exists() and output_path.is_dir()
	name_index = cl.defaultdict(int)
	
	if dir_mode:
		# Per-name folder strategy
		with open(args.input) as fin:
			for raw in fin:
				if not raw.strip() or raw.startswith("#"):
					continue
				fields = raw.strip().split()
				if len(fields) < 4:
					raise ValueError(f"Expected >=4 cols: chrom start end name; got: {raw.strip()}")
					
				chr_raw, start0, end0, name = fields[0], int(fields[1]), int(fields[2]), fields[3]
				# Select reference and chr alias
				try:
					query, chr_name = select_query_and_chr(chr_raw, queries, singleref)
				except ValueError:
					continue  # skip if mapping missing
				
				# --- NEW: clamp/validate against .fai ---
				ok, start0_c, end0_c = clamp_to_fai(query, chr_name, start0, end0)
				if not ok:
					continue
				start0, end0 = start0_c, end0_c
				
				seq = faidx_get_seq(query, chr_name, start0, end0)
				if not seq:
					continue  # skip failed extraction
				
				abs_targets = parse_target_absolute_col5(fields)
				if abs_targets:
					rel = clip_abs_to_rel(abs_targets, start0, end0)
					seq = uppercase_ranges(seq.lower(), rel)
					
				outdir = output_path / name
				outdir.mkdir(parents=True, exist_ok=True)
				out_fa = outdir / f"{name}.fa"
				
				name_index[name] += 1
				header_text = f"{name}_{name_index[name]}\t{raw.strip()}"
				
				with open(out_fa, "a") as fout:
					fout.write(f">{header_text}\n")
					fout.write(wrap80(seq))
	else:
		# Single multi-FASTA file
		output_path.parent.mkdir(parents=True, exist_ok=True)
		with open(args.input) as fin, open(output_path, "w") as fout:
			for raw in fin:
				if not raw.strip() or raw.startswith("#"):
					continue
				fields = raw.strip().split()
				if len(fields) < 4:
					raise ValueError(f"Expected >=4 cols: chrom start end name; got: {raw.strip()}")
					
				chr_raw, start0, end0, name = fields[0], int(fields[1]), int(fields[2]), fields[3]
				try:
					query, chr_name = select_query_and_chr(chr_raw, queries, singleref)
				except ValueError:
					continue  # skip if mapping missing
				
				# --- NEW: clamp/validate against .fai ---
				ok, start0_c, end0_c = clamp_to_fai(query, chr_name, start0, end0)
				if not ok:
					continue
				start0, end0 = start0_c, end0_c
				
				seq = faidx_get_seq(query, chr_name, start0, end0)
				if not seq:
					continue  # skip failed extraction
				
				abs_targets = parse_target_absolute_col5(fields)
				if abs_targets:
					rel = clip_abs_to_rel(abs_targets, start0, end0)
					seq = uppercase_ranges(seq.lower(), rel)
					
				name_index[name] += 1
				header = f"{name}_{name_index[name]}\t{raw.strip()}"
				fout.write(f">{header}\n")
				fout.write(wrap80(seq))
				
if __name__ == "__main__":
	main()
