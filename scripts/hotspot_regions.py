#!/usr/bin/env python3

"""
Integrated hotspot pipeline.

Steps per block directory
-------------------------
1. Build hotspot FASTA from   <block>/hotspots/<sample>_hotspot.txt
2. Run BLASTN (genome + exon DB) → PAF → filtered PAF
3. Generate breakpoint table   via getbreakpoints.py
4. Append extracted regions to <block>/<block>_samples.fasta
    (thread-safe with file locking)

Use --touch-align / --touch-region for Snakemake sentinels.
"""

from __future__ import annotations

import argparse
import collections as cl
import pathlib
import subprocess
import sys
import time
import uuid
from typing import Dict, List, Tuple
import os

script_dir = pathlib.Path(__file__).resolve().parent
# --------------------------------------------------------------------------- #
# Tunables
# --------------------------------------------------------------------------- #
BLAST_IDENTITY = 90
BLAST_THREADS  = 4
OVERLAP_GAP    = 20_000   # bp allowed between merged loci
SENTINEL_EMPTY = 5        # file size threshold in bytes

# --------------------------------------------------------------------------- #
# Generic helpers
# --------------------------------------------------------------------------- #
def sh(cmd: str, *, check: bool = False) -> None:
    """Run a shell command with I/O passthrough."""
    subprocess.run(cmd, shell=True, check=check)
    
    
def read_fasta(path: pathlib.Path) -> Dict[str, str]:
    """Load a single-line FASTA into {header: sequence}."""
    seqs: Dict[str, str] = {}
    with path.open() as fh:
        head, buf = "", []
        for ln in fh:
            ln = ln.rstrip()
            if ln.startswith(">"):
                if head:
                    seqs[head] = "".join(buf)
                head, buf = ln[1:].split()[0], []
            else:
                buf.append(ln)
        if head:
            seqs[head] = "".join(buf)
    return seqs

def merge_intervals(spans: List[Tuple[int, int]], gap: int = OVERLAP_GAP) -> List[Tuple[int, int]]:
    """Greedy merge of (start,end) pairs allowing ≤ *gap* bp separation."""
    if not spans:
        return []
    
    spans.sort()
    merged = [list(spans[0])]
    for s, e in spans[1:]:
        last = merged[-1]
        if s - last[1] <= gap:
            last[1] = max(last[1], e)
        else:
            merged.append([s, e])
            
    return [(max(0, s - gap // 2), e + gap // 2) for s, e in merged]


def run_blast(fasta: pathlib.Path, db: pathlib.Path) -> None:
    """Run BLASTN against genome DB and exon DB, appending to *.out."""
    out = f"{fasta}_blast.out"
    cmd = f"runblastn -query {fasta} -db {db}.fa_db -perc_identity {BLAST_IDENTITY}  -outfmt 17 -num_threads {BLAST_THREADS} -out {out}"    
    os.system(
        f"bash {script_dir}/runblastn -query {fasta} -db {db}.fa_db "
        f"-perc_identity {BLAST_IDENTITY}  "
        f"-outfmt 17 -num_threads {BLAST_THREADS} -out - > {out}"
    )
    os.system(
        f"bash {script_dir}/runblastn -query {fasta} -db {db}_exons.fa_db "
        f"-perc_identity {BLAST_IDENTITY} -outfmt 17 -num_threads {BLAST_THREADS} "
        f"-out - >> {out}"
    )
    
# --------------------------------------------------------------------------- #
# Locked-file helper
# --------------------------------------------------------------------------- #
    
class LockedFile:
    """Robust file lock using exclusive creation of a shared .lock file."""
    def __init__(self, filepath: pathlib.Path, mode: str = "a", max_age: int = 60):
        self.filepath = filepath
        self.mode     = mode
        self.max_age  = max_age
        self.lockfile = filepath.with_suffix(".lock")
        self._fh      = None
        
    def _cleanup(self) -> None:
        try:
            if self.lockfile.exists() and (time.time() - self.lockfile.stat().st_mtime > self.max_age):
                self.lockfile.unlink()
        except FileNotFoundError:
            pass
            
    def __enter__(self):
        self._cleanup()
        while True:
            try:
                # atomic creation
                fd = os.open(self.lockfile, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
                with os.fdopen(fd, 'w') as f:
                    f.write(str(time.time()))
                break
            except FileExistsError:
                time.sleep(0.5)
        if self.mode == "r+" and not self.filepath.exists():
            self.filepath.touch()
        self._fh = self.filepath.open(self.mode)
        return self._fh
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        if self._fh:
            self._fh.close()
        try:
            self.lockfile.unlink()
        except FileNotFoundError:
            pass
            
            
            
# --------------------------------------------------------------------------- #
# Pipeline helpers
# --------------------------------------------------------------------------- #
def build_hotspot_fasta(bed: pathlib.Path, fasta_out: pathlib.Path, query: Dict[str, str], extendsize: int) -> None:
    """Convert <sample>_hotspot.txt → FASTA of merged loci."""
    if not bed.exists():
        return
    
    groups: Dict[str, List[Tuple[int, int]]] = cl.defaultdict(list)
    
    with bed.open("r") as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split('\t')
            contig = parts[-3].split()[0]

            if contig not in query or 'chrM' in contig:
                continue

            beg = max(0, int(parts[-2]) - extendsize)
            end = min(len(query[contig]) , int(parts[-1])+extendsize)
            groups[contig].append((beg, end))
            
            
    with fasta_out.open("w") as out:
        for contig, spans in groups.items():
            if contig not in query or 'chrM' in contig:
                continue
            for beg, end in merge_intervals(spans):
                out.write(f">{contig}_{beg}_{end}\t{contig}:{beg}-{end}\n")
                out.write(query[contig][beg:end] + "\n")
                
                
def run_breakpoints_and_append(
    txt_out: pathlib.Path,
    fasta_fn: pathlib.Path,
    title_fn: pathlib.Path,
    sample: str,
    query: Dict[str, str],
    ) -> None:
    """Clean (if needed) and append breakpoint regions into *fasta_fn* (thread-safe)."""
    if not txt_out.exists() or txt_out.stat().st_size < SENTINEL_EMPTY:
        return
    
    with LockedFile(fasta_fn, "r+") as dst:
        # Step 1: read all existing lines
        try:
            lines = dst.readlines()
        except Exception as e:
            print(f"Error reading {fasta_fn}: {e}")
            return
        
        # Step 2: check if any header lines match this sample
        needs_cleaning = any(
            line.startswith(">") and sample in line.split('\t')[0]
            for line in lines
        )
        
        if needs_cleaning:
            # Step 3: clean and rewrite filtered lines
            dst.seek(0)
            dst.truncate(0)
            keep = True
            for line in lines:
                if line.startswith(">"):
                    keep = sample not in line.split('\t')[0]
                if keep:
                    dst.write(line)
        else:
            # If no cleanup needed, seek to end for appending
            dst.seek(0, os.SEEK_END)
            
        # Step 4: append new breakpoint regions
        prefix = fasta_fn.stem.split("_")[0]
        with txt_out.open() as src:
            counter = 1
            for ln in src:
                parts = ln.split()
                if len(parts) < 4:
                    continue  # malformed line
                title, chrom, beg, end = parts[:4]
                exon  = parts[6] if len(parts) > 6 else "."
                header = f"{prefix}_{sample}_{counter}\t{chrom}:{beg}-{end}\t{exon}\t{title}"
                #seq = query[chrom][(int(beg)):int(end)]
                #dst.write(f">{header}\n{seq}\n")
                dst.write(f">{header}\n")
                counter += 1                
                

def run_merge_breakpoints_and_append(
    txt_out: pathlib.Path,
    fasta_fn: pathlib.Path,
    title_fn: pathlib.Path,  # unused, kept for signature consistency
    sample: str,
    query: Dict[str, str],
    max_gap: int = 30_000
    ) -> None:
    """
    Merge breakpoint intervals in-memory and append merged sequences to FASTA (thread-safe).
    
    FASTA header format:
        >{prefix}_{sample}_merge{i}   chrom:start-end   Exon/Intron   title1;title2;title3
    """
    if not txt_out.exists() or txt_out.stat().st_size < SENTINEL_EMPTY:
        return
    
    try:
        # Step 1: parse lines into per-chromosome groups
        lines_by_chrom: Dict[str, List[Tuple[int, int, str, str]]] = cl.defaultdict(list)
        with txt_out.open() as f:
            for ln in f:
                parts = ln.strip().split()
                if len(parts) < 4:
                    continue
                title, chrom, beg, end = parts[:4]
                exon = parts[6] if len(parts) > 6 else "."
                beg, end = int(beg), int(end)
                if chrom in query:
                    lines_by_chrom[chrom].append((beg, end, exon, title))
                    
        if not lines_by_chrom:
            return
        
        # Step 2: merge logic
        merged_records = []  # (chrom, beg, end, Exon/Intron, [titles])
        for chrom, spans in lines_by_chrom.items():
            spans.sort()
            group = [spans[0]]
            for entry in spans[1:]:
                if entry[0] - group[-1][1] <= max_gap:
                    group.append(entry)
                else:
                    beg = max(0, group[0][0])
                    end = group[-1][1]
                    exon_status = "Exon" if any(x[2] != "." for x in group) else "Intron"
                    titles = ";".join([x[3] for x in group])+";"
                    merged_records.append((chrom, beg, end, exon_status, titles))
                    group = [entry]
            beg = max(0, group[0][0])
            end = group[-1][1]
            exon_status = "Exon" if any(x[2] != "." for x in group) else "Intron"
            titles = ";".join([x[3] for x in group])+";"
            merged_records.append((chrom, beg, end, exon_status, titles))
            
        # Step 3: safe FASTA append with cleaning
        with LockedFile(fasta_fn, "r+") as dst:
            try:
                lines = dst.readlines()
            except Exception as e:
                print(f"Error reading {fasta_fn}: {e}")
                return
            
            # Check for previous merged entries of this sample
            needs_cleaning = any(
                line.startswith(">") and sample in line.split('\t')[0]
                for line in lines
            )
            
            if needs_cleaning:
                dst.seek(0)
                dst.truncate(0)
                keep = True
                for line in lines:
                    if line.startswith(">"):
                        keep = sample not in line.split('\t')[0]
                    if keep:
                        dst.write(line)
            else:
                dst.seek(0, os.SEEK_END)
                
            # Write new merged entries
            prefix = fasta_fn.stem.split("_")[0]
            for i, (chrom, beg, end, exon_status, titles) in enumerate(merged_records, start=1):
                header = f"{prefix}_{sample}_{i}\t{chrom}:{beg}-{end}\t{exon_status}\t{titles}"
                seq = query[chrom][beg:end]
                dst.write(f">{header}\n{seq}\n")
                
    except Exception as e:
        print(f"Error in run_merge_breakpoints_and_append for {txt_out}: {e}")


def process_block(
    block: pathlib.Path,
    sample: str,
    query: Dict[str, str],
    undo: bool,
    script_dir: pathlib.Path,
    ) -> None:
    print("running: "+ block.name)
    
    hotspot_txt = block / "hotspots" / f"{sample}_hotspot.txt"
    hotspot_fa  = hotspot_txt.with_suffix(".txt.fa")
    db_prefix   = block / block.name
    
    paf        = hotspot_fa.with_suffix(".fa.paf")
    paf_filt   = hotspot_fa.with_suffix(".fa.paf_filter.PAF")
    title_txt  = hotspot_fa.with_suffix(".fa_title.txt")
    blast_out  = hotspot_fa.with_suffix(".fa_blast.out")
    txt_out    = paf_filt.with_suffix(".PAF.txt")
    query_fn   = block / f"{block.name}_origin.fa"
    fasta_title = txt_out.with_suffix(".txt.titles")
    # 1  Build FASTA
    
    # OPTIONAL: skip empty blocks to save time (remove if you prefer)
    # if not hotspot_fa.exists() or hotspot_fa.stat().st_size < SENTINEL_EMPTY:
    #     return
    # 2  BLAST → PAF → filter
    if undo or not (blast_out.exists() and blast_out.stat().st_size > SENTINEL_EMPTY):
        extendsize = query_fn.stat().st_size
        build_hotspot_fasta(hotspot_txt, hotspot_fa, query, extendsize)
        run_blast(hotspot_fa, db_prefix)
        
        
    os.system(f"grep '^>' {hotspot_fa} > {title_txt} || true")
    os.system(f"python {script_dir}/blastntopaf.py -i {blast_out} -q {title_txt} > {paf}")
    os.system(f"python {script_dir}/filter_alignment.py -i {paf} -r {db_prefix}.fa -o {paf_filt}")
    
    # 3  Breakpoints
    if paf_filt.exists() and paf_filt.stat().st_size >= SENTINEL_EMPTY:
        os.system(f"python {script_dir}/getbreakpoints.py -r {query_fn} -i {paf_filt} -o {txt_out}")
        
    # 4  Append regions to per-block FASTA
    breaks_fa = block / f"{block.name}_breaks.list"
    run_breakpoints_and_append(txt_out, breaks_fa, fasta_title, sample, query)

    samples_fa = block / f"{block.name}_samples.fasta"
    run_merge_breakpoints_and_append(txt_out, samples_fa, fasta_title, sample, query)

    hotspot_fa.unlink(missing_ok=True)   # tidy
    
    
# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def cli() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Hotspot → BLAST → region-FASTA pipeline")
    p.add_argument("-i", "--input",   required=True, help="Root folder with block dirs")
    p.add_argument("-q", "--query",   required=True, help="Reference FASTA for this sample")
    p.add_argument("-s", "--sample",  required=True, help="Sample prefix (matches hotspot files)")
    p.add_argument("-u", "--undo", action="store_true",help="Force re-BLAST even if previous results exist")
    p.add_argument("--touch-align",  help="Sentinel file to touch after align phase")
    p.add_argument("--touch-region", help="Sentinel file to touch after region phase")
    return p.parse_args()


def main() -> None:
    args       = cli()
    root       = pathlib.Path(args.input).expanduser().resolve()
    script_dir = pathlib.Path(__file__).resolve().parent
    query_seqs = read_fasta(pathlib.Path(args.query))
    with open(root) as f:
        for line in f:
            blk = pathlib.Path(line.strip())
            if blk.is_dir():
                process_block(blk, args.sample, query_seqs, args.undo, script_dir)
            
            
            
if __name__ == "__main__":
    main() 
