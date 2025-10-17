#!/usr/bin/env python3

import argparse
import collections as cl
import shutil
import os
from pathlib import Path
from multiprocessing import Pool
import re

script_folder = Path(__file__).resolve().parent

class checkforoverlap:
    def __init__(self, thereffile):  # <-- Fix: missing `self`
        self.refregions = cl.defaultdict(list)
        with open(thereffile, mode='r') as f:
            for line in f:
                if line.startswith(">"):
                    line = line.split()
                    contig, start, end, name = line[1], line[2], line[3], line[-2]
                    self.refregions[contig].append([int(start), int(end), name])
                    
    def checkfile(self, file):
        
        overlaps = []
        with open(file, mode='r') as f:
            for line in f:
                if line.startswith(">"):
                    line = line.strip().split()
                    region = line[1]
                    if region[-1] in ['-', '+']:
                        contig, coordi = region[:-1].split(":")
                    else:
                        contig, coordi = region.split(":")
                    if contig in self.refregions:
                        coordi = [int(x) for x in coordi.split('-')]
                        for region in self.refregions[contig]:
                            if max( region[0], region[1], coordi[1], coordi[0]) - min(region[0], region[1], coordi[1], coordi[0]) - (region[1] - region[0]) - (coordi[1] - coordi[0]) < -1000 :
                                
                                overlaps.append(f"{region[2]}:{contig}:{region[0]}_{region[1]}")
                                
                                
        return ";".join(overlaps)
    
def run_filter(apart_args):
    apart, similar, filtercutoff = apart_args
    cmd = (
        f"python {script_folder}/querylist_filterbykmer.py "
        f"-i {apart} -o {apart}_filtered.fa "
        f"-r {similar} -c {filtercutoff} -p 1 >> {apart}_info.txt"
    )
    os.system(cmd)
    
def nameedit(fastafile, outputfile, partname):
    
    with open(fastafile, mode = 'r') as f, open(outputfile, mode = 'w') as w:
        
        for line in f:
            
            if line.startswith(">"):
                line = line.split('\t')
                line[0] = line[0].split("_")
                line[0] = "_".join([line[0][0]+partname] + line[0][1:])
                line = "\t".join(line)
                
            w.write(line)
            
            
def main(args):
    inputfile = Path(args.input)
    outputfolder = Path(args.output)
    kmerfile = Path(args.kmer)
    nthreads = args.nthreads
    cutoff = args.cutoff
    filtercutoff = args.filter
    similar = args.similar
    wkfolder = inputfile.parent
    thereffile = args.query
    if len(thereffile) == 0:
        thereffile = wkfolder / f"{wkfolder.name}_origin.fa"
        
        
    if outputfolder.exists():
        pass
        shutil.rmtree(outputfolder)
    outputfolder.mkdir(parents=True, exist_ok=True)
    
    os.system(f"{script_folder}/kmerpartition -i {inputfile} -k {kmerfile} -o {outputfolder}/ -n {nthreads} -c {cutoff} -s {similar}")
    os.system(f"mv {outputfolder}/_filtered.fa {outputfolder}/p0.fa")
    refregion = checkforoverlap(thereffile)
    
    allparts = []
    filtered = []
    
    overlaps = cl.defaultdict(str)
    for x in os.listdir(outputfolder):
        fa_path = outputfolder / x
        kmer_path = outputfolder / f"{x}_kmer.list"
        
        if x.endswith(".fa") and fa_path.is_file() and fa_path.stat().st_size > 100 and \
            kmer_path.is_file() and kmer_path.stat().st_size > 100:
                
            overlap = refregion.checkfile(fa_path)
        
            if len(overlap):
                if re.search(r"p\d+\.fa$", x):
                    match = re.search(r"p(\d+)", x)
                    octal = oct(int(match.group(1)))[2:]
                    octal = octal
                    name = x.replace(match.group(0), f"p{octal}")[:-3]  # strip .fa
                else:
                    name = x[:-3]
                    
                
                prefix, rest = inputfile.stem.split("_", 1)
                newname = f"{prefix}{name}_samples.fasta"
                
                partfolder = outputfolder / name
                partfolder.mkdir(exist_ok=True)
                
                # Define destination paths
                new_fa_path = partfolder / newname
                new_kmer_path = partfolder / (newname  + "_kmer.list")
                
                # Edit FASTA file → write modified version into partfolder
                nameedit(fa_path, new_fa_path, name)
                
                # Move the kmer file as-is
                shutil.move(str(kmer_path), new_kmer_path)
                
                # Track overlap
                overlaps[new_fa_path] = overlap
                allparts.append(new_fa_path)
                
                # Delete the original .fa file only if you want to remove it (optional)
                fa_path.unlink()
                
            else:
                
                filtered.append(fa_path)
                if fa_path.exists():
                    fa_path.unlink()
                if kmer_path.exists():
                    kmer_path.unlink()
                    
        else:
            filtered.append(fa_path)
            if fa_path.exists():
                fa_path.unlink()
            if kmer_path.exists():
                kmer_path.unlink()
                
    with open(f"{inputfile}_partitions.list", mode='w') as f:
        for apart in allparts:
            f.write(f"{apart}\t{overlaps[apart]}\n")
            
    with open(f"{inputfile}_filtered.list", mode='w') as f:
        for apart in filtered:
            f.write(f"{apart}\tNA\n")
            
def run():
    parser = argparse.ArgumentParser(description="Program to determine pseudogene")
    parser.add_argument("-i", "--input", required=True, type=str)
    parser.add_argument("-o", "--output", type=str, default="")
    parser.add_argument("-k", "--kmer", type=str, default="")
    parser.add_argument("-n", "--nthreads", type=int, default=1)
    parser.add_argument("-c", "--cutoff", type=int, default=1000)
    parser.add_argument("-f", "--filter", type=int, default=1000)
    parser.add_argument("-s", "--similar", type=float, default=0.1)
    parser.add_argument("-q", "--query", type=str, default="")
    parser.set_defaults(func=main)
    args = parser.parse_args()
    args.func(args)
    
if __name__ == "__main__":
    run()
