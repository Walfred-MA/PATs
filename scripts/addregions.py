#!/usr/bin/env python3


import argparse
import collections as cl
import pathlib
import subprocess
import sys
import time
import uuid
import os

def process_locus2(name, contig, start, end, haplopath, outputfile):
    
    if "#" not in contig: 
        haplo = "CHM13_h1" if "NC_0609" in contig else "HG38_h1"
    else:
        haplo = "_h".join(contig.split("#")[:2])
        
    path = haplopath[haplo]
    
    start = max(0,int(start))
    header = ">{}_{}_0\t{}:{}-{}".format(name, haplo, contig, start, end)
    
    cmd1 = f'echo "{header}" >> {outputfile} && samtools faidx {path} {contig}:{start+1}-{end} | tail -n +2 >> {outputfile}'
    os.system(cmd1)
    
    
def makenewfasta(regions, queryfile, outputfile):
    
    haplopath = dict()
    with open(queryfile, mode = 'r') as f:
        for line in f:
            line = line.strip().split()
            haplopath[line[0]] = line[1]
            
    os.system("rm {} || true ".format(outputfile))
    
    alloutputs = []
    index = 0
    for contig,regions in regions.items():
        for region in regions:
            index += 1
            start, end = region
            name = "seq"+str(index)
            process_locus2( name, contig, start, end, haplopath, outputfile)
            
            
def OverlapChr(refs, queries):
    
    results = [[] for x in refs]
    lr = len(refs)
    length = len(refs) + len(queries)
    lgenes = len(refs)
    
    coordinates = [x for y in refs+queries for x in y[:2]] 
    
    coordinates_sortindex = sorted( range(len(coordinates)), key = lambda x: coordinates[x])
    
    overlaps = []
    cover = 0
    for index in coordinates_sortindex:
        
        index2 = index // 2
        coordi = coordinates[index]
        
        if index % 2 == 0:
            cover += 1
            if cover == 1:
                overlaps.append([coordi,coordi])
        else:
            cover -= 1
            if cover == 0:
                overlaps[-1][1] = coordi
                
    return overlaps



def main(args):
    
    inputfile = args.input
    addfile = args.add
    outputfile  = args.out
    queryfile = args.query
    inputregions = cl.defaultdict(list)
    with open(inputfile, mode='r') as f:
        for line in f:
            if line.startswith(">"):
                name, region = line.split()[:2]
                name = name[1:]
                contig, coordi = region.split(":")
                if coordi.endswith("-") or coordi.endswith("+"):
                    coordi = coordi[:-1] 
                start, end = coordi.split('-')
                start,end = int(start), int(end)
                inputregions[contig].append([start,end])
                
    addregions = cl.defaultdict(list)
    with open(addfile, mode='r') as f:
        for line in f:
            if line.startswith(">"):
                name, region = line.split()[:2]
                name = name[1:]
                contig, coordi = region.split(":")
                if coordi.endswith("-") or coordi.endswith("+"):
                    coordi = coordi[:-1] 
                start, end = coordi.split('-')
                start,end = int(start), int(end)
                addregions[contig].append([start,end])
                
    for chrom, regions in addregions.items():
        
        inputregions[chrom] = OverlapChr(inputregions[chrom] ,regions)

    ifoverlap = 0
    for chrom, add_regs in addregions.items():
        before_len = sum(e - s for s, e in inputregions[chrom]) + sum(e - s for s,e in add_regs)
        new_regions = OverlapChr(inputregions[chrom], add_regs)
        after_len  = sum(e - s for s, e in new_regions)
        if after_len < before_len:
            ifoverlap = 1
            inputregions[chrom] = new_regions

    # --- No overlap means total coverage didn’t increase ---
    if ifoverlap:
        app_path = outputfile + "._app_"
        os.system(f"rm -f {outputfile} {app_path}")
        os.link(inputfile, outputfile)
        os.link(addfile, app_path)
        return

    makenewfasta(inputregions, queryfile, outputfile)
    
def run():
    """
        Parse arguments and run
    """
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-i", "--input", help="path to input data file",dest="input", type=str, required=True)
    parser.add_argument("-a", "--add", help="path to input data file",dest="add", type=str, required=True)
    parser.add_argument("-o", "--output", help="path to input data file",dest="out", type=str, required=True)
    parser.add_argument("-q", "--query", help="path to input data file",dest="query", type=str, required=True)
    
    parser.set_defaults(func=main)
    args = parser.parse_args()
    args.func(args)
    
    
if __name__ == "__main__":
    run()
