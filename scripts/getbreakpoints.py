#!/usr/bin/env python3

import pandas as pd
import os
import argparse
import sys
import collections as cl

def readfile(querypath):
    
    read = ""
    header = ""
    query = {}
    with open(querypath, mode = 'r') as f:
        for line in f:
            if len(line) == 0:
                continue
            if line[0]==">":
                query[header] = read
                read = ""
                header = line[1:].split()[0]
            else:
                read+=line.strip()
                
    query[header] = read
    read = ""

    return query

def fixovermerge(aligns, exonposi, cutoff = 10000, samegene_cutoff = 20000):
    
    all_coordinates = [x for y in exonposi for x in y[:2]]
    all_names = [y[2].split("_")[-1] for y in aligns]
    
    
    sort_index = sorted(range(len(all_coordinates)), key = lambda x: all_coordinates[x])
    
    gene_ends = cl.defaultdict(lambda: -samegene_cutoff - 1)
    for index in sort_index:
        
        coordinate = all_coordinates[index]
        genename = all_names[index//2]
        gene_ends[genename] = coordinate
        
    lastgene_coordinate = cl.defaultdict(lambda: -samegene_cutoff - 1)
    lastgene_coordinate[""] = -samegene_cutoff - 1
    last_coordinate = -cutoff - 1
    
    numchunks = 0
    allmidpoints = [ ]
    
    for index in sort_index:
        
        coordinate = all_coordinates[index]
        genename = all_names[index//2]
        
        if index %2 == 0 :
            
            
            if numchunks == 0 and coordinate - last_coordinate> cutoff and coordinate - max(list(lastgene_coordinate.values())) > samegene_cutoff:
                
                midpoint =  ( coordinate + last_coordinate ) // 2 
                
                allmidpoints.append( midpoint )
                
            numchunks += 1
            
        else:
            
            last_coordinate = coordinate
            
            if coordinate >= gene_ends[genename]:
                lastgene_coordinate[genename] = -samegene_cutoff - 1
            else:
                lastgene_coordinate[genename] = coordinate 
                
            numchunks -= 1
            
            
    all_coordinates = [x for y in aligns for x in y[:2]] + allmidpoints
    
    lintron = len(aligns) * 2
    
    sort_index = sorted(range(len(all_coordinates)), key = lambda x: all_coordinates[x])
    
    allchunks = [ set([]) for x in allmidpoints]
    
    
    chunkindex = -1
    for index in sort_index:
        
        if index >= lintron :
            chunkindex += 1 
            
        else:
            allchunks[chunkindex].add( index//2 )
            
    return [list(x) for x in allchunks], allmidpoints


def querycontig_overlap(allaligns, allow_gap = 20000):
    
    all_coordinates = [x for y in allaligns for x in y[:2]]
    
    sort_index = sorted(range(len(all_coordinates)), key = lambda x: all_coordinates[x])
    
    allgroups = []
    current_group = set([])
    current_group_size = 0 
    current_group_start = 0
    allgroups_sizes = []
    
    number_curr_seq = 0
    coordinate = 0
    
    last_coordinate = -allow_gap-1
    for index in sort_index:
        
        coordinate = all_coordinates[index]
        if index %2 == 0 :
            
            number_curr_seq += 1
            
            if number_curr_seq == 1:
                
                current_group_start = coordinate
                
                if coordinate - last_coordinate> allow_gap :
                    
                    allgroups_sizes.append(current_group_size)
                    allgroups.append(current_group)
                    
                    current_group = set([])
                    current_group_size = 0
                    
            current_group.add(index//2)
            
        else:
            
            number_curr_seq -= 1
            
        last_coordinate = coordinate
        
    allgroups.append(current_group)
    
    allgroups = allgroups[1:]
    
    return allgroups

def connectblocks(allgroups, scores):
    
    allgroups = sorted(allgroups,key = lambda x: x[2])
    
    merged = cl.defaultdict(list)
    gene_past = dict()
    
    block_name = [group[0] for index, group in enumerate(allgroups)]
    
    for index, group in enumerate(allgroups):
        
        score =  group[1]
        
        start, end = group[2], group[3]
        size = end - start
        
        if group[0] in gene_past and start - gene_past[group[0]][2] <= max(10000, gene_past[group[0]][1]):
            gene_past[group[0]][2] = end
            gene_past[group[0]][1] += size - (start - gene_past[group[0]][2])
            
            merged[group[0]][-1][1] = index
            merged[group[0]][-1][2] += size - (start - gene_past[group[0]][2])
            for i in range(merged[group[0]][-1][0],index):
                block_name[i] = group[0]
                
        else:
            
            gene_past[group[0]] = [score, size, end]
            merged[group[0]].append([index,index,score])
            
        if size > 10000:
            gene_past = {name:info for name, info in gene_past.items() if name == group[0]}
        else:
            gene_past = {name:info for name, info in gene_past.items() if info[1] >= size} 
    
    startindex = 0
    lastname = ""
    newgroups = []
    for genename, merges in merged.items():
        
        for merge in merges :
            if len(merge):
                newgroups.append([genename,allgroups[merge[0]][2], allgroups[merge[1]][3], allgroups[merge[1]][3] - allgroups[merge[0]][2], merge[2] ])
            
            startindex = index
                    
        
    newgroups =sorted(newgroups, key = lambda x:x[1])
    
    return newgroups


def querycontig_overlap2(allaligns):
    
    all_coordinates = [x for y in allaligns for x in y[:2]]
    
    scores = [6*x[4] - 5*x[5] for x in allaligns]
    
    sort_index = sorted(range(len(all_coordinates)), key = lambda x: all_coordinates[x])
    
    allgroups = []
    current_group = []
    current_highest = -1
    current_start = 0
        
    coordinate = 0
    
    last_coordinate = 0
    for index in sort_index:
        
        coordinate = all_coordinates[index]
        score = scores[index//2]
        if index %2 == 0 :
            
            current_group.append((score,index//2))
            current_group = sorted(current_group, reverse = 1)
            
            if current_highest < 0 or score > scores[current_highest]:
                
                if current_highest >= 0:
                    allgroups.append([current_highest, current_start, coordinate])
                
                current_start = coordinate
                current_highest = index//2
            
        else:
            
            current_group.remove((score,index//2))
            if index//2 == current_highest:
                
                allgroups.append([current_highest, current_start, coordinate])
                current_highest = current_group[0][1] if len(current_group) else -1
                current_start = coordinate

    
    allgroups = [[allaligns[x[0]][3], scores[x[0]]]+x[1:]+[x[2]-x[1]] for x in allgroups if x[2] - x[1] > 100]
    
    allgroups = connectblocks(allgroups, scores)
    
    return allgroups

def querycontig_overlap3(allaligns,blocks, allow_gap = 0):
    all_coordinates = [x for y in blocks+allaligns for x in y[:2]]
    all_names = [x[3] for x in allaligns]
    
    sort_index = sorted(range(len(all_coordinates)), key = lambda x: all_coordinates[x])
    

    lblocks = len(blocks)
    allgroups = [set([]) for x in blocks]
    
    current_aligns = []
    current_group = -1
    number_curr_seq = 0
    coordinate = 0
    
    last_coordinate = -allow_gap-1
    for index in sort_index:
                
        if index < 2*lblocks:
            
            if index % 2 == 0:
                
                number_curr_seq += 1
                                
                current_group = index//2 
                allgroups[current_group].update(current_aligns)
                

            else:
            
                number_curr_seq -= 1
                
        else:
            if index % 2 == 0:
                current_aligns.append(index//2)
            
                if number_curr_seq:
                    allgroups[current_group].add(index//2)
            else:
                current_aligns.remove(index//2)
    
    allgroups = [[all_names[y-lblocks] for y in x] for x in allgroups]
 
    return allgroups


def main(args):
    
    inputpath = args.input

    refs = args.refs.split(",")
    
    allrefs = dict()
    for ref in refs:
        
        refseqs = readfile(ref)
        allrefs.update(refseqs)

    
    table = pd.read_csv(inputpath, header = None , sep = "\t")
    
    eachcontig = cl.defaultdict(list)
    
    for index, row in table.iterrows():
        
        contig = row.iat[0]
        
        gene = row.iat[5]
        gstart,gend = int(row.iat[7]), int(row.iat[8])
        similar = float(row.iat[12])       
        length = row.iat[6]
        totalsize = gend - gstart

        if not gene.startswith("Trans_"):
            unmasksize = len([a for a in allrefs[gene][gstart:gend] if a.isupper()])

        if not  gene.startswith("Trans_") and (similar < 100*args.similar or unmasksize < 1000 ):

            continue
 
        eachcontig[(contig)].append(index)
        
    allqueryloci = []
    
    for (contig), indexes in eachcontig.items():
        
        subdata = table.iloc[indexes]
        
        maxgenesize = max(list(subdata[6]))
        mingenesize = min(list(subdata[6]))
        
        usecols = subdata[[2,3,4,5,9,10]].values.tolist()
        
        usecols_exon = [x for x in usecols if x[3].startswith("Trans_")]
        
        usecols_gene = [x for x in usecols if not x[3].startswith("Trans_")]
               
        group_gene = querycontig_overlap2(usecols_gene)
        group_exon = querycontig_overlap3(usecols_exon, [x[1:3] for x in group_gene], 0)
       
        left = []
        loci = []
        for geneinfo, exoninfo in zip(group_gene, group_exon):

            exon = "Exon" if len(exoninfo) else "Intron"
            if len(exoninfo) and geneinfo[-1] > 3000 :
                loci.append(geneinfo+[exon])
            elif geneinfo[-1] > 5000:
                loci.append(geneinfo+[exon])
               

        loci = [x[:1]+[contig]+x[1:] for x in loci]

        allqueryloci.extend(loci)
        


    allqueryloci = sorted(allqueryloci)
    
    with open( args.output, mode = "w") as f:
    
        for locus in allqueryloci:
        
            f.write("\t".join([str(x) for x in locus])+"\n")
        
        
        
        
        
        
def run():
    """
        Parse arguments and run
    """
    parser = argparse.ArgumentParser(description="program distract overlap genes and alignments on contigs")
    parser.add_argument("-i", "--input", help="path to input data file",dest="input", type=str, required=True)
    parser.add_argument("-r", "--ref", help="path to input data file",dest="refs", type=str, required=True)
    parser.add_argument("-o", "--output", help="path to input data file",dest="output", type=str, required=True)
    parser.add_argument("-s", "--similar", help="path to input data file",dest="similar", type=float, default = 0.9)
    parser.set_defaults(func=main)
    args = parser.parse_args()
    args.func(args)
    
    
if __name__ == "__main__":
    run()
