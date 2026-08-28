#!/usr/bin/env python3

import os
import argparse
import collections as cl
from bitarray import bitarray
import re
import sys
import numpy as np
import pyximport
pyximport.install(
    language_level=3,
    setup_args={"include_dirs": [np.get_include()]}
)
from hashmap import KmerDB,create_kmerdb  # hashmap.pyx must be in the current folder or PYTHONPATH
from matrixformat import INDEX_VERSION_LINE



mainchr = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM']

reciprocals = [200, 196, 192, 189, 185, 182, 179, 175, 172, 169, 167, 164, 161, 159, 156, 154, 152, 149, 147, 145, 143, 141, 139, 137, 135, 133, 132, 130, 128, 127, 125, 123, 122, 120, 119, 118, 116, 115, 114, 112, 111, 110, 109, 108, 106, 105, 104, 103, 102, 101, 100, 99, 98, 97, 96, 95, 94, 93, 93, 92, 91, 90, 89, 88, 88, 87, 86, 85, 85, 84, 83, 83, 82, 81, 81, 80, 79, 79, 78, 78, 77, 76, 76, 75, 75, 74, 74, 73, 72, 72, 71, 71, 70, 70, 69, 69, 68, 68, 68, 67, 67]
reciprocals_uindex = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 57, 58, 59, 60, 61, 62, 62, 63, 64, 65, 65, 66, 67, 67, 68, 69, 69, 70, 71, 71, 72, 72, 73, 74, 74, 75, 75, 76, 76, 77, 78, 78, 79, 79, 80, 80, 81, 81, 82, 82, 82, 83, 83]

anchorsize = 3000
fixcol = 15


def pack_numbers(first, second, third):
    
    result = (first << 63) | (second << 32) | third
    return result

def unpack_numbers(packed):
    third = packed & ((1 << 32) - 1)           # bits 0–31
    second = (packed >> 32) & ((1 << 31) - 1)  # bits 32–62
    first = (packed >> 63) & 1                 # bit 63
    return first, second, third

def intdecode(value):
    v = 0
    for c in value:
        v *= 64
        v += ord(c) - ord("0")
    return v

def intencode(value):
    
    code = ""
    
    code = chr( ord('0') + (value % 64) ) + code
    value = value//64
    
    while value:
        
        code = chr( ord('0') + (value % 64) ) + code
        value = value//64
        
    return code

def intencode_single(value, start = '0'):
    
    code = chr( ord(start) + value ) 
    
    return code

def kmerencode(kmer):
    
    kmerint = 0
    for base in kmer:
        
        kmerint *= 4
        if base=='a' or base =='A':
            
            kmerint += 0
            
        elif base=='c' or base =='C':
            
            kmerint += 1
            
        elif base=='g' or base =='G':
            
            kmerint += 2
        elif base=='t' or base =='T':
            
            kmerint += 3
            
    return kmerint

def kmerdecompress(kmerstr, size = 31):
    
    kmerint = 0
    for char in kmerstr:
        
        kmerint *= 64
        kmerint += ord(char) - ord('0')
        
    return kmerint

def parse_encoded(text):
    def intdecode(value):
        v = 0
        for c in value:
            v *= 64
            v += ord(c) - ord("0")
        return v
    def eachelement(text):
        if "*" in text:
            d = intdecode(text.split("*")[1])
            text = text.split("*")[0]
        else:
            d = 1
            
            
        if "~" in text:
            start, size = text.split("~")
            start = intdecode(start)
            size = intdecode(size)
            c = [x for y in range(start, start + size + 1) for x in [y] * d]
        else:
            c = [intdecode(text)] * d
            
        return c
    
    
    return [x for y in text.split(",")[:-1] for x in eachelement(y)]


def getmax(text):
    
    if "*" in text:
        text = text.split("*")[0]
        
    if "~" in text:
        c = intdecode(text.split('~')[1]) + intdecode(text.split('~')[0])
    else:
        c = intdecode(text)
        
    return c

def adderrorcorr(inputfile, errorfile ):
    
    kmdatabase = create_kmerdb()
    kmdatabase.load(errorfile)
    biasbins = [0]*101
    meanbins =  [0]*101
    varbins = [0]*101
    with open(inputfile, mode = 'r') as r, open(inputfile+"_", mode = 'w') as w:
        
        line = r.readline()
        while line:
            if len(line):
                if line[0] == "&" or line[0] == '^':
                    
                    line = line.split('\t')
                    
                    """
                    if len(line[2]) < 3:
                        line[2] = "0"*(3-len(line[2])) + line[2]
                        
                    if len(line[3]) < 5:
                        line[3] = "0"*(5-len(line[3])) + line[3]

                    if len(line[4]) < 11:
                        line[4] = "0"*(11-len(line[4])) + line[4]
                    """
                    
                    kmerint = kmerdecompress(line[4])
                    
                    packed = kmdatabase.find(kmerint)
                    
                    if packed is not None:
                        valid, biasflag, kmerratio, kmervar = packed
                    else:
                        valid, biasflag, kmerratio, kmervar = 0, 0, 50,100
                        
                    biasbins[biasflag] += 1
                    meanbins[kmerratio] += 1
                    varbins[kmervar] += 1
                    line[1] = intencode_single(biasflag,chr(33))+intencode_single(kmervar - 10,chr(33))+intencode_single(reciprocals_uindex[kmerratio ] , chr(33))
                    line = "\t".join(line)
                w.write(line)
            line = r.readline()
            
    os.system("mv {} {}".format(inputfile+"_", inputfile))
    
    print("Frequecy summary of values:\n")
    for i in range(101):
        print(f"{i}\t{biasbins[i]}\t{meanbins[i]}\t{varbins[i]}")
        
        
        
def getdupsindex(inputfile):
    kmer_array = np.empty(2**32 - 1, dtype=np.uint64)  # ~32GB
    lindex = 0
    
    with open(inputfile, "r") as r:
        for line in r:
            if line and line[0] in "&^*.":
                fields = line.rstrip("\n").split("\t")
                if len(fields) > 4:
                    kmer_array[lindex] = kmerdecompress(fields[4])
                    lindex += 1
                    
    ifdup = bitarray(lindex)
    ifdup.setall(0)
    
    order = np.argsort(kmer_array[:lindex], kind="mergesort")
    lastkmer = np.uint64(2**64 - 1)
    
    for idx in order:
        kmer = kmer_array[idx]
        if kmer == lastkmer:
            ifdup[idx] = 1
        lastkmer = kmer
        
    return ifdup

            
def clearreds(inputfile):
    
    ifdup = getdupsindex(inputfile)
    
    tmp = inputfile + "_"
    lindex = 0
    
    with open(inputfile, "r") as r, open(tmp, "w") as w:
        for line in r:
            if line and line[0] in "&^*.":
                fields = line.split("\t")
                
                if ifdup[lindex]:
                    if fields[0].startswith("&"):
                        fields[0] = "^" + fields[0][1:]
                    elif fields[0].startswith("*"):
                        fields[0] = "." + fields[0][1:]
                        
                w.write("\t".join(fields))
                lindex += 1
            else:
                w.write(line)
                
    os.replace(tmp, inputfile)
    
    
def matrixindex(inputfile, genedbfile):
    
    if len(genedbfile):
        genecode = GenodeDB(genedbfile)
        
        
    genes = []
    allrefregions = []
    regions = cl.defaultdict(list)
    refregions = ""
    
    names = [] 
    locations = []
    gene_nums = []
    kmer_nums = []
    eles_nums = []
    extd_nums = []
    extd_kmers = []
    passedkmers = set()
    with open(inputfile, mode = 'r') as r:
        line = r.readline()
        kmer_num = 0
        gene_num = 0
        eles_num = 0
        extd_num = 0
        extd_num = 0
        extd_kmer = 0
        lindex = 0 
        while line:
            if len(line):
                if line.startswith('#'):
                    locations.append(r.tell() - len(line))
                    gene_nums.append(gene_num)
                    kmer_nums.append(kmer_num)
                    eles_nums.append(eles_num)
                    names.append(line.split()[0]) 
                    extd_nums.append(extd_num)
                    extd_kmers.append(extd_kmer)
                    gene_num = 0
                    kmer_num = 0
                    extd_kmer = 0
                    eles_num = 0
                    extd_num = 0
                    if len(genedbfile):
                        genes.append(genecode.Overlap(regions))
                    else:
                        genes.append("")
                        
                    allrefregions.append("".join(refregions))
                    refregions = []
                    regions = cl.defaultdict(list)
                
                if line.startswith(('&', '^', '*', '.')):
                    
                    if line.startswith(('&','^')):
                        kmer_num += 1
                    else:
                        extd_kmer += 1
                    line = line.split('\t')
                    
                    ranges = 0
                    if line[0][1] in ['-','+']:
                        eles = line[5].split(" ")
                        eles = [parse_encoded(x) for x in  eles]
                        ranges = len(eles[0])
                        extd1 = len(eles[1]) if len(eles) > 1 else 0
                        extd1 += len(eles[2] ) if len(eles) > 2 else 0
                    eles_num += ranges+ fixcol
                    extd_num += (extd1 + 2)
                elif line.startswith( '>'):
                    gene_num += 1
                    for x in line.split(";"):
                        region = x.split('\t')[1]
                        contig,start,end = addcoordi(region, regions)
                        if "#" not in contig:
                            refregions.append(f"{contig}:{max(1,start-3000)}-{end+3000}{region[-1]};")
            line = r.readline()           
            lindex += 1
            
            
        locations.append(r.tell())
        gene_nums.append(gene_num)
        kmer_nums.append(kmer_num)
        eles_nums.append(eles_num)
        extd_nums.append(extd_num)
        extd_kmers.append(extd_kmer)
        if len(genedbfile):
            genes.append(genecode.Overlap(regions))
        else:
            genes.append("")
        allrefregions.append("".join(refregions))
        
    with open(inputfile+".index", mode = 'w') as f:

        f.write(INDEX_VERSION_LINE + "\n")
        
        last_location = 0
        
        for name, location,kmer_num, eles_num, extd_num, extd_kmer, genum, gene, refregion in zip(names, locations[1:], kmer_nums[1:], eles_nums[1:], extd_nums[1:],extd_kmers[1:], gene_nums[1:], genes[1:], allrefregions[1:]):
            
            if name.startswith("#HLA"):
                refregion += "HLA;"
                
            f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(name, last_location, location, kmer_num, eles_num, extd_num, extd_kmer, genum,refregion,"".join([ x+";" for x in gene if x.startswith(name[1:].split("_")[0]) or name.startswith("#group")]) ))
            last_location = location
            
    return names, locations, kmer_nums, eles_nums, allrefregions



class GenodeDB:
    
    def __init__(self, path):
        
        self.pattern = re.compile(r'(?:gene_name=)([^;]+)')
        paths = [x for x in path.split(",") if x]
        self.gene = cl.defaultdict(list)
        
        for path in paths:
            self.load_gff(path)

    def load_gff(self, path):

        with open(path, mode = 'r') as r:
            
            for line in r:
                if len(line) == 0 or line[0]=="#":
                    continue
                
                row = line.split("\t")
                
                
                if row[2] != "gene":
                    continue
                
                gene_name = re.findall( self.pattern ,row[8])
                
                if len(gene_name) == 1:
                    
                    geneinfo = [int(row[3])-anchorsize, int(row[4])+anchorsize] + gene_name
                    
                    self.gene[row[0]].append(geneinfo)
                    
    def OverlapChr(self, chr , queries):
        
        results = [[] for x in queries]
        
        refs = self.gene[chr]
        
        length = len(refs) + len(queries)
        lgenes = len(refs)
        
        coordinates = [x for y in refs + queries for x in y[:2]] 
        names = [x[2] for x in refs]
        
        coordinates_sortindex = sorted( range(length * 2), key = lambda x: coordinates[x])
        
        current_refs = []
        current_queries = []
        for index in coordinates_sortindex:
            
            index2 = index // 2
            
            if index % 2 == 0:
                
                if index2 < lgenes:
                    
                    current_refs.append(index2)
                    
                    for query in current_queries:
                        
                        results[query].append(index2)
                        
                else:
                    
                    current_queries.append(index2 - lgenes)
                    
                    results[index2 - lgenes].extend(current_refs)
                    
            else:
                
                if index2 < lgenes:
                    
                    current_refs.remove(index2)
                    
                else:
                    
                    current_queries.remove(index2 - lgenes)
                    
                    
        outputs = cl.defaultdict(list)
        for i,result in enumerate(results):
            
            queryname = queries[i][2]
            
            outputs[queryname].extend([refs[refindex][2] for refindex in result])
            
        return outputs
    
    def Overlap(self, locations):
        
        results = set()
        for chr, locations_bychr in locations.items():
            
            result_chr = self.OverlapChr(chr, locations_bychr)
            
            for queryname,genes in result_chr.items():
                
                results.update(genes)
                
        return sorted(results)
    
def addcoordi(text,regions):
    
    contig, coordi = text.split(":")
    if coordi[-1] in ['+','-','\r']:
        coordi = coordi[:-1]
    start, end = coordi.split('-')
    
    regions[contig].append([int(start), int(end), contig])
    return contig, int(start), int(end)


def main(args):
    
    
    if args.clear:
        clearreds(args.input)
        
    if len(args.error):
        adderrorcorr(args.input,args.error )
        
    matrixindex(args.input, args.gene)
    
    
def run():
    """
        Parse arguments and run
    """
    parser = argparse.ArgumentParser(description="program run partition")
    parser.add_argument("-i", "--input", help="path to input data file", dest="input", type=str, required=True)
    parser.add_argument("-e", "--error", help="path to input data file", dest="error", type=str, default= "")
    parser.add_argument("-r", "--clearredun", help="path to input data file", dest="clear", type=int, default= 1)
    parser.add_argument("-g", "--gene", help="path to input data file", dest="gene", type=str, default="")
    parser.set_defaults(func=main)
    args = parser.parse_args()
    args.func(args)
    
    
if __name__ == "__main__":
    run()
