#!/usr/bin/env python3

import os
import argparse
import collections as cl
from bitarray import bitarray

def makereverse(sequence):
    
    tran=str.maketrans('ATCGatcg', 'TAGCtagc')
    
    return sequence[::-1].translate(tran)

def kmerencode(kmer):
    
    kmerint = 0
    for base in kmer:
        
        kmerint *= 4
        if base=='a' or base =='A':
            
            kmerint += 0
            
        elif base=='t' or base =='T':
            
            kmerint += 3
            
        elif base=='c' or base =='C':
            
            kmerint += 1
        elif base=='g' or base =='G':
            
            kmerint += 2
            
    return kmerint

def kmerdecode(kmerint, size = 31):
    
    index = 0
    text = ['A' for a in range(size)]
    
    while kmerint:
        
        text[index] = "ACGT"[kmerint % 4]
        kmerint //= 4
        index += 1
    
    
    
    return "".join(text[::-1])



class KmerReader:
    
    def __init__(self,  kmersize = 31):
        
        self.kmersize = kmersize
        self.kmerindexes = cl.defaultdict(int)
        
        
    def loadsamples(self, infile):
        
        self.samplesizes = [0]
        self.samples = []
        with open(infile, mode = 'r') as f:
            
            for line in f:
                
                if line.startswith(">"):
                    
                    name = line.split()[0][1:]
                    self.samples.append(name)
                    self.samplesizes.append(0)
                else:
                    self.samplesizes[-1] += len(line.strip())
                    
        self.samplesizes = self.samplesizes[1:]
        self.samplenum = len(self.samples)
        return self
    
    def loadkmers(self, kmerfile, ifmask = 1):
        
        maskcount = 0
        kmercount = 0
        with open(kmerfile, mode = 'r') as f:
            
            for line in f:
                
                if line.startswith(">") == False and (ifmask == 1 or line[:1].isupper()):
                    
                    kmerint = kmerencode(line.strip().split()[0])
                    
                    if line[:1].isupper():
                        kmercount += 1
                        self.kmerindexes[kmerint] = kmercount
                    else:
                        maskcount -= 1
                        self.kmerindexes[kmerint] = maskcount
                
                    
        self.kmernum = kmercount + 1
        self.masknum = abs(maskcount) + 1
        
        return self
    
    def read(self, file):
        
        self.kmercols = bitarray (self.samplenum *  self.kmernum)
        self.kmercounter = [0] * (self.samplenum + 1)
        
        self.maskcols = bitarray (self.samplenum *  self.masknum)
        
        intoperator = 4**(self.kmersize-1)
        kmersize = self.kmersize
        
        groupnum = self.samplenum 
        current_k = 0
        reverse_k = 0
        current_size = 0
        groupindex =-1 
        
        currentline = ""
        kmercounter = 1
        with open(file,mode = 'r') as r:
            
            for line in r:
                
                line = line.strip()
                
                if line.startswith(">"):
                    
                    current_k=0
                    reverse_k=0
                    current_size = 0
                    
                    groupindex += 1
                    continue
                
                for char in line.upper():
                    
                    if char not in ['A','C','G','T']:
                        
                        current_k = 0
                        reverse_k = 0
                        current_size = 0
                        continue 
                    
                    kmercode = ['A','C','G','T'].index(char)
                    
                    if current_size >= self.kmersize:
                        
                        current_k %= intoperator
                        current_k <<= 2
                        current_k += kmercode
                        
                        reverse_k >>= 2
                        reverse_k += (3-kmercode) * intoperator
                        
                    else:
                        
                        
                        current_k <<= 2
                        current_k += kmercode
                        reverse_k += (3-kmercode) * ( 1 << (2*current_size))
                        
                        current_size += 1
                        
                    kmer = max(current_k, reverse_k)
                    
                    if current_size == kmersize:
                        
                        kindex = self.kmerindexes.get(kmer, 0)
                        
                        if kindex > 0:
                            self.kmercols[kindex* groupnum + groupindex] = 1
                            self.kmercounter[groupindex] += 1
                        elif kindex < 0:
                            self.maskcols[abs(kindex)* groupnum + groupindex] = 1
        return self
    
    
    
def selectsamples(inputfile,outputfile,selectlist):
    
    
    iffilter = 0
    with open(inputfile, mode = 'r') as r, open(outputfile, mode ='w') as w:
        
        for line in r:
            
            if line.startswith(">"):
                
                name = line.split()[0][1:]
                
                iffilter = 1
                if name in selectlist:
                    iffilter = 0
                    
            if not iffilter:
                w.write(line)
                
                
def selectkmers(reader, outputfile , filterlist, notusekmer):
    
    masktable = reader.maskcols
    kmertable = reader.kmercols
    samplenum = reader.samplenum 
    
    with open(outputfile, mode = 'w') as w:
        
        for kmer,kindex in reader.kmerindexes.items():
            
            if notusekmer[kindex] == 1:
                continue
            
            if kindex > 0:
                index_start = kindex*samplenum 
                iffilter = any(kmertable[index_start + groupindex] for groupindex in filterlist)
            else:
                index_start = (-kindex)*samplenum 
                iffilter = any(masktable[index_start + groupindex] for groupindex in filterlist)
                
            if not iffilter :
                w.write(">\n{}\n".format(kmerdecode(kmer) if kindex > 0 else kmerdecode(kmer).lower()) )
                

                
def recountkmers(reader, outputfile , filterlist, notusekmer):
    
    kmertable = reader.kmercols
    samplenum = reader.samplenum 
    samplekmer = reader.kmercounter
    
    with open(outputfile, mode = 'w') as w:
        
        for kmer,kindex in reader.kmerindexes.items():
            
            if kindex < 0 or notusekmer[kindex] == 1:
                continue
            
            index_start = kindex*samplenum 
            
            ifany = 0
            for groupindex in filterlist:
                
                if kmertable[index_start + groupindex]:
                    
                    ifany = 1
                    samplekmer[groupindex] += 1
                    
            if not ifany:
                notusekmer[kindex] == 1
                
                
                
def main(args):
    
    kmertype = args.size
    
    kmerfile = args.kmer
    
    if len(kmerfile) == 0:
        kmerfile = args.input + "_kmer.list"
        
    infile = args.input
    
    outfile = args.output
    
    kmeroutfile = args.kmerlist
    if len(kmeroutfile) == 0:
        kmeroutfile = outfile+"_kmer.list"
        
    reader = KmerReader().loadsamples(infile).loadkmers(kmerfile, args.mask).read(infile)
    
    icircle = 0
    finished = 0
    filterlist = set()
    notusekmer = [0] * (reader.kmernum + 1)
    while finished == 0:
        
        icircle += 1
        print("circle: ", icircle)
        
        lastfilternum = len(filterlist)
        outlist = set()
        for i, sample in enumerate(reader.samples):
            
            if i in filterlist:
                continue
            
            size = reader.samplesizes[i]
            count = reader.kmercounter[i]
            
            print(sample+"\t" + str(size)+ "\t" +str(count))
            
            if count> args.cutoff or count > ( size - 2 * args.anchor) * args.rcutoff :
                outlist.add(i)
                
            else:
                filterlist.add(i)
                
        selectsamples(infile,outfile,set([reader.samples[i] for i in outlist]))
        
        selectsamples(infile,outfile+"_exclude",set([reader.samples[i] for i in filterlist]))
        
        if args.repeat == 0 or len(filterlist) == lastfilternum:
            
            selectkmers(reader, kmeroutfile, filterlist, notusekmer)
            selectkmers(reader, kmeroutfile, filterlist, notusekmer)
            break
        
        recountkmers(reader, kmeroutfile, filterlist, notusekmer)
        
        
def run():
    """
        Parse arguments and run
    """
    parser = argparse.ArgumentParser(description="program determine psuedogene")
    parser.add_argument("-i", "--input", help="path to input data file",dest="input", type=str, required=True)
    parser.add_argument("-k", "--kmer", help="path to output file", dest="kmer",type=str, default="")
    parser.add_argument("-s", "--size", help="kmer size", dest="size",type=int, default = 31)
    parser.add_argument("-o", "--output", help="path to output file", dest="output",type=str, required=True)
    parser.add_argument("-a", "--anchor", help="path to output file", dest="anchor",type=int, default=0)
    parser.add_argument("-l", "--list", help="path to kmer output file", dest="kmerlist",type=str, default="")
    parser.add_argument("-c", "--cutoff", help="path to kmer output file", dest="cutoff",type=int, default=1000)
    parser.add_argument("-r", "--rcutoff", help="path to kmer output file", dest="rcutoff",type=float, default=0.1)
    parser.add_argument("-p", "--repeat", help="if repeat", dest="repeat",type=int, default=0)
    parser.add_argument("-m", "--mask", help="if mask", dest="mask",type=bool, default=1)
    parser.set_defaults(func=main)
    args = parser.parse_args()
    args.func(args)
    
    
if __name__ == "__main__":
    run()
