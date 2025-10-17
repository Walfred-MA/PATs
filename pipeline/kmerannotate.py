#!/usr/bin/env python3

import collections as cl
import re
import argparse

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

def kmerdecode(kmerint, size = 31):
    
    index = 0
    text = ['A' for a in range(size)]
    
    while kmerint:
        
        text[index] = "ACGT"[kmerint % 4]
        kmerint //= 4
        index += 1
        
    return "".join(text[::-1])

class CIGAR:
    
    def __init__(self, fullcigar, allpaths, allpathindex, pathsizes):
        
        self.fullcigar = fullcigar
        
        self.cigars = re.findall("[><][^>^<]+",fullcigar)
        self.allpaths =  re.findall(r'[><][^>^<]+', allpaths)
        
        self.rposi = 0
        self.qposi = 0
        self.strd_t = 0
        self.path = 0
        self.segmentindex = 0
        
        self.pathindex = 0
        self.nextqposi = 0
        
        qposi = 0
        self.segments = []
        for i,(cigar,path) in enumerate(zip(self.cigars,self.allpaths)):
            
            cigar = cigar[0] + cigar[1:].split(":")[1]
            
            strd = path[0]
            pathindex = allpathindex.get(path[1:],0)
            rstart = 0 if strd == ">" else pathsizes[i]-1
            strd = 1 if strd == ">" else -1
            
            cigar_segs = re.findall(r'\d+[A-Z=]',cigar[1:])
            cigar_segs = [(int(x[:-1]), x[-1]) for x in cigar_segs]
            
            rposi = rstart
            for (cigar_s, cigar_t) in cigar_segs:
                
                if cigar_t == 'I':
                    self.segments.append((qposi,rposi - 16 * strd,0,pathindex))
                else:
                    self.segments.append((qposi,rposi - 16 * strd ,strd,pathindex))
                    
                if cigar_t in ['H', 'D']:
                    rposi += strd * cigar_s
                elif cigar_t == 'I':
                    qposi += cigar_s
                else:
                    rposi += strd * cigar_s 
                    qposi += cigar_s
                    
        self.segments.append((-1,0,0,0))
        self.segments.append((1000000000000,0,0,0))
        
    def pop(self):
        
        
        while self.qposi == self.nextqposi:
            
            self.qposi,self.rposi,self.strd_t,self.path = self.segments[self.segmentindex]
            self.segmentindex += 1
            
            self.nextqposi = self.segments[self.segmentindex][0]
            
            
            
        self.qposi += 1
        self.rposi += self.strd_t
        

        return self.path, self.rposi
    
class KmerData:
    
    def __init__(self, kmersize = 31):
        
        self.kmersize = kmersize

        self.mask = set()
        
        self.kmerslist = {}
        
        self.samplesizes = []   
        
        self.genekmercounts = []
        
        self.kmerflag=cl.defaultdict(lambda: ('0', '0') )
        
        self.sampleloc = dict()

        self.multicounter = cl.defaultdict(set)
        
    def LoadSeqs(self, seqfile):
        
        index = 0
        samples = []
        with open(seqfile,mode = 'r') as r:
            for line in r:
                
                line = line.strip()
                
                if len(line):
                    if line[0]==">":
                        self.samplesizes.append(0)
                        name = line.split()[0][1:]
                        samples.append(name)
                        self.sampleloc[name] = line.split()[1]
                        self.genekmercounts.append( 0 )
                        index += 1
                    else:
                        self.samplesizes[-1] += len(line.strip())
                        
                        
        self.sampleslist = list(samples)
        
        self.kmerlocation = [dict() for x in samples]
        
        
    def LoadKmers(self, kmerfile):
        
        self.kmerindex = 0
        with open(kmerfile, mode = 'r') as f:
            
            for line in f:
                
                if len(line) and line[0] != ">":
                    
                    self.kmerslist[kmerencode(line.strip().split()[0])] = self.kmerindex
                    self.kmerindex += 1
                    if not line.isupper():
                        self.mask.add(self.kmerindex - 1)

                    
        self.kmerlocations = [(0,0) for x in range(self.kmerindex+2)]
        
    def addflag(self, flagfile):
        
        with open(flagfile, mode = 'r') as f:
            
            for line in f:
                
                if len(line) and line[0] != ">":
                    
                    line = line.strip().split()
                    
                    kmer = kmerencode(line[0])
                    
                    if kmer in self.kmerslist:
                        
                        kmerflag = line[1]
                        kmerratio = round ( float(line[2])/max(1.0, float(line[3])) , 6)
                        
                        if int(kmerflag) > int(self.kmerflag[kmer][0]):
                            self.kmerflag[kmer] =  ( kmerflag , str( kmerratio )  )
                            
                            
    def LoadCigars(self, graphfile, pathindex, pathsizes):
        
        self.allcigars = {}
        
        with open(graphfile, mode = 'r') as f:
            
            for line in f:
                
                if len(line) == 0:
                    continue
                
                seqname, path, fullcigar, rranges, qranges = line.split('\t')
               

                sizes = [pathsizes[x]  for x in path.replace("<",">").split(">")[1:]]
                self.allcigars[seqname] = CIGAR(fullcigar,path, pathindex, sizes)
                
    def AddKmer(self,kmer,pathinfo):
        
        kmerindex = self.kmerslist.get(kmer, -1) 
        
        if kmerindex >= 0:
            
            if self.kmerlocations[kmerindex] != (0,1):
                
                pathname, pathloc = self.kmerlocations[kmerindex]
                
                if pathinfo[0] != pathname:
                    if pathname == 0:
                        self.kmerlocations[kmerindex] = (0,2)

                elif pathinfo[1] != pathloc:

                    self.multicounter[kmerindex].add(pathinfo[1])

                        
            else:
                self.kmerlocations[kmerindex] = pathinfo[:2]



                
        return kmerindex
    
    
    def ReadKmers(self, seqfile):
        
        intoperator = 4**(self.kmersize-1)
        
        current_k = 0
        reverse_k = 0
        current_size = 0
        
        currentline = ""
        
        with open(seqfile,mode = 'r') as r:
            
            index = 0
            for line in r:
                
                if len(line) ==0:
                    continue
                
                if line[0]==">":
                    index += 1
                    
                    name = line.split()[0][1:]
                    cigars = self.allcigars[name]
                    
                    current_k=0
                    reverse_k=0
                    current_size = 0
                    
                    continue
                
                for posi,char in enumerate(line.upper()):
                    
                    if char =="\n":
                        continue
                    
                    pathinfo = cigars.pop()
                    if char not in ['A','T','C','G']:
                        
                        current_k = 0
                        reverse_k = 0
                        current_size = 0
                        
                        continue
                    
                    if current_size >= self.kmersize:
                        
                        current_k %= intoperator
                        current_k <<= 2
                        current_k += ['A','C','G','T'].index(char)
                        
                        reverse_k >>= 2
                        reverse_k += (3-['A','C','G','T'].index(char)) * intoperator
                        
                    else:
                        
                        
                        current_k <<= 2
                        current_k += ['A','C','G','T'].index(char)
                        reverse_k += (3-['A','C','G','T'].index(char)) * ( 1 << (2*current_size))
                        
                        current_size += 1
                        
                    if current_size >= self.kmersize :  
                        kmer = max(current_k, reverse_k)
                        findindex = self.AddKmer(kmer,(pathinfo[0],pathinfo[1]))
                        
                        
                        
                        
                        
                        
def kmerannotate(align, graph, seq, kmer,output, errorlist):
    
    KmerReader = KmerData()
    KmerReader.LoadKmers(kmer)
   
    headers = []
    fullheaders = cl.defaultdict(str)
    if len(graph):
        pathname = ""
        pathsizes = cl.defaultdict(int)
        unmassizes = cl.defaultdict(int)
        headers = []
        with open(graph, mode = 'r') as f:
            
            for line in f:
                
                if len(line) and line[0] == '>':
                    pathname = line.split()[0][1:]
                    fullheaders[pathname] = "+"+"\t".join(line.strip().split()[:2])[1:]
                    headers.append(pathname)
                elif len(line):
                    unmassizes[pathname] += sum(1 for c in line.strip() if c.isupper())
                    pathsizes [pathname] = len(line.strip())

        headers = [x for x in headers if unmassizes[x] > 300]
        header_index = {header:i+1 for i,header in enumerate(headers)}

        KmerReader.kmerlocations = [(0,1) for x in range(KmerReader.kmerindex+2)]

        KmerReader.LoadSeqs(seq)
        KmerReader.LoadCigars(align, header_index, pathsizes)
        KmerReader.ReadKmers(seq)
        if len(errorlist):
            KmerReader.addflag(errorlist)
        
        
    with open(output, mode ='w') as f:
        
        for header in headers:
            
            f.write(fullheaders[header] + "\n")
            
        for kmerint, index, in KmerReader.kmerslist.items():
            
            pathindex, location = KmerReader.kmerlocations[index]
            
            if pathindex <= 0:
                pathindex =0 
                location = 0
        
            if index in KmerReader.multicounter:
                #print(kmerdecode(kmerint), pathindex, location,KmerReader.multicounter[index])
                location = min(14,len(KmerReader.multicounter[index])+1)

            kmer = kmerdecode(kmerint) if index not in KmerReader.mask else kmerdecode(kmerint).lower()

            f.write(">\n{}\t{}\t{}\t{}\n".format( kmer, pathindex, location , "\t".join(list(KmerReader.kmerflag.get(kmerint,["0","1.0"]))) ) )
            
def main(args):
    
    kmerannotate(args.align, args.graph, args.seq, args.kmer,args.output,args.error)
    
def run():
    """
        Parse arguments and run
    """
    parser = argparse.ArgumentParser(description="program distract overlap genes and alignments on contigs")
    
    parser.add_argument("-s", "--seq", help="path to output file", dest="seq", type=str,default = "")
    parser.add_argument("-k", "--kmer", help="path to output file", dest="kmer", type=str,required=True)
    parser.add_argument("-g", "--graph", help="path to output file", dest="graph", type=str,default = "")
    parser.add_argument("-a", "--align", help="path to output file", dest="align", type=str,default = "")
    parser.add_argument("-o", "--output", help="path to output file", dest="output", type=str,required=True)
    parser.add_argument("-e", "--error", help="path to output file", dest="error", type=str,default = "")
    
    parser.set_defaults(func=main)
    args = parser.parse_args()
    args.func(args)
    
    
if __name__ == "__main__":
    run()
