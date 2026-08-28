#!/usr/bin/env python3

import multiprocessing as mul
import os
import argparse
import re
import collections as cl
pathminsize = 50

script_folder = os.path.dirname(os.path.abspath(__file__))
lock1 = mul.Lock()


def makereverse(seq, strd = "-"):
    
    if strd == "+":
        
        return seq
    
    tran=str.maketrans('ATCGatcg', 'TAGCtagc')
    
    return seq[::-1].translate(tran)

def selfblastn(input,nthreads):
    
    #cmd = "makeblastdb -in {}  -dbtype nucl -parse_seqids -out {}_db ".format(input, input)
   
    #dbcmd = "{}/makefastdb.sh  {} {}_db {}".format(script_folder,nthreads, input, input) 
    #os.system(dbcmd)
    
    #cmd = "blastn -task megablast  -query {} -db {}_db -gapopen 10 -gapextend 2 -word_size 30  -perc_identity 95 -evalue 1e-200 -outfmt 17 -out {}_selfblast.out  -num_threads {} -max_target_seqs 100".format( input , input, input,nthreads)
    
    #cmd = "{}/runfast.sh {} {}_db {} 0.95 300 > {}_selfblast.out".format(script_folder,  input,input, nthreads, input)
    #os.system(cmd)
    
    return "{}_selfblast.out".format(input)

def selfalign(input,output,nthreads):
    
    cmd = "{}/runwinnowmap.sh {} {} {} 0.95 150 > {}_selfblast.out".format(script_folder,  input,input, nthreads, input)
    os.system(cmd)
        
    

def getsegments(seq, header, segments, output, anchor = pathminsize):
    
    seq_s = 0
    header = header.strip('>').split()[0]
    
    if header.count("_") > 2:
        
        header_info = header.split("_")
        
        if header_info[-1].replace("-","").isnumeric() and header_info[-2].replace("-","").isnumeric():
            
            seq_s = max( 0, int(header_info[-2]) )
            seq_e = int(header_info[-1])
            
            header = "_".join(header_info[:-2])
            
    #else:
        #header += "_{}_{}".format(0, len(seq))
            
    segments[-1][1] = len(seq)
    
    
    headers_set = set()
    with open(output, mode = 'a') as f:
        
        for seg in segments:
            
            start = max(0,  seg[0] - anchor)
            
            end = min( len(seq),  seg[1] + anchor)
            
            size = end - start
            
            outseq = seq[start:end]
            unmasked =   sum(1 for c in outseq if c.isupper())
            masked = len(outseq) - unmasked
            if unmasked + 0.1*masked < pathminsize:
                continue
            
            if size < 151:
                
                if start == 0:
                    
                    end += 151 - size
                    
                else:
                    
                    start -= min(start, 151 - size)
                    
                    
            if (header, seq_s+start , seq_s+end) in headers_set:
                continue
            
            
            headers_set.add((header, seq_s+start , seq_s+end))
            
                
            f.write(">{}_{}_{}\n{}\n".format(header,  seq_s+start , seq_s+end , seq[start:end]))
    
            
class kmeraligns:
    
    def __init__(self, cigars, ksize = 30, cutoff = 50):
        
        self.cigars = re.findall(r'\d+[a-zA-Z=]', cigars)
        self.ksize = ksize
        self.cutoff = cutoff
        self.min = -cutoff
        self.max = cutoff*3
        self.consensus_f = []
        self.consensus_r = []
        
    def getalignsegs(self):
        
        lastzero = 0
        self.consensus_f = [[0,1]]
        rposi = 0
        qposi = 0 
        lastvar = 0
        lastscore = 0
        score = 0 
        currseg = 1
        for index, cigar in enumerate(self.cigars):
            
            lastscore = score
            thesize = int(cigar[:-1])
            thetype = cigar[-1]
            if thetype == "=" or thetype == "M":
                
                qposi += thesize
                rposi += thesize
                
                score += thesize - self.ksize
                score = min(self.max,max(self.min, score))
                
            elif thetype == "X":
                
                qposi += thesize
                rposi += thesize
                
                score -= ( min(self.ksize, qposi - lastvar) + thesize - 1 )
                score = min(self.max,max(self.min, score))
                
                lastvar = qposi
                
            elif thetype == "H":
                
                rposi += thesize
                continue
            
            elif thetype == "I" :
                
                rposi += thesize
                
                score -= ( min(self.ksize, qposi - lastvar)  )
                score = min(self.max,max(self.min, score))
                
                lastvar = qposi
                
            elif thetype == "D":
                
                qposi += thesize
                
                score -= ( min(self.ksize, qposi - lastvar) + thesize - 1 )
                score = min(self.max,max(self.min, score))
                
                lastvar = qposi
                
            if lastscore * score < 0 or score == 0:
                
                lastzero = index
                
            if currseg <= 0 and score >= self.cutoff:
                
                self.consensus_f.append([lastzero,lastzero])
                currseg = 1
                
            if currseg > 0 and score <= 0:
                
                
                self.consensus_f[-1][1] = index 
                currseg = -1
                
        if currseg > 0:
            if score >= lastscore:
                self.consensus_f[-1][1] = index + 1
            else:
                self.consensus_f[-1][1] = index 
                
        self.consensus_f = [x for x in self.consensus_f if x[1] > x[0]]
        
        lastzero = len(self.cigars)
        self.consensus_r = [[lastzero-1,lastzero]]
        rposi = 0
        qposi = 0 
        lastvar = 0
        lastscore = 0
        score = 0 
        currseg = 1
        for index, cigar in enumerate(self.cigars[::-1]):
            
            index = len(self.cigars) - index
            lastscore = score
            thesize = int(cigar[:-1])
            thetype = cigar[-1]
            
            if thetype == "=" or thetype == "M":
                
                qposi += thesize
                rposi += thesize
                
                score += thesize - self.ksize
                score = min(self.max,max(self.min, score))
                
            elif thetype == "X":
                
                qposi += thesize
                rposi += thesize
                
                if score > 0 and thesize > self.cutoff:
                    score = 0
                    
                score -= ( min(self.ksize, qposi - lastvar) + thesize - 1 )
                score = min(self.max,max(self.min, score))
                
                lastvar = qposi
                
            elif thetype == "H":
                
                rposi += thesize
                continue
            
            elif thetype == "I" :
                
                rposi += thesize
                
                score -= ( min(self.ksize, qposi - lastvar)  )
                score = min(self.max,max(self.min, score))
                
                lastvar = qposi
                
            elif thetype == "D":
                
                qposi += thesize
                
                if score > 0 and thesize > self.cutoff:
                    score = 0
                    
                score -= ( min(self.ksize, qposi - lastvar) + thesize - 1 )
                score = min(self.max,max(self.min, score))
                
                lastvar = qposi
                
            if lastscore * score < 0 or score == 0:
                
                lastzero = index
                
            if currseg <= 0 and score >= self.cutoff:
                
                self.consensus_r.append([lastzero,lastzero])
                currseg = 1
                
            if currseg > 0 and score <= 0:
                
                self.consensus_r[-1][0] = index 
                currseg = -1
                
        if currseg > 0:
            if score >= lastscore:
                self.consensus_r[-1][0] = 0
            else:
                self.consensus_r[-1][0] = 1
                
        self.consensus_r = [x for x in self.consensus_r if x[1] > x[0]]
        
        
        return self
    
    def combineseg(self):
        
        consensus = self.consensus_f + self.consensus_r
        consensus_pos = [a for b in consensus for a in b]
        
        sortindexes = sorted(range(len(consensus_pos)), key = lambda x: consensus_pos[x])
        
        num_curr = 0
        self.consensus_comb = []
        for index, sortindex in enumerate(sortindexes):
            
            if sortindex % 2 == 0:
                num_curr += 1
                
                if num_curr == 1:
                    
                    gapsize = 10000000
                    if len(self.consensus_comb):
                        lastend = self.consensus_comb[-1][-1]
                        gapsize = sum([int(self.cigars[x][:-1]) for x in range(lastend,sortindex) if self.cigars[x][-1]])
                        
                    if gapsize > self.cutoff:
                        self.consensus_comb.append([consensus_pos[sortindex], consensus_pos[sortindex]])
                        
            else:
                num_curr -= 1
                
                if num_curr == 0:
                    
                    self.consensus_comb[-1][-1] = consensus_pos[sortindex]
                    
        return self
    
    def trim(self):
        
        self.consensus_trimed = []
        for segment in self.consensus_comb:
            
            lastvar = 0
            score = -1
            lastzero = -1
            rposi = 0
            qposi = 0 
            lastvar = 0
            lindex = segment[1]-segment[0]
            rindex = segment[1]-segment[0]
            lpass = 0
            for index in range(segment[0],segment[1]):
                
                lastscore = score
                
                cigar = self.cigars[index] 
                thesize = int(cigar[:-1])
                thetype = cigar[-1]
                
                if thetype == "=" or thetype == "M":
                    
                    if score < 0 and thesize > self.cutoff:
                        score = 0
                        
                    score += thesize - self.ksize
                    score = min(self.max,max(self.min, score))
                    
                    if lastscore < 0 and score >= 0:
                        
                        lastzero = index
                        
                        
                    if score > self.ksize:
                        lindex = lastzero if lastscore >= 0 else index
                        lpass = 1
                        break
                    
                elif thetype == "X":
                    
                    if score > 0 and thesize > self.cutoff:
                        score = 0
                        
                    score -= ( min(self.ksize, qposi - lastvar) + thesize - 1 )
                    score = min(self.max,max(self.min, score))
                    
                    lastvar = qposi
                    
                elif thetype == "H":
                    
                    rposi += thesize
                    
                    
                elif thetype == "I" :
                    
                    rposi += thesize
                    
                    score -= ( min(self.ksize, qposi - lastvar)  )
                    score = min(self.max,max(self.min, score))
                    
                    lastvar = qposi
                    
                elif thetype == "D":
                    
                    qposi += thesize
                    
                    score -= ( min(self.ksize, qposi - lastvar) + thesize - 1 )
                    score = min(self.max,max(self.min, score))
                    
                    lastvar = qposi
                    
            score = -1
            lastzero = -1
            rposi = 0
            qposi = 0 
            lastvar = 0
            rpass = 0
            for index in list(range(segment[0],segment[1]))[::-1]:
                
                lastscore = score
                cigar = self.cigars[index] 
                thesize = int(cigar[:-1])
                thetype = cigar[-1]
                
                if thetype == "=" or thetype == "M":
                    
                    if score < 0 and thesize > self.cutoff:
                        score = 0
                        
                    score += thesize - self.ksize
                    score = min(self.max,max(self.min, score))
                    
                    if lastscore < 0 and score >= 0:
                        
                        lastzero = index
                        
                    if score >= self.ksize:
                        rindex = lastzero if lastscore >= 0 else index
                        rpass = 1
                        break
                    
                elif thetype == "X":
                    
                    score -= ( min(self.ksize, qposi - lastvar) + thesize - 1 )
                    score = min(self.max,max(self.min, score))
                    
                    lastvar = qposi
                    
                elif thetype == "H":
                    
                    rposi += thesize
                    
                    
                elif thetype == "I" :
                    
                    rposi += thesize
                    
                    score -= ( min(self.ksize, qposi - lastvar)  )
                    score = min(self.max,max(self.min, score))
                    
                    lastvar = qposi
                    
                elif thetype == "D":
                    
                    qposi += thesize
                    
                    score -= ( min(self.ksize, qposi - lastvar) + thesize - 1 )
                    score = min(self.max,max(self.min, score))
                    
                    lastvar = qposi
                    
            if rindex >= lindex and lpass and rpass:
                self.consensus_trimed.append([lindex, rindex+1])
                
        return self
    def coordinate(self):
        
        self.segment_coordinates = [[] for x in self.consensus_comb]
        
        allstartends = {}
        for index,segment in enumerate(self.consensus_trimed):
            
            allstartends[2*segment[0]] = index
            allstartends[2*segment[1]-1] = index
            
        peakscore = 0
        score = 0
        rposi = 0
        qposi = 0 
        lastvar = 0
        for index, cigar in enumerate(self.cigars):
            
            if 2*index in allstartends:
                self.segment_coordinates[allstartends[2*index]].append((qposi,rposi))
                score = 0
                peakscore = 0
                
            thesize = int(cigar[:-1])
            thetype = cigar[-1]
            
            if thetype == "=" or thetype == "M":
                
                qposi += thesize
                rposi += thesize
                score += thesize - self.ksize
                score = max(self.min, score)
                peakscore = max(peakscore, score)
                
            elif thetype == "X":
                
                qposi += thesize
                rposi += thesize
                score -= ( min(self.ksize, qposi - lastvar) + thesize - 1 )
                score = max(self.min, score)
                
                lastvar = qposi
                
            elif thetype == "H":
                
                rposi += thesize
                
            elif thetype == "I" :
                
                rposi += thesize
                
                score -= ( min(self.ksize, qposi - lastvar)  )
                score = max(self.min, score)
                
                lastvar = qposi
                
            elif thetype == "D":
                
                qposi += thesize
                
                score -= ( min(self.ksize, qposi - lastvar)  )
                score = max(self.min, score)
                
                lastvar = qposi
                
            if 2*index + 1 in allstartends:
                
                
                segindex = allstartends[2*index + 1]
                self.segment_coordinates[segindex].append((qposi,rposi))
                self.segment_coordinates[segindex].append(peakscore)
                #cigar_segs = "".join([self.cigars[x] for x in range(*self.consensus_trimed[segindex])])
                
                #self.segment_coordinates[segindex].append(cigar_segs)
                
                
        return self
    
    
def findkmeraligns(cigars, ksize = 30, cutoff = 300):
    
    consensus = kmeraligns(cigars, ksize, cutoff).getalignsegs().combineseg().trim().coordinate()
    
    output = [x for x in consensus.segment_coordinates if len(x) and x[1][0] - x[0][0] > cutoff and x[1][1] - x[0][1] > cutoff and x[2] > cutoff]
    
    return  output  




def findrepeats(queries_index,alignfile,qindex = 1,filterlength = 300):
    
    
    allaligns = []
    with open(alignfile, mode = 'r') as f:
        
        for line in f:
            
            line = line.strip()
            
            if len(line) == 0 or line[0]=="@":
                continue
            
            elements = line.strip().split()

            qname = elements[2]
            thisindex = queries_index[qname] if qname in queries_index else int(qname.split("_")[-1]) -1 
            
            if thisindex != qindex:
                continue
            
            rindex = queries_index[elements[0]]
            
            if rindex > qindex:
                continue
            
            
            #name, qsize, qstart, qend, path, rsize, rstart, rend = elements[:8]
            
            
            #match  = [int(x[:-1]) for x in re.findall(r'\d+M',cigar )] + [0]
            #largestmatch = max(match)
            #allmatch = sum(match)
            #mismatch = int( [x for i,x in enumerate(elements) if x[:5] == "NM:i:"][0][5:] )
            #if  largestmatch <= 30 or allmatch < 130 or mismatch >= 0.1 * allmatch  or elements[4] != '255':
            
            identity =float( [x for i,x in enumerate(elements) if x[:5] == "PI:f:"][0][5:])
            match = int( [x for i,x in enumerate(elements) if x[:5] == "AS:i:"][0][5:] ) 
            
            if identity < 90 or match < 100 or elements[4] != '255':
                continue
            
            
            strand, qstart, cigar = elements[1], int(elements[3])-1, elements[5]
            
            strand = 1 if strand == '0' else -1
            
            rstart = [int(x[:-1]) for x in re.findall(r'^\d+H',cigar)]
            if len(rstart):
                rstart = rstart[0]
            else:
                rstart = 0
                
            size_withH =  sum([int(x[:-1]) for x in re.findall(r'\d+[HMXI]',cigar)])
            rsize = sum([int(x[:-1]) for x in re.findall(r'\d+[=MXI]',cigar)])
            qsize = sum([int(x[:-1]) for x in re.findall(r'\d+[=MXD]',cigar)])
            
            #rend = rstart + rsize
            
            #if strand == -1:
            #rstart = size_withH  - rstart
            #rend = size_withH - (rend)
            
            #aligns = findkmeraligns(cigar)
            aligns = [[(0, rstart), (qsize,  rstart+rsize), match]]

            for align in aligns:
                
                (qstart_a,rstart_a),(qend_a,rend_a),score = align
                qstart_a += qstart
                qend_a += qstart
                
                if strand == -1:
                    rstart_a = size_withH  - rstart_a
                    rend_a = size_withH - rend_a
                    
                    
                if (qindex, min(qstart_a,qend_a)) > (rindex, min(rend_a, rstart_a)+filterlength):
                    if qindex == rindex:
                        cstart = max(min(qstart_a,qend_a), max(rend_a, rstart_a))
                        cend = max(qstart_a,qend_a)
                        align = [cstart, cend, rstart_a, rend_a]
                    else:
                        align = [qstart_a, qend_a, rstart_a, rend_a]
                    
                    allaligns.append(align)
                    
    return allaligns


def mergeregion(allaligns, allow_gap = pathminsize):
    
    all_coordinates = [x for y in allaligns for x in y[:2]]
    
    sort_index = sorted(range(len(all_coordinates)), key = lambda x: all_coordinates[x])
    
    allgroups = []
    current_group_coordi = []
    
    number_curr_seq = 0
    last_coordinate = -allow_gap - 1
    for index in sort_index:
        
        coordinate = all_coordinates[index]
        
        if index %2 == 0 :
            
            if number_curr_seq == 0 and coordinate - last_coordinate>allow_gap :
                
                allgroups.append([coordinate, coordinate])
                
            number_curr_seq += 1
            
        else:
            allgroups[-1][-1] = coordinate
            
            number_curr_seq -= 1
            
        last_coordinate = coordinate
        
        
    return allgroups

def process_eachquery(queries_index, alignfile, index, seqs, qname, header, output):
   
    allaligns = findrepeats(queries_index,alignfile,index)
    allaligns = mergeregion(allaligns)
    
    allaligns = [0] + [x for y in allaligns for x in y] + [10000000000000000]
    
    allaligns = [[allaligns[2*i], allaligns[2*i+1]] for i in range(len(allaligns)//2)]
   
 
    lock1.acquire()
    getsegments(seqs[qname], header, allaligns, output)
    lock1.release()

    return sum([x[1]-x[0] for x in allaligns] + [0])

def cleanrepeats(inputfile, output, nthreads = 1):
   
    alignfile = "{}_selfblast.out".format(inputfile) 
    seqs = {}
    headers = {}
    queries = []
    totalbefore = 0
    with open(inputfile, mode = 'r') as f:
        for line in f:
            if len(line) and line[0] == '>':
                qname = line[1:].split()[0]
                header = line.strip()
                headers[qname] = header
                seqs[qname] = ""
                queries.append(qname)
            else:   
                seqs[qname] += line.strip()
                totalbefore += len(seqs[qname])
    queries_index = {x:i for i,x in enumerate(queries)}
    try:
        os.remove(output)
    except:
        pass
    
    manager = mul.Manager()
    seqs = manager.dict(seqs)
    p=mul.Pool(processes=nthreads)
    async_results = []
    for index,qname in enumerate(queries):
      
        #process_eachquery(queries_index, alignfile, index, seqs, qname, headers[qname], output) 
        ar = p.apply_async(process_eachquery, (queries_index, alignfile, index, seqs, qname, headers[qname], output))
        async_results.append(ar)
    p.close()
    p.join()

    total = 0
    for ar in async_results:
        r = ar.get()          # if worker crashed, this will raise (good)
        total += (r or 0)     # None -> 0

    return totalbefore-total

        
def main(args):
    
    total = 0
    cleaned = 100000000000
    cycle = 0
    theinput = args.input
    theoutput = args.output
    thetemp  = args.output + ".temp"
    while cleaned > 300 and  cycle < args.cycle:
        if cycle > 0:
            theinput = theoutput
            theoutput = args.output if theoutput  == thetemp else thetemp

        selfalign(theinput, theoutput, args.threads)
        cleaned = cleanrepeats(theinput, theoutput, args.threads)
        cycle += 1
        total += cleaned
        print(f"[cycle] {args.input} cycle {cycle} of {args.cycle} clean = {cleaned}")

    print(f"[cleanrepeats] {args.input} total clean = {total}") 
    
def run():
    """
        Parse arguments and run
    """
    parser = argparse.ArgumentParser(description="program distract overlap genes and alignments on contigs")
    
    parser.add_argument("-i", "--input", help="path to output file", dest="input", type=str,required=True)
    parser.add_argument("-o", "--output", help="path to output file", dest="output", type=str,required=True)
    parser.add_argument("-t", "--threads", help="path to output file", dest="threads", type=int,default = 1)
    parser.add_argument("-c", "--cycle", help="path to output file", dest="cycle", type=int,default = 1)
    parser.set_defaults(func=main)
    args = parser.parse_args()
    args.func(args)
    
    
if __name__ == "__main__":
    run()
