#!/usr/bin/env python3

import os
import argparse
import re
import collections as cl
import datetime
import multiprocessing as mul
from statistics import mode
import sys

sys.setrecursionlimit(1000000)

pathminsize = 50


def cigar_cutrrange(cigars, rstart , rend, strd =1):
    
    cigars = re.findall(r'\d+[HSMIDX=]', cigars)
    if strd == -1:
        cigars = cigars[::-1]
        
    qstart = -1
    qend = -1
    
    rposi = 0
    qposi = 0
    newcigars = []
    
    newref = ""
    for cigar in cigars:
        
        lastrposi = rposi
        lastqposi = qposi
        
        thesize = int(cigar[:-1])
        thetype = cigar[-1]
        
        if thetype == "=" or thetype == "M":
            
            qposi += thesize
            rposi += thesize
            
            
        elif thetype == "X" or thetype == "S":
            
            qposi += thesize
            rposi += thesize
            
        elif thetype == "D" or thetype == "H":
            
            rposi += thesize
            
        elif thetype == "I":
            
            qposi += thesize
            
        if rposi <= rstart:
            continue
        elif lastrposi > rend:
            break
        
        segmentsize =  min(rposi, rend)- max(lastrposi, rstart)
        if segmentsize:
            newcigars.append("{}{}".format(segmentsize,thetype))
        elif thetype == 'I':
            newcigars.append(cigar)
            
        if qstart == -1 and rstart >= lastrposi:
            
            if thetype in ['M','=', 'X','S', 'I']:
                
                qstart = lastqposi + rstart - lastrposi
            else:
                qstart = lastqposi
                
        if qend == -1 and rend < rposi :
            
            if thetype in ['M','=', 'X','S', 'I']:
                
                qend = lastqposi + rend - lastrposi
            else:
                qend = lastqposi
                
    if rend >= rposi:
        qend = qposi
        if rend-rposi:
            newcigars.append('{}H'.format(rend-rposi))
        
                    
    return  qstart, qend, "".join(newcigars)


def makereverse(seq):
    
    if type(seq) == type(""):
        tran=str.maketrans('ATCGatcg', 'TAGCtagc')
        return seq[::-1].translate(tran)
    else:

        seq[-1] = -1 * seq[-1]
        return seq




def rename_rep(allpathes):
    
    index_count = cl.defaultdict(int)
    
    allpathes_new = []
    for index in allpathes:
        
        index = index.split('_')[0]
        index_count[index] += 1
        if index_count[index] > 1:
            
            index +="_" + str(index_count[index])
            
        allpathes_new.append(index)
        
    return allpathes_new

def getchildren(vertices_children, node):
    
    children = vertices_children[node]
    
    if len(children) == 0:
        return []
    
    allchildren = sorted([[child]+getchildren(vertices_children, child) for child in children], key = lambda x: len(x))
    
    allchildren = sum(allchildren, [])
    
    return allchildren



def getlinearpath(allchunks_onquery):
    
    links = cl.defaultdict(list)
    reverselink = cl.defaultdict(list)
    
    allindex = set(['0'])   
    for query, chunks_thisquery in allchunks_onquery.items():
        
        last_index = '0'
        for query, (chunk_index, breakinfo) in enumerate(chunks_thisquery):
            
            links[last_index].append(chunk_index)
            last_index = chunk_index
            
        links[last_index] = links[last_index]
        
        
    links_rev = cl.defaultdict(list)
    for name, link in links.items():
        for key in link:
            links_rev[key].append(name)
            
            
    vertices_all = list(links.keys()) 
    vertices_sort = []
    
    for nextname,link in links_rev.items():
        
        counter = sorted([(x,v) for v,x in cl.Counter(link ).items()], reverse = 1)
        
        for (count,name) in counter:
            vertices_sort.append( ( 1.0*count/(len(link) ) , count, name, nextname ) )
            
            
    vertices_sort = sorted(vertices_sort, reverse = 1)
    vertices_sort = [x[-2:] for x in vertices_sort]
    
    vertices_linknodes = cl.defaultdict(list)
    #vertices_parent = {'0':None}
    
    for (lastname, nextname) in vertices_sort:
        
        vertices_linknodes[lastname].append(nextname)
        #vertices_parent[nextname] = lastname
        
        
    allnames = sorted(list(links.keys()))
    allpathes = ['0']
    
    vertices_used = set(['0'])
    vertices_children = cl.defaultdict(list)
    pastnodes = ['0']
    currentnode = '0'
    
    for i in range(2*len(allnames)):
        
        nextnodes = [x for x in vertices_linknodes[currentnode] if x not in vertices_used]
        if len(nextnodes) > 0:
            nextnode = nextnodes[0]
            vertices_children[currentnode].append(nextnode)
            vertices_used.add(nextnode)
            pastnodes.append(currentnode)
        else:
            nextnode = pastnodes.pop()
            
        currentnode = nextnode 
        
        if len(vertices_used) == len(allnames):
            break
        
        
    allpathes = getchildren(vertices_children, '0') 
    
    allpathes = rename_rep(allpathes)
    
    if 'e' in allpathes:
        allpathes.remove('e')
        
    if '0' in allpathes:
        allpathes.remove('0')
        
        
        
    allpathes_order = {x:i for i,x in enumerate(allpathes)}
    
    
    
    return allpathes_order


def getchunkstrand(chunks, chunks_onquery, breaks):
    
    chunks_strd = cl.defaultdict(list)
    
    all_strands = [y[4] for y in breaks for x in y[:2]]
    
    chunks_strd = cl.defaultdict(list)
    for query, query_chunks in chunks_onquery.items():
        
        for chunkindex, coordinateindex in query_chunks:
            
            
            chunks_strd[chunkindex].append(all_strands[coordinateindex])
            
    for chunkindex, strands in chunks_strd.items():
        
        if len(strands) == 0:
            strands += [1]
            
        chunks_strd[chunkindex] = mode(strands)
        
    return chunks_strd


def getchunksoncontig(breaks):
    
    breaks = sorted(breaks)
    
    allchunks = []
    chunkstarts = []
    numchunks = 0
    
    for break_index, chunk_index in breaks:
        
        if break_index % 2 ==0 :
            
            chunkstarts.append(chunk_index)
            
        else:
            
            chunkstart = chunkstarts.pop()
            
            allchunks.append([chunkstart, chunk_index, break_index])
            
    return sum( [ [ [x, chunk[2]] for x in list(range(chunk[0], chunk[1]))] for chunk in allchunks], [])


def combinebreaks(breaks):
    
    all_coordinates = [x for y in breaks for x in y[:2]]
    
    sort_index = sorted(range(len(all_coordinates)), key = lambda x: all_coordinates[x])
    
    break_groups = []
    break_group = []
    
    break_group_index = []
    break_group_indexs = []
    
    last_coordinate = 0
    for index in sort_index:
        
        coordinate = all_coordinates[index]
        
        if coordinate - last_coordinate < pathminsize:
            
            break_group.append(coordinate)
            break_group_index.append(index)
            
        else:
            
            break_groups.append(break_group)
            break_group_indexs.append(break_group_index)
            
            break_group = [coordinate]
            break_group_index = [index]
            
        last_coordinate = coordinate
        
    break_groups.append(break_group)
    break_group_indexs.append(break_group_index)
    
    
    return break_groups, break_group_indexs

#get all chunks on each gene, and each chunk information is index of chunk, index of breaks among breaks of all genes.
def getchunks_frombreak(break_groups, break_group_indexs, breaks):
    
    all_queries = [(y[2],y[3]) for y in breaks for x in y[:2]]
    
    chunk_starts = []
    breaks_onquery = cl.defaultdict(list)
    
    chunk_index = 0 
    for (break_group, break_group_index) in zip(break_groups, break_group_indexs):
        
        if break_group == []:
            continue
        
        themode = mode(break_group)
        
        chunk_starts.append(themode)
        
        
        for coodinate_index in break_group_index:
            
            #for each query, a list information on all chunks it has, each chunk information is index of break, index of chunk after merging breaks
            breaks_onquery[all_queries[coodinate_index]].append(( coodinate_index, chunk_index))
            
            
        chunk_index += 1
        
    chunks = [ [first, second] for first, second in zip(chunk_starts,chunk_starts[1:])]
    
    chunks_onquery = {query: getchunksoncontig(breaks) for query, breaks in breaks_onquery.items()}
    
    return chunks, chunks_onquery


def getrefsizes(reffile):
    
    with open(reffile, mode = 'r') as f:
        
        reads = [read.splitlines() for read in f.read().split(">")[1:]]
        sizes = {read[0].split()[0]: len("".join(read[1:])) for read in reads}
        
    return sizes

def readalignments(inputfile):
    
    index = 0
    allpathes =  []
    lines = []
    allrows = []
    allqueries = []
    with open(inputfile, mode = 'r') as f:
        
        for line in f:
            
            line = line.strip().split('\t')
            if len(line) < 5:
                continue
            query, path, cigar, pathranges, qranges = line
            pathranges = pathranges.split(";")
            strands = [1 if x =='>' else -1 for x in re.findall(r'[><]', path)]
            path = path.replace("<",">").split(">")[1:]
            cigars = cigar.replace("<",">").split(">")[1:]
            qranges = list(map(int,qranges.split(";")))
            
            
            allalignments = [(p,r,index,i,s) for i, (p, r, s,c) in enumerate(zip(path,pathranges,strands,cigars)) ]
            
            allqueries.append(query)
            allpathes.append(allalignments)
            
            index += 1
            
            
    return allpathes,allqueries


def getchunks(allpathaligns):
    
    allchunks = []
    
    
    #get all chunks by merging close breaks, each chunk get is [contig, start, end, strand]
    #also get chunks information on all genes, each gene has a list of chunks information, each piece of information is index of chunk, the alignment associated to that chunk on gene
    allchunks_onpath = cl.defaultdict(list)
    for (refname), pathaligns in allpathaligns.items():
        
        
        #each break breaks is path_start, path_end, queryindex, alignindex, strand
        
        #combine the overlap chunks
        break_groups, break_group_indexs = combinebreaks(pathaligns)
        
        #get querychunkr and chunk alignments on all genes
        #chunks here is [start, end]
        chunks,  chunks_onquery = getchunks_frombreak(break_groups, break_group_indexs, pathaligns)
        
        #get most common strand of chunks
        chunks_strd = getchunkstrand(chunks, chunks_onquery, pathaligns)
        
        index_offsite = len(allchunks)
        
        allchunks.extend([[refname]+x+[chunks_strd[i]] for i,x in enumerate(chunks)])
                
        for (query,pathindex), chunks_onquerythis in chunks_onquery.items():
            
            for (chunkindex, pathalignindex) in chunks_onquerythis:
                
                allchunks_onpath[(query,pathindex)].append([ chunkindex + index_offsite, pathaligns[pathalignindex//2] ])
                
    sorted_keys = sorted(allchunks_onpath.keys())
    
    chunkonquery_counter = cl.defaultdict(int)
    allchunks_onquery = cl.defaultdict(list)
    #assignning uniq names for duplicated chunks 
    for (query,pathindex) in sorted_keys:
        
        chunks_onquerythis = allchunks_onpath[(query,pathindex)]
        
        chunk_indexes_uniq = []
        for (chunkindex, originpathinfo) in chunks_onquerythis:
            
            #originbreakinfo is start, end, queryindex, alignindex, strand
            
            chunkonquery_counter[(chunkindex,query)] += 1
            
            if chunkonquery_counter[(chunkindex,query)] > 1:
                index_uniq = str(chunkindex+1)+"_"+str(chunkonquery_counter[(chunkindex,query)])
                chunk_indexes_uniq.append( (index_uniq , originpathinfo) )
            else:
                index_uniq = str(chunkindex+1)
                chunk_indexes_uniq.append( (index_uniq , originpathinfo) )
                
        allchunks_onquery[query].extend(chunk_indexes_uniq)
    
    
    return allchunks, allchunks_onquery


def getpathaligns(allpathes, refsizes):
    
    allpathaligns = cl.defaultdict(list)
    
    for path in allpathes:
                
        for refname, region ,queryindex, alignindex, strand in path:
            
            
            
            region = list(map(int, region.split("_")))
            
            if strand == -1 :
                region = [refsizes[refname]-region[1], refsizes[refname]-region[0]]
                
            allpathaligns[refname].append( region + [queryindex, alignindex, strand] )
            
    
    return allpathaligns

def getchunkseq(reffile, allchunks):
    
    with open(reffile, mode = 'r') as f:
        
        reads = [read.splitlines() for read in f.read().split(">")[1:]]
        reads = {read[0].split()[0]:"".join(read[1:]) for read in reads}
        
    chunk_reads = []
    
    for index, chunk in enumerate(allchunks):
        
        
        read = reads[chunk[0]][chunk[1]:chunk[2]]
        
        if chunk[-1] == -1:
            
            read = makereverse(read)
            
        chunk_reads.append(read)
        
    return chunk_reads

def getchunkpath(reffile, allchunks):
        
    return chunk_reads


def putchunks(output, allchunks_withorders ):
    
    with open(output, mode = 'w') as w:
        
        for (querychunk, chunk_read, chunkorders, chunkindex) in allchunks_withorders:
            
            orders = ";".join(list(map(str, chunkorders))) + ";"
            
            if querychunk[-1] == -1:
                temp = querychunk[2] 
                querychunk[2] = querychunk[1]
                querychunk[1] = temp
                
            querychunk = "_".join(list(map(str,querychunk[:3])))
            
            w.write("S\t{}\t{}\t{}\tLN:i:{}\tSN:Z:{}\n".format("s_"+str(chunkindex+1), chunk_read, orders, len(chunk_read),querychunk))
            
def putalignments(output, allquery_cigars, query_length, query_order):
    
    with open(output, mode = 'a') as w:
        
        for queryname in query_order:
            
            (leftgap, rightgap, cigar, chunknames) = allquery_cigars.get(queryname, [0,0,"",""])
            
            w.write("L\t{}\t+\t{}\t+\t{}\tSR:i:{}\tL1:i:{}\tL2:i:{}\n".format(queryname, chunknames,cigar,  rightgap - leftgap, leftgap , query_length[queryname] - rightgap ))
            
            
def getqueryseq(queryfile):
    
    with open(queryfile, mode = 'r') as f:
        
        reads = [read.splitlines() for read in f.read().split(">")[1:]]
        queryorder = [read[0].split()[0] for read in reads]
        reads = {read[0].split()[0]:"".join(read[1:]) for read in reads}
        querylength = {name:len(read) for name,read in reads.items()}
        
    return reads, querylength, queryorder


def getquerycigar(queryfile, inputfile):
    
    current_name = ""
    querylength = cl.defaultdict(int)
    queryorder = []
    with open(queryfile, mode = 'r') as f:
        
        for line in f:
            
            if line.startswith(">"):
                current_name = line[1:].split()[0]
                queryorder.append(current_name)
                querylength [current_name] = 0
            else:
                querylength [current_name] += len(line.strip())
    
    reads = cl.defaultdict(list)
    pattern = re.compile(r'[<>][^<>]+')
    with open(inputfile, mode = 'r') as f:
    
        for line in f:
            
            line = line.strip().split()
            if len(line) < 3:
                continue
            name = line[0]
            reads[name] = list( pattern.findall(line[2]))
    
    return reads, querylength, queryorder


def getgappedcodi(seq, posi):
    if posi==[]:
        return []
    
    piece_coordi = [a.span() for a in re.finditer(r'[A-Za-z]+',seq)]
    
    lastend = 0
    index = 0
    posi_offsite = 0
    
    newposi = [a if i%2 == 0 else a - 1 for i,a in enumerate(posi)]
    for i, (start, end) in enumerate(piece_coordi):
        
        posi_offsite += start - lastend 
        
        while index<len(newposi) and newposi[index]<end - posi_offsite:
            
            newposi[index] += posi_offsite
            
            index += 1
            
        lastend = end
        
    for i in range(index,len(newposi)):
        
        newposi[i] +=  posi_offsite
        
    newposi = [a if i%2 == 0 else a + 1 for i,a in enumerate(newposi)]
    
    
    return newposi

def getungappedcodi(seq, posi):
    if posi==[]:
        return []
    
    piece_coordi = [a.span() for a in re.finditer(r'[A-Za-z]+',seq)]
    
    lastend = 0
    index = 0
    posi_offsite = 0
    
    newposi = [a if i%2 == 0 else a - 1 for i,a in enumerate(posi)]
    
    start = piece_coordi[0][0] if len(piece_coordi) else 0
    
    while index<len(newposi) and newposi[index] < start:
        
        newposi[index] = 0
        index += 1
        
    for i, (start, end) in enumerate(piece_coordi):
        
        posi_offsite += start - lastend 
        
        while index<len(newposi) and newposi[index] < end:
            
            newposi[index] -= posi_offsite
            index += 1
            
        lastend = end
        
    posi_offsite += len(seq) - lastend
    
    for i in range(index,len(newposi)):
        newposi[index] -= posi_offsite
        
    newposi = [a if i%2 == 0 else a + 1 for i,a in enumerate(newposi)]
    
    
    return newposi


def runstrecher(query, ref, name):
    
    from graphmaskstretcher import maskfile
    
    tempfile_q = name+"_q.fa"
    tempfile_r = name+"_r.fa"
    tempfile_a = name+"_o.txt"
    
    with open(tempfile_q, mode = 'w') as f:
        f.write(">query\n{}".format(query))
        
    with open(tempfile_r, mode = 'w') as f:
        f.write(">ref\n{}".format(ref))
        
    cm='stretcher {:s} {:s} -snucleotide2  -gapopen 16  -gapextend 4  {:s}'.format(tempfile_q,tempfile_r,tempfile_a)
    
    #print(cm)
    os.system(cm)
    
    maskfile(tempfile_a)
    
    with open(tempfile_a, mode = 'r') as f:
        read = f.read()
        
    #os.remove(tempfile_r)
    #os.remove(tempfile_q)
    #os.remove(tempfile_a)
        
    alllines=[a.strip() for a in read.split('\n')]
    
    alllinesindex=[a for a in range(len(alllines)) if len(alllines[a])>0 and alllines[a][0]!='#' and alllines[a][0]!=':' and  re.match(r'\S+\s\S+',alllines[a])!=None]
    
    eachline=alllinesindex[::2]
    qlines=[a for a in eachline]
    rlines=[a+2 for a in eachline]
    alines=[a+1 for a in eachline]
    
    
    query=''.join([re.search(r'\s([A-Za-z]|-)+\s*',alllines[a]).group().strip() for a in qlines]).replace("\n","")
    ref=''.join([re.search(r'\s([A-Za-z]|-)+\s*',alllines[a]).group().strip() for a in rlines]).replace("\n","")
    
    return query, ref

def aligntocigar(query, ref):
    
    align = ['']*len(query)
    for i,(q,r) in enumerate(zip(query, ref)):
        
        if q.upper()==r.upper():
            align[i] = 'M'
        elif q == '-':
            align[i] = 'D'
        elif r == '-':
            align[i] = 'I'
        else:
            align[i] = 'X'
            
            
    aligns = [x[0]+x[1] for x in re.findall(r"(.)(\1*)", "".join(align))]
    rposi = 0
    qposi = 0
    cigars = []
    for thealign in aligns:
        
        thetype = thealign[0]
        
        size = len(thealign)
        
        cigars.append("{}{}".format(size, thetype))
        
        if thetype == 'M':
            rposi += size
            qposi += size
        elif thetype == 'X':
            cigars.append("{}{}".format(query[qposi:(qposi+size)],ref[rposi:(rposi+size)]))
            #cigars.append("{}".format(query[qposi:(qposi+size)]))
            rposi += size
            qposi += size
        elif thetype == 'D':
            cigars.append("{}".format(ref[rposi:(rposi+size)]))
            rposi += size
            qposi += size
        elif thetype == 'I':
            cigars.append("{}{}".format(ref[rposi:(rposi+2)],query[qposi:(qposi+size)]))
            #cigars.append("{}".format(query[qposi:(qposi+size)]))
            rposi += size
            qposi += size
            
            
    return "".join(cigars)

def findglobalalign(query, ref, rranges, queryindex):
    
    
    
    query_align, ref_align = runstrecher(query, ref, str(queryindex))
    
    rranges_onalign = getgappedcodi(ref_align, rranges)
    
    qranges = getungappedcodi(query_align,rranges_onalign)
    
    rranges_onalign = [[rranges_onalign[2*i], rranges_onalign[2*i+1]] for i in range(len(rranges_onalign)//2)]
    aligns = [(query_align[x[0]:x[1]], ref_align[x[0]:x[1]]) for x in rranges_onalign]
    
    laligns = len(aligns)
    gaps = []
    cigars = []
    for i,(align_q, align_r ) in enumerate(aligns):
        
        
        cigar = aligntocigar(align_q, align_r)
        
        if i < laligns - 1:
            gap = query[qranges[2*i+1]:qranges[2*i+2]]
            
            if len(gap):
                cigar += "{}I{}".format(len(gap), gap)
                
        cigars.append(cigar)
        
    return qranges[0],qranges[-1],cigars


        



def connectchunks(queryaligns, chunks):
    
    pposi = 0
    qposi = 0
    pposis = []
    qpathstarts = []
    for align in queryaligns:
        
        qpathstarts.append(qposi)
        pposis.append(0)
        qposi += sum([ int(x[:-1]) for x in re.findall(r'\d+[M=XI]',align)])
        #pposi += sum([ int(x[:-1]) for x in re.findall(r'\d+[M=XI=HD]',align)])
    
    allchunkcigars = []
    lastqposi = 0
    for chunk in chunks:
                        
        pathname, pstart, pend, strd, alignindex, chunksize, chunkindex = chunk
        if strd == -1:
            pathsize = int(pathname.split('_')[-1]) - int(pathname.split('_')[-2])
            pstart, pend = pathsize - pend, pathsize - pstart
        
        align = queryaligns[alignindex]
        
        pgap = pstart - pposis[alignindex]
        
        qstart, qend, aligncigar = cigar_cutrrange(align, max(pstart, pposis[alignindex]), pend)
        
        
        qexpect = qpathstarts[alignindex] + qstart
        
        qgap = qexpect  - lastqposi
            
        
        pposis[alignindex] = pend
        lastqposi = (qend - qstart + qexpect)
        
        
        gap = ""
        if pgap < 0:
            gap +="{}H".format(abs(pgap))
        if qgap > 0:
            qgap ="{}I".format(qgap) 
            if len(allchunkcigars):
                allchunkcigars[-1] += qgap
            else:
                gap += qgap
                
        allchunkcigars.append(align[0]+gap+aligncigar)
            
    return allchunkcigars, lastqposi

    

def corr_cigar(allcigars, queryread, chunkread):
    
    
    rposi = 0
    qposi = 0 
    newcigars = []
    for allcigar in allcigars:
        
        cigars = re.findall('\d+[M=XDIH]',allcigar)
        
        newcigar = []
        for cigar in cigars:
            
            thetype = cigar[-1]
            
            lastrposi = rposi
            lastqposi = qposi
            size = int(cigar[:-1])
            
            if thetype == 'M' or thetype =='=':
                rposi += size
                qposi += size
            elif thetype == 'X':
                rposi += size
                qposi += size
            elif thetype == 'D' or thetype =='H':
                rposi += size
                
            elif thetype == 'I':
                qposi += size
                
            if thetype != 'M' and thetype != "=" :
                newcigar.append(cigar)
                continue
        
            qchunk = queryread[lastqposi:qposi].upper()
            rchunk = chunkread[lastrposi:rposi].upper()
            
            
            compare = [int(q==r) for (q,r) in zip(qchunk, rchunk)]
            segments = []
            ifmatch = 1
            lastposi = 0
            for i,ifmatch_new in enumerate(compare):
                
                if ifmatch_new != ifmatch:
                    
                    if i - lastposi:
                        segments.append("{}{}".format(i - lastposi,"XM"[ifmatch]))
                    lastposi = i
                    ifmatch = ifmatch_new
            
            
            segments.append("{}{}".format(i+1 - lastposi,"XM"[ifmatch]))
            
            newcigar.append("".join(segments))
        
        newcigars.append(allcigar[0]+"".join(newcigar))
        
    return newcigars
        
    

def findglobalcigar(queryaligns, chunks, querylength,queryread, chunkread):
    
    
    cigars , qend = connectchunks(queryaligns, chunks)
    
    
    newcigar = "{}I".format(querylength - qend)
    
    if len(cigars) and querylength > qend:
        cigars[-1]+=newcigar
    
    #cigars = corr_cigar(cigars, queryread, chunkread)
    

    return 0, 0, cigars


def getlinearalign(allqueries, linearpath, allchunks, allchunks_onquery, query_reads, chunk_reads, thread, tempfolder):
    
    numchunks = len(linearpath)
    
    chunks_null = ["" for x in range(numchunks)]
    
    for chunkname, order in linearpath.items():
        
        chunkindex = int(chunkname.split("_")[0])-1
        
        chunks_null[order] = ">{}:{}H".format(chunkname,allchunks[chunkindex][2]-allchunks[chunkindex][1])
        
        
    try:
        os.mkdir(tempfolder)
    except:
        pass
        
    p=mul.Pool(processes=thread)
    
    queries = []
    additioninfo = []
    results = []
    allquery_cigars = dict()
    for query, chunks_onquerythis in allchunks_onquery.items():
        
        queryname = allqueries[query]
        queries.append(query)
        
        
        refseqs = []
        pathorders = []
        for index0,querychunk in enumerate(chunks_onquerythis):
            
            chunkindex = int(querychunk[0].split("_")[0]) - 1
            
            refchunk = allchunks[chunkindex]
            
            refsize = refchunk[2] - refchunk[1]
            
            refstrand = refchunk[-1]
            querystrand = querychunk[-1][-1]
            
            strd = 1
            chunkread = chunk_reads[chunkindex]
            if refstrand != querystrand:
                chunkread = makereverse(chunkread)
                strd = -1
                
            refseqs.append(  ( querychunk[1][-2], (index0+1)*querystrand ,chunkindex, refsize, chunkread ) )
            
            pathorders.append( ( querychunk[0], strd ) )
            
        additioninfo.append(pathorders)
        
        refseq = "".join([x[-1] for x in sorted(refseqs)])
        
        refposi = 0
        rranges = [0]
        for ref0 in refseqs:
            refposi += ref0[-2]
            
            rranges.append(refposi)
            rranges.append(refposi)
        rranges = rranges[:-1]
        
        #findglobalalign(query_reads[queryname], refseq, rranges, tempfolder+str(query))
        
        results.append(p.apply_async(findglobalalign, (query_reads[queryname], refseq, rranges, tempfolder+str(query))))
        
        #results.append(findglobalalign(query_reads[queryname], refseq, rranges, tempfolder+str(query)))
        
        
    p.close()
    p.join()
    
    for query,info, result in zip(queries,additioninfo, results):
        
        queryname = allqueries[query]
        
        chunks_onquerythis = allchunks_onquery[query]
        
        start, end, cigars = result.get()
        
        chunkcigars = [x for x in chunks_null]
        
        chunknames = []
        for querychunk, cigar, ( chunkname, strd) in zip(chunks_onquerythis, cigars, info):
            
            strd_s = ">" if strd == 1 else "<"
            
            chunkorder = linearpath[chunkname] 
            
            chunkcigars[chunkorder] = strd_s+ chunkname + ":" + cigar
            
            chunknames.append(strd_s+chunkname)
            
        
        #chunk_cigars = [cigar for cigar in chunk_cigars if sum([int(x[:-1]) for x in re.findall(r'\d+[S=MXD]',cigar ) ]) > 150 ]
        allquery_cigars[queryname] = [start, end, "".join(chunkcigars), "".join(chunknames)]
        
        
    try:
        os.system( "rm -rf {}".format(tempfolder) )
    except:
        pass
        
        
    return allquery_cigars

def getlinearcigar(allqueries, linearpath, allchunks, allchunks_onquery, cigar_reads,query_reads,query_lengths,  chunk_reads, thread, tempfolder):
    
    numchunks = len(linearpath)
    
    chunks_null = ["" for x in range(numchunks)]
    
    for chunkname, order in linearpath.items():
        
        chunkindex = int(chunkname.split("_")[0])-1
        
        chunks_null[order] = ">{}:{}H".format(chunkname,allchunks[chunkindex][2]-allchunks[chunkindex][1])
    
    
    for query, chunks_onquerythis in allchunks_onquery.items():
        
        queryname = allqueries[query]
        
    p=mul.Pool(processes=thread)
    
    queries = []
    additioninfo = []
    results = []
    allquery_cigars = dict()
    for query, chunks_onquerythis in allchunks_onquery.items():
        
        queryname = allqueries[query]
        queries.append(query)
        
        refseqs = []
        chunks = []
        pathorders = []
        for index0,querychunk in enumerate(chunks_onquerythis):
            
            #querychunk is path_start, path_end, queryindex, alignindex, path_strand
            
            chunkindex = int(querychunk[0].split("_")[0]) - 1
            
            
            refchunk = allchunks[chunkindex]
                        
            refsize = refchunk[2] - refchunk[1]
            
            refstrand = refchunk[-1]
            querystrand = querychunk[-1][-1]
            
            
            strd = 1
            chunkread = chunk_reads[chunkindex]
            if refstrand != querystrand:
                chunkread = makereverse(chunkread)
                strd = -1
                
                            
            refseqs.append(  ( querychunk[1][-2], (index0+1)*querystrand ,chunkindex, refsize, chunkread ) )
            
            chunks.append(  refchunk[:-1] + [ querystrand, querychunk[1][3],  refsize, querychunk[0].split("_")[0]])
            
            pathorders.append((querychunk[0], strd))
        
        
        sortindex = sorted(range(len(refseqs)), key = lambda x: refseqs[x])
        
        refseq = "".join([refseqs[index][-1] for index in sortindex])
        chunks = [chunks[index] for index in sortindex ]
        pathorders = [pathorders[index] for index in sortindex]
        
        additioninfo.append(pathorders)
         
        results.append(p.apply_async(findglobalcigar, (cigar_reads[queryname], chunks,query_lengths[queryname], query_reads[queryname], refseq) ) )   
        #results.append(findglobalcigar(cigar_reads[queryname], chunks,query_lengths[queryname], query_reads[queryname], refseq))
    
    p.close()
    p.join()        
            
    for query,info, result in zip(queries,additioninfo, results):
        
        queryname = allqueries[query]
        
        chunks_onquerythis = allchunks_onquery[query]
        
        start, end, cigars = result.get()
        
        
        chunkcigars = [x for x in chunks_null]
        
        chunknames = []
        for querychunk, cigar, ( chunkname, strd) in zip(chunks_onquerythis, cigars, info):
            
            strd_s = ">" if strd == 1 else "<"
            
            chunkorder = linearpath[chunkname] 
            
            chunkcigars[chunkorder] = strd_s+ chunkname + ":" + cigar[1:]
            
            chunknames.append(strd_s+chunkname)
            
            
        #chunk_cigars = [cigar for cigar in chunk_cigars if sum([int(x[:-1]) for x in re.findall(r'\d+[S=MXD]',cigar ) ]) > 150 ]
        allquery_cigars[queryname] = [start, end, "".join(chunkcigars), "".join(chunknames)]
        
                
        
    return allquery_cigars

def graphtolinear(inputfile, queryfile, reffile, output, thread = 8, tempfolder = "", polish = 0):
    
    allpathes, allqueries = readalignments(inputfile)
    
    
    refsizes = getrefsizes(reffile)
    
    #allbreaks: a dict of list of break information, each piece of break information: [breakstart, breakend, index of query, index of alignment, strand]
    allpathaligns = getpathaligns(allpathes, refsizes)
    
    #distract chunks of the linear path from breaks informatoin
    allchunks, allchunks_onquery = getchunks(allpathaligns)
    
    
    #get a global linear path based on chunk information of all genes.
    linearpath = getlinearpath(allchunks_onquery)
        
    chunk_reads = getchunkseq(reffile, allchunks)
    
    if polish:
        
        query_reads, query_lengths, query_order = getqueryseq(queryfile)
        
        if len(tempfolder) == 0:
            tempfolder = output + datetime.datetime.now().strftime("_stretchertemp_%y%m%d_%H%M%S/")

        allquery_cigars = getlinearalign(allqueries, linearpath, allchunks, allchunks_onquery, query_reads, chunk_reads, thread, tempfolder)
    else:
        
        cigar_reads, query_lengths, query_order = getquerycigar(queryfile , inputfile)
        
        query_reads = cl.defaultdict(str)
        allquery_cigars = getlinearcigar(allqueries, linearpath, allchunks, allchunks_onquery, cigar_reads, query_reads,query_lengths,chunk_reads, thread, tempfolder)
        
    
    allchunks_orders = [[] for x in range(len(allchunks))]
    for chunkname, chunkorder in linearpath.items():
        allchunks_orders[int(chunkname.split('_')[0]) - 1].append(chunkorder)
        
    allchunks_withorders = [(chunk, chunk_read, order, index) for index, (chunk, chunk_read, order) in enumerate( zip(allchunks, chunk_reads, allchunks_orders) )]
    
    allchunks_withorders = sorted(allchunks_withorders, key = lambda x: min(x[2]) if len(x[2]) else len(allchunks) + 1)
    
    putchunks(output, allchunks_withorders)
    
    putalignments(output, allquery_cigars, query_lengths, query_order)
    
    
    
def main(args):
    
    graphtolinear(args.input, args.query, args.ref, args.output, args.threads, args.tempfolder, args.polish)
    
def run():
    """
        Parse arguments and run
    """
    parser = argparse.ArgumentParser(description="program distract overlap genes and alignments on contigs")
    
    parser.add_argument("-i", "--input", help="path to output file", dest="input", type=str,required=True)
    parser.add_argument("-q", "--query", help="path to output file", dest="query", type=str,required=True)
    parser.add_argument("-r", "--ref", help="path to output file", dest="ref", type=str,required=True)
    parser.add_argument("-o", "--output", help="path to output file", dest="output", type=str,required=True)
    parser.add_argument("-t", "--threads", help="path to output file", dest="threads", type=int,default = 1)
    parser.add_argument("-f", "--folder", help="path to output file", dest="tempfolder", type=str,default = "")
    parser.add_argument("-p", "--polish", help="path to output file", dest="polish", type=int ,default = 0)
    
    parser.set_defaults(func=main)
    args = parser.parse_args()
    args.func(args)
    
    
if __name__ == "__main__":
    run()
