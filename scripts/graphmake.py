#!/usr/bin/env python3

import multiprocessing as mul
import argparse
import os
import subprocess
import datetime
from graphfilesplit import filesplit, record_filename
from graphinserts import insertdistract
from graphcigar import graphcigar
from graphcleanrepeats import cleanrepeats
from graphtolinear import graphtolinear
import re

script_folder = os.path.dirname(os.path.abspath(__file__))
lock1 = mul.Lock()

ifblast_g = 1

def addgraphglobalcoords(inputfile, graphfile):
    
    def parse_input_coord(coord_text):
        coord_text = coord_text.strip()
        strand = coord_text[-1] if len(coord_text) and coord_text[-1] in "+-" else "+"
        core = coord_text[:-1] if len(coord_text) and coord_text[-1] in "+-" else coord_text
        
        m = re.match(r"^(.*):(\d+)-(\d+)$", core)
        if m is None:
            raise ValueError(f"Cannot parse input coordinate field: {coord_text}")
        
        contig = m.group(1)
        globalstart = int(m.group(2))
        globalend = int(m.group(3))
        return contig, globalstart, globalend, strand

    def parse_graph_header_name(header_name):
        m = re.match(r"^(.*)_(\d+)_(\d+)$", header_name)
        if m is None:
            raise ValueError(f"Cannot parse graph header: {header_name}")
        
        sequencename = m.group(1)
        localstart = int(m.group(2))
        localend = int(m.group(3))
        return sequencename, localstart, localend

    def resolve_seqname(graph_seqname, inputinfo):
        if graph_seqname in inputinfo:
            return graph_seqname
        
        parts = graph_seqname.split("_")
        
        # graph header may truncate one or more leading "_" fields
        for i in range(1, len(parts)):
            candidate = "_".join(parts[i:])
            if candidate in inputinfo:
                return candidate
        
        return None

    inputinfo = {}
    with open(inputfile, mode="r") as r:
        for line in r:
            if len(line) == 0 or line[0] != ">":
                continue
            
            fields = line.strip().split()
            if len(fields) < 2:
                continue
            
            seqname = fields[0][1:]
            contig, globalstart, globalend, strand = parse_input_coord(fields[1])
            inputinfo[seqname] = (contig, globalstart, globalend, strand)

    tempfile = graphfile + ".addglobalcoords.temp"
    
    with open(graphfile, mode="r") as r:
        with open(tempfile, mode="w") as w:
            for line in r:
                if len(line) == 0 or line[0] != ">":
                    w.write(line)
                    continue
                
                fields = line.strip().split()
                header_name = fields[0][1:]
                
                graph_seqname, localstart, localend = parse_graph_header_name(header_name)
                matched_name = resolve_seqname(graph_seqname, inputinfo)
                
                if matched_name is None:
                    raise KeyError(f"Cannot find matching input header for graph sequence: {graph_seqname}")
                
                contig, globalstart, globalend, strand = inputinfo[matched_name]
                seqlen = globalend - globalstart
                
                if localstart < 0 or localend < localstart or localend > seqlen:
                    raise ValueError(
                        f"Local coordinates out of range for {header_name}: "
                        f"{localstart}-{localend} not within 0-{seqlen}"
                    )
                
                if strand == "+":
                    graph_globalstart = globalstart + localstart
                    graph_globalend = globalstart + localend
                else:
                    graph_globalstart = globalend - localend
                    graph_globalend = globalend - localstart
                
                graph_globalcoord = f"{contig}:{graph_globalstart}-{graph_globalend}{strand}"
                
                newfields = [fields[0], graph_globalcoord] + fields[1:]
                w.write("\t".join(newfields) + "\n")
    
    os.replace(tempfile, graphfile)
    return graphfile

def checkmask(fasta_file):
    
    bash_cmd = f'''
    grep -v '^>' "{fasta_file}" | tr -d '\\n' | tr -cd 'a-z' | wc -c
    '''
    result = subprocess.run(["bash", "-c", bash_cmd], capture_output=True, text=True)
    lowercase_count = int(result.stdout.strip())
    
    return lowercase_count

def fileordered(allfiles, priors=None):
    priors = priors or set()
    reference_files = []
    query_files = []
    seen = set()
    for path in allfiles:
        if path in seen:
            continue
        seen.add(path)
        fn = os.path.basename(path)
        stem = fn[:-3] if fn.endswith(".fa") else os.path.splitext(fn)[0]
        (reference_files if stem in priors else query_files).append(path)
    return reference_files + query_files


def input_file_order(input_fasta, folder):
    ordered = []
    with open(input_fasta, mode="r") as handle:
        for line in handle:
            if line.startswith(">"):
                path = os.path.join(folder, record_filename(line[1:].strip()))
                if os.path.isfile(path):
                    ordered.append(path)
    return ordered

def selfalign(input,output,nthreads, ifblast = 0):
    
    if ifblast:
        dbcmd = "bash {}/runmakeblastdb -in {} -out {}_db ".format(script_folder,input, input)
        os.system(dbcmd)
        cmd = "bash {}/runblastn  -task megablast  -query {} -db {}_db -gapopen 10 -gapextend 2 -word_size 30  -perc_identity 95 -dust yes -lcase_masking -evalue 1e-200 -outfmt 17 -out {}_selfblast.out  -num_threads {} -max_target_seqs 100 ".format( script_folder,input , input, input,nthreads)
        os.system(cmd)
        totalclean = cleanrepeats(input, input+"_", nthreads)
        
        cmd = "{}/runwinnowmap.sh {} {} {} 0.95 150  > {}_selfblast.out".format(script_folder,  input+"_",input+"_", nthreads, input+"_")
        os.system(cmd)
        totalclean += cleanrepeats(input+"_", output, nthreads)
    else:  
        cmd = "{}/runwinnowmap.sh {} {} {} 0.95 150 > {}_selfblast.out".format(script_folder,  input,input, nthreads, input)
        os.system(cmd)
        
        totalclean = cleanrepeats(input, output, nthreads)

    return totalclean

def callblastn(query, graphfile,  output,nthreads, simi = 0.90, ifrepeat = 1):
    
    if ifrepeat:
        repeatopts = "-dust yes -lcase_masking"
    else:
        repeatopts = ""
        
    cmd = "bash {}/runblastn -task megablast -query {} -db {} -gapopen 10 -gapextend 2 -word_size 30   -perc_identity {} {} -evalue 1e-200  -outfmt 17 -out {}  -num_threads {} -max_target_seqs 100 ".format( script_folder, query,graphfile+"_db" , int(simi * 100),repeatopts,output,nthreads)
    os.system(cmd)
    
    
    
def callalign(query, graphfile,  output,nthreads, simi = 0.90,  ifhighqual = 1, ifrepeat = 1):
    if ifhighqual:
        callblastn(query, graphfile,  output,nthreads, simi, ifrepeat)
        
        
        #if ifrepeat:
            #repeatopts = "-dust yes -lcase_masking"
        #else:
            #repeatopts = ""
        
        #cmd = "bash {}/runblastn -task megablast -query {} -db {} -gapopen 10 -gapextend 2 -word_size 30   -perc_identity {} {} -evalue 1e-200  -outfmt 17 -out {}  -num_threads {} -max_target_seqs 100 ".format( script_folder, query,graphfile+"_db" , int(simi * 100),repeatopts,output,nthreads)
        #os.system(cmd)
        
        if (not os.path.exists(output) or os.path.getsize(output) <= 10):
            
            print(f"WARNING: alignment fail :{output}, realigning\n")
            
            callblastn(query, graphfile,  output,nthreads, simi, ifrepeat)
            
            #fallback_cmd = "bash {}/runblastn -task megablast -query {} -db {} -gapopen 10 -gapextend 2 -word_size 30 -perc_identity {} {} -evalue 1e-200 -outfmt 17 -out {} -num_threads {} -max_target_seqs 100".format(script_folder, query, graphfile+"_db", int(simi * 100), repeatopts, output, nthreads)
            
            #os.system(fallback_cmd)
            
            
        if ifrepeat:
            
            cmd = "{}/runwinnowmap.sh {} {} {} {} 150 > {}".format(script_folder,  query,graphfile, nthreads, simi, output+"_wm")
            os.system(cmd)
            if (not os.path.exists(output+"_wm") or os.path.getsize(output+"_wm") <= 10):
                os.system(cmd)
            os.system(f" ( cat  {output}_wm  >> {output} &&  rm {output}_wm  ) || true ")
            
    else:
        #cmd = "{}/runfast.sh {} {} {} 0.95 300 -c > {}".format(script_folder,  query,dbpath, nthreads,output)
        cmd = "{}/runwinnowmap.sh {} {} {} {} 150 > {}".format(script_folder,  query,graphfile, nthreads, simi, output)
        os.system(cmd)
        
        
def addnewseq(folder, graphfile, addfile,nthreads, simi, ifblast = 1):
    
    filealign = addfile + "_blastn.out"
    fileinserts = addfile + "_inserts.fasta"
    
    try:
        os.remove(fileinserts)
    except:
        pass
        
    if ifblast:
        callblastn( addfile,  folder+graphfile.split("/")[-1], filealign,nthreads, simi, 1)
    else:
        callalign( addfile,  folder+graphfile.split("/")[-1], filealign,nthreads, simi, 0)
        
    insertfiles = insertdistract( filealign, addfile, fileinserts, 1, 150)
    
    if len(insertfiles) and os.path.isfile(fileinserts) != False and os.path.getsize(fileinserts) > 10:
        os.system("cat {} >> {}".format(fileinserts, graphfile))
        
    else:
        return ""
    
    return insertfiles

def runcleanrepeats(folder, graphfile, newgraphfile,nthreads, ifhighqual = 1):


    size = 100000000
    
    ncycle = 1
    if ifhighqual:
        ncycle = 10
   
    theinput = graphfile
    theoutput = newgraphfile
    thetemp  = newgraphfile + ".temp"

    icycle = 0
    while icycle < ncycle and size > 300:
        if icycle > 0:
            ifhighqual = 0
            theinput = theoutput
            theoutput = newgraphfile if theoutput  == thetemp else thetemp
        size = selfalign(theinput, theoutput, nthreads, ifhighqual)
        icycle += 1

    if theoutput == thetemp:
        os.replace(thetemp, newgraphfile)
    elif icycle > 1:
        os.remove(thetemp)


    return newgraphfile

def editname(thefile, graphfile):
    
    size = 0
    with open(thefile, mode = 'r') as r:
        for line in r:
            if len(line) and line[0] != '>':
                size += len(line.strip())
                
                
    with open(thefile, mode = 'r') as r:
        with open(graphfile, mode = 'w') as w:
            for line in r:
                if len(line) and line[0] == '>':
                    
                    header = line.strip().split()
                    header[0] += "_{}_{}".format(0, size)
                    
                    if len(header[0]) > 50 -2:
                        header[0] = ">"+"_".join(header[0].split("_")[2:])
                    w.write("\t".join(header)+"\n")
                else:
                    w.write(line)
                    
def creategraph(folder, graphfile,nthreads, ifhighqual = 1, prior = "", ordered_files = None):
    candidates = ordered_files
    if candidates is None:
        candidates = [folder + "/" + file for file in os.listdir(folder) if file.endswith(".fa")]
    allfiles = fileordered(candidates, set([x for x in prior.split(",") if len(x)]))

    if len(allfiles) == 0:
            return
    
    editname(allfiles[0],graphfile) 
    
    graphfile_filename = graphfile.split("/")[-1]
    
    dbpath = folder+graphfile_filename+"_db"
    
    if ifhighqual:
        dbcmd = "bash {}/runmakeblastdb -in {} -out {}".format(script_folder,graphfile, dbpath)
    else:
        dbcmd = "ln -f {} {}".format(graphfile, dbpath)
        
    os.system(dbcmd)
    
    for addfile in allfiles[1:]:
        
        insertfiles = addnewseq(folder, graphfile, addfile,nthreads, 0.90, ifhighqual)
        if len(insertfiles):
            os.system(dbcmd)
            
            
def run_align(graphfile, thefile, graphalign, simi, ifhighqual = 1):
    
    lowcasecount = checkmask(thefile)
    
    ifrepeat = 0
    if lowcasecount > 5000:
        ifrepeat = 1
        
    output = thefile + "_align.out"
    
    callalign(thefile,graphfile, output, 1, simi, ifhighqual, ifrepeat)
    queryname, fullpath, fullcigar, pathranges,qranges = graphcigar(graphfile, thefile , output)
    
    lock1.acquire()
    with open(graphalign, mode = 'a') as w:
        w.write("{}\t{}\t{}\t{}\t{}\n".format(queryname, fullpath, fullcigar, pathranges,qranges))
    lock1.release()
    
    
def run_cigar(graphfile, thefile, dbpath):
    
    output = thefile + "_align.out"
    
    queryname, fullpath, fullcigar, pathranges,qranges = graphcigar(graphfile, thefile , output)
    
    return [queryname, fullpath, fullcigar, pathranges,qranges]

def aligngraph(folder, graphfile, graphalign,nthreads,ifblast = 0, ifunfinish = 0, ordered_files = None):
    
     
    finished = set()
    if os.path.isfile(graphalign):
    
        if ifunfinish:
            iferror = 0
            with open(graphalign, mode = 'r') as f:
                for line in f:
                    line = line.split()
                    if len(line) == 5:
                        name = line[0]
                        finished.add(name+".fa")
                    elif len(line) != 5:
                        iferror = 1
            if iferror:
                with open(graphalign, mode = 'r') as f, open(graphalign+".temp", mode = 'w') as w:
                    for line in f:
                        if len(line.split()) != 5:
                            continue
                        else:
                            w.write(line)
                os.system("mv {} {}".format(graphalign+".temp", graphalign))
        else:
            os.remove(graphalign)

    graphfile_filename = graphfile.split("/")[-1]
    
    dbpath = folder+graphfile_filename+"_db"

    if ifblast:
        dbcmd = "bash {}/runmakeblastdb -in {} -out {} ".format(script_folder,graphfile, dbpath)
        os.system(dbcmd)
    
    dbcmd = "ln -f {} {}".format(graphfile, folder+graphfile_filename)
    os.system(dbcmd)
    
    candidates = ordered_files
    if candidates is None:
        candidates = [folder + "/" + file for file in os.listdir(folder) if file.endswith(".fa")]
    allfiles = [
        path for path in candidates
        if os.path.isfile(path) and os.path.basename(path) not in finished
    ]
    
    
    with mul.Pool(processes=nthreads) as pool:
        
        pool.starmap(run_align, [(folder+graphfile_filename, thefile, graphalign, 0.9,ifblast) for thefile in allfiles])
        
        
        
        
def main(args):
    
    if len(args.dir)==0:
        folder = args.input +  datetime.datetime.now().strftime("%y%m%d_%H%M%S/")
        #folder = "./tempx/"
        
        folder = folder.replace("%","_")
        folder2 = folder + "/graphtemp/"
        
        try:
            os.system("rm -rf {}".format(folder))
            os.mkdir(folder)
        except:
            pass
            
        try:
            os.system("rm -rf {}".format(folder2))
            os.mkdir(folder2)
        except:
            pass
            
    else:
        folder = args.dir
        folder2 = folder + "/graphtemp/"
        try:
            os.system("rm -rf {}".format(folder2))
            os.mkdir(folder2)
        except:
            pass
            
    if 1 == 1:   
        inputfile_name = args.input.split("/")[-1]
        
        graphfile_raw = folder + inputfile_name+"_graphwithrep.FA"
        graphfile = args.input +"_graph.FA"
        graphalign = args.input + "_allgraphalign.out"
        graphlinear = args.input + "_lineargraph.gaf"
        
        if args.split:
            ordered_files = filesplit(args.input, folder, 0)
        else:
            ordered_files = input_file_order(args.input, folder)

        ordered_files = fileordered(
            ordered_files,
            set([name for name in args.prior.split(",") if name]),
        )
            
        if args.create and (args.unfinish == 0 or not os.path.isfile(graphfile_raw)):
            creategraph(folder, graphfile_raw,args.thread, args.ifhighqual, args.prior, ordered_files)
            
        if args.fine and (args.unfinish == 0 or not os.path.isfile(graphfile)):
            runcleanrepeats(folder2, graphfile_raw, graphfile,args.thread, args.ifhighqual)
            addgraphglobalcoords(args.input, graphfile)

        if args.align:
            aligngraph(folder, graphfile, graphalign,args.thread, args.ifhighqual, args.unfinish, ordered_files)
            
        if args.linear:
            graphtolinear(graphalign, args.input, graphfile, graphlinear, args.thread, folder+"stretcherouts/", args.globalalign)
            
        if args.dir == "":
            pass
            #os.system("rm -rf {} || true".format(folder))
            
            
def run():
    """
            Parse arguments and run
    """
    parser = argparse.ArgumentParser(description="program distract overlap genes and alignments on contigs")
    
    parser.add_argument("-i", "--input", help="path to output file", dest="input", type=str,required=True)
    parser.add_argument("-c", "--create", help="if generate fine graph", dest="create", type=int,default = 1)
    parser.add_argument("-f", "--fine", help="if generate fine graph", dest="fine", type=int,default = 1)
    parser.add_argument("-a", "--align", help="if align to graph", dest="align", type=int,default = 1)
    parser.add_argument("-l", "--linear", help="if align to graph", dest="linear", type=int,default = 1)
    parser.add_argument("-t", "--thread", help="if align to graph", dest="thread", type=int,default = 1)
    parser.add_argument("-s", "--split", help="if align to graph", dest="split", type=int,default = 1)
    parser.add_argument("-d", "--dir", help="if align to graph", dest="dir", type=str,default = "")
    parser.add_argument("-g", "--globalalign", help="if align to graph", dest="globalalign", type=int,default = 0)
    parser.add_argument("-q", "--highqual", help="if align to graph", dest="ifhighqual", type=int,default = 1)
    parser.add_argument("-u", "--unfinish", help="if restart", dest="unfinish", type=int,default = 1)
    parser.add_argument("-p", "--prior", help="comma-separated reference record names to place first", dest="prior", type=str,default = "")
    parser.set_defaults(func=main)
    args = parser.parse_args()
    args.func(args)
    
    
if __name__ == "__main__":
    run()
