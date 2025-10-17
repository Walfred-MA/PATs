#!/usr/bin/env python3

import os
import subprocess
import sys
from pathlib import Path

def run_or_fail(command):
    try:
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Command failed: {command}")
        print(f"[ERROR] Exit code: {e.returncode}")
        if e.returncode == -9:
            print("[OOM] Likely killed by out-of-memory (SIGKILL).")
        sys.exit(1)
        
def main():
    rerun = 1
    inputfile = sys.argv[1]
    rerun = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    kmerfile = inputfile + "_kmer.list"
    normfile = inputfile + "_norm.gz"
    treefile = inputfile + "_tree.fa"
    newkmerfile = inputfile + "_tree.fa_kmer.list"
    matrixfile = inputfile + "_kmatrix.txt"
    
    print("running: " + inputfile)
    
    folder = os.path.dirname(os.path.abspath(inputfile))
    script_folder = os.path.dirname(os.path.abspath(__file__))
   
    folders = str(folder).split("/")

    graphfile = "{}/{}{}{}_samples.fasta_fixed.fa_loci.txt.fasta".format("/".join(folders[:-2]),folders[-5].split("_")[0], folders[-3], "_".join(folders[-5].split("_")[1:]))
    graphinfo = "{}/{}{}{}_samples.fasta_fixed.fa_loci.txt".format("/".join(folders[:-2]),folders[-5].split("_")[0], folders[-3], "_".join(folders[-5].split("_")[1:]))


    # Step 1: skip if matrix file exists and is large enough
    if not rerun and (os.path.isfile(matrixfile) and os.path.getsize(matrixfile) > 100):
        return
    
    # Step 2: run kmernorm
    if rerun or not (os.path.isfile(normfile) and os.path.getsize(normfile) > 100):
        run_or_fail(f"python {script_folder}/querylist_filterbykmer.py -c 1000 -r 0.10 -i {inputfile} -k {kmerfile} -o {inputfile}_filed.fa -l {inputfile}_filed.fa_kmer.list -p 1 > {inputfile}_filed.fa_info.txt")
        run_or_fail(f"{script_folder}/kmernorm -i {inputfile}_filed.fa -k {inputfile}_filed.fa_kmer.list -o {inputfile}_filed.fa_norm.gz  -w 1 -h 1")
        
        inputfile = inputfile + "_filed.fa"
        kmerfile = inputfile + "_kmer.list"
        normfile = inputfile + "_norm.gz"
        
    # Step 4: run kmertree
    if rerun or not (os.path.isfile(treefile) and os.path.getsize(treefile) > 100):
        run_or_fail(f"{script_folder}/kmertree -i {inputfile} -n {normfile} -o {treefile}")
    
    # Step 5: run kmerannotate
    if rerun or not (os.path.isfile(newkmerfile) and os.path.getsize(newkmerfile) > 100):
        if os.path.isfile(graphfile+"_allgraphalign.out") and os.path.getsize(graphfile+"_allgraphalign.out") > 100:
            print("find graphfile: "+graphfile)
            run_or_fail(f"python {script_folder}/kmerannotate.py -k {kmerfile} -s {graphfile} -g {graphfile}_graph.FA -a {graphfile}_allgraphalign.out -o {newkmerfile}")
            run_or_fail(f"bash {script_folder}/headeraddinfor.sh {newkmerfile} {graphinfo}")

        else:
            run_or_fail(f"python {script_folder}/kmerannotate.py -k {kmerfile} -o {newkmerfile}")
            
    # Step 6: compile matrix
    if rerun or not (os.path.isfile(matrixfile) and os.path.getsize(matrixfile) > 100):
        run_or_fail(f"python {script_folder}/matrixcompile.py -s {treefile} -k {newkmerfile} -o {matrixfile}")
        
if __name__ == "__main__":
    main()
