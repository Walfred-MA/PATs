PAT — Pangenome-Allele Annotation Toolkit
========================================

This manual describes how to build an allele-based pangenome database for
both functional annotation and NGS genotyping. The pipeline consumes (1)
haplotype-resolved assemblies and reference genomes, and (2) a BED file of
target genes (optionally with exon targets). Major stages are: input
preparation, masking, target extraction, and running the Snakemake workflow.


0) Requirements (high level)
----------------------------
- Linux environment with bash, Python 3, samtools, and Snakemake
- Compilers/tools for C++ utilities in `snakemake/` (see that folder’s README)
- PAT repository layout with subfolders:
  masking/   tools/   snakemake/   scripts/


1) Prepare your input files
---------------------------

A. Create the assemblies/reference path list (TSV)
   ------------------------------------------------
   Make a **two-column, no-header** TSV file listing all samples and references.
   Column 1: `sampleName_h{1|2}`   Column 2: absolute or relative FASTA path

   Example (`query_pathes_withrefs.txt`):
     apr001_h1    Assemblies/apr001.1.polished.fa
     apr001_h2    Assemblies/apr001.2.polished.fa
     CHM13_h1     References/CHM13.v2.0.fa
     HG38_h1      References/hg38.fa

   Notes:
   - Use concise sample names (≤ ~10 letters recommended).
   - `_h1` and `_h2` denote haplotypes; when known, prefer `_h1`=paternal, `_h2`=maternal.
   - For **references**, use only `_h1` (e.g., `CHM13_h1`, `HG38_h1`).
   - If you need to distinguish contigs by reference (recommended when mixing
     CHM13 and hg38), use a prefix convention in contig names, e.g.:
       - CHM13: `NC_0609…` (e.g., `NC_060925.1` for chr1)
       - hg38 : `chr…`    (e.g., `chr1`)
     Make sure these prefixes match your FASTA headers.

B. Repeat-mask all references and assemblies
   -----------------------------------------
   Use the scripts in `masking/` (see `masking/readme.txt` for details).
   This step produces masked FASTAs needed by downstream steps.

C. Define your target regions (genes) as a BED
   -------------------------------------------
   Use the helper script to export per-gene BED from a GFF3:

     python tools/gff_toGeneBed.py -i <genes.gff3> -o genes.bed

   Optional arguments:
     -g X      Select only specific genes (file with one per line or
               comma-list). Supports prefixes like `SMN*`.
     -e        Add merged exon intervals as `Target=...` in column 5
               (absolute, 0-based, half-open).
     -a N      Add N bp anchors upstream/downstream to gene BED start/end.

   Help:
     python tools/gff_toGeneBed.py -h

   BED formats accepted by the pipeline:
     4-col (gene-level)
       chrom   start   end     gene_or_group

     5-col (gene with exon targets; recommended for gene-based studies)
       chrom   start   end     gene_or_group   Target=s1-e1,s2-e2,...
     Example:
       chr5_GL339449v2_alt  456848  485731  SMN2  Target=456848-457609,457774-458352

D. Run the Snakemake workflow
-----------------------------
- See `pipeline/` for the workflow, required C++ utilities, and detailed
  run instructions. Compile the provided C++ tools as directed in that
  folder’s README.

Typical invocation (example):
  snakemake --cores 32 --use-conda -s snakemake/Snakefile

SLURM example (adjust to your site):
  snakemake --profile slurm -s snakemake/Snakefile

(If using the JSON config described below, the Snakemake profile or
wrapper scripts may consume it to set paths and partitioning.)


3) JSON configuration (example + field descriptions)
-------------------------------------------------------------
Many wrappers expect a JSON file to centralize common settings:

{
  "slurm": " --account=<acct> --time 50:00:00 --partition=qcb ",
  "QueryPath": "query_pathes_withrefs.txt",
  "ScriptFolder": "/path/to/PAT/scripts/",
  "TargetFolder": "groups",
  "TempFolder": "snaketemp",
  "genelist": "./genes.bed",
  "ReferencePrefix": "NC_0609",
  "NumPartitions": 1
}

Field descriptions:
- slurm            SLURM submission arguments passed to jobs (account, time, partition, etc.).
- QueryPath        Path to the two-column TSV listing all assemblies/references (see 1A).
- ScriptFolder     Directory containing your pipeline scripts.
- TargetFolder     Root output directory for pipeline results.
- TempFolder       Directory used for intermediate/temporary files.
- genelist         BED file of target genes (4- or 5-column as described in 1C).
- ReferencePrefix  Prefix used to identify reference contigs, e.g. `"chr"` (hg38) or `"NC_0609"` (CHM13).
                   Must be a prefix of chromosome names in your reference FASTA headers.
- NumPartitions    Number of partitions to split the run into. As a guideline, aim for ~1000 genes
                   per partition (e.g., 3000 genes → 3 partitions).


4) Tips & gotchas
-----------------
- Ensure your FASTA indices (`.fai`) are present and correspond to the masked FASTAs.
- Keep naming consistent across:
    (a) TSV sample IDs (e.g., `apr001_h1`),
    (b) reference IDs (e.g., `CHM13_h1`),
    (c) contig/chromosome prefixes in FASTA headers (`NC_0609…`, `chr…`).
- For exon-centric analyses, prefer the 5-column BED with `Target=` exon blocks.
- If mixing references (CHM13 + hg38), using distinct prefixes prevents accidental
  cross-mapping of contig names.

Questions / Issues
------------------
Please see the READMEs inside `masking/` and `snakemake/`, and run:
  python tools/gff_toGeneBed.py -h
for BED generation options. If problems persist, check input paths and
FASTA indices first.
