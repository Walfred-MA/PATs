PAT — Pangenome-Allele Annotation Toolkit
========================================
This manual describes how to build an allele-based pangenome database for
both functional annotation and NGS genotyping. The pipeline consumes:
(1) haplotype-resolved assemblies and reference genomes
(2) a BED file of target genes (optionally with exon targets)

Major stages: input preparation, masking, target extraction, and running
the Snakemake workflow.

Requirements (high level)
----------------------------
- Linux environment with bash, Python 3.10, samtools, and Snakemake 6.5.0
- Compilers/tools for C++ utilities in `snakemake/` (see that folder’s README)
- PAT repository layout with subfolders:
  masking/  tools/  snakemake/  scripts/

A. Installing via install/install.py
-----------------------------------------------------

We provide an automatic script to help set up the required binaries and Python dependencies for the PAT toolkit.

Requirements:
- Python 3.10
- Snakemake 6.5.0 with PuLP 2.7.0 (installed automatically from `requirements.txt`)
- A valid Conda environment (recommended)
- C++ compiler (g++ ≥ 8)
- Internet access to install dependencies via conda

How to Use:

1. Create and activate a Python 3.10 environment (if not already done):
   conda create -n patenv python=3.10 -y
   conda activate patenv

2. Install required dependencies:
   conda install -y htslib=1.21 eigen zlib

3. Run the install script from the root of the PAT repository:
   python install/install.py --prefix /path/to/PAT/scripts/

   Replace /path/to/PAT/scripts/ with the absolute path to your ScriptFolder. This should match what you set in the Snakemake config JSON.

   The script will:
   - Compile all C++ binaries (kmercounter8, kmer_selector, etc.)
   - Use the active conda environment's include/ and lib/ directories for htslib/sam.h and zlib.h
   - Copy compiled binaries into the provided --prefix folder

4. You should see output like:
   [BUILD] kmercounter8 not found, attempting to build from /.../src/kmercounter
   [DONE] Copied kmercounter8 → /.../scripts/kmercounter8

Troubleshooting:

- If you see errors like:
     fatal error: zlib.h: No such file or directory
  or
     fatal error: htslib/sam.h: No such file or directory

  then your environment may be missing required headers. Try:
     conda install -y htslib=1.21 zlib

B. Create the assemblies/reference path list (TSV)
--------------------------------------------------
Make a two-column, no-header TSV listing all samples and references.
Column 1: sampleName_h{1|2}
Column 2: absolute or relative FASTA path

Example (`query_pathes_withrefs.txt`):
```
apr001_h1    Assemblies/apr001.1.polished.fa
apr001_h2    Assemblies/apr001.2.polished.fa
CHM13_h1     References/CHM13.v2.0.fa
HG38_h1      References/hg38.fa
```
Notes:
- Use concise sample names (≤10 letters recommended).
- `_h1` and `_h2` denote haplotypes; `_h1` is typically paternal.
- For references, only use `_h1` (e.g., `CHM13_h1`, `HG38_h1`).
- If mixing CHM13 and hg38, distinguish contigs using FASTA header prefixes:
    - CHM13: NC_0609... (e.g., NC_060925.1 for chr1)
    - hg38 : chr... (e.g., chr1)

Reference-name resolution during graph construction:
- A contig name containing `#` uses its explicit `sample#haplotype#contig`
  identity. For example, `CHM13#1#chr1` resolves to query-table sample
  `CHM13_h1`; the `chr1` suffix does not make it hg38.
- Implicit CHM13/hg38 inference is used only when the contig name contains no
  `#` characters.
- Without `#`, a contig containing `NC_0609` resolves to `CHM13_h1`.
  Therefore a `CHM13_h1` reference using implicit inference must retain its
  `NC_0609...` chromosome names.
- Without `#`, all other contig names—including `chr1`, `chr2`, and other
  `chr*` names—resolve to `HG38_h1`.

Graph construction order:
- Reference records identified by `graphmake.py -p ref1,ref2,...` are placed
  first, retaining their relative order from the input multi-FASTA.
- Every non-reference record follows in its original multi-FASTA order. There
  is no sorting by assembly name, sequence size, masking, or filename.
- The PATs workflow supplies `-p` automatically from the first query-table
  reference. When running `graphmake.py` directly without `-p`, put the
  reference records first in the input multi-FASTA.

Make sure these prefixes match your FASTA headers.
Example file can be seen in pipeline/query_pathes.txt_example.txt. 

C. Repeat-mask all references and assemblies
--------------------------------------------
Use the scripts in `masking/` (see `masking/readme.txt` for details).
This step generates the masked FASTAs used downstream. You may skip if you have already done so

D. Define your target regions (genes) as a BED file
---------------------------------------------------
Use the helper script to extract per-gene BED from a GFF3 (gene annotation file, for example, some gff3 files can be found: https://www.gencodegenes.org/):
```
python tools/gff_toGeneBed.py -i <genes.gff3> -o genes.bed

Optional arguments:
  -g X   Select specific genes (file or comma-list). Supports wildcards like SMN*.
  -e     Add merged exon intervals as `Target=...` in column 5.
  -a N   Add N bp flanking anchors to each gene region.
```
Help:
```
python tools/gff_toGeneBed.py -h
```
Accepted BED formats:
  4-column (gene-level):
  ```
    chrom  start  end  gene_or_group
  ```
  5-column (recommended):
  ```
    chrom  start  end  gene_or_group  Target=s1-e1,s2-e2,...
  ```

Example:
```
chr5_GL339449v2_alt 456848 485731 SMN2 Target=456848-457609,457774-458352
```
Example output format can be seen in pipeline/SMNexample.bed. 

E. Run the Snakemake workflow
-----------------------------
See `pipeline/` for workflow logic, C++ utilities, and instructions.
Compile all C++ tools as directed in that folder’s README.

There are two snakemake files included, Snakefile and Snakefile_light. Snakefile is used for a large number of genes, while Snakefile_light is a lighter version that only works for a small number of genes and requires fewer resources. 

Run example:
```
snakemake -k  --cluster "sbatch --account=mchaisso_100 --partition=qcb --time=500:00:00 {resources.slurm_extra}" --default-resources "mem_mb=3000" --jobs 500  --rerun-incomplete  --notemp --latency-wait 100 --resources mem_gb=1000
```

It requires a JSON config, the Snakemake profile or wrappers will load it.

2) JSON configuration (example + field descriptions)
-----------------------------------------------------
```
Example JSON:
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
```

Descriptions:
```
- slurm:         SLURM job submission args
- QueryPath:     TSV of assemblies/references, example file "query_pathes.txt_example.txt", included.
- ScriptFolder:  Path to scripts folder
- TargetFolder:  Output folder for results
- TempFolder:    Temp file directory, default: ./snaketemp/
- genelist:      BED file for genes (4- or 5-col), example file "SMNexample.bed", included.
- blocksize:     the max sequence length used for pangenome alleles, larger will be splitted, default: 80,000, 
- ReferencePrefix: Prefix used in reference contigs (e.g., "NC_0609" or "chr")
- NumPartitions: Number of partitions to divide genes (aim ~1000 genes/partition)
```

F. Find, merge and index compiled k-mer matrices
-----------------------------
After Snakemake completes, run:

```
cat tempfolder/*_bfixpartitions{small,large}.list | sed 's/$/_kmatrix.txt/' | sort -u > allpartitions.list
```

Then each line will be each compiled matrix

Then concatenate all matrices:
```
xargs -a allpartitions.list -I {} bash -c '[[ -f "{}" ]] && cat "{}"' > All_matrix.txt
```
Index the matrix:
```
python scripts/matrixindex.py -i All_matrix.txt -g <genes.gff3> -r 1
```
This creates:
```
- All_matrix.txt         → merged matrix
- All_matrix.txt.index   → index file
```

Matrix format v2.0.1
--------------------

`matrixcompile.py` writes fixed-width k-mer metadata. A k-mer row has six
tab-separated fields:

```
row marker | tag(3) | path+strand(3) + qindex(3) | size(4) + qpos(4) + rpos(4) | k-mer(11) | allele indexes
```

The three-character path field stores `(path_index << 1) | strand`, where
strand is `0` for `+` and `1` for `-`. The k-mer begins at byte offset 27 and
the allele-index field begins at byte offset 39. Internal `|` separators are
not used. All numeric fields use the existing 64-character integer encoding.

Limits are path index `0..32767`, query index `0..262143`, and
size/query-position/reference-position `0..16777215`. Compilation stops with
an error if a value is outside its field range; paths are never combined into
a virtual overflow path.

`matrixindex.py` writes this as the first index line:

```
@v2.0.1,support:v1.2.0
```

Note: You might need a background kmer file, which can be obtained from Ctyper's GitHub site, and you may name it as All_matrix.txt.bgd for your convenient

4) Tips & Gotchas
------------------
- Ensure `.fai` index files match the masked FASTAs.
- Be consistent across:
    - TSV sample IDs (e.g., apr001_h1)
    - Reference IDs (e.g., CHM13_h1)
    - FASTA header contig prefixes (e.g., chr1 vs NC_0609...)
- Prefer 5-col BED for exon-based annotation (`Target=`).
- Mixing CHM13 and hg38? Use distinct prefixes to avoid accidental matches.

Questions / Issues
------------------
- See `masking/` and `snakemake/` READMEs
- For BED help: python tools/gff_toGeneBed.py -h
- If problems arise, double-check FASTA paths and index files.
