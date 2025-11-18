PATs: Pangenome Allele Typing Snakemake Pipeline
------------------------------------------------

This pipeline performs pangenome-based allele analysis and k-mer genotyping using 
partitioned sequence regions from pangenome assemblies and short-read queries.

Snakefile: the main pipeline file used for large number of genes
Snakefile_light: the alternative lighter pipeline file used for only a few genes


📁 Repository Root:
    Place this pipeline script and the `config.json` file at the root.
    Compiled binaries and Python scripts are expected in the folder specified 
    by the "ScriptFolder" variable in config.json.

🧰 Required Tools:
    - BLAST (blastn, makeblastdb)
    - bedtools, samtools, minimap2, winnowmap
    - zlib (for gz support)

🔧 Binaries Automatically Built:
    - kmercounter8
    - kmerselector
    - kmerpartition
    - kmernorm
    - KmerStrd
    - kmertree
  If not found in ScriptFolder, they will be compiled from ScriptFolder/src/<toolname>.

📂 Inputs:
    - `config.json`: pipeline configuration (see below)
    - `genelist`: BED-like file of target blocks
    - `QueryPath`: TSV file mapping sample name to query FASTA + index

🛠️ Config Fields (example config.json):
{
    "slurm": " --account=mchaisso_100 --time 50:00:00 --partition=qcb ",
    "QueryPath": "query_pathes.txt_withrefs.txt",
    "ScriptFolder": "../scripts/",
    "TargetFolder": "groups/",
    "TempFolder": "snaketemp",
    "genelist": "regions.bed",
    "blocksize": 80000,
    "NumPartitions": 10,
    "ReferencePrefix": "NC_0609"
}

▶️ Example the pipeline running command:
    snakemake -k  --cluster "sbatch --account=mchaisso_100 --partition=qcb --time=500:00:00 {resources.slurm_extra}" --default-resources "mem_mb=3000" --jobs 500  --rerun-incomplete  --notemp --latency-wait 100 --resources mem_gb=1000  &

The pipeline handles:  
    - block preparation and partitioning  
    - hotspot detection and alignment  
    - k-mer selection and normalization  
    - fixed k-mer reanalysis and matrix compilation

This is designed for high-throughput, parallel processing of pangenome haplotype alignments 
and short-read query integration.

🧪 Author: Wangfei "Walfred" Ma, Chaisson Lab, USC
