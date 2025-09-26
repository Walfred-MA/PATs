Used for masking and prepare assemblies for PATs

1. Soft link all your assemblies in the Assemblies/ folder, for example

for f in $orginal_assemblies/; do ln -s "$f" "Assemblies/$(basename "$f")"; done

2. create a path tsv file for all novel assemblies. Tsv file has two columns without header: 1. samplename_h[1/2] 2. filepath. 

for example:
apr001_h1       Assemblies/apr001.1.polished.fa
apr001_h2       Assemblies/apr001.2.polished.fa

apr001 is the sample name (recommend do not use more than 10 letters), _h1 and _h2 mean two haplotypes, if possible, use _h1 as paternal and _h2 as maternal. For CHM13 and HG38, use only _h1

save this as query_pathes_withrefs.txt

3. run snakemake

cd masking

$account: you slurm account, if your slurm does not require, can skip
$partition: you slurm partition, if your slurm does not require, can skip

snakemake -k  --cluster "sbatch --account=$account --partition=$partition --time=500:00:00 {resources.slurm_extra}" --default-resources "mem_mb=3000" --jobs 500  --rerun-incomplete  --notemp --latency-wait 100 --resources mem_gb=1000  &

4. Replace the masked files to the soft links

mv output/*.fa ../Assemblies/
mv output/*.fasta ../Assemblies/
