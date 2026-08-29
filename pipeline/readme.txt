PATs: Pangenome Allele Typing
================================

This document describes the redesigned, target-driven PATs workflow.

Supported entry point:

    python PATs-main/pats.py [target] -r REFERENCE.fa -q QUERY_PATHS.tsv -o PREFIX

The root PATs-main/Snakefile is launched automatically by pats.py. The older
pipeline/Snakefile, pipeline/Snakefile_lite, and pipeline/config.json are kept
only as legacy examples. They are not required by the redesigned workflow.


1. WHAT PATs PRODUCES
---------------------

PATs builds a pangenome-allele FASTA and an indexed matrix for one or more
target genes, genomic regions, or FASTA sequences.

For output prefix PREFIX, the final products are:

    PREFIX.fa
    PREFIX.matrix.txt
    PREFIX.matrix.txt.index

PREFIX.fa is a primary output, not a disposable intermediate. It contains the
final gfixbreaks.py allele/locus sequences from every target group, combined in
the same group order used to combine the matrix shards. Sequences temporarily
reintroduced by addregion.py are not included in PREFIX.fa.

Intermediate and resumable files are stored in:

    PREFIX.work/

Do not delete PREFIX.work/ if you want Snakemake to resume an interrupted run.


2. REQUIREMENTS AND INSTALLATION
--------------------------------

Supported environment:

    Linux cluster
    Python 3.10
    Snakemake 6.5.0
    PuLP 2.7.0
    C/C++ compiler with OpenMP support

Snakemake 6.5.0 requires the older pulp.list_solvers() API. PATs therefore
pins PuLP to 2.7.0. Do not upgrade PuLP independently in this environment.

Create and activate an environment:

    conda create -n PATs python=3.10 -y
    conda activate PATs

Install from the repository root:

    cd PATs-main
    python install.py --jobs 4

The installer installs Python dependencies and required conda/bioconda tools,
then builds the included C++ tools. Important external tools include samtools,
minimap2, blastn, makeblastdb, and winnowmap.

Included compiled PATs tools include:

    kmer_convertor8
    kmer_searcher
    kmerpartition
    kmer_selector
    kmernorm
    kmerstrd
    kmertree

OpenMP is supported and expected for cluster execution.

Check an installation without rebuilding it:

    python install.py --check

Check only Python package versions:

    python install.py --check --skip-external --skip-build

Expected versions can also be checked directly:

    snakemake --version
    python -c 'import pulp; print(pulp.__version__, hasattr(pulp, "list_solvers"))'

Expected output includes:

    6.5.0
    2.7.0 True


3. REFERENCE AND QUERY TABLE
----------------------------

The -r/--reference argument is one reference FASTA file.

The -q/--queries argument is a whitespace-delimited table with no header:

    SAMPLE_NAME    /absolute/path/to/assembly.fa    YES|NO

Column 1: sample or haplotype name
Column 2: assembly FASTA path
Column 3: whitelist status; optional, default NO

Example:

    CHM13_h1     /data/references/chm13v2.fa       YES
    HG38_h1      /data/references/hg38.fa          YES
    sample01_h1  /data/assemblies/sample01.1.fa    NO
    sample01_h2  /data/assemblies/sample01.2.fa    NO
    sample02_h1  /data/assemblies/sample02.1.fa

Rules for the query table:

    * The reference supplied with -r must be the FASTA in the FIRST data row.
    * Blank lines and lines beginning with # are ignored.
    * Sample names must be unique and cannot contain / or backslash.
    * Relative FASTA paths are resolved relative to the query table.
    * YES, Y, TRUE, and 1 are accepted as whitelist values.
    * NO, N, FALSE, 0, a missing value, and an empty value mean non-whitelist.

The first reference row is always treated as whitelist YES internally. If its
input value is missing or NO, PATs reports that it is forcing the generated
working manifest to YES.

INPUT SAFETY

The user-supplied query table is read-only. PATs does not edit, copy, hard-link,
or symlink it. PATs parses it and creates two independent internal manifests:

    PREFIX.work/inputs/queries.tsv
    PREFIX.work/inputs/queries.search.tsv

The first contains normalized absolute paths and whitelist values. The second
is a two-column manifest for the compiled k-mer tools. PATs has an explicit
guard that aborts rather than allowing either generated path to overwrite the
input table.

Assembly FASTA content is also read-only. samtools faidx may create or update
an adjacent FASTA index named ASSEMBLY.fa.fai when it is absent or older than
the FASTA.


4. SAMPLE AND CONTIG NAMING
---------------------------

PATs first uses explicit sample information from record names and descriptions.
For pangenome-style contig names containing #, the first two components encode
the sample and haplotype:

    SAMPLE#1#CONTIG  -> SAMPLE_h1
    SAMPLE#2#CONTIG  -> SAMPLE_h2

Examples:

    HG38#1#chr22       -> HG38_h1
    sample01#2#chr7    -> sample01_h2

The implicit historical fallback is used only when # is absent:

    contig containing NC_0609 -> CHM13_h1
    any other contig          -> HG38_h1

Therefore, an unencoded CHM13_h1 reference should use RefSeq chromosome
accessions such as NC_060946.1. An unencoded HG38_h1 reference normally uses
chr1, chr2, ..., chrX, chrY.

Some derived files append an internal numeric locus suffix. PATs can resolve:

    NC_060946.1_1 -> NC_060946.1

This happens only when the suffixed name is absent and the base name exists in
the relevant FASTA index. If NC_060946.1_1 is an actual indexed FASTA contig,
the exact name is preserved.


IMPORTANT GFF3/REFERENCE NAME WARNING

The sequence names in GFF3 column 1 must match the reference FASTA index. A
GFF3 using chr1, chr2, ..., chrX/chrY generally will not work directly with the
RefSeq CHM13v2 FASTA used here, whose primary chromosomes are named
NC_060925.1 through NC_060948.1.

PATs includes a helper for CAT/LiftOff GFF3 files using chr* primary names:

    python PATs-main/tools/fixcatLiftOffGenes.py \
      -i catLiftOffGenesV1.gff3 \
      -o catLiftOffGenesV1_namefix.gff3

The helper writes a NEW GFF3 and refuses to overwrite its input. It converts
chr1..chr22, chrX, and chrY to the corresponding primary CHM13v2 RefSeq
accessions. Other sequence names are retained unchanged.

Check the result before running PATs:

    awk -F '\t' '!/^#/ {print $1}' catLiftOffGenesV1_namefix.gff3 \
      | sort -u | head
    grep '^NC_060946.1' reference.fa.fai

Use the converted file with --gff3 when -r points to the NC_0609* CHM13 FASTA.
Do not use this conversion merely because coordinates came from hg38: changing
names does not lift coordinates between assemblies.


5. CHOOSING TARGETS
-------------------

Exactly one target mode is required:

    -g/--gene     gene-name prefix or prefixes
    -b/--bed      target BED file
    -f/--fasta    target multi-FASTA file

All three modes support multiple targets.


5.1 Gene mode

Gene mode requires a GFF3 file:

    -g PREFIX --gff3 genes.gff3

Gene matching uses the beginning of the GFF3 gene_name attribute. Therefore:

    -g CYP2D

can match CYP2D6, CYP2D7, and other names beginning with CYP2D. PATs prints the
complete gene names found for every requested prefix.

Multiple prefixes can be supplied as:

    -g CYP2D -g CYP2A
    -g CYP2D,CYP2A
    -g gene_prefixes.txt

The prefix file may contain whitespace-separated or comma-separated names.
Blank lines and # comments are ignored.

Gene coordinates are expanded upstream and downstream by 5,000 bp by default:

    --gene-extension 5000

Nearby selected gene intervals separated by at most 10,000 bp are merged by
default:

    --gene-merge-distance 10000

Exons are selected automatically from exon features carrying the matched
gene_name values in the same GFF3. Do not add --exon in gene mode.


5.2 BED mode

BED mode accepts a standard zero-based, half-open BED file:

    chromosome    start    end    optional_name

Run with:

    -b targets.bed

An optional pooled exon BED can be supplied:

    -b targets.bed --exon exons.bed

Target and exon contig names must exist in the reference FASTA index.

BED mode does not automatically merge nearby rows by genomic distance. Each
row begins as a separate target record. The later k-mer partitioner may place
similar records in the same group, but that is not equivalent to defining one
continuous locus interval.


5.3 FASTA mode

FASTA mode accepts a multi-FASTA containing one or more target sequences:

    -f targets.fa

FASTA record names must be unique. An optional pooled exon FASTA can be used:

    -f targets.fa --exon exons.fa

When genomic coordinates are unavailable, PATs uses find_consistent_units.py
to find consistent initial blocks among the hit reference loci before building
local graphs.


5.4 Exon search behavior

When exon input is available, PATs pools the exon sequences and searches using
only pooled exon k-mers. The k-mers do not need to retain their exon of origin.

Default hit criteria are:

    target search: 1,000 k-mers in a 2,000-bp window
    exon search:     300 k-mers in a 2,000-bp window

When no exon sequences are available, target k-mers and the target threshold
are used.


6. QUICK-START COMMANDS
-----------------------

Run commands from the directory containing PATs-main and the query table.


6.1 CYP2D gene family with four threads

    python PATs-main/pats.py \
      -g CYP2D \
      --gff3 catLiftOffGenesV1_namefix.gff3 \
      -r /data/references/GCF_009914755.1_T2T-CHM13v2.0_genomic.fa \
      -q query_paths.txt \
      -o CYP2D \
      --threads 4 \
      --sample-threads 1


6.2 Multiple gene families

    python PATs-main/pats.py \
      -g CYP2D,CYP2A \
      --gff3 genes.gff3 \
      -r reference.fa \
      -q query_paths.txt \
      -o CYP_panel \
      --threads 16 \
      --sample-threads 2


6.3 BED targets with exon BED

    python PATs-main/pats.py \
      -b targets.bed \
      --exon exons.bed \
      -r reference.fa \
      -q query_paths.txt \
      -o BED_panel \
      --threads 16


6.4 FASTA targets with exon FASTA

    python PATs-main/pats.py \
      -f targets.fa \
      --exon exons.fa \
      -r reference.fa \
      -q query_paths.txt \
      -o FASTA_panel \
      --threads 16


6.5 Validate and inspect before running

Validate inputs and generate the working configuration without launching
Snakemake:

    python PATs-main/pats.py [normal arguments] --prepare-only

Build the Snakemake DAG without executing jobs:

    python PATs-main/pats.py [normal arguments] \
      --dry-run \
      --snakemake-arg=--reason


7. IMPORTANT OPTIONS AND DEFAULTS
---------------------------------

    Option                        Default    Meaning
    ----------------------------  ---------  -----------------------------------
    --gene-extension              5000       Gene upstream/downstream extension
    --gene-merge-distance         10000      Maximum gap for merging gene targets
    --anchor-size                 5000       Added to both ends of merged hits
    --edge-distance               20000      Non-whitelist edge filter distance
    --kmer-size                   31         K-mer length
    --partition-min-kmers         500        Minimum shared-k-mer requirement
    --partition-similarity        0.20       Shared-k-mer fraction requirement
    --hit-kmers                   1000       Target k-mers required per window
    --exon-hit-kmers              300        Exon k-mers required per window
    --hit-window                  2000       KmerSearcher window size
    --threads                     CPU count  Total Snakemake/group threads
    --sample-threads              1          KmerSearcher threads per FASTA
    --scripts-dir                 scripts/   Compiled PATs binary directory
    --reference-prefix            inferred   Override gfixbreaks ref prefix
    --profile                     empty      Optional Snakemake profile

The partitioner groups multiple targets by k-mer similarity. The default group
criterion uses at least max(500 shared k-mers, 0.20 times the target k-mers).
Targets not assigned by the compiled partitioner are retained as singleton
groups. Final identifiers are deterministic: group1, group2, ...

Record prefixes include the output name and group, for example:

    CYP2Dgroup1_

Additional Snakemake arguments can be passed by repeating:

    --snakemake-arg=ARGUMENT

Use the equals sign when ARGUMENT begins with --.


8. PIPELINE STAGES
------------------

1. Validate the reference and query table and ensure FASTA indexes exist.

2. Prepare targets:
       gene mode  -> match gene_name prefixes and extract reference regions
       BED mode   -> extract BED regions from the reference
       FASTA mode -> normalize the supplied target sequences

3. Prepare pooled exon sequences when exon input is available.

4. Convert target sequences to k-mers and partition multiple targets using
   k-mer similarity. Each partition receives a groupN identifier.

5. Run KmerSearcher against all query assemblies. Exon k-mers and the lower
   exon threshold are used when exon sequences are present.

6. Merge hit gaps within two times the anchor size, then extend both ends by
   the anchor size. With defaults, gaps up to 10 kb are merged and 5 kb is
   added at each end.

7. Classify regions containing N gaps as filtered. Non-whitelist regions less
   than 20 kb from a contig edge are also filtered. Filtered regions are kept
   for later reinclusion rather than discarded permanently.

8. Build a local graph for every group with gfixbreaks.py and graphmake.py.
   gfixbreaks caching is disabled. BED/gene coordinates provide initial regions;
   FASTA mode obtains initial reference blocks with find_consistent_units.py.

9. Temporarily reinclude filtered regions with addregion.py so valid exclusive
   k-mers are not lost solely because of sequence filtration. Run kmer_selector
   across all query assemblies. Its scratch files are isolated under
   PREFIX.work/kmers/ and removed after selection; the augmented FASTA is also
   removed after k-mer selection completes.

10. Run packedrun.py and matrixcompile.py against the final fixed.fa sequences,
    not the temporary augmented FASTA. Combine fixed FASTA shards into PREFIX.fa,
    combine group matrix shards, and create the index with matrixindex.py.


9. OUTPUT AND WORK DIRECTORY
----------------------------

Final files:

    PREFIX.fa
        Combined final allele/locus FASTA produced by gfixbreaks.py. This is a
        primary PATs result and contains no temporary addregion.py sequences.

    PREFIX.matrix.txt
        Combined PATs allele/k-mer matrix.

    PREFIX.matrix.txt.index
        Seek/index metadata used by Ctyper.

Important intermediate files:

    PREFIX.work/config.json
        Normalized configuration generated by pats.py.

    PREFIX.work/inputs/queries.tsv
        Internal three-column query/whitelist manifest.

    PREFIX.work/inputs/queries.search.tsv
        Internal two-column KmerSearcher manifest.

    PREFIX.work/targets/targets.bed
    PREFIX.work/targets/targets.fa
    PREFIX.work/targets/exons.fa
    PREFIX.work/targets/report.json
        Normalized targets, pooled exons, and gene-match report.

    PREFIX.work/targets/groups.tsv
        Target-group manifest generated by the partition checkpoint.

    PREFIX.work/search/
        KmerSearcher outputs.

    PREFIX.work/groups/groupN/hits.bed
    PREFIX.work/groups/groupN/hits.report.json
    PREFIX.work/groups/groupN/fixed.fa
    PREFIX.work/groups/groupN/matrix.txt
        Per-group hit, final graph FASTA, and matrix outputs.

Temporary selector files:

    PREFIX.work/groups/groupN/augmented.fa
        addregion.py output used only by kmer_selector. Snakemake removes it
        after successful exclusive-k-mer selection.

    PREFIX.work/kmers/groupN.exclusive.kmers.txt
        Exclusive k-mers consumed by packedrun.py and then removed as a
        Snakemake temporary output. Per-assembly kmer_selector scratch files use
        the same kmers/ directory and are removed immediately after selection.


10. MATRIX FORMAT AND CTYPER
----------------------------

The current PATs writer emits matrix encoding v2.0.1. The matrix index begins
with version information equivalent to:

    @v2.0.1,support:v1.2.0

The v2.0.1 k-mer metadata uses fixed-width base-64 fields:

    path index + strand    3 characters
    query index            3 characters
    size                   4 characters
    query position         4 characters
    reference position     4 characters

The path index is limited to 15 bits, with a maximum index of 32,767. PATs
raises an error instead of silently creating a virtual appended path when a
value exceeds a supported field limit.

The shared encoding constants are in:

    scripts/matrixformat.py

This format is intended for the updated Ctyper 1.2.0 reader with v2.0.1 format
information and automatic detection for matrices lacking an explicit version.
An older, unmodified Ctyper binary may not decode a v2.0.1 matrix correctly.


10.1 PROFILING WITH CTYPER FOR FAST TARGET GENOTYPING
------------------------------------------------------

PATs builds the matrix database. Ctyper profiling is a separate, usually
one-time step that examines aligned reads to locate the reference intervals
from which informative or mismapped target reads originate. The resulting BED
file can then be supplied with -B so later Ctyper runs read only those regions
instead of scanning every aligned read.

This workflow is intended for indexed BAM or CRAM files. It has two phases:

    1. Profile one or more representative aligned samples to make a target BED.
    2. Reuse that BED with -B for fast genotyping of other samples.

Required files
~~~~~~~~~~~~~~

Keep the PATs matrix and its index together:

    CYP2D.matrix.txt
    CYP2D.matrix.txt.index

The Ctyper binary must include support for PATs matrix encoding v2.0.1. BAM
files should have a .bai index and CRAM files should have a .crai index. When
reading CRAM, provide the exact decoding reference with -T unless REF_CACHE and
REF_PATH have already been configured.

The profiling BED is tied to the reference coordinate system used by the
aligned reads, not necessarily the reference used by PATs to construct the
matrix. For example, reads aligned to hg38 require an hg38 profiling BED, even
if the PATs matrix was anchored on CHM13. Do not reuse an hg38 profiling BED
for CHM13-aligned reads, or the reverse.

Record the alignment-reference MD5 in the BED filename:

    md5sum aligned_reference.fa

For example:

    TargetRegions.<reference_md5>.bed

Profiling one representative sample
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Run profiling without -B. Ctyper must scan the aligned file to discover all
useful source intervals:

    ctyper \
      -m CYP2D.matrix.txt \
      -p representative.cram \
      -o TargetRegions.<reference_md5>.bed \
      -T aligned_reference.fa \
      -d 1 \
      -N 4 \
      2> CYP2D.profile.log

In the current Ctyper 1.2.0 command-line implementation, -d 1 is needed to
satisfy coverage-input validation during profiling. Profiling does not use
that placeholder value for genotyping.

For a CYP2D-only PATs matrix, no -g option is necessary because every matrix
group is already a CYP2D target. For a database containing many gene families,
add a quoted prefix or an exact gene/matrix name, for example:

    -g 'CYP2D*'
    -g CYP2D6
    -g '#CYP2Dgroup1'

Use the exact names present in CYP2D.matrix.txt.index. A trailing * selects a
prefix and should be quoted so the shell does not expand it.

Profiling several samples
~~~~~~~~~~~~~~~~~~~~~~~~~

Several representative samples can find mapping locations that are absent
from one individual. Make one input list and one output list, with one path per
line and the same number and order of lines:

    profile.inputs.txt                 profile.outputs.txt
    ------------------                 -------------------
    sample1.cram                       profiles/sample1.bed
    sample2.cram                       profiles/sample2.bed
    sample3.cram                       profiles/sample3.bed

Then run:

    ctyper \
      -m CYP2D.matrix.txt \
      -P profile.inputs.txt \
      -O profile.outputs.txt \
      -o TargetRegions.<reference_md5>.bed \
      -T aligned_reference.fa \
      -d 1 \
      -n 3 \
      -N 2 \
      2> CYP2D.profile.log

Keep -P before -O and keep the merged-summary -o after -O. The -O paths receive
the per-sample profiling BEDs and -o receives the merged BED. Use a new summary
pathname rather than mixing regions from an unrelated reference or an older
profiling experiment.

Fast genotyping of one sample
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After profiling, reuse the BED with -B:

    ctyper \
      -m CYP2D.matrix.txt \
      -i new_sample.cram \
      -o new_sample.CYP2D.ctyper.txt \
      -B TargetRegions.<reference_md5>.bed \
      -T aligned_reference.fa \
      -d 24 \
      -N 4 \
      2> new_sample.CYP2D.log

The example -d 24 is the expected 31-mer depth for approximately 30x coverage
with 150-bp reads:

    31-mer depth = (1 - 30/read_length) * sequencing depth

Replace 24 with the appropriate value for the sample. Do not use -d together
with -b. If CYP2D.matrix.txt.bgd exists, Ctyper uses it automatically and -d
can be omitted. A depth file can instead be provided with -D for a cohort.

Ctyper appends normal genotyping output when an -o path already exists. Use a
new output pathname for each run unless appending is intentional.

Fast genotyping of a cohort
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Create matching files with one input path and one output path per line:

    cohort.inputs.txt                  cohort.outputs.txt
    -----------------                  ------------------
    sampleA.cram                       results/sampleA.ctyper.txt
    sampleB.cram                       results/sampleB.ctyper.txt
    sampleC.cram                       results/sampleC.ctyper.txt

For samples with individual depth estimates, also create cohort.depths.txt
with one 31-mer depth value per line in input order. Then run:

    ctyper \
      -m CYP2D.matrix.txt \
      -I cohort.inputs.txt \
      -O cohort.outputs.txt \
      -D cohort.depths.txt \
      -B TargetRegions.<reference_md5>.bed \
      -T aligned_reference.fa \
      -n 8 \
      -N 2 \
      2> CYP2D.cohort.log

Here -n is the number of samples processed in parallel and -N is the number of
threads used within each sample. Approximate simultaneous CPU use is -n times
-N, so select both values to fit the cluster allocation. Parallelizing across
samples is usually the most efficient choice for a cohort; 1-4 threads per
sample is a practical starting point when storage I/O is limiting.

For a multi-family matrix, restrict both profiling and genotyping consistently
with -g 'CYP2D*' or a matching -G gene-list file. When -g/-G and -B are used
together, Ctyper uses only BED records whose names match the selected targets.

If no profiling BED is available
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Ctyper can use matrix reference intervals as a fallback for a selected target:

    ctyper \
      -m multi_family.matrix.txt \
      -i sample.cram \
      -o sample.CYP2D.ctyper.txt \
      -g 'CYP2D*' \
      -r gene \
      -T aligned_reference.fa \
      -d 24 \
      -N 4 \
      2> sample.CYP2D.log

This fallback can miss useful reads that map outside the expected locus and is
less robust to mismapping and reference bias. A BED learned by profiling
representative aligned samples is preferred for repeated fast target runs.


11. RESUMING, DRY RUNS, AND UNLOCKING
-------------------------------------

PATs automatically passes these options to Snakemake:

    --rerun-incomplete
    --printshellcmds

To resume after a failed job, run exactly the same pats.py command again. Do
not delete PREFIX.work/. Successfully completed, up-to-date outputs are reused.
Generated manifests and config are rewritten only when their content changes,
preserving timestamps for Snakemake resume behavior.

The target partition step is a Snakemake checkpoint. Before that checkpoint
finishes, the initial job table cannot display every later group-specific job.
This is normal and is not evidence that outputs are being deleted.

See why jobs need to run:

    python PATs-main/pats.py [normal arguments] \
      --dry-run \
      --snakemake-arg=--reason

Unlock only after confirming that no PATs/Snakemake process is running:

    jobs -l

Then rerun the original command with:

    --snakemake-arg=--unlock

Remove --unlock before the real run. An equivalent direct unlock is:

    snakemake \
      --snakefile PATs-main/Snakefile \
      --configfile PREFIX.work/config.json \
      --cores 1 \
      --unlock

For a background run with a persistent log:

    nohup python PATs-main/pats.py [normal arguments] \
      > PREFIX.log 2>&1 &


12. FREQUENTLY ASKED QUESTIONS
------------------------------

Q1. I have three genes. Should they be three BED rows or one merged interval?

It depends on whether the genes should be represented as independent targets
or as one structural locus.

Use separate BED rows when:

    * each gene should be typed and partitioned independently;
    * the genes are far apart; or
    * structural alleles are not expected to span more than one gene.

Use one merged BED interval when:

    * the genes are adjacent members of one gene family;
    * CNVs, deletions, duplications, conversions, or rearrangements can span
      multiple genes; and
    * the database should preserve haplotype context across the whole locus.

This choice matters. In BED mode, separate rows remain separate target records
and are not automatically proximity-merged. The k-mer partitioner may group
similar rows, but that does not make them one continuous reference interval.
One merged interval builds the target from the entire continuous sequence and
can preserve variants spanning gene boundaries.

Do not merge unrelated distant genes into one very large interval. That adds
irrelevant k-mers, increases graph size, and may reduce locus specificity.

For a nearby family such as CYP2D, gene mode is often easier:

    -g CYP2D --gff3 matching_reference.gff3

Gene mode selects every matching gene and merges selected intervals within
--gene-merge-distance (10 kb by default), with --gene-extension (5 kb by
default) added around each merged locus.


Q2. Do gene coordinates depend on the reference system?

Yes. Target coordinates and contig names must match the exact FASTA supplied
with -r.

    T2T/CHM13 reference -> use T2T/CHM13 coordinates and matching contig names
    hg38 reference      -> use hg38 coordinates and matching contig names

BED coordinates from hg38 cannot be used directly with a T2T FASTA, and T2T
coordinates cannot be used directly with hg38. Likewise, renaming GFF3 column
1 from chr22 to NC_060946.1 changes only the sequence name; it does not convert
hg38 coordinates into T2T coordinates.

BED coordinates are zero-based and half-open. GFF3 coordinates are interpreted
according to the one-based, closed GFF3 convention and converted internally.

The reference FASTA must also be the first query-table row. Other query
assemblies do not need to share the reference coordinate system because PATs
finds their loci using target/exon k-mers and builds local graphs.


Q3. Can I build a database for only one family, such as CYP2D?

Yes. PATs and the updated Ctyper do not require multiple gene families in one
matrix. A CYP2D-only matrix and index are valid, and Ctyper can use them.

Such a database will type only the loci and alleles represented in that matrix.
It cannot report gene families that were not included. Accuracy still depends
on the completeness and diversity of the query assemblies used to build the
database, the target boundaries, and the quality of the selected exclusive
k-mers. A one-family database is often smaller, faster to build, and faster to
use than a broad multi-family database.


13. TROUBLESHOOTING
-------------------

Problem: module 'pulp' has no attribute 'list_solvers'

Cause: Snakemake 6.5.0 is running with PuLP 2.8 or newer.

Fix:

    python -m pip install --force-reinstall \
      'snakemake==6.5.0' 'PuLP==2.7.0'


Problem: Python version rejected by install.py

Use Python 3.10. The pinned legacy Snakemake stack is not supported by this
PATs installer under Python 3.11 or newer.


Problem: samtools faidx reports that a contig is absent

Check exact FASTA index names:

    grep '^NC_060946' reference.fa.fai
    samtools faidx reference.fa 'NC_060946.1:START-END' | head

GFF3, BED, and unencoded FASTA contig names must resolve to names in the
reference index. Review the GFF3 warning and conversion helper in Section 4.


Problem: Snakemake reports a locked working directory

Confirm no process is active with jobs -l or ps. Then use the unlock command
in Section 11. Never unlock a directory used by an active run.


Problem: the initial job table shows only a small number of rules

partition_targets is a checkpoint. Snakemake expands remaining group jobs after
the checkpoint produces PREFIX.work/targets/groups.tsv.


Problem: a failed job appears again after restart

This is expected. A failed rule and missing downstream outputs must run. Use
--dry-run with --reason to see why each job is scheduled.


Problem: a compiled PATs executable cannot be found

Run:

    python PATs-main/install.py --check

By default, pats.py searches PATs-main/scripts/. A different binary directory
can be supplied with --scripts-dir.


14. GETTING HELP
----------------

Show all command-line options:

    python PATs-main/pats.py --help

For a reproducible problem report, include:

    * complete pats.py command
    * first query-table row, with sensitive paths shortened if necessary
    * snakemake --version
    * python --version
    * failing rule name
    * complete traceback and external-command stderr


Author: Wangfei "Walfred" Ma, Chaisson Lab, USC
