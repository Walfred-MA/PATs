# Experimental PAT transcripts run

This opt-in command builds WES or RNA matrices from existing ParaTyper results.
It is intentionally not listed in the main README or the ordinary `pats.py`
help. The genomic Snakemake workflow, shared scripts and their defaults are
unchanged. This command does not run ParaTyper or infer new isoforms.

## Invocation

Create a headered, tab-separated manifest:

```text
sample	assembly	calls	fragments
S1_h1	assemblies/S1_h1.fa	ParaTyper/S1_h1.transcript_calls.tsv
S2_h1	assemblies/S2_h1.fa	ParaTyper/S2_h1.transcript_calls.tsv	ParaTyper/S2_h1.pseudofragments.tsv
```

Paths are relative to the manifest, including paths with spaces. Samples must
be unique; contig names may repeat across samples. The optional `fragments`
column explicitly includes the corresponding ParaTyper fragment table. An
empty calls table is accepted when other samples or its fragment table supply
models. Do not list the same parent in both call tables.

```bash
python pats.py transcripts run --mode wes --manifest samples.tsv -o results/exome
python pats.py transcripts run --mode rna --manifest samples.tsv -o results/rna
```

Python 3.10+, NumPy and SciPy are required (already PATs dependencies). No graph
or compiled PATs executables are required for this command. Assemblies must be
uncompressed, regularly wrapped FASTA. A current `.fai` is reused when present;
otherwise the file is scanned once to build an in-memory contig index. Neither
the assemblies nor their indexes are modified. Call tables may be plain TSV or
gzip compressed TSV.

Options:

- `--reference-sample NAME`: populate native matrix index reference intervals
  using this sample's extracted regions. Without it the reference-region field
  is empty. Supplied regions use the assembly's coordinates; for reads aligned
  to another reference, obtain a matching Ctyper profiling BED.
- `--scratch-dir DIR`: disk location for the temporary SQLite catalog and k-mer
  postings. Allow space for candidate sequences, counts and database indexes.
- `--max-partition-candidates N`: default 5,000 distinct retained count profiles
  per partition. Tree distances use quadratic memory; the native Ctyper tree
  index limits the maximum to 32,768 leaves. Exceeding a limit fails explicitly.
- `--force`: replace outputs only after the build and format checks succeed.

## Sequence and partition rules

Input is the native ParaTyper 3.9.2 headered schema, parsed by column names.
Extra columns are retained. Exons use target `exon_query_coordinates`, which
must refer to the supplied assembly; source-reference coordinates are not used
for extraction. Arrays are paired within each alternative, shared scalar
values are broadcast, and the original rows are preserved in annotations.

**WES:** one candidate per sample, `GENE_index` and parent gene alternative.
Collect the parent's exon intervals, its owned `merged_interval_coordinates`
and inserted exon blocks, and its matching children's exon intervals. Expand
each interval by 30 bp at both ends, clip at contig boundaries, then union
overlapping/touching intervals within that locus. Count k-mers independently
in each resulting interval and sum them into the locus column. Disconnected
intervals are separate FASTA records; they are never concatenated. Alternative
isoforms do not multiply genomic exon counts. WES requires full-gene parent
rows; transcript-only ParaTyper calls cannot reliably assign genomic copies.

**RNA:** one candidate per transcript alternative. Assemble target exons in
transcript orientation, including reverse complementation on the minus strand.
No anchors are added. Join adjacent reference exon numbers only when both
exons report 100% coverage. Missing exon numbers or incomplete alignment
coverage split a model into separate segments: no k-mer crosses that break.
Partial models are retained and labeled, including contiguous partial models.
`call_status=complete` is ParaTyper's completeness among eligible exons, not a
guarantee of a full biological transcript, correct splice ends or expression.

Both modes normalize lowercase bases to uppercase, preserve other alphabetic
IUPAC bases in FASTA, and break k-mer counting at any non-ACGT base. Exon,
sample, gene, transcript, parent, strand, ambiguity, source and fragment
information are retained. Fragment candidates are labeled background evidence,
not automatically interpreted as intact genes.

Partitions begin with gene IDs pooled across all samples, copies and isoforms.
ParaTyper `GA&GB` shared-exon units remain one unit; their child isoforms retain
the original gene IDs. Competing `GA;GB` alternatives remain separately
annotated hypotheses, linked by sample/locus/tie IDs. They do not assert extra
physical copies. Same-evidence alternatives collapse when they occur in the
same final partition; consumers must respect ambiguity links if different
profiles remain in separate columns or partitions.

An edge joins original gene units that share **more than 100 distinct
canonical 31-mers**. PATs uses the numerically larger of a k-mer and its reverse
complement. Final partitions are connected components of those edges, not
repeated threshold comparisons against growing unions. WES uses anchored
interval k-mers; RNA uses segment/splice k-mers, so the modes can partition
differently. Repeated occurrences do not increase partition overlap counts.

K-mers still shared across final partitions are excluded from every matrix
partition. Within a partition, real occurrence multiplicities are retained.
Identical retained count profiles share one matrix column, preserving every
alias in the title and annotation table. Candidates without informative
31-mers remain in FASTA/annotations with `matrix=null`; all-empty builds fail.

Specificity is relative to the supplied candidate catalog. This version does
not perform a whole-genome off-target screen or calibrate WES capture/RNA
expression bias. Supply relevant fragment annotations when available.

## Ctyper compatibility contract

The matrix uses PATs/Ctyper **v2.0.1** and the index starts with
`@v2.0.1,support:v1.2.0`. It reuses `scripts/matrixformat.py`'s native field
encoders. Matrix block headers, candidate titles, counts, trees, allele groups,
k-mer membership/multiplicity and tree annotations use existing encodings.

- K-mer rows use `&+`, six tab-separated fields, the 11-character k-mer at byte
  27, and the allele-index field at byte 39.
- Graph path, query index and coordinate fields are zero (no graph). Titles
  retain `NA:NA` for absent graph alignments; coordinates are in the sidecar.
- Tags use Ctyper's post-correction neutral values: bias 0, variance 100 and
  reciprocal multiplier 100, encoded as `!{S`. This avoids treating raw,
  uncalibrated compiler tags as low-variance evidence. No empirical calibration
  is implied.
- Allele indexes are zero-based and preserve repeated occurrences. Native
  index offsets are byte offsets into the final matrix, with correct row and
  allocation sizes. Candidate groups are singletons after equivalence collapse.
- Trees use average linkage on weighted Jaccard distances between retained
  k-mer count profiles. They are sequence-similarity trees, affected by splice
  structure, with native preorder node-delta annotations. Tree leaf order and
  matrix column order are identical.

No Ctyper source change is included. Its matrix reader can load these files;
WES depth/capture normalization, RNA abundance/isoform inference, partial model
handling and physical-locus ambiguity constraints belong to the consumer's
assay-specific logic. A genomic `.bgd` is not generated or copied: it would not
establish appropriate WES/RNA depth calibration.

## Files and annotation mapping

For prefix `results/exome`, outputs are:

```text
results/exome.fa
results/exome.matrix.txt
results/exome.matrix.txt.index
results/exome.annotations.jsonl
results/exome.partitions.tsv
results/exome.merges.tsv
results/exome.trees.tsv
results/exome.report.json
```

Each JSONL record describes one candidate with its original ParaTyper source
rows, sample, gene unit, locus, isoform and segment coordinates (0-based,
half-open). `fasta_records` lists its sequence records in segment order.
`matrix.partition`, `matrix.column` (zero-based) and `matrix.representative`
map it to the native matrix. Different aliases may map to the same column.
The sample plus `locus_id` identifies the physical parent across isoforms and
alternatives. WES segment coordinates are genomic-order intervals; RNA
segments and their exon lists are in transcript order.

`partitions.tsv` records every gene unit, even if no matrix evidence survived.
`merges.tsv` records the threshold-crossing edges and exact shared counts.
`trees.tsv` contains the emitted matrix trees. The report records filtering,
equivalence collapse, partial/segmented models, omitted candidates and the
specificity scope. Inputs are processed through a temporary disk-backed
catalog, one matrix partition is loaded at a time, and scratch data is removed
on completion or failure. Whole-pangenome runtime/memory has not been benchmarked.

## Regression and native-reader checks

```bash
python -m unittest discover -s tests -p 'test_transcripts.py' -v
```

The optional integration harness loads the output with **unchanged Ctyper**
`PriorData`, `KtableReader` and tree code. It checks native row membership,
multiplicity, tags, counts, singleton/multileaf trees and multiple partition
offsets. To build it, point the include/library paths at your Ctyper, Eigen and
HTSlib installations (shown here with a conda dependency environment):

```bash
c++ -std=c++17 -O1 -g \
  -I/path/to/Ctyper/src -I"$CONDA_PREFIX/include" \
  -I"$CONDA_PREFIX/include/eigen3" \
  tests/ctyper_transcripts_reader.cpp \
  /path/to/Ctyper/src/PriorData.cpp \
  /path/to/Ctyper/src/KtableReader.cpp \
  /path/to/Ctyper/src/FileReader.cpp \
  /path/to/Ctyper/src/TreeRound.cpp \
  -L"$CONDA_PREFIX/lib" -lhts -o /tmp/ctyper-transcripts-reader

CTYPER_READER=/tmp/ctyper-transcripts-reader \
  python -m unittest discover -s tests -p 'test_transcripts.py' -v
```

Integration checked against Ctyper commit
`f22433742162cd843ef45225c5b3e741e8b31840`, with address/undefined-behavior
sanitizers. This verifies serialization/reading, not biological genotype accuracy.

`tests/data/paratyper_smn.transcript_calls.tsv.gz` contains unmodified output generated
by ParaTyper commit `6a4e0c00bec4fe52ec8db260b12537f12e1c22ad` from its frozen
`SMN_CHM13_gene_filter` regression alignments, GFF3 and eligible exon table:

```bash
python /path/to/ParaTyper/scripts/call_genes_from_exon_alignments_v3.py \
  -i /path/to/ParaTyper/test/regression/data/SMN_CHM13_gene_filter.alignments.tsv \
  -g /path/to/ParaTyper/test/regression/data/SMN_CHM13_gene_filter.gff3 \
  --eligible-exon-info /path/to/ParaTyper/test/regression/data/SMN_CHM13_gene_filter.exon_ids.tsv \
  --query-coordinate-mode local --threads 1 \
  -o tests/data/paratyper_smn.transcript_calls.tsv
gzip -n tests/data/paratyper_smn.transcript_calls.tsv
```
