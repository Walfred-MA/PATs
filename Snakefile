from pathlib import Path
import csv
import json
import shlex
import subprocess
import sys


# Keep the workflow portable: its Python helpers always come from the
# scripts/ directory beside this Snakefile.  Only compiled-tool lookup is
# configurable through tools_dir.
RUNTIME = Path(workflow.basedir).resolve()
WORK = Path(config["work_dir"])
SCRIPTS = RUNTIME / "scripts"
TOOLS = Path(config["tools_dir"])
REFERENCE = config["reference"]
QUERY_TABLE = config["query_table"]
SEARCH_TABLE = config["search_table"]
MODE = config["mode"]
THREADS = int(config["threads"])

TARGET_DIR = WORK / "targets"
TARGET_BED = TARGET_DIR / "targets.bed"
TARGET_FASTA = TARGET_DIR / "targets.fa"
EXON_FASTA = TARGET_DIR / "exons.fa"
TARGET_REPORT = TARGET_DIR / "report.json"
TARGET_KMERS = TARGET_DIR / "targets.kmers.txt"
EXON_KMERS = TARGET_DIR / "exons.kmers.txt"
GROUP_DIR = TARGET_DIR / "groups"
GROUP_MANIFEST = TARGET_DIR / "groups.tsv"
INDEX_SENTINEL = WORK / "inputs" / "query_indexes.tsv"
SEARCH_DIR = WORK / "search"
SEARCH_DONE = SEARCH_DIR / ".done"
KMERS_DIR = WORK / "kmers"
MATRIX = Path(config["matrix"])
MATRIX_INDEX = Path(config["matrix_index"])
FINAL_FASTA = Path(config.get("fasta", f"{config['output_prefix']}.fa"))

INDEX_QUERIES_SCRIPT = shlex.quote(str(SCRIPTS / "index_queries.py"))
PREPARE_TARGETS_SCRIPT = shlex.quote(str(SCRIPTS / "prepare_targets.py"))
CONVERT_KMERS_SCRIPT = shlex.quote(str(SCRIPTS / "convert_kmers.py"))
PARTITION_TARGETS_SCRIPT = shlex.quote(str(SCRIPTS / "partition_targets.py"))
RUN_SEARCH_SCRIPT = shlex.quote(str(SCRIPTS / "run_kmer_search.py"))
BUILD_HITS_SCRIPT = shlex.quote(str(SCRIPTS / "build_hit_fastas.py"))
CONSISTENT_SCRIPT = shlex.quote(str(SCRIPTS / "find_consistent_units.py"))
RUN_GFIXBREAKS_SCRIPT = shlex.quote(str(SCRIPTS / "run_gfixbreaks.py"))
ADDREGION_SCRIPT = shlex.quote(str(SCRIPTS / "addregion.py"))
SELECT_KMERS_SCRIPT = shlex.quote(str(SCRIPTS / "select_exclusive_kmers.py"))
PACKEDRUN_SCRIPT = shlex.quote(str(SCRIPTS / "packedrun.py"))
COMBINE_FASTAS_SCRIPT = str(SCRIPTS / "combine_fastas.py")

Q_QUERY_TABLE = shlex.quote(str(QUERY_TABLE))
Q_SEARCH_TABLE = shlex.quote(str(SEARCH_TABLE))
Q_REFERENCE = shlex.quote(str(REFERENCE))
Q_MODE = shlex.quote(str(MODE))
Q_SCRIPTS = shlex.quote(str(SCRIPTS))
Q_TOOLS = shlex.quote(str(TOOLS))
Q_INDEX_SENTINEL = shlex.quote(str(INDEX_SENTINEL))
Q_TARGET_BED = shlex.quote(str(TARGET_BED))
Q_TARGET_FASTA = shlex.quote(str(TARGET_FASTA))
Q_EXON_FASTA = shlex.quote(str(EXON_FASTA))
Q_TARGET_REPORT = shlex.quote(str(TARGET_REPORT))
Q_TARGET_KMERS = shlex.quote(str(TARGET_KMERS))
Q_EXON_KMERS = shlex.quote(str(EXON_KMERS))
Q_SEARCH_DIR = shlex.quote(str(SEARCH_DIR))
PYTHON = shlex.quote(sys.executable)


def q(value):
    return shlex.quote(str(value))


def target_sources(_wildcards):
    paths = [REFERENCE]
    for key in ("target_input", "gff3", "exon"):
        value = config.get(key, "")
        if value:
            paths.append(value)
    return paths


def query_sources(_wildcards):
    paths = [QUERY_TABLE]
    with Path(QUERY_TABLE).open() as handle:
        for raw in handle:
            if raw.strip() and not raw.lstrip().startswith("#"):
                fields = raw.split()
                if len(fields) >= 2:
                    paths.append(fields[1])
    return paths


def checkpoint_manifest(_wildcards):
    # This checkpoint has no wildcards of its own.  Downstream group rules do,
    # so forwarding their ``group`` wildcard to ``get`` makes Snakemake reject
    # it as an unexpected checkpoint wildcard.
    return checkpoints.partition_targets.get().output.manifest


def read_group_rows(wildcards):
    manifest = Path(checkpoint_manifest(wildcards))
    with manifest.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"target group manifest is empty: {manifest}")
    return rows


def group_value(field):
    def lookup(wildcards):
        for row in read_group_rows(wildcards):
            if row["group"] == wildcards.group:
                return row[field]
        raise ValueError(f"unknown target group: {wildcards.group}")
    return lookup


def matrix_shards(wildcards):
    return [
        str(WORK / "groups" / row["group"] / "matrix.txt")
        for row in read_group_rows(wildcards)
    ]


def fasta_shards(wildcards):
    return [
        str(WORK / "groups" / row["group"] / "fixed.fa")
        for row in read_group_rows(wildcards)
    ]


rule all:
    input:
        str(MATRIX_INDEX),
        str(FINAL_FASTA),


rule index_queries:
    input:
        query_sources
    output:
        str(INDEX_SENTINEL)
    shell:
        "{PYTHON} {INDEX_QUERIES_SCRIPT} -q {Q_QUERY_TABLE} -o {Q_INDEX_SENTINEL}"


rule prepare_targets:
    input:
        sources=target_sources,
        indexes=str(INDEX_SENTINEL),
    output:
        bed=str(TARGET_BED),
        fasta=str(TARGET_FASTA),
        exons=str(EXON_FASTA),
        report=str(TARGET_REPORT),
    params:
        target_input=(
            "--input " + q(config["target_input"])
            if config.get("target_input") else ""
        ),
        genes=json.dumps(config.get("gene_prefixes", [])),
        gff3=("--gff3 " + q(config["gff3"]) if config.get("gff3") else ""),
        exon=("--exon " + q(config["exon"]) if config.get("exon") else ""),
    shell:
        "{PYTHON} {PREPARE_TARGETS_SCRIPT} "
        "--mode {Q_MODE} {params.target_input} --genes-json {params.genes:q} "
        "{params.gff3} {params.exon} --reference {Q_REFERENCE} "
        "--gene-extension {config[gene_extension]} "
        "--gene-merge-distance {config[gene_merge_distance]} "
        "--targets-bed {Q_TARGET_BED} --targets-fasta {Q_TARGET_FASTA} "
        "--exons-fasta {Q_EXON_FASTA} --report {Q_TARGET_REPORT}"


rule target_kmers:
    input:
        str(TARGET_FASTA)
    output:
        str(TARGET_KMERS)
    shell:
        "{PYTHON} {CONVERT_KMERS_SCRIPT} -i {Q_TARGET_FASTA} -o {Q_TARGET_KMERS} "
        "-k {config[kmer_size]} --scripts-dir {Q_TOOLS}"


rule exon_kmers:
    input:
        str(EXON_FASTA)
    output:
        str(EXON_KMERS)
    shell:
        "{PYTHON} {CONVERT_KMERS_SCRIPT} -i {Q_EXON_FASTA} -o {Q_EXON_KMERS} "
        "-k {config[kmer_size]} --scripts-dir {Q_TOOLS} --allow-empty"


checkpoint partition_targets:
    input:
        fasta=str(TARGET_FASTA),
        kmers=str(TARGET_KMERS),
    output:
        groups=directory(str(GROUP_DIR)),
        manifest=str(GROUP_MANIFEST),
    threads:
        THREADS
    shell:
        "{PYTHON} {PARTITION_TARGETS_SCRIPT} -i {input.fasta:q} -k {input.kmers:q} "
        "-o {output.groups:q} --manifest {output.manifest:q} --mode {Q_MODE} "
        "--prefix {config[name_prefix]} --scripts-dir {Q_TOOLS} "
        "--kmer-size {config[kmer_size]} --min-kmers {config[partition_min_kmers]} "
        "--similarity {config[partition_similarity]} --threads {threads}"


rule search_kmers:
    input:
        manifest=checkpoint_manifest,
        exons=str(EXON_KMERS),
        indexes=str(INDEX_SENTINEL),
    output:
        str(SEARCH_DONE)
    threads:
        THREADS
    shell:
        "{PYTHON} {RUN_SEARCH_SCRIPT} --manifest {input.manifest:q} "
        "--queries {Q_SEARCH_TABLE} --exon-kmers {input.exons:q} --output {Q_SEARCH_DIR} "
        "--scripts-dir {Q_TOOLS} --kmer-size {config[kmer_size]} "
        "--hit-kmers {config[hit_kmers]} --exon-hit-kmers {config[exon_hit_kmers]} "
        "--window {config[hit_window]} --threads {threads} "
        "--sample-threads {config[sample_threads]}"


rule build_hit_fastas:
    input:
        search=str(SEARCH_DONE),
        indexes=str(INDEX_SENTINEL),
        manifest=checkpoint_manifest,
        target=group_value("fasta"),
        initial=group_value("initial_bed"),
    output:
        passing=str(WORK / "groups" / "{group}" / "hits.pass.fa"),
        filtered=str(WORK / "groups" / "{group}" / "hits.filtered.fa"),
        bed=str(WORK / "groups" / "{group}" / "hits.bed"),
        report=str(WORK / "groups" / "{group}" / "hits.report.json"),
    threads:
        THREADS
    params:
        name_prefix=config["name_prefix"],
        reference_sample=config["reference_sample"],
    shell:
        "{PYTHON} {BUILD_HITS_SCRIPT} --group {wildcards.group} "
        "--name-prefix {params.name_prefix:q} --mode {Q_MODE} --search-dir {Q_SEARCH_DIR} "
        "--queries {Q_QUERY_TABLE} --reference-sample {params.reference_sample:q} "
        "--target-fasta {input.target:q} --initial-bed {input.initial:q} "
        "--consistent-script {CONSISTENT_SCRIPT} "
        "--anchor-size {config[anchor_size]} --edge-distance {config[edge_distance]} "
        "--threads {threads} --passing {output.passing:q} --filtered {output.filtered:q} "
        "--regions-bed {output.bed:q} --report {output.report:q}"


rule gfixbreaks:
    input:
        fasta=str(WORK / "groups" / "{group}" / "hits.pass.fa"),
        target=group_value("fasta"),
    output:
        fixed=str(WORK / "groups" / "{group}" / "fixed.fa"),
        loci=str(WORK / "groups" / "{group}" / "fixed.fa_loci.txt"),
        graph_base=str(WORK / "groups" / "{group}" / "fixed.fa_loci.txt.fasta"),
        graph=str(WORK / "groups" / "{group}" / "fixed.fa_loci.txt.fasta_graph.FA"),
        alignment=str(WORK / "groups" / "{group}" / "fixed.fa_loci.txt.fasta_allgraphalign.out"),
    threads:
        THREADS
    params:
        reference_prefix=(
            "--reference-prefix " + q(config["reference_prefix"])
            if config.get("reference_prefix") else ""
        ),
    shell:
        "{PYTHON} {RUN_GFIXBREAKS_SCRIPT} -i {input.fasta:q} -o {output.fixed:q} "
        "-q {Q_QUERY_TABLE} -r {Q_REFERENCE} --target-fasta {input.target:q} "
        "{params.reference_prefix} "
        "--group-prefix {config[name_prefix]}{wildcards.group} "
        "--anchor-size {config[anchor_size]} --threads {threads} "
        "--scripts-dir {Q_SCRIPTS} --tools-dir {Q_TOOLS}"


rule add_filtered_regions:
    input:
        fixed=str(WORK / "groups" / "{group}" / "fixed.fa"),
        filtered=str(WORK / "groups" / "{group}" / "hits.filtered.fa"),
    output:
        temp(str(WORK / "groups" / "{group}" / "augmented.fa"))
    shell:
        "{PYTHON} {ADDREGION_SCRIPT} -i {input.fixed:q} -a {input.filtered:q} "
        "-o {output:q} -q {Q_QUERY_TABLE} "
        "--prefix {config[name_prefix]}{wildcards.group}"


rule exclusive_kmers:
    input:
        fasta=str(WORK / "groups" / "{group}" / "augmented.fa")
    output:
        temp(str(KMERS_DIR / "{group}.exclusive.kmers.txt"))
    threads:
        THREADS
    shell:
        "{PYTHON} {SELECT_KMERS_SCRIPT} -i {input.fasta:q} "
        "-q {Q_SEARCH_TABLE} -o {output:q} --scripts-dir {Q_TOOLS} --threads {threads}"


rule packedrun:
    input:
        fasta=str(WORK / "groups" / "{group}" / "fixed.fa"),
        kmers=str(KMERS_DIR / "{group}.exclusive.kmers.txt"),
        reference=group_value("fasta"),
        graph_base=str(WORK / "groups" / "{group}" / "fixed.fa_loci.txt.fasta"),
        graph=str(WORK / "groups" / "{group}" / "fixed.fa_loci.txt.fasta_graph.FA"),
        alignment=str(WORK / "groups" / "{group}" / "fixed.fa_loci.txt.fasta_allgraphalign.out"),
    output:
        str(WORK / "groups" / "{group}" / "matrix.txt")
    params:
        reference_sample=config["reference_sample"],
    shell:
        "{PYTHON} {PACKEDRUN_SCRIPT} -i {input.fasta:q} -k {input.kmers:q} "
        "-r {input.reference:q} -g {input.graph_base:q} -a {input.alignment:q} -o {output:q} "
        "--scripts-dir {Q_SCRIPTS} --tools-dir {Q_TOOLS} "
        "--reference-sample {params.reference_sample:q}"


rule combine_fastas:
    input:
        fasta_shards
    output:
        str(FINAL_FASTA)
    run:
        subprocess.run(
            [sys.executable, COMBINE_FASTAS_SCRIPT, "-o", output[0], *input],
            check=True,
        )


rule combine_matrices:
    input:
        matrix_shards
    output:
        str(MATRIX)
    run:
        subprocess.run(
            [sys.executable, str(SCRIPTS / "combine_matrices.py"), "-o", output[0], *input],
            check=True,
        )


rule matrix_index:
    input:
        str(MATRIX)
    output:
        str(MATRIX_INDEX)
    params:
        gff=config.get("gff3", "")
    run:
        command = [sys.executable, str(SCRIPTS / "matrixindex.py"), "-i", input[0], "-r", "1"]
        if params.gff:
            command.extend(["--gene", params.gff])
        subprocess.run(command, check=True)
