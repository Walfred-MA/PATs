from __future__ import annotations

from collections import Counter
import csv
import gzip
import json
import os
from pathlib import Path
import random
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from patslib.transcripts.catalog import candidates, count_kmers, create_catalog, partition_genes, reverse_complement
from patslib.transcripts.input import Fasta, Sample, alternatives, read_calls, read_manifest
from patslib.core import read_fasta


# Native ParaTyper 3.9.2 schema; column order may change without conversion.
FIELDS = """GENE_index transcript_id gene_id gene_name transcript_type ifmane_transcript
weighted_score raw_exon_score total_aligned_bases mean_identity mean_exon_coverage
expected_exons found_exons missing_exons fraction_expected_found found_exon_numbers
query_id query_contig query_start query_end strand exon_query_coordinates
exon_reference_coordinates exon_ids exon_coverages exon_identities alignment_AS
call_status tie_group_id tie_count assignment_status pipeline_version
merged_interval_coordinates merged_interval_gene_scores merged_interval_score_exon_ids
merged_interval_score_exon_lengths transcript_exon_length model_type inserted_exons
insertion_penalty inserted_exon_numbers inserted_exon_query_coordinates inserted_query_blocks
insertion_run_unique_exons""".split()


def decode(text):
    result = 0
    for char in text:
        result = 64 * result + ord(char) - ord("0")
    return result


def decode_list(text):
    result = []
    for entry in text.split(","):
        if not entry:
            continue
        value, *multiplicity = entry.split("*")
        start, *length = value.split("~")
        for index in range(decode(start), decode(start) + (decode(length[0]) if length else 0) + 1):
            result.extend([index] * (decode(multiplicity[0]) if multiplicity else 1))
    return result


def naive_counts(parts):
    result = Counter()
    for part in parts:
        for i in range(len(part) - 30):
            word = part[i:i+31].upper()
            if set(word) <= set("ACGT"):
                word = max(word, reverse_complement(word))
                result[int("".join(str("ACGT".index(b)) for b in word), 4)] += 1
    return result


class TranscriptTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.directory = Path(self.temp.name)
        self.sequence = "".join(random.Random(7).choices("ACGT", k=1600))
        self.fasta = self.directory / "assembly.fa"
        self.write_fasta(self.sequence)
        self.calls = self.directory / "sample.transcript_calls.tsv"

    def write_fasta(self, sequence):
        self.fasta.write_text(">chr1\n" + "\n".join(sequence[i:i+61] for i in range(0, len(sequence), 61)) + "\n")

    def row(self, gene="GA", transcript="TA", exons=((100, 200), (350, 450)), **changes):
        parent = transcript == gene
        row = dict.fromkeys(FIELDS, "")
        row.update(GENE_index="GENE_1", transcript_id=transcript, gene_id=gene, gene_name=gene,
                   model_type="full_gene" if parent else "transcript", query_contig="chr1", query_id="chr1",
                   query_start=str(min(s for s, _ in exons)), query_end=str(max(e for _, e in exons)), strand="+",
                   exon_query_coordinates=",".join(f"{s}-{e}" for s, e in exons),
                   found_exon_numbers=",".join(str(i+1) for i in range(len(exons))),
                   exon_coverages=",".join("100" for _ in exons), call_status="complete",
                   assignment_status="unique", tie_count="1", pipeline_version="3.9.2")
        row.update(changes)
        return row

    def write_calls(self, rows, path=None):
        path = path or self.calls
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(reversed(FIELDS)), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        return path

    def sample(self):
        return Sample("S1", self.fasta, self.calls, None)

    def run_cli(self, mode, rows, expected=0, extra=()):
        self.write_calls(rows)
        manifest = self.directory / "manifest.tsv"
        manifest.write_text("sample\tassembly\tcalls\nS1\tassembly.fa\tsample.transcript_calls.tsv\n")
        output = self.directory / mode
        result = subprocess.run([sys.executable, str(ROOT / "pats.py"), "transcripts", "run", "--mode", mode,
                                 "--manifest", str(manifest), "-o", str(output), *extra], capture_output=True, text=True)
        self.assertEqual(result.returncode, expected, result.stdout + result.stderr)
        return output, result

    def test_fasta_seek_wrapping_crlf_and_bounds(self):
        self.fasta.write_bytes(b">chr1\r\n" + b"\r\n".join(self.sequence[i:i+61].encode() for i in range(0, 1600, 61)))
        fasta = Fasta(self.fasta)
        self.addCleanup(fasta.close)
        self.assertEqual(fasta.fetch("chr1", 57, 143), self.sequence[57:143])
        self.assertEqual(fasta.fetch("chr1", 0, 1600), self.sequence)
        with self.assertRaises(ValueError):
            fasta.fetch("chr1", 1500, 1700)
        with self.assertRaises(ValueError):
            fasta.fetch("missing", 0, 1)

    def test_fasta_external_index(self):
        Path(str(self.fasta) + ".fai").write_text("chr1\t1600\t6\t61\t62\n")
        fasta = Fasta(self.fasta)
        self.addCleanup(fasta.close)
        self.assertEqual(fasta.fetch("chr1", 57, 143), self.sequence[57:143])

    def test_native_schema_ties_shared_scalars_and_shared_units(self):
        row = self.row(gene="GA&GB;GC", transcript="TA;TC",
                       exon_query_coordinates="100-200,350-450;100-200,350-450",
                       found_exon_numbers="1,2;1,2", exon_coverages="100,100;100,100",
                       assignment_status="tied", tie_count="2")
        self.write_calls([row])
        expanded = list(read_calls(self.calls))[0][2]
        self.assertEqual([r["gene_id"] for r in expanded], ["GA&GB", "GC"])
        self.assertEqual(expanded[1]["strand"], "+")
        row["exon_query_coordinates"] = "100-200,350-450"
        with self.assertRaises(ValueError):
            alternatives(row)

    def test_wes_anchor_union_and_separate_intervals(self):
        parent = self.row(transcript="GA", exons=((5, 90), (115, 200), (1500, 1590)))
        child = self.row(exons=((5, 90), (115, 200), (1500, 1590)))
        duplicate = self.row(transcript="TB", exons=((5, 90), (120, 190)))
        self.write_calls([parent, child, duplicate])
        models = list(candidates(self.sample(), "wes", Counter()))
        self.assertEqual(len(models), 1)
        parts, metadata = models[0]
        self.assertEqual(metadata["segments"], [[(0, 230)], [(1470, 1600)]])
        self.assertEqual(parts, [self.sequence[:230], self.sequence[1470:]])
        self.assertEqual(count_kmers(parts), naive_counts(parts))
        self.assertNotEqual(count_kmers(parts), naive_counts(["".join(parts)]))

    def test_rna_minus_strand_and_missing_exon_break(self):
        row = self.row(strand="-", found_exon_numbers="2,1")
        self.write_calls([self.row(transcript="GA", strand="-", found_exon_numbers="2,1"), row])
        parts, metadata = list(candidates(self.sample(), "rna", Counter()))[0]
        self.assertEqual(parts, [reverse_complement(self.sequence[100:200] + self.sequence[350:450])])
        self.assertEqual(metadata["segments"], [[(350, 450), (100, 200)]])
        row.update(found_exon_numbers="3,1", call_status="partial")
        self.write_calls([self.row(transcript="GA", strand="-"), row])
        parts, metadata = list(candidates(self.sample(), "rna", Counter()))[0]
        self.assertEqual(len(parts), 2)
        self.assertEqual(metadata["reconstruction"], "segmented")

    def test_partial_alignment_breaks_junction_even_when_call_complete(self):
        self.write_calls([self.row(GENE_index="", exon_coverages="98,100")])
        parts, _ = list(candidates(self.sample(), "rna", Counter()))[0]
        self.assertEqual(parts, [self.sequence[100:200], self.sequence[350:450]])
        with self.assertRaisesRegex(ValueError, "WES requires"):
            list(candidates(self.sample(), "wes", Counter()))

    def test_shared_unit_children_pool_in_wes_and_remain_isoforms_in_rna(self):
        self.write_calls([self.row(gene="GA&GB", transcript="GA&GB"), self.row(), self.row(gene="GB", transcript="TB")])
        wes = list(candidates(self.sample(), "wes", Counter()))
        rna = list(candidates(self.sample(), "rna", Counter()))
        self.assertEqual(len(wes), 1)
        self.assertEqual(len(rna), 2)
        self.assertTrue(all(m["gene_unit"] == "GA&GB" for _, m in rna))

    def test_wes_preserves_owned_exon_blocks_beyond_representative(self):
        parent = self.row(transcript="GA", merged_interval_coordinates="100-225,350-450")
        self.write_calls([parent, self.row()])
        _, metadata = list(candidates(self.sample(), "wes", Counter()))[0]
        self.assertEqual(metadata["segments"], [[(70, 255)], [(320, 480)]])

    def test_partition_strict_threshold_transitivity_and_reverse_complement(self):
        db = create_catalog(self.directory / "partition.sqlite")
        self.addCleanup(db.close)
        # A-B share 101; B-C share 101; A-C share none. D-A share exactly 100.
        sets = [set(range(101)), set(range(202)), set(range(101,202)), set(range(100))]
        for i, values in enumerate(sets, 1):
            db.execute("INSERT INTO genes(id,name) VALUES (?,?)", (i, chr(64+i)))
            db.executemany("INSERT INTO gene_kmers VALUES (?,?)", ((k,i) for k in values))
        groups = partition_genes(db)
        self.assertEqual(groups, [[1,2,3], [4]])
        self.assertEqual(db.execute("SELECT count(*) FROM shared_kmers").fetchone()[0], 100)
        sequence = self.sequence[20:250]
        self.assertEqual(count_kmers([sequence]), count_kmers([reverse_complement(sequence)]))
        self.assertEqual(count_kmers(["a"*50 + "N" + "t"*50]), naive_counts(["A"*50 + "N" + "T"*50]))

    def check_matrix(self, output):
        raw = Path(str(output) + ".matrix.txt").read_bytes()
        indexes = Path(str(output) + ".matrix.txt.index").read_text().splitlines()
        self.assertEqual(indexes[0], "@v2.0.1,support:v1.2.0")
        sequences = {r.name: r.sequence for r in read_fasta(Path(str(output) + ".fa"))}
        annotations = [json.loads(line) for line in Path(str(output) + ".annotations.jsonl").read_text().splitlines()]
        native_expected = []
        for line in indexes[1:]:
            fields = line.split("\t")
            block = raw[int(fields[1]):int(fields[2])].decode().splitlines()
            self.assertTrue(block[0].startswith(fields[0] + "\t"))
            columns = int(fields[7])
            native_expected.append(f"{fields[0]}\t{columns}\t{fields[3]}")
            profiles = [Counter() for _ in range(columns)]
            matrix_size = phylo_size = 0
            row_count = 0
            # Use PATs' existing tree parser as an independent annotation oracle.
            sys.path.insert(0, str(ROOT / "scripts"))
            from matrixcompile import phylotree
            phylo = phylotree(block[2][1:])
            if columns == 1:
                phylo.tree[0].index = 0  # legacy parser does not assign a singleton's index
            for row in block:
                if not row.startswith("&"):
                    continue
                row_count += 1
                self.assertEqual(row[27:38], row.split("\t")[4])
                self.assertEqual(row[39:], row.split("\t")[5])
                self.assertEqual(row[3:6], "!{S")
                self.assertEqual(row[7:13], "000000")
                allele, positive, negative = row[39:].split(" ")
                indexes_row = decode_list(allele)
                positive, negative = decode_list(positive), decode_list(negative)
                self.assertEqual([i+1 for i in positive] + [-i-1 for i in negative],
                                 sorted(phylo.treeencodes(indexes_row), key=lambda x: (x < 0, abs(x))))
                native_expected.append("ROW" + "".join(f"\t{i}" for i in indexes_row))
                matrix_size += len(indexes_row) + 15
                phylo_size += len(positive) + len(negative) + 2
                for i in indexes_row:
                    profiles[i][decode(row[27:38])] += 1
            self.assertEqual(row_count, int(fields[3]))
            self.assertEqual(matrix_size, int(fields[4]))
            self.assertEqual(phylo_size, int(fields[5]))
            self.assertEqual([sum(p.values()) for p in profiles], list(map(int, block[1][1:].split())))
            retained = set().union(*profiles)
            for record in annotations:
                assignment = record["matrix"]
                if assignment and assignment["partition"] == fields[0][1:]:
                    expected = naive_counts([sequences[key] for key in record["fasta_records"]])
                    self.assertEqual(profiles[assignment["column"]], {k:n for k,n in expected.items() if k in retained})
        binary = os.environ.get("CTYPER_READER")
        if binary:
            result = subprocess.run([binary, str(output) + ".matrix.txt"], text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(result.stdout.splitlines(), native_expected)
        return raw

    def test_end_to_end_native_matrix_two_partitions_and_equivalent_isoforms(self):
        rows = [self.row(transcript="GA"), self.row(), self.row(transcript="TB"),
                self.row(gene="GC", transcript="GC", exons=((800,1100),), GENE_index="GENE_2"),
                self.row(gene="GC", transcript="TC", exons=((800,1100),), GENE_index="GENE_2")]
        for mode in ("wes", "rna"):
            output, _ = self.run_cli(mode, rows, extra=("--reference-sample", "S1"))
            self.check_matrix(output)
            report = json.loads(Path(str(output) + ".report.json").read_text())
            self.assertEqual(len(report["partitions"]), 2)
            annotation = [json.loads(line) for line in Path(str(output) + ".annotations.jsonl").read_text().splitlines()]
            self.assertEqual(len(annotation), 2 if mode == "wes" else 3)
            if mode == "rna":
                self.assertEqual(annotation[0]["matrix"], annotation[1]["matrix"])

    def test_multiple_tree_leaves_and_repeat_multiplicity(self):
        self.write_fasta("A"*200 + self.sequence[200:])
        rows = [self.row(transcript="GA", exons=((10,180),(350,500))),
                self.row(exons=((10,180),(350,500))),
                self.row(transcript="TB", exons=((10,100),(350,500))),
                self.row(transcript="TC", exons=((10,180),(350,400)))]
        output, _ = self.run_cli("rna", rows)
        self.check_matrix(output)
        text = Path(str(output) + ".matrix.txt").read_text()
        self.assertIn("*", text)
        self.assertEqual(int(text.splitlines()[0].split("\t")[2]), 3)

    def test_cross_partition_kmers_removed_from_both_matrices(self):
        sequence = self.sequence[:500] + self.sequence[:80] + ("A" if self.sequence[80] != "A" else "C") + self.sequence[581:]
        self.write_fasta(sequence)
        rows = [self.row(transcript="GA", exons=((0,250),)), self.row(exons=((0,250),)),
                self.row(gene="GB", transcript="GB", exons=((500,750),), GENE_index="GENE_2"),
                self.row(gene="GB", transcript="TB", exons=((500,750),), GENE_index="GENE_2")]
        output, _ = self.run_cli("rna", rows)
        self.check_matrix(output)
        report = json.loads(Path(str(output) + ".report.json").read_text())
        self.assertEqual(len(report["partitions"]), 2)
        self.assertGreater(report["shared_partition_kmers_excluded"], 0)
        kmers = [line[27:38] for line in Path(str(output) + ".matrix.txt").read_text().splitlines() if line.startswith("&")]
        self.assertEqual(len(kmers), len(set(kmers)))

    def test_max_partition_candidates_and_failed_force_preserve_outputs(self):
        rows = [self.row(transcript="GA"), self.row(), self.row(transcript="TB", exons=((100,190),(350,450)))]
        output, _ = self.run_cli("rna", rows)
        before = Path(str(output) + ".matrix.txt").read_bytes()
        _, result = self.run_cli("rna", rows, expected=2, extra=("--force", "--max-partition-candidates", "1"))
        self.assertIn("exceed the partition limit", result.stderr)
        self.assertEqual(before, Path(str(output) + ".matrix.txt").read_bytes())

    def test_failure_does_not_publish_and_existing_output_is_preserved(self):
        output, result = self.run_cli("rna", [self.row(GENE_index="", query_end="1700", exon_query_coordinates="1500-1700", found_exon_numbers="1", exon_coverages="100")], expected=2)
        self.assertIn("outside assembly", result.stderr)
        self.assertFalse(Path(str(output) + ".matrix.txt").exists())
        output, _ = self.run_cli("rna", [self.row(GENE_index="")])
        before = Path(str(output) + ".matrix.txt").read_bytes()
        _, result = self.run_cli("rna", [self.row(GENE_index="")], expected=2)
        self.assertIn("output already exists", result.stderr)
        self.assertEqual(before, Path(str(output) + ".matrix.txt").read_bytes())

    def test_empty_table_gzip_manifest_and_duplicate_sample(self):
        self.write_calls([])
        self.assertEqual(list(read_calls(self.calls)), [])
        self.write_calls([self.row(GENE_index="")])
        compressed = self.directory / "calls.tsv.gz"
        with gzip.open(compressed, "wb") as handle:
            handle.write(self.calls.read_bytes())
        self.assertEqual(len(list(read_calls(compressed))), 1)
        manifest = self.directory / "manifest.tsv"
        manifest.write_text("sample\tassembly\tcalls\nS1\tassembly.fa\tcalls.tsv.gz\nS1\tassembly.fa\tcalls.tsv.gz\n")
        with self.assertRaisesRegex(ValueError, "duplicate"):
            read_manifest(manifest)

    def test_tied_parent_alternatives_keep_physical_locus_and_collapse_evidence(self):
        changes = dict(gene="GA;GB", transcript="GA;GB", gene_name="GA;GB", assignment_status="tied",
                       tie_count="2", tie_group_id="TIE_1", exon_query_coordinates="100-200,350-450;100-200,350-450",
                       exon_coverages="100,100;100,100", found_exon_numbers="1,2;1,2")
        parent = self.row(**changes)
        parent["model_type"] = "full_gene;full_gene"
        child = self.row(**{**changes, "transcript": "TA;TB", "model_type": "transcript;transcript"})
        for mode in ("wes", "rna"):
            output, _ = self.run_cli(mode, [parent, child])
            self.check_matrix(output)
            annotations = [json.loads(x) for x in Path(str(output) + ".annotations.jsonl").read_text().splitlines()]
            self.assertEqual(len(annotations), 2)
            self.assertEqual(annotations[0]["locus_id"], annotations[1]["locus_id"])
            self.assertEqual(annotations[0]["matrix"], annotations[1]["matrix"])
            self.assertEqual(annotations[0]["assignment_status"], "tied")

    def test_two_samples_same_contig_and_explicit_fragment_table(self):
        self.write_calls([self.row(transcript="GA"), self.row()])
        self.write_calls([self.row(gene="GC", transcript="GC", exons=((800,1100),),
                                                GENE_index="GENE_9", call_status="partial"),
                                     self.row(gene="GC", transcript="TC", exons=((800,1100),),
                                              GENE_index="GENE_9", call_status="partial")], self.directory / "fragments.tsv")
        manifest = self.directory / "multi.tsv"
        manifest.write_text("sample\tassembly\tcalls\tfragments\n"
                            "S1\tassembly.fa\tsample.transcript_calls.tsv\tfragments.tsv\n"
                            "S2\tassembly.fa\tsample.transcript_calls.tsv\t\n")
        output = self.directory / "multi"
        result = subprocess.run([sys.executable, str(ROOT / "pats.py"), "transcripts", "run", "--mode", "wes",
                                 "--manifest", str(manifest), "-o", str(output)], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.check_matrix(output)
        records = [json.loads(x) for x in Path(str(output) + ".annotations.jsonl").read_text().splitlines()]
        self.assertEqual(len(records), 3)
        self.assertEqual(records[0]["matrix"], records[2]["matrix"])
        self.assertTrue(records[1]["candidate_fragment"])

    def test_native_paratyper_smn_output(self):
        rows = list(read_calls(ROOT / "tests/data/paratyper_smn.transcript_calls.tsv.gz"))
        self.assertEqual(len(rows), 10)
        parents = [row for _, _, alts in rows for row in alts if row["model_type"] == "full_gene"]
        self.assertEqual({(p["gene_name"], p["strand"]) for p in parents}, {("SMN1", "+"), ("SMN2", "-")})

    def test_existing_cli_help_is_unchanged(self):
        result = subprocess.run([sys.executable, str(ROOT / "pats.py"), "--help"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0)
        self.assertIn("--bed", result.stdout)
        self.assertNotIn("transcripts", result.stdout)


if __name__ == "__main__":
    unittest.main()
