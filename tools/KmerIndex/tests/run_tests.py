#!/usr/bin/env python3

import pathlib
import struct
import subprocess
import sys
import tempfile


UINT32_MAX = (1 << 32) - 1


def canonical(sequence):
    table = {"A": 0, "C": 1, "G": 2, "T": 3}
    forward = 0
    for base in sequence:
        forward = (forward << 2) | table[base]
    reverse = 0
    for base in reversed(sequence):
        reverse = (reverse << 2) | (3 - table[base])
    return max(forward, reverse)


def write_fasta_and_fai(directory, contigs, width=17):
    fasta = directory / "input.fa"
    fai = directory / "input.fa.fai"
    rows = []
    with fasta.open("wb") as output:
        for name, sequence in contigs:
            output.write(f">{name}\n".encode())
            offset = output.tell()
            for start in range(0, len(sequence), width):
                output.write(sequence[start : start + width].encode() + b"\n")
            rows.append((name, len(sequence), offset, width, width + 1))
    with fai.open("w", encoding="ascii") as output:
        for row in rows:
            output.write("\t".join(map(str, row)) + "\n")
    return fasta, fai


def parse_index(path, maxcover):
    data = path.read_bytes()
    records = []
    offset = 0
    marker_floor = UINT32_MAX - maxcover
    while offset < len(data):
        if offset + 12 > len(data):
            raise AssertionError("truncated record")
        kmer, value = struct.unpack_from("<QI", data, offset)
        offset += 12
        if marker_floor <= value <= UINT32_MAX - 2:
            count = UINT32_MAX - value
            positions = list(struct.unpack_from(f"<{count}I", data, offset))
            offset += 4 * count
        else:
            positions = [value]
        records.append((kmer, positions))
    return records


def run(binary, *arguments, ok=True):
    result = subprocess.run(
        [binary, *map(str, arguments)], text=True, capture_output=True
    )
    if ok and result.returncode != 0:
        raise AssertionError(result.stderr)
    if not ok and result.returncode == 0:
        raise AssertionError("command unexpectedly succeeded")
    return result


def main():
    binary = str(pathlib.Path(sys.argv[1]).resolve())
    with tempfile.TemporaryDirectory() as temporary:
        directory = pathlib.Path(temporary)
        repeated = "A" * 31
        unique = "ACGTTGCATGTCAGTACGATCGTACCTAGCA"
        contigs = [
            ("one", repeated + "N" + unique),
            ("two", repeated),
            ("short", "ACGT"),
        ]
        fasta, fai = write_fasta_and_fai(directory, contigs)

        index_one = directory / "one.bin"
        index_four = directory / "four.bin"
        run(binary, "--fasta", fasta, "--fai", fai, "--output", index_one,
            "--threads", 1, "--maxcover", 10)
        run(binary, "--fasta", fasta, "--fai", fai, "--output", index_four,
            "--threads", 4, "--maxcover", 10)
        if index_one.read_bytes() != index_four.read_bytes():
            raise AssertionError("thread counts produced different output")

        records = dict(parse_index(index_one, 10))
        repeated_key = canonical(repeated)
        unique_key = canonical(unique)
        expected_second_contig_position = len(contigs[0][1])
        if records[repeated_key] != [0, expected_second_contig_position]:
            raise AssertionError(f"wrong repeated positions: {records[repeated_key]}")
        if records[unique_key] != [32]:
            raise AssertionError(f"wrong unique position: {records[unique_key]}")

        covered = directory / "covered.bin"
        run(binary, "--fasta", fasta, "--output", covered,
            "--maxcover", 1)
        covered_records = dict(parse_index(covered, 1))
        if repeated_key in covered_records or covered_records[unique_key] != [32]:
            raise AssertionError("maxcover filtering failed")

        spanned = directory / "spanned.bin"
        run(binary, "--fasta", fasta, "--output", spanned,
            "--maxcover", 10, "--maxspan", 0)
        if spanned.read_bytes():
            raise AssertionError("maxspan filtering failed")

        invalid = run(binary, "--fasta", fasta, "--output", directory / "bad.bin",
                      "--maxcover", 65536, ok=False)
        if "65535" not in invalid.stderr:
            raise AssertionError("maxcover validation message is missing")

    print("all tests passed")


if __name__ == "__main__":
    main()
