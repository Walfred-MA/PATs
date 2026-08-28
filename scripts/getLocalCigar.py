#!/usr/bin/env python3

from __future__ import annotations

import argparse
import os
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Tuple


# The legacy parser deliberately ignores optional sequence payloads in CIGARs.
LEGACY_CIGAR_RE = r"(\d+)([=MXDIHS])"

# The combined-file mode preserves optional base payloads while clipping.
CIGAR_RE = re.compile(r"(\d+)([M=XDISH])([A-Za-z]*)")
SPAN_RE = re.compile(r"\d+_\d+(?:;\d+_\d+)*")
QUERY_OPS = frozenset("M=XIS")
REFERENCE_OPS = frozenset("MXHD=")


@dataclass(frozen=True)
class Coordinate:
    contig: str
    start: int
    end: int
    strand: int


@dataclass(frozen=True)
class CigarOp:
    length: int
    code: str
    sequence: str = ""


@dataclass(frozen=True)
class CigarSegment:
    path: str
    ops: Tuple[CigarOp, ...]


@dataclass(frozen=True)
class LocalSegment:
    segment: CigarSegment
    query_start: int
    query_end: int


def parse_coord_token(token: str) -> Coordinate:
    token = token.strip()
    strand = 1
    if token.endswith(("+", "-")):
        strand = 1 if token[-1] == "+" else -1
        token = token[:-1]

    try:
        contig, interval = token.rsplit(":", 1)
    except ValueError as exc:
        raise ValueError(f"cannot parse coordinate token: {token}") from exc

    match = re.fullmatch(r"(-?\d+)-(-?\d+)", interval)
    if match is None:
        raise ValueError(f"cannot parse coordinate token: {token}")
    start, end = (int(value) for value in match.groups())
    if start < 0 or end < start:
        raise ValueError(f"invalid 0-based, right-open coordinate: {token}")
    return Coordinate(contig, start, end, strand)


# ---------------------------------------------------------------------------
# Original -i/-s/-a mode
# ---------------------------------------------------------------------------

def legacy_cigar_cutqrange(fullcigar, qranges):
    def finish_currpath(path, cigar_tokens, right_clip, strand=1):
        if len(cigar_tokens) == 0:
            return ""

        toks = list(cigar_tokens)
        if right_clip > 0:
            toks.append(f"{right_clip}H")

        if strand == -1:
            path = path.replace("<", "`").replace(">", "<").replace("`", ">")
            reversed_toks = []
            for tok in reversed(toks):
                match = re.fullmatch(LEGACY_CIGAR_RE, tok)
                if match:
                    reversed_toks.append(f"{match.group(1)}{match.group(2)}")
            toks = reversed_toks

        return path + ":" + "".join(toks)

    allcigars = re.findall(r"[><][^><]+", fullcigar)
    total_q = len(qranges)
    if total_q == 0:
        return []

    results = [[] for _ in range(total_q)]
    results_strd = [value[2] if len(value) > 2 else 1 for value in qranges]

    curr_qindex = 0
    curr_qregion = qranges[0]
    curr_results = results[0]
    qstart, qend = curr_qregion[0], curr_qregion[1]
    qstrand = curr_qregion[2] if len(curr_qregion) > 2 else 1
    qposi = 0

    for segment_text in allcigars:
        path, cigar_text = segment_text.split(":", 1)
        cigars = [
            (int(length), code)
            for length, code in re.findall(LEGACY_CIGAR_RE, cigar_text)
        ]

        path_len = sum(
            length
            for length, code in cigars
            if code in {"M", "X", "H", "D", "="}
        )
        segment_right_clip = (
            cigars[-1][0] if cigars and cigars[-1][1] == "H" else 0
        )
        rposi = 0
        newcigars = []

        for size, code in cigars:
            lastrposi = rposi
            lastqposi = qposi

            if code in {"=", "M", "X"}:
                qposi += size
                rposi += size
            elif code in {"D", "H"}:
                rposi += size
            elif code in {"I", "S"}:
                qposi += size

            while qstart < qposi:
                if len(newcigars) == 0:
                    if code not in {"I", "D", "H", "S"} and qstart > lastqposi:
                        newrposi = lastrposi + (qstart - lastqposi)
                    else:
                        newrposi = lastrposi
                    if newrposi > 0:
                        newcigars.append(f"{newrposi}H")

                overlap_start = max(lastqposi, qstart)
                overlap_end = min(qposi, qend)
                segment_size = overlap_end - overlap_start

                if code in {"=", "M", "X", "I", "S"}:
                    if segment_size > 0:
                        newcigars.append(f"{segment_size}{code}")
                elif code == "D":
                    if qstart < qposi and qend > lastqposi:
                        newcigars.append(f"{size}D")

                if qposi >= qend:
                    if code in {"=", "M", "X"}:
                        used_rposi = lastrposi + max(
                            0, min(qend, qposi) - lastqposi
                        )
                    else:
                        used_rposi = lastrposi

                    right_clip = path_len - used_rposi
                    curr_results.append(
                        finish_currpath(path, newcigars, right_clip, qstrand)
                    )

                    curr_qindex += 1
                    if curr_qindex >= total_q:
                        for index in range(total_q):
                            if results_strd[index] == -1:
                                results[index] = "".join(results[index][::-1])
                            else:
                                results[index] = "".join(results[index])
                        return results

                    newcigars = []
                    curr_qregion = qranges[curr_qindex]
                    curr_results = results[curr_qindex]
                    qstart, qend = curr_qregion[0], curr_qregion[1]
                    qstrand = curr_qregion[2] if len(curr_qregion) > 2 else 1
                else:
                    break

        if len(newcigars):
            curr_results.append(
                finish_currpath(path, newcigars, segment_right_clip, qstrand)
            )

    for index in range(total_q):
        if results_strd[index] == -1:
            results[index] = "".join(results[index][::-1])
        else:
            results[index] = "".join(results[index])
    return results


def legacy_getoverlaps(columns, input_regions):
    merge_coordinate = parse_coord_token(columns[1])
    overlap_names = [
        name
        for name, region in input_regions.items()
        if region["contig"] == merge_coordinate.contig
        and max(merge_coordinate.end, region["end"])
        - min(merge_coordinate.start, region["start"])
        + 10
        < (region["end"] - region["start"])
        + (merge_coordinate.end - merge_coordinate.start)
    ]

    outputs = []
    for input_name in overlap_names:
        region = input_regions[input_name]
        input_start = region["start"]
        input_end = region["end"]
        input_strand = region["strand"]

        if merge_coordinate.strand == 1:
            q0 = input_start - merge_coordinate.start
            q1 = input_end - merge_coordinate.start
        else:
            q0 = merge_coordinate.end - input_end
            q1 = merge_coordinate.end - input_start

        merge_length = merge_coordinate.end - merge_coordinate.start
        if q1 <= 0 or q0 >= merge_length:
            continue

        q0 = max(0, q0)
        q1 = min(merge_length, q1)
        strand = 1 if input_strand == merge_coordinate.strand else -1
        overlap_score = (
            (region["end"] - region["start"])
            + (merge_coordinate.end - merge_coordinate.start)
            - max(merge_coordinate.end, region["end"])
            + min(merge_coordinate.start, region["start"])
        )
        outputs.append([q0, q1, strand, input_name, overlap_score])
    return outputs


def legacy_getpathinfo(sub_full):
    segments = re.findall(r"[><][^><]+", sub_full)
    if not segments:
        return []

    path_list, cigar_list, reference_spans, query_spans = [], [], [], []
    query_position = 0

    for segment_text in segments:
        path, cigar = segment_text.split(":", 1)
        path_list.append(path)
        cigar_list.append(cigar)
        ops = [
            (int(length), code)
            for length, code in re.findall(LEGACY_CIGAR_RE, cigar)
        ]
        path_length = sum(
            length
            for length, code in ops
            if code in {"M", "X", "H", "D", "="}
        )
        query_size = sum(
            length
            for length, code in ops
            if code in {"M", "X", "I", "=", "S"}
        )
        left_hard_clip = ops[0][0] if ops and ops[0][1] == "H" else 0
        right_hard_clip = ops[-1][0] if ops and ops[-1][1] == "H" else 0
        reference_spans.append(
            f"{left_hard_clip}_{path_length - right_hard_clip}"
        )
        query_spans.append(f"{query_position}_{query_position + query_size}")
        query_position += query_size

    return (
        "".join(path_list),
        "".join(cigar_list),
        ";".join(reference_spans),
        ";".join(query_spans),
    )


def legacy_readinput(inputfile):
    input_names = []
    input_regions = {}
    with open(inputfile) as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            columns = line[1:].strip().split()
            name, coordinate_text = columns[0], columns[1]
            coordinate = parse_coord_token(coordinate_text)
            input_names.append(name)
            input_regions[name] = {
                "ref": coordinate.contig,
                "start": coordinate.start,
                "end": coordinate.end,
                "strand": coordinate.strand,
                "contig": coordinate.contig,
            }
    return input_names, input_regions


def legacy_readaligns(alignfile):
    alignments = {}
    with open(alignfile) as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 3:
                alignments[parts[0]] = parts[2]
    return alignments


def run_legacy_mode(inputfile: str, seqfile: str, alignfile: str, outputfile: str):
    alignments = legacy_readaligns(alignfile)
    input_names, input_regions = legacy_readinput(inputfile)
    overlap_sizes = {}
    output_cigars = {}

    with open(seqfile) as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            columns = line[1:].strip().split()
            if len(columns) < 2 or columns[0] not in alignments:
                continue

            overlaps = sorted(legacy_getoverlaps(columns, input_regions))
            clipped_cigars = legacy_cigar_cutqrange(
                alignments[columns[0]], overlaps
            )
            for overlap, cigar in zip(overlaps, clipped_cigars):
                if (
                    overlap[3] in overlap_sizes
                    and overlap[4] <= overlap_sizes[overlap[3]]
                ):
                    continue
                overlap_sizes[overlap[3]] = overlap[4]
                output_cigars[overlap[3]] = cigar

    with open(outputfile, "w") as output:
        for name in input_names:
            if name not in output_cigars:
                continue
            full_cigar = output_cigars[name]
            paths, _cigar, reference_spans, query_spans = legacy_getpathinfo(
                full_cigar
            )
            output.write(
                "\t".join(
                    [name, paths, full_cigar, reference_spans, query_spans]
                )
                + "\n"
            )


# ---------------------------------------------------------------------------
# Combined alignment-row mode: -i/-o only
# ---------------------------------------------------------------------------

def parse_cigar_ops(text: str) -> Tuple[CigarOp, ...]:
    ops: List[CigarOp] = []
    position = 0
    for match in CIGAR_RE.finditer(text):
        if match.start() != position:
            raise ValueError(f"cannot parse CIGAR near: {text[position:]}")
        length = int(match.group(1))
        if length:
            ops.append(CigarOp(length, match.group(2), match.group(3)))
        position = match.end()
    if position != len(text) or not ops:
        raise ValueError(f"cannot parse CIGAR: {text}")
    return tuple(ops)


def parse_graph_cigar(full_cigar: str) -> Tuple[CigarSegment, ...]:
    segments: List[CigarSegment] = []
    position = 0
    for match in re.finditer(r"[><][^><]+", full_cigar):
        if match.start() != position:
            raise ValueError(
                f"cannot parse graphical CIGAR near: {full_cigar[position:]}"
            )
        segment_text = match.group(0)
        try:
            path, cigar = segment_text.split(":", 1)
        except ValueError as exc:
            raise ValueError(
                f"graphical CIGAR segment has no colon: {segment_text}"
            ) from exc
        segments.append(CigarSegment(path, parse_cigar_ops(cigar)))
        position = match.end()
    if position != len(full_cigar) or not segments:
        raise ValueError(f"cannot parse graphical CIGAR: {full_cigar}")
    return tuple(segments)


def parse_spans(
    text: str, require_sorted: bool = False
) -> Tuple[Tuple[int, int], ...]:
    if not SPAN_RE.fullmatch(text):
        raise ValueError(f"invalid span list: {text}")
    spans = []
    previous_end = -1
    for item in text.split(";"):
        start_text, end_text = item.split("_", 1)
        start, end = int(start_text), int(end_text)
        if (
            start < 0
            or end < start
            or (require_sorted and start < previous_end)
        ):
            raise ValueError(f"invalid or unsorted span list: {text}")
        spans.append((start, end))
        previous_end = end
    return tuple(spans)


def op_query_length(op: CigarOp) -> int:
    return op.length if op.code in QUERY_OPS else 0


def op_reference_length(op: CigarOp) -> int:
    return op.length if op.code in REFERENCE_OPS else 0


def segment_query_length(segment: CigarSegment) -> int:
    return sum(op_query_length(op) for op in segment.ops)


def segment_reference_length(segment: CigarSegment) -> int:
    return sum(op_reference_length(op) for op in segment.ops)


def slice_op_sequence(op: CigarOp, offset: int, length: int) -> str:
    if not op.sequence:
        return ""
    if len(op.sequence) == op.length:
        return op.sequence[offset : offset + length]
    return op.sequence


def merge_ops(ops: Iterable[CigarOp]) -> Tuple[CigarOp, ...]:
    merged: List[CigarOp] = []
    for op in ops:
        if op.length <= 0:
            continue
        if (
            merged
            and merged[-1].code == op.code
            and not merged[-1].sequence
            and not op.sequence
        ):
            previous = merged[-1]
            merged[-1] = CigarOp(previous.length + op.length, op.code)
        else:
            merged.append(op)
    return tuple(merged)


def reverse_complement(text: str) -> str:
    return text.translate(str.maketrans("ACGTacgtNn", "TGCAtgcaNn"))[::-1]


def reverse_segment(segment: CigarSegment) -> CigarSegment:
    path = ("<" if segment.path.startswith(">") else ">") + segment.path[1:]
    ops = (
        CigarOp(
            op.length,
            op.code,
            reverse_complement(op.sequence) if op.sequence else "",
        )
        for op in reversed(segment.ops)
    )
    return CigarSegment(path, merge_ops(ops))


def clip_segment(
    segment: CigarSegment,
    segment_query_start: int,
    query_start: int,
    query_end: int,
) -> Optional[CigarSegment]:
    query_position = segment_query_start
    prefix_reference = 0
    suffix_reference = 0
    selected_query = 0
    selected_ops: List[CigarOp] = []

    for op in segment.ops:
        query_length = op_query_length(op)
        reference_length = op_reference_length(op)
        if query_length:
            op_start = query_position
            op_end = op_start + query_length
            before = min(op.length, max(0, query_start - op_start))
            selected_end = min(op.length, max(0, query_end - op_start))
            selected = max(0, selected_end - before)
            after = op.length - before - selected

            if reference_length:
                prefix_reference += before
                suffix_reference += after
            if selected:
                selected_ops.append(
                    CigarOp(
                        selected,
                        op.code,
                        slice_op_sequence(op, before, selected),
                    )
                )
                selected_query += selected
            query_position = op_end
            continue

        if query_position <= query_start:
            prefix_reference += reference_length
        elif query_position >= query_end:
            suffix_reference += reference_length
        else:
            selected_ops.append(op)

    if selected_query == 0:
        return None

    output_ops: List[CigarOp] = []
    if prefix_reference:
        output_ops.append(CigarOp(prefix_reference, "H"))
    output_ops.extend(selected_ops)
    if suffix_reference:
        output_ops.append(CigarOp(suffix_reference, "H"))

    clipped = CigarSegment(segment.path, merge_ops(output_ops))
    if segment_reference_length(clipped) != segment_reference_length(segment):
        raise AssertionError(f"reference length changed while clipping {segment.path}")
    return clipped


def clip_graph_cigar_with_spans(
    full_cigar: str,
    query_spans: Sequence[Tuple[int, int]],
    query_start: int,
    query_end: int,
    reverse: bool,
) -> Tuple[LocalSegment, ...]:
    segments = parse_graph_cigar(full_cigar)
    if len(segments) != len(query_spans):
        raise ValueError(
            f"CIGAR has {len(segments)} graph segments but query-span column "
            f"has {len(query_spans)} spans"
        )
    if query_start < 0 or query_end <= query_start:
        raise ValueError(f"invalid local query interval [{query_start}, {query_end})")

    localized: List[LocalSegment] = []
    for segment, (span_start, span_end) in zip(segments, query_spans):
        overlap_start = max(query_start, span_start)
        overlap_end = min(query_end, span_end)
        if overlap_start >= overlap_end:
            continue

        clipped = clip_segment(
            segment, span_start, overlap_start, overlap_end
        )
        if clipped is None:
            continue
        if reverse:
            local_start = query_end - overlap_end
            local_end = query_end - overlap_start
        else:
            local_start = overlap_start - query_start
            local_end = overlap_end - query_start
        localized.append(LocalSegment(clipped, local_start, local_end))

    if not localized:
        raise ValueError(
            f"no graph CIGAR segment overlaps query interval "
            f"[{query_start}, {query_end})"
        )

    if reverse:
        localized = [
            LocalSegment(
                reverse_segment(item.segment),
                item.query_start,
                item.query_end,
            )
            for item in reversed(localized)
        ]

    for item in localized:
        expected = item.query_end - item.query_start
        observed = segment_query_length(item.segment)
        if observed != expected:
            raise AssertionError(
                f"clipped {item.segment.path} consumes {observed} query bases; "
                f"expected {expected}"
            )
    return tuple(localized)


def leading_hard_clip(ops: Sequence[CigarOp]) -> int:
    total = 0
    for op in ops:
        if op.code != "H":
            break
        total += op.length
    return total


def trailing_hard_clip(ops: Sequence[CigarOp]) -> int:
    total = 0
    for op in reversed(ops):
        if op.code != "H":
            break
        total += op.length
    return total


def format_op(op: CigarOp) -> str:
    return f"{op.length}{op.code}{op.sequence}"


def format_local_alignment(
    name: str, localized: Sequence[LocalSegment]
) -> str:
    paths = "".join(item.segment.path for item in localized)
    full_cigar = "".join(
        f"{item.segment.path}:"
        f"{''.join(format_op(op) for op in item.segment.ops)}"
        for item in localized
    )
    reference_spans = []
    query_spans = []
    for item in localized:
        segment = item.segment
        reference_length = segment_reference_length(segment)
        reference_start = leading_hard_clip(segment.ops)
        reference_end = reference_length - trailing_hard_clip(segment.ops)
        reference_spans.append(f"{reference_start}_{reference_end}")
        query_spans.append(f"{item.query_start}_{item.query_end}")
    return "\t".join(
        [
            name,
            paths,
            full_cigar,
            ";".join(reference_spans),
            ";".join(query_spans),
        ]
    )


def find_graph_columns(
    fields: Sequence[str],
) -> Tuple[str, str, str, str, str, int]:
    for index in range(6, len(fields) - 4):
        name = fields[index].strip().lstrip(">")
        paths = fields[index + 1].strip()
        full_cigar = fields[index + 2].strip()
        reference_spans = fields[index + 3].strip()
        query_spans = fields[index + 4].strip()
        if (
            name
            and not name.startswith(("<", ">"))
            and paths.startswith(("<", ">"))
            and full_cigar.startswith(("<", ">"))
            and ":" in full_cigar
            and SPAN_RE.fullmatch(reference_spans)
            and SPAN_RE.fullmatch(query_spans)
        ):
            return (
                name,
                paths,
                full_cigar,
                reference_spans,
                query_spans,
                index,
            )
    raise ValueError("cannot locate name/path/CIGAR/reference-span/query-span columns")


def find_full_query_header(
    fields: Sequence[str],
    graph_name_index: int,
    target_contig: str,
    fallback_name: str,
) -> Tuple[str, Coordinate]:
    packed_name = None
    candidates = []
    for field in fields[6:graph_name_index]:
        for token in field.split():
            if token.startswith(">") and len(token) > 1 and packed_name is None:
                packed_name = token[1:]
            try:
                coordinate = parse_coord_token(token.lstrip(">"))
            except ValueError:
                continue
            if coordinate.contig == target_contig:
                candidates.append(coordinate)
    if not candidates:
        raise ValueError(
            f"cannot find a full-query coordinate for contig {target_contig}"
        )
    return packed_name or fallback_name, candidates[-1]


def split_combined_fields(line: str) -> List[str]:
    fields = [field.strip() for field in line.rstrip("\n").split("\t")]
    if len(fields) >= 11:
        return fields
    return line.strip().split()


def localize_combined_row(line: str) -> str:
    fields = split_combined_fields(line)
    if len(fields) < 11:
        raise ValueError(f"expected at least 11 columns, found {len(fields)}")

    try:
        target = Coordinate(
            fields[2],
            int(fields[4]),
            int(fields[5]),
            1 if fields[3] == "+" else -1 if fields[3] == "-" else 0,
        )
    except (IndexError, ValueError) as exc:
        raise ValueError("invalid contig/strand/start/end columns 3-6") from exc
    if target.strand == 0:
        raise ValueError(f"invalid target strand: {fields[3]}")
    if target.start < 0 or target.end <= target.start:
        raise ValueError(
            f"invalid target interval: {target.contig}:{target.start}-{target.end}"
        )

    (
        graph_name,
        paths,
        full_cigar,
        reference_span_text,
        query_span_text,
        graph_name_index,
    ) = find_graph_columns(fields)
    name, full_coordinate = find_full_query_header(
        fields,
        graph_name_index,
        target.contig,
        graph_name,
    )

    if (
        target.start < full_coordinate.start
        or target.end > full_coordinate.end
    ):
        raise ValueError(
            f"target {target.contig}:{target.start}-{target.end} is outside "
            f"full query {full_coordinate.contig}:"
            f"{full_coordinate.start}-{full_coordinate.end}"
        )

    if full_coordinate.strand == 1:
        query_start = target.start - full_coordinate.start
        query_end = target.end - full_coordinate.start
    else:
        query_start = full_coordinate.end - target.end
        query_end = full_coordinate.end - target.start

    segments = parse_graph_cigar(full_cigar)
    parsed_paths = "".join(segment.path for segment in segments)
    if paths != parsed_paths:
        raise ValueError(
            f"path column does not match graphical CIGAR paths: "
            f"{paths} != {parsed_paths}"
        )

    reference_spans = parse_spans(reference_span_text)
    query_spans = parse_spans(query_span_text, require_sorted=True)
    if len(reference_spans) != len(segments):
        raise ValueError(
            f"CIGAR has {len(segments)} graph segments but reference-span "
            f"column has {len(reference_spans)} spans"
        )

    localized = clip_graph_cigar_with_spans(
        full_cigar,
        query_spans,
        query_start,
        query_end,
        target.strand != full_coordinate.strand,
    )
    return format_local_alignment(name, localized)


def run_combined_mode(
    inputfile: str,
    outputfile: str,
) -> Tuple[int, int, int]:
    input_path = Path(inputfile)
    output_path = Path(outputfile)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(f"{output_path}.tmp.{os.getpid()}")
    read_count = 0
    written_count = 0
    skipped_count = 0

    try:
        with input_path.open() as source, temporary.open("w") as output:
            for line_number, line in enumerate(source, 1):
                if not line.strip() or line.startswith("#"):
                    continue
                read_count += 1
                try:
                    output.write(localize_combined_row(line) + "\n")
                except (ValueError, AssertionError) as exc:
                    skipped_count += 1
                    print(
                        f"[getLocalCigar.py] warning: skipped "
                        f"{input_path}:{line_number}: {exc}",
                        file=sys.stderr,
                        flush=True,
                    )
                    continue
                written_count += 1
        os.replace(temporary, output_path)
    finally:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass
    return read_count, written_count, skipped_count


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Extract local graphical CIGARs. Use -i/-o for corrected combined "
            "*_fasta.txt_align.out.txt rows, or -i/-s/-a/-o for the original "
            "three-input mode."
        )
    )
    parser.add_argument("-i", "--input", required=True, help="input file")
    parser.add_argument("-o", "--output", required=True, help="output file")
    parser.add_argument(
        "-s",
        "--seq",
        help="legacy mode: FASTA containing the full query headers",
    )
    parser.add_argument(
        "-a",
        "--align",
        help="legacy mode: five-column full graphical alignment file",
    )
    args = parser.parse_args()

    if args.seq is None and args.align is None:
        read_count, written_count, skipped_count = run_combined_mode(
            args.input,
            args.output,
        )
        print(
            f"[getLocalCigar.py] combined rows={read_count} "
            f"localized={written_count} skipped={skipped_count} "
            f"output={args.output}",
            flush=True,
        )
        return

    if args.seq is not None and args.align is not None:
        run_legacy_mode(args.input, args.seq, args.align, args.output)
        return

    parser.error("-s/--seq and -a/--align must be supplied together")


if __name__ == "__main__":
    main()

