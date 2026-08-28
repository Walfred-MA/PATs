#!/usr/bin/env python3

import argparse
import collections as cl
import heapq
import os
import shutil
import subprocess
import sys
import tempfile


def read_fasta(path):
    seqs = {}
    name = None
    chunks = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if name is not None:
                    seqs[name] = ''.join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        seqs[name] = ''.join(chunks)
    return seqs


def run_minimap2(seq_a, seq_b, threads, preset, max_secondary, minimap2, paf_path):
    # seqA contains the reference/unit sequences we want to place onto seqB.
    # Therefore seqB is the minimap2 target and seqA is the minimap2 query.
    cmd = [
        minimap2,
        '-x', preset,
        '-t', str(threads),
        '-c',
        '--eqx',
        '-N', str(max_secondary),
        seq_b,
        seq_a,
    ]
    with open(paf_path, 'w') as out:
        proc = subprocess.run(cmd, stdout=out, stderr=sys.stderr)
    if proc.returncode != 0:
        raise RuntimeError('minimap2 failed with exit code %d' % proc.returncode)


def parse_paf(paf_path, seq_a_seqs, min_identity=0.90, min_unmasked=1000):
    """Return alignments grouped by seqB target contig.

    Each alignment is a dict with seqB coordinates and seqA unit identity.
    Trans_* queries are retained without the ordinary identity/unmasked filters,
    matching the behavior of the original script.
    """
    by_contig = cl.defaultdict(list)

    with open(paf_path) as f:
        for line_no, line in enumerate(f, 1):
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 12:
                raise ValueError('Malformed PAF line %d: expected >=12 fields' % line_no)

            qname = fields[0]          # seqA record/unit
            qlen = int(fields[1])
            qstart = int(fields[2])
            qend = int(fields[3])
            strand = fields[4]
            tname = fields[5]          # seqB contig
            tstart = int(fields[7])
            tend = int(fields[8])
            nmatch = int(fields[9])
            alnlen = int(fields[10])
            mapq = int(fields[11])

            if alnlen <= 0 or tend <= tstart:
                continue

            identity = nmatch / alnlen
            mismatches = max(0, alnlen - nmatch)
            raw_score = 6 * nmatch - 5 * mismatches
            is_trans = qname.startswith('Trans_')

            if not is_trans:
                seq = seq_a_seqs.get(qname)
                if seq is None:
                    raise KeyError('PAF query %r is absent from seqA FASTA' % qname)
                unmasked = sum(1 for b in seq[qstart:qend] if b.isupper())
                if identity < min_identity or unmasked < min_unmasked:
                    continue

            by_contig[tname].append({
                'start': tstart,
                'end': tend,
                'name': qname,
                'strand': strand,
                'identity': identity,
                'nmatch': nmatch,
                'alnlen': alnlen,
                'score': raw_score,
                'mapq': mapq,
                'is_trans': is_trans,
                'qstart': qstart,
                'qend': qend,
                'qlen': qlen,
            })

    return by_contig


def best_supported_segments(alignments, min_segment=100):
    """Partition covered seqB sequence by the highest-scoring active alignment."""
    if not alignments:
        return []

    events = cl.defaultdict(lambda: {'start': [], 'end': []})
    for i, aln in enumerate(alignments):
        events[aln['start']]['start'].append(i)
        events[aln['end']]['end'].append(i)

    positions = sorted(events)
    active = set()
    heap = []  # (-score, -mapq, -alnlen, index)
    segments = []
    prev = None

    for pos in positions:
        if prev is not None and pos > prev and active:
            while heap and heap[0][3] not in active:
                heapq.heappop(heap)
            if heap:
                idx = heap[0][3]
                if pos - prev >= min_segment:
                    a = alignments[idx]
                    segments.append({
                        'start': prev,
                        'end': pos,
                        'name': a['name'],
                        'strand': a['strand'],
                        'score': a['score'],
                        'identity': a['identity'],
                        'mapq': a['mapq'],
                        'source_idx': idx,
                    })

        # End events first: PAF/BED intervals are half-open [start,end).
        for idx in events[pos]['end']:
            active.discard(idx)
        for idx in events[pos]['start']:
            active.add(idx)
            a = alignments[idx]
            heapq.heappush(heap, (-a['score'], -a['mapq'], -a['alnlen'], idx))
        prev = pos

    return segments


def merge_same_units(segments, merge_gap=10000):
    """Merge nearby segments assigned to the same seqA unit and strand."""
    if not segments:
        return []

    segments = sorted(segments, key=lambda x: (x['start'], x['end']))
    merged = []

    for seg in segments:
        if (merged and
                seg['name'] == merged[-1]['name'] and
                seg['strand'] == merged[-1]['strand'] and
                seg['start'] - merged[-1]['end'] <= merge_gap):
            m = merged[-1]
            old_len = m['covered_bp']
            add_len = seg['end'] - seg['start']
            total = old_len + add_len
            m['identity'] = ((m['identity'] * old_len) + (seg['identity'] * add_len)) / total
            m['covered_bp'] = total
            m['end'] = max(m['end'], seg['end'])
            m['score'] += seg['score']
            m['mapq'] = max(m['mapq'], seg['mapq'])
        else:
            x = dict(seg)
            x['covered_bp'] = seg['end'] - seg['start']
            merged.append(x)

    return merged


def mark_trans_overlap(units, trans_alignments):
    """Mark a resolved unit Exon if any Trans_* alignment overlaps it."""
    trans = sorted(trans_alignments, key=lambda x: x['start'])
    j = 0
    active = []
    out = []

    for unit in sorted(units, key=lambda x: x['start']):
        while j < len(trans) and trans[j]['start'] < unit['end']:
            active.append(trans[j])
            j += 1
        active = [x for x in active if x['end'] > unit['start']]
        x = dict(unit)
        x['class'] = 'Exon' if any(t['start'] < unit['end'] and t['end'] > unit['start'] for t in active) else 'Intron'
        out.append(x)
    return out


def bed_score(identity):
    return max(0, min(1000, int(round(identity * 1000))))


def main(args):
    if shutil.which(args.minimap2) is None and not os.path.exists(args.minimap2):
        sys.exit('ERROR: minimap2 was not found. Put it in PATH or provide --minimap2 /path/to/minimap2')

    seq_a_seqs = read_fasta(args.seqA)

    tmp_paf = None
    if args.paf_out:
        paf_path = args.paf_out
    else:
        fd, paf_path = tempfile.mkstemp(prefix='consistent_units.', suffix='.paf')
        os.close(fd)
        tmp_paf = paf_path

    try:
        run_minimap2(
            args.seqA, args.seqB, args.threads, args.preset,
            args.max_secondary, args.minimap2, paf_path
        )

        by_contig = parse_paf(
            paf_path,
            seq_a_seqs,
            min_identity=args.similar,
            min_unmasked=args.min_unmasked,
        )

        output_rows = []
        for contig, all_alns in by_contig.items():
            ordinary = [x for x in all_alns if not x['is_trans']]
            trans = [x for x in all_alns if x['is_trans']]

            segments = best_supported_segments(ordinary, args.min_segment)
            units = merge_same_units(segments, args.merge_gap)
            units = mark_trans_overlap(units, trans)

            for u in units:
                span = u['end'] - u['start']
                min_size = args.min_exon_size if u['class'] == 'Exon' else args.min_intron_size
                if span < min_size:
                    continue

                name = u['name']
                if args.annotate_class:
                    name = '%s|%s' % (name, u['class'])

                output_rows.append((
                    contig,
                    u['start'],
                    u['end'],
                    name,
                    bed_score(u['identity']),
                    u['strand'],
                ))

        output_rows.sort(key=lambda x: (x[0], x[1], x[2], x[3]))
        with open(args.output, 'w') as out:
            for row in output_rows:
                out.write('\t'.join(map(str, row)) + '\n')

    finally:
        if tmp_paf is not None:
            try:
                os.unlink(tmp_paf)
            except FileNotFoundError:
                pass


def run():
    p = argparse.ArgumentParser(
        description=(
            'Align seqA reference/unit FASTA records directly to seqB with minimap2, '
            'resolve overlapping mappings into locally consistent units on seqB, and output BED6.'
        )
    )
    p.add_argument('-a', '--seqA', required=True,
                   help='seqA FASTA: reference/unit sequences to place onto seqB')
    p.add_argument('-b', '--seqB', required=True,
                   help='seqB FASTA: target assembly on which BED regions are reported')
    p.add_argument('-o', '--output', required=True,
                   help='output BED6 file on seqB coordinates')
    p.add_argument('-t', '--threads', type=int, default=8,
                   help='minimap2 threads [8]')
    p.add_argument('-x', '--preset', default='asm10',
                   help='minimap2 preset [asm10]')
    p.add_argument('-N', '--max-secondary', type=int, default=100,
                   help='minimap2 -N value [100]')
    p.add_argument('-s', '--similar', type=float, default=0.90,
                   help='minimum PAF identity for non-Trans_ alignments [0.90]')
    p.add_argument('--min-unmasked', type=int, default=1000,
                   help='minimum uppercase aligned bases in seqA for non-Trans_ records [1000]')
    p.add_argument('--min-segment', type=int, default=100,
                   help='minimum locally winning segment length [100]')
    p.add_argument('--merge-gap', type=int, default=10000,
                   help='merge nearby pieces assigned to the same seqA unit [10000]')
    p.add_argument('--min-exon-size', type=int, default=3000,
                   help='minimum output span for units overlapping Trans_* records [3000]')
    p.add_argument('--min-intron-size', type=int, default=5000,
                   help='minimum output span for other units [5000]')
    p.add_argument('--annotate-class', action='store_true',
                   help='append |Exon or |Intron to BED name')
    p.add_argument('--paf-out', default=None,
                   help='optional path to retain raw minimap2 PAF')
    p.add_argument('--minimap2', default='minimap2',
                   help='minimap2 executable or path [minimap2]')
    args = p.parse_args()
    main(args)


if __name__ == '__main__':
    run()
