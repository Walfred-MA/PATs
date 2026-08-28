# KmerIndex

`kmerindex` builds a binary location index for canonical 31-mers in an
uncompressed FASTA assembly. It reads contig lengths and random-access geometry
from the FASTA `.fai`, processes contigs in parallel, sorts all 12-byte
`(uint64 kmer, uint32 position)` occurrences, and then writes the filtered
index.

The 2-bit encoding matches KmerCounter: `A=0`, `C=1`, `G=2`, `T=3`. The stored
key is the larger of the forward 31-mer and its reverse complement. Any base
other than A/C/G/T breaks the current k-mer.

## Build and run

```sh
make
./kmerindex --fasta assembly.fa --output assembly.k31.idx --threads 8
```

The default FASTA index is `assembly.fa.fai`. Create it, if necessary, with
`samtools faidx assembly.fa`, or pass another path with `--fai`.

Options:

- `--maxcover N`: omit a 31-mer when its occurrence count is greater than N.
  Default 255; allowed range 1 through 65535.
- `--maxspan N`: omit a 31-mer when the union of its +/-10,000-base occurrence
  intervals is greater than N. Default 1,000,000.
- `--threads N`: number of contigs processed concurrently. Default is the
  machine's reported CPU count.

Intervals use half-open coordinates and are clipped to their contig. Intervals
from different contigs never merge.

## Coordinates and binary format

Locations are zero-based starts in a cumulative assembly coordinate system.
Contig zero starts at 0; each later contig starts after the preceding contig's
`.fai` length. The original `.fai` therefore maps every stored location back to
a contig and local coordinate.

All integers are little-endian and records are ordered by ascending k-mer key,
then ascending location:

- One location: `uint64 kmer`, `uint32 location` (12 bytes).
- Multiple locations: `uint64 kmer`, `uint32 (UINT32_MAX - count)`, followed by
  `count` `uint32` locations.

There is no file header; parsing ends at EOF. To keep count markers distinct
from coordinates, the program rejects an assembly whose possible 31-mer start
coordinates overlap the marker range selected by `--maxcover`.

Run `make test` for the integration tests.
