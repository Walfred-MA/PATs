#!/usr/bin/env python3

"""Shared constants and encoders for the PATs/Ctyper matrix format."""

MATRIX_VERSION = "v2.0.1"
SUPPORTED_CTYPER_VERSION = "v1.2.0"
INDEX_VERSION_LINE = f"@{MATRIX_VERSION},support:{SUPPORTED_CTYPER_VERSION}"

PATH_STRAND_WIDTH = 3
QUERY_INDEX_WIDTH = 3
SIZE_WIDTH = 4
QUERY_POSITION_WIDTH = 4
REFERENCE_POSITION_WIDTH = 4

MAX_PATH_INDEX = (1 << 15) - 1


def encode_fixed(value, width, field_name, maximum=None):
    """Encode a non-negative integer as exactly ``width`` base-64 characters."""

    value = int(value)
    capacity = (64**width) - 1
    if maximum is None:
        maximum = capacity
    else:
        maximum = min(int(maximum), capacity)

    if value < 0 or value > maximum:
        raise ValueError(
            f"{field_name}={value} is outside the supported range 0..{maximum} "
            f"for the {width}-character {MATRIX_VERSION} matrix field"
        )

    encoded = ""
    remaining = value
    for _ in range(width):
        encoded = chr(ord("0") + (remaining % 64)) + encoded
        remaining //= 64

    if remaining:
        raise ValueError(f"{field_name}={value} does not fit in {width} characters")
    return encoded


def encode_tag(kmer_flag, kmer_ratio):
    """Encode the existing one-character flag and two-character ratio fields."""

    ratio = min(1000, int(100 * float(kmer_ratio) + 0.5))
    return (
        encode_fixed(max(0, int(kmer_flag)), 1, "k-mer flag")
        + encode_fixed(max(0, ratio), 2, "k-mer ratio")
    )


def encode_path_query(path_index, negative_strand, query_index):
    """Encode a 15-bit path plus strand bit, followed by the query index."""

    path_index = int(path_index)
    if path_index < 0 or path_index > MAX_PATH_INDEX:
        raise ValueError(
            f"path index={path_index} is outside the supported range "
            f"0..{MAX_PATH_INDEX} for the {MATRIX_VERSION} matrix format"
        )

    packed_path = (path_index << 1) | int(bool(negative_strand))
    return encode_fixed(
        packed_path, PATH_STRAND_WIDTH, "packed path+strand"
    ) + encode_fixed(query_index, QUERY_INDEX_WIDTH, "query index")


def encode_coordinates(size, query_position, reference_position):
    """Encode size, query position, and reference/graph position."""

    return (
        encode_fixed(size, SIZE_WIDTH, "size")
        + encode_fixed(query_position, QUERY_POSITION_WIDTH, "query position")
        + encode_fixed(
            reference_position, REFERENCE_POSITION_WIDTH, "reference position"
        )
    )


def encode_metadata(path_index, location, size, query_index, query_position,
                    kmer_flag, kmer_ratio):
    """Return the three fixed-width metadata fields used in a k-mer row."""

    location = int(location)
    return (
        encode_tag(kmer_flag, kmer_ratio),
        encode_path_query(path_index, location < 0, query_index),
        encode_coordinates(size, query_position, abs(location)),
    )
