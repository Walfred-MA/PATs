#!/bin/bash
# Usage: ./add_abs_coords_inplace.sh file1.txt file2.txt
# Updates file1.txt in place, appending contig:abs_start-abs_end<strd>

file1="$1"
file2="$2"
tmpfile="$(mktemp)"

awk -v ref="$file2" '
BEGIN {
    # Read file2 into memory
    while ((getline < ref) > 0) lines[++n] = $0
    close(ref)
}
{
    # Split by underscores: e.g. foo_2_100_200
    split($0, parts, "_")
    idx = parts[2]
    start_rel = parts[3]
    end_rel   = parts[4]

    # If invalid index, skip
    if (!(idx ~ /^[0-9]+$/) || !(idx in lines)) {
        print $0
        next
    }

    # Get reference line and parse
    split(lines[idx], a, " ")
    contig = a[1]
    strd   = a[2]
    start_ref = a[3]
    end_ref   = a[4]

    # Compute absolute coords
    if (strd == "+") {
        abs_start = start_ref + start_rel
        abs_end   = start_ref + end_rel
    } else if (strd == "-") {
        abs_start = end_ref - end_rel
        abs_end   = end_ref - start_rel
    } else {
        abs_start = abs_end = "NA"
    }

    # Append formatted info
    printf "%s\t%s:%d-%d%s\n", $0, contig, abs_start, abs_end, strd
}' "$file1" > "$tmpfile" && mv "$tmpfile" "$file1"

