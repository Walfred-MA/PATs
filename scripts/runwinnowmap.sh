#!/bin/bash

query=$1
database=$2
threads="${3:-4}"
simi="${4:-90}"
size="${5:-150}"
mode="${6:-nonmasking}" # 'nonmask' (default) or 'mask'

# Set base options common to both modes
base_opts=" -t $threads -cx asm5 --secondary=yes "

# Decide on profile
if [[ "$mode" == "masking" ]]; then

    # Parameters for non-masked sensitive profile
    repfile="${database}_repetitive_k13_H0.1.txt"
    [[ ! -f "$repfile" ]] && {
        winnowmap index -k 13 -H 0.1 "$database" -o "$repfile"
    }

    # Additional profile-specific options
    extra_opts="-k 15 -w 5 -m $size -p 0.000 --no-end-flt -W $repfile"

else

    # Parameters for masked profile
    extra_opts="-k 15 -w 5 -m 300 -p 0.01 "
fi

winnowmap $extra_opts $base_opts "$database" "$query" 2>/dev/null |

awk -v cutoff=$size -v simicutoff=$simi '
BEGIN {
    OFS="\t"
} 
$1 ~ /^@/ {print; next} 
{
    
    # Determine strand: if leading H, assume reverse
    qname=$1
    qsize=$2
    qstart=$3
    qend=$4
    strand=($5 == "+") ? 1 : -1
    rname=$6
    rsize=$7
    rstart=$8
    rend=$9
    qual=$12
    nm=$13
    as=$15
    
    for (i = 13; i <= NF; i++) {
        if ($i ~ /^cg:Z:/) {
            split($i, cigartag, ":")
            cigar = cigartag[3]
            $i = ""
        }
        if ($i ~ /^NM:i:/) {
            split($i, a, ":"); nm = a[3]
            nmindex = i
        }
        if ($i ~ /^AS:i:/) {
            asindex = i
        }
        if ($i ~ /^ms:i:/) {
            split($i, b, ":"); len_aln = b[3]
        }
        if ($i ~ /^de:f:/) {
            split($i, b, ":"); df = b[3]
        }
    }

    nm = int(df * len_aln + 0.5)
    
    gsub(/I/, "_", cigar)
    gsub(/D/, "I", cigar)
    gsub(/_/, "D", cigar)
    
    oldcigar = cigar
    len_aln = 0
    len_nm = 0
        cindex=0 
        while (match(cigar, /[0-9]+[M=XIDH]/)) {
            n = substr(cigar, RSTART, RLENGTH - 1)  # Number
            op = substr(cigar, RSTART + length(n), 1)  # Operation
           
            # Update alignment length
            if (op ~ /[M=]/) len_aln += n 
            if (op ~ /[X]/) len_nm += n 
            if (op ~ /[I]/)  len_nm += 5 
            if (op ~ /[D]/) len_nm += 5

            if (strand == -1){
                cindex++
                tok[cindex] = n op
            }
                # Remove parsed part from cigar
            cigar = substr(cigar, RSTART + RLENGTH)
        }
       
    if (strand == -1){ 
        cigar = ""
        for (i = cindex; i >= 1; i--) {
            cigar = cigar tok[i]
        }
    } else {cigar = oldcigar}
    
    hr_1 = rstart
    hr_2 = rsize - rend  # safer calculation
    
    hr_1_s = (hr_1 > 0) ? hr_1 "H" : ""
    hr_2_s = (hr_2 > 0) ? hr_2 "H" : ""
    
    if (strand == -1) {
        cigar = hr_2_s cigar hr_1_s
    } else {
        cigar = hr_1_s cigar hr_2_s
    }
    
    $1=rname
    $2= (strand == 1) ? 0 : 16
    $3=qname
    $4=qstart + 1
    $5 = 255
    $6 = cigar
    $7="*"
    $8="0"
    $9="0"
    $10="*"
    $11="*"
    idscore = 100*(1 - len_nm / ( len_aln + len_nm ) )
    $12 = sprintf("PI:f:%.2f", idscore) # Set the computed PI
    # Parse NM:i: (mismatches/indels)
    $asindex ="AS:i:" len_aln 
    $nmindex ="NM:i:" len_nm


    if (len_aln > cutoff ) {
        # Calculate identity score (based on mismatches/indels and alignment length)
        # If identity score passes threshold
        if (idscore > simicutoff) {
            print len_aln, $0 
            # Print the updated SAM line
        }
    }
}' | sort -k1,1nr | cut -f2- 
