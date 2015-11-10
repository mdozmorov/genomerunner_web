#!/bin/bash

tawk() {
    awk -v OFS='	' "$@"
}

# NOTE: tabix -B slurps the entire region BED file into memory
#       and output is already sorted and uniq'ed.
#       Could be better to compile an alternate tabix binary
#       to incrementally parse a SORTED region bed file, 
#       then sort/uniq externallly.

# Count the number of unique intervals in $1 (the index)
# that overlap with at least one interval in $2 (the track)
overlapCount() {
    case $1 in
        *.bb|*.bigBed) 
            toBED $2 | while IFS=$'\t' read -a fields; do
                bigBedToBed -chrom="${fields[0]}" \
                    -start="${fields[1]}" -end="${fields[2]}" \
                    $1 stdout
            done | sort -u ;;
        # $1 is an index
        *.bed.gz) 
            tabix -B $1 <(toBED $2)
            ;;
        # $1 is not an index, so build one 
        *)
            #bedtools intersect -u -a $1 -b $2 
            tmp=$(mktemp).bed.gz
            toBED $1 | sort -k1,1b -k2,2n | bgzip > $tmp
            tabix $tmp
            tabix -B $tmp <(toBED $2)
            rm $tmp ${tmp}.tbi
        ;;
    esac | wc -l
}
export -f overlapCount

# Count the number of intervals in $1
# FIXME: ignore comment characters?
intervalCount() {
    case "$1" in
        *.bb|*.bigBed) bigBedInfo $1 \
            | awk '$1==itemCount {print $2}' | tr -d ',' ;;
        *) lineCount $1 | cut -f1;;
    esac
}
export -f intervalCount

# If $1 is a BGZF file, ensure that the tabix index is built.
# Otherwise, do nothing.
# In either case, return the path on stdout.
ensureTabix() {
    path=$1

    if [[ ${path} != *.gz ]]; then
        echo $path
        return
    fi

    if [ ! -f ${path}.tbi ]; then
        tabix -p bed ${path} 2> /dev/null || {
            cat <<EOF 1>&2
ERROR: Couldn't locate or create tabix index for ${path}.
* Is tabix on your PATH?
* Is the file a sorted BGZF file?
EOF
            exit 1
        }
    fi
    echo ${path}
}

# FIXME: This may have endian issues
isBGZF() {
    test "$( hexdump -n 4 $1 | head -1 | cut -d' ' -f2- | tr -d ' ' )" == "8b1f0408"
}

basesCovered() {
    case $1 in
        *.bigBed|*.bb)
            bigBedInfo $1 | awk '$1=="basesCovered:" {print $2}' | tr -d ','
            ;;
        *) zcat -f $1 \
                | bedtools sort | bedtools merge \
                | awk 'BEGIN {n=0} {n+=$3-$2} END {print n}'
            ;;
    esac
}
