#!/bin/bash

# TAB SEPARATED

OPTIND=1
while getopts ":hi:o:" optname; do
    case "$optname" in
        "i") input=${OPTARG} ;;
        "o") output=${OPTARG} ;;
        "h") usage ;;
        ":") usage ;;
    esac
done

if [[ -z "$input" ]];then
    echo "ERROR: no input file specified"
    exit 1
fi

if [[ ! -e "$input" ]];then
    echo "ERROR: input file ($input) does not exist"
    exit 1
fi

if [[ -z "$output" ]];then
    echo "ERROR: no output file specified"
    exit 1
fi

tmpfile=$(mktemp /tmp/compress_table.XXXXXX)

echo "##fileformat=VCFv4.1" > "$tmpfile"
echo "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" >> "$tmpfile"
cat "$input" | sort -k1,1n -k2,2n | awk 'BEGIN{FS="\t";OFS="\t";}{print $1,$2,".",$3,$4,".",".",".";}' >> "$tmpfile"
cat "$tmpfile" | bgzip > "$output"
tabix -p vcf "$output"

rm "$tmpfile"



   
