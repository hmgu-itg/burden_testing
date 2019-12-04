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

cat "$input" | sort -n -k1,1 -n -k2,2 | bgzip > "$output"
tabix -c C -b 2 -e 2 -s 1 "$output"



   
