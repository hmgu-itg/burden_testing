#!/bin/bash

function usage {
    echo ""
    echo "SMMAT to SAIGE group list conversion"
    echo "Usage:" $(basename $0) "-i <input.file> -o <output.file> { -c <chr> -w }"
    echo " -h : this help message"
    echo " -i : input SMMAT group file"
    echo " -o : output SAIGE group file"
    echo " -c : optional; select records from only one chromosome; default: all"
    echo " -w : optional; output weights; default: false"

    exit 0
}

OPTIND=1

infile=""
outfile=""
chr=""
weights=0

while getopts "i:o:wc:h" optname; do
    case "$optname" in
        "i" ) infile="${OPTARG}";;
        "o" ) outfile="${OPTARG}";;
        "w" ) weights=1;;
        "c" ) chr="${OPTARG}";;
        "h" ) usage ;;
        "?" ) usage ;;
        *) usage ;;
    esac;
done

if [[ $# -eq 0 ]];then
    usage
    exit 0
fi

if [[ -z "$infile" ]];then
    echo "ERROR: input file not specified"
    exit 1
fi

if [[ ! -f "$infile" ]];then
    echo "ERROR: input file does not exist"
    exit 1
fi

cmd="cat"
if [[ ! -z "$chr" ]];then
    cmd="awk -v c=$chr 'BEGIN{FS=\"\t\";OFS=\"\t\";}\$2==c{print \$0;}'"
fi

eval $cmd "$infile" | sort -k2,2 -k3,3n | perl -slne 'BEGIN{$H={};$,=" ";}{@a=split(/\t/);$str=$a[1].":".$a[2]."_".$a[3]."/".$a[4]; if ($w=="1"){$str.=";".$a[5];} push @{$H{$a[0]}},$str;}END{foreach $x (sort keys %H){print $x,join(" ",@{$H{$x}});}}' -- -w="$weights" > "$outfile"
exit 0
