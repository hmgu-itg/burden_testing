#!/bin/bash

function usage {
    echo ""
    echo "Usage: $0 -o <outdir>"
    echo ""
    exit 0
}

outdir=""
OPTIND=1
while getopts "o:h" optname; do
    case "$optname" in
        "o" ) outdir="${OPTARG}" ;;
        "h" ) usage ;;
        "?" ) usage ;;
        *) usage ;;
    esac;
done

if [[ -z ${outdir} ]];then
    echo "Error: output directory not specified"
    usage
    exit 1
fi

if [[ ! -e ${outdir} ]];then
    echo "Info: creating ${outdir}"
    mkdir -p ${outdir}
    if [[ $? -ne 0 ]];then
	echo "Error: could not create ${outdir}"
	exit 1
    fi    
fi

cd ${outdir}

wget https://krishna.gs.washington.edu/download/CADD/v1.5/GRCh38/whole_genome_SNVs.tsv.gz
wget https://krishna.gs.washington.edu/download/CADD/v1.5/GRCh38/whole_genome_SNVs.tsv.gz.tbi
for i in $(seq 1 22);do
    wget https://xioniti01.u.hpc.mssm.edu/v1.1/Eigen_hg19_noncoding_annot_chr"$i".tab.bgz
    wget https://xioniti01.u.hpc.mssm.edu/v1.1/Eigen_hg19_noncoding_annot_chr"$i".tab.bgz.tbi
done
wget https://xioniti01.u.hpc.mssm.edu/v1.1/Eigen_hg19_coding_annot_04092016.tab.bgz
wget https://xioniti01.u.hpc.mssm.edu/v1.1/Eigen_hg19_coding_annot_04092016.tab.bgz.tbi
