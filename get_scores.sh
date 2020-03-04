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

#for i in $(seq 1 22);do
#    echo "Info: downloading noncoding Eigen scores, chromosome $i"
#    wget https://xioniti01.u.hpc.mssm.edu/v1.1/Eigen_hg19_noncoding_annot_chr"$i".tab.bgz
#    wget https://xioniti01.u.hpc.mssm.edu/v1.1/Eigen_hg19_noncoding_annot_chr"$i".tab.bgz.tbi
#done

echo "Info: downloading coding Eigen scores"
wget https://xioniti01.u.hpc.mssm.edu/v1.1/Eigen_hg19_coding_annot_04092016.tab.bgz
wget https://xioniti01.u.hpc.mssm.edu/v1.1/Eigen_hg19_coding_annot_04092016.tab.bgz.tbi

for i in $(seq 1 22);do
    echo "Info: thinning chr $i"
    thin $i 5000
done

echo "Info: creating all.thin"
cat *.eigen.v11 | sort -k1,1n -k2,2n > all.thin

echo "Info: splitting coding Eigen scores"

for i in $(seq 1 22);do
    tabix Eigen_hg19_coding_annot_04092016.tab.bgz $i > $i.coding
done

for i in $(seq 1 22);do
    echo "Info: creating EigenPhred for coding chrom $i"
    cat $i.coding | cut -f1-4,18 | Rscript --vanilla -e 'a=read.table(file("stdin"),header=F,sep="\t");b=read.table("all.thin",header=F,sep="\t");Eigen=ecdf(b$V5);scaleEigen = function(x) (-10*log(1-Eigen(x), base=10)); a$V6=scaleEigen(a$V5);colnames(a)=NULL;write.table(a, file="'$i'.coding.phred", row.names=F, col.names=F, quote=F);'
done


#for i in $(seq 21 21);do
#    cat eigen.noncoding.tmp.$i >> temp
#    rm eigen.noncoding.tmp.$i
#done

#sort -t ' ' -k1,1n -k2,2n temp | sed 's/ /\t/g' > eigen.noncoding.tab
#bgzip eigen.noncoding.tab
#tabix -b 2 -e 2 -s 1 eigen.noncoding.tab.gz
#rm temp

#wget https://krishna.gs.washington.edu/download/CADD/v1.5/GRCh38/whole_genome_SNVs.tsv.gz
#wget https://krishna.gs.washington.edu/download/CADD/v1.5/GRCh38/whole_genome_SNVs.tsv.gz.tbi
