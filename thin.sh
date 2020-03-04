#!/bin/bash
zcat Eigen_hg19_noncoding_annot_chr"$1".tab.bgz | awk '(NR % '$2') == 0' | cut -f1-4,30> "$1".eigen.v11
zcat Eigen_hg19_coding_annot_04092016.tab.bgz | cut -f1-4,18 | awk -v j="$1" '$1==j && (NR % '$2') == 0' >>  "$1".eigen.v11
sort -k2,2n $1.eigen.v11 | sponge $1.eigen.v11
