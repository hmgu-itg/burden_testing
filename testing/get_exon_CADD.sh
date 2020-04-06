#!/bin/bash

gene=$1
gencode=$2
lfeatures=$3

# get gene coordinates from GENCODE
chr="NA"
start="NA"
end="NA"
read -r chr start end <<<$(zcat gencode.basic.annotation.tsv.gz| while read line;do read -r c s e <<<$(echo $line|awk -v g=$gene '{if ($5==g){print $1,$2,$3;}else{print "NA","NA","NA";}}'); if [[ $c != "NA" ]];then echo $c $s $e;break;fi;done)

tabix $lfeatures $chr:$start-$end | grep "\"gene_ID\":\"$gene\""| cut -f 5| perl -MJSON `%h = %{decode_json($line)};if ($h{"source"} eq "GENCODE" && $h{"class"} eq "exon"){$,="\t";print $h{"chr"},$h{"start"}-50,$h{"end"}+50;}`


