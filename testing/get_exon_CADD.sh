#!/bin/bash

gene=$1
gencode=$2
lfeatures=$3
ilist=$4
caddfile=$5
extend=$6

if [[ -z ${extend} ]];then
    extend=50
fi

# get gene coordinates from GENCODE
chr="NA"
start="NA"
end="NA"
read -r chr start end <<<$(zcat ${gencode} | while read line;do read -r c s e <<<$(echo $line|awk -v g=$gene '{if ($5==g){print $1,$2,$3;}else{print "NA","NA","NA";}}'); if [[ $c != "NA" ]];then echo $c $s $e;break;fi;done)

# select lines from lfeatures that correspond to the gene, select exon records only; then merge the BED records
export extend=$extend
tabix $lfeatures $chr:$start-$end | grep "\"gene_ID\":\"$gene\""| cut -f 5| perl -MJSON -lne 'BEGIN{$extend=$ENV{extend};}{%h = %{decode_json($_)};if ($h{source} eq "GENCODE" && $h{class} eq "exon"){$,="\t";print $h{"chr"},$h{"start"}-$extend,$h{"end"}+$extend;}}' | sort -k1,1 -k2,2n > 01.regions.bed

mergeBed -i 01.regions.bed > 02.regions.merged.bed

# selecting variants from ilist
tabixstr=""
cat 02.regions.merged.bed| perl -lne '@a=split(/\t/);$s=$a[1]+1;print $a[0].":".$s."-".$a[2];'| while read line;do tabixstr=$tabixstr" $line";done

echo "TABIXSTR: $tabixstr"

tabix $ilist $tabixstr > 03.variants.txt

# get CADD scores
cut -f 1-5 03.variants.txt | perl -lne '@a=split(/\t/);$,="\t";if (length($a[3])==1 && length($a[4])==1){print $a[0],$a[1],$a[3],$a[4];}'|while read c p r a;do score=$(tabix $caddfile $c":"$p"-"$p| awk -v a=$a 'BEGIN{FS="\t";}$4==a{print $6;}');echo $gene $c $p "." $r $a $score;done | tr ' ' '\t' > 04.output.list.txt








