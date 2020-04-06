#!/bin/bash

#gene=$1
gencode=$1
lfeatures=$2
ilist=$3
caddfile=$4
extend=$5

if [[ -z ${extend} ]];then
    extend=50
fi
export extend=$extend

zcat $gencode|while read chr start end gname ID;do
    tabix $lfeatures $chr:$start-$end | grep "\"gene_ID\":\"$ID\""| cut -f 5| perl -MJSON -lne 'BEGIN{$extend=$ENV{extend};}{%h = %{decode_json($_)};if ($h{source} eq "GENCODE" && $h{class} eq "exon"){$,="\t";print $h{"chr"},$h{"start"}-$extend,$h{"end"}+$extend;}}' | sort -k1,1 -k2,2n > 01.regions.bed

    mergeBed -i 01.regions.bed > 02.regions.merged.bed

    # selecting variants from ilist
    tabixstr=$(cat 02.regions.merged.bed| perl -lne '@a=split(/\t/);$s=$a[1]+1;print $a[0].":".$s."-".$a[2];'| tr '\n' ' ')
    tabix $ilist $tabixstr > 03.variants.txt

    # get CADD scores, omitting indels
    cat 03.variants.txt | perl -lne '@a=split(/\t/);$,="\t";if (length($a[3])==1 && length($a[4])==1){print $a[0],$a[1],$a[3],$a[4];}'|while read c p r a;do score=$(tabix $caddfile $c":"$p"-"$p| awk -v a=$a 'BEGIN{FS="\t";}$4==a{print $6;}');echo $ID $c $p "." $r $a $score;done | tr ' ' '\t' >> 04.output.list.txt
    echo $ID
done










