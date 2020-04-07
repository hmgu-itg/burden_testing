#!/bin/bash

gencode=$1
lfeatures=$2
ilist=$3
caddfile=$4
extend=$5
outfile=$6

if [[ -z ${extend} ]];then
    extend=50
fi
export extend=$extend

zcat $gencode|while read chr start end gname ID;do
#    echo "Tabix linked features"
    tabix $lfeatures $chr:$start-$end | grep "\"gene_ID\":\"$ID\""| cut -f 5| perl -MJSON -lne 'BEGIN{$extend=$ENV{extend};}{%h = %{decode_json($_)};if ($h{source} eq "GENCODE" && $h{class} eq "exon"){$,="\t";print $h{"chr"},$h{"start"}-$extend,$h{"end"}+$extend;}}' | sort -k1,1 -k2,2n > 01.regions.CADD.bed

#    echo "mergeBed"
    mergeBed -i 01.regions.CADD.bed > 02.regions.merged.CADD.bed

    # selecting variants from ilist
    # list is 1-based, BEDs are 0-based
    tabixstr=$(cat 02.regions.merged.CADD.bed| perl -lne '@a=split(/\t/);$s=$a[1]+1;print $a[0].":".$s."-".$a[2];'| tr '\n' ' ')
    tabix $ilist $tabixstr > 03.variants.CADD.txt

    # get CADD scores, omitting indels#
#    echo "CADD scores"
    cat 03.variants.CADD.txt | perl -lne '@a=split(/\t/);$,="\t";if (length($a[3])==1 && length($a[4])==1){print $a[0],$a[1],$a[3],$a[4];}'|while read c p r a;do score=$(tabix $caddfile $c":"$p"-"$p| awk -v a=$a 'BEGIN{FS="\t";}$4==a{print $6;}');echo $ID $c $p "." $r $a $score;done | tr ' ' '\t' >> $outfile
    echo $ID
done










