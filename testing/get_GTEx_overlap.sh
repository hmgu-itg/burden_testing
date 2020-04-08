#!/bin/bash

gencode=$1
lfeatures=$2
ilist=$3
eigenfile=$4
outfile=$5

chainfile="/usr/local/bin/burden_testing/hg38ToHg19.over.chain"

# Equivalent of: -e promoter,enhancer,TF_bind -l promoter,enhancer,TF_bind -s EigenPhred

zcat $gencode|while read chr start end gname ID;do
    start=$((start-1)) # linked features is 0-based
    echo "tabix $lfeatures $chr:$start-$end"
    tabix $lfeatures $chr:$start-$end | grep "\"gene_ID\":\"$ID\""| cut -f 5| perl -MJSON -lne '%h = %{decode_json($_)};if ( ($h{source} eq "GTEx" && ( $h{class} eq "enhancer" || $h{class} eq "promoter" || $h{class} eq "TF_binding_site")) || ($h{source} eq "overlap" && ( $h{class} eq "enhancer" || $h{class} eq "promoter" || $h{class} eq "TF_binding_site")) ){$,="\t";print $h{"chr"},$h{"start"},$h{"end"};}' | sort -k1,1 -k2,2n > 01.regions.GTEx.overlap.bed

    n=$(cat 01.regions.GTEx.overlap.bed | wc -l)
    if [[ $n -eq 0 ]];then
	continue
    fi    
    
    mergeBed -i 01.regions.GTEx.overlap.bed > 02.regions.merged.GTEx.overlap.bed

    n=$(cat 02.regions.merged.GTEx.overlap.bed | wc -l)
    if [[ $n -eq 0 ]];then
	continue
    fi    

    # selecting variants from ilist
    # list is 1-based, BEDs are 0-based
    tabixstr=$(cat 02.regions.merged.GTEx.overlap.bed| perl -lne '@a=split(/\t/);$s=$a[1]+1;print $a[0].":".$s."-".$a[2];'| tr '\n' ' ')
    tabix $ilist $tabixstr > 03.variants.GTEx.overlap.txt

    n=$(cat 03.variants.GTEx.overlap.txt | wc -l)
    if [[ $n -eq 0 ]];then
	continue
    fi    
    
    # remove indels, output 0-based bed
    cat 03.variants.GTEx.overlap.txt| perl -lne '@a=split(/\t/);$,="\t";if (length($a[3])==1 && length($a[4])==1){print "chr".$a[0],$a[1]-1,$a[1],$a[0]."_".$a[1]."_".$a[3]."_".$a[4];}' > 04.variants.GTEx.overlap.noindels.bed

    n=$(cat 04.variants.GTEx.overlap.noindels.bed | wc -l)
    if [[ $n -eq 0 ]];then
	continue
    fi    
    
    # liftOver
    liftOver 04.variants.GTEx.overlap.noindels.bed $chainfile 05.liftover.out.GTEx.overlap.bed 05.unmapped.GTEx.overlap.bed
    
    n=$(cat 05.liftover.out.GTEx.overlap.bed | wc -l)
    if [[ $n -eq 0 ]];then
	continue
    fi    

    # getEigenPhred scores

    cat 05.liftover.out.GTEx.overlap.bed | while read c s e d;do c=${c/#chr};read -r x z r a <<<$(echo $d| tr '_' ' '); echo "tabix $eigenfile $c:$e-$e";lines=$(tabix $eigenfile $c":"$e"-"$e);score=$(echo "$lines" | awk -v r=$r -v a=$a 'BEGIN{FS="\t";b="NA";}($3==r && $4==a) || ($3==a && $4==r){b=$7;}END{print b;}');if [[ -z $score ]];then score="NA";fi;echo $ID $x $z "." $r $a $score | tr ' ' '\t' >> $outfile;done 
    echo $ID
    
done










