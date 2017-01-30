#!/usr/local/bin/bash

## Project: HELIC 15X burden testing
## Job: The aim is to create a script that process the output and generates
### as comprehensive output as possible.

## Version: 1.0
## Date: 2017.01.18

## Analysis steps:
# 1. The script takes the working directory of the input folders of all runs.
# 2. Summarizes the p-values from the results files.
# 3. Cretes Manhattan plots using the mid points of the genes.
# 4. Reads all hits from all phenotypes.
# 5. Summarizes the hits into a table.

# Print out help:
if [[ "${1}" == "-h" ]]; then
    echo ""
    echo "Usage:"
    echo "   $0 <workingDir> <sourceDir>"
    echo ""
    echo "<workingDir> - folder in which the summaries will be created."
    echo "<sourceDir> - folder from which the results of the burden runs will be read."
    echo ""
    echo ""

    exit 1
fi


# Setting up environment:
workingDir="$1"
sourceDir="$2"
#workingDir=/lustre/scratch115/realdata/mdt0/projects/t144_helic_15x/analysis/HA/burdentesting/2017.01.17_noweight_severe_variants
scriptDir=/nfs/team144/ds26/scripts/burden_testing

# If no working and source dir are given the current dir will be used:
if [[ -z "${workingDir}" ]]; then workingDir=$(pwd); fi
if [[ -z "${sourceDir}" ]]; then sourceDir=$(pwd); fi

# Create a summary folder in the workingDir:
mkdir -p ${workingDir}/summaries
mkdir -p ${workingDir}/temp_for_hits

# Get a list of all relevant folders:
echo -n "[Info] Summarizing results files."
find $sourceDir -maxdepth 1 -type d | grep Pheno | while read folder; do
    echo -n "."
    trait=$( echo $folder | perl -lane '$_ =~ /Pheno\.(.+)/; print $1')

    # Summarizing the results files from the given path:
    cat <(echo -e "GeneName\tp-value\tSNP_count") \
        <(find ${folder} -name results | xargs -n1 grep -v "GeneName") \
        | gzip > ${workingDir}/summaries/Pooled_results.Pheno.${trait}.tsv.gz

    # Calling Manhattan plot creating script.
    #. ...
done
echo " done."

# Read hits files and extact relevant information:
echo -n "[Info] Summarizing hits."
echo -e "trait\tgene\tburden_pval\tSNP_count\thighest_b_pval\tsingle_point\tlowest_b_pval\tsingle_point"  \
    > ${workingDir}/summaries/hits_summary.tsv
for file in $( find ${sourceDir} | grep hits | grep gz) ; do
    echo -n "."

    # Get trait and gene: Pheno.RDWPC/hits/DCHS1.tar.gz
    read trait gene <<< $(echo $file | perl -lane '$_ =~ /Pheno\.(.+?)\/hits\/(.+?)\.tar\.gz/; print "$1 $2"')

    # Get burden p-value: # Is not working now, as the the p-values are not the correct ones in thefile...
    read bgene bpval bsnps <<< $(zgrep -w "${gene}" \
        ${workingDir}/summaries/Pooled_results.Pheno.${trait}.tsv )

    # Get highest and lowest p-value from the tests:
    # extract file from achive:
    tar -xf ${file} -C ${workingDir}/temp_for_hits/ ${gene}/${gene}.pvalues.txt

    # Get least significant p-value:
    read HSP HBP <<< $( tail -n+2 ${workingDir}/temp_for_hits/${gene}/${gene}.pvalues.txt \
             | sort -k4,4g | cut -f3,4 | tail -n1 )
    # Get the most significant p-value:
    read LSP LBP <<< $( tail -n+2 ${workingDir}/temp_for_hits/${gene}/${gene}.pvalues.txt \
             | sort -k4,4g | cut -f3,4 | head -n1 )

    # Print report:
    echo -e "${trait}\t${gene}\t${bpval}\t${bsnps}\t${HBP}\t${HSP}\t${LBP}\t${LSP}"  >> ${workingDir}/summaries/hits_summary.tsv

    rm -rf ${workingDir}/temp_for_hits/${gene}
done
echo " done."
