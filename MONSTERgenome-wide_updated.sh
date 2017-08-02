#!/usr/local/bin/bash

# A wrapper script to automate genome-wide burden testing using MONSTER.
# For more information on the applied method see: http://www.stat.uchicago.edu/~mcpeek/software/MONSTER/

# In this version when the script founds a significant gene, it tests if the association is
## driven by a single variant or more by repeating the test with removing a variant
## each time.

# For more information on how to use the wrapper, please see the README.md file or check the repository:
# https://github.com/wtsi-team144/burden_testing

# Todo for the next generation of the wrapper:
## An older verision of MONSTER had an issue why we had to wrap genes individually for 

version="v11 Last modified: 2017.08.01" # This version mainly changes in the documentation.

today=$(date "+%Y.%m.%d") # Get the date

# The variant selector script, that generates snp and genotype input for MONSTER:
regionSelector=burden_get_regions.pl

# Folder with all the scripts:
export scriptDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Folder with all the phenotypes: # No longer hardwired. Accept as a command line parameter.
# phenotypeDir=/lustre/scratch115/projects/t144_helic_15x/analysis/HA/phenotypes
# phenotypeDir=/lustre/scratch115/projects/t144_helic_15x/analysis/HA/phenotypes/correct_names.andmissing

# Default f4ile with all the gene names (only on autosomes):
geneListFile="${scriptDir}/gene_list.lst" 

# Kinship matrix file: # No longer hardwired, accepted as command line parameter. 
# kinshipMatrix=/nfs/team144/ds26/burden_testing/kinship/2016.10.20_fix_diagonal/kinship.fixdiag.txt

# Single point p-values are read from here: # Single point values are no longer needed.
# singlePointDir=/lustre/scratch115/projects/t144_helic_15x/analysis/HA/single_point/output
# singlePointDir=/lustre/scratch115/projects/t144_helic_15x/analysis/HA/single_point/output/missing_chunks

# Path to MONSTER executable: # No longer hardwired, expected to be in the path.<- not sure about this though...
# MONSTER=/nfs/team144/software/MONSTER_v1.3/MONSTER
MONSTER=$(which MONSTER)
missing_cutoff=1 # Above the set missingness threshold, variants will be excluded. Below the missing genotype will be imputed.
imputation_method='' # The default imputation method is BLUP, slowest, but the most accurate. For other options, see documentation. 

# Testing if MONSTER is in the path, exiting if not:
if [[ -z "${MONSTER}" ]]; then echo "[Error] MONSTER is not in the path. Exiting."; exit; fi

##
## Defining default values:
##
export signifThreshold=1e-5 # By default the hit threshold is 1e-5
export keep_temp="No" # By default we are not keeping temporary files:
export chunkCount=1 # By default we process all genes at one chunk.
export chunkNo=1 # By default we are processing the first chunk.
if [[ $LSB_JOBINDEX > 0 ]]; then chunkNo=$LSB_JOBINDEX; fi # If jobindex is available, we set that as chunk number, it might be overridden by -c

# --- print out help message and exit:
function display_help() {
    echo "$1"
    echo ""
    echo "Genome-wide Monster wrapper!"
    echo "version: ${version}"
    echo ""
    echo "This script was written to run MONSTER genome wide. This script takes a series of arguments
based on which it calls downstream helper scripts, and generates specific directory for the output files.
This is supposed to be called by bsub, where the chunk number is read from LSB_JOBINDEX.
It pools results together within one chunk."
    echo ""
    echo "Usage: $0 <parameters>"
    echo ""
    echo "Variant filter options:"
    echo "     -g  - list of gencode features."
    echo "     -e  - list of linked GTEx featuress"
    echo "     -l  - list of linked overlapping features."
    echo "     -m  - upper maf thresholds"
    echo "     -x  - extend genomic regions (bp)."
    echo "     -o  - exclude all non loss-of-function variants from the test (less than missense in severity)."
    echo "     -f  - include only HC and LC loftee variants."
    echo ""
    echo "Parameters to set up scores for variants:"
    echo "     -s  - turn weights on. Arguments: CADD, Eigen, EigenPC, EigenPhred, EigenPCPhred,"
    echo "           Linsight, Mixed"
    echo "     -t  - the value with which the scores will be shifted"
    echo "     -k  - below the specified cutoff value, the variants will be excluded"
    echo ""
    echo "Gene list and chunking:"
    echo "     -L  - list file with the gene names (default: ${geneListFile})."
    echo "     -d  - chunk count: how many chunks the gene list should be split into (default: 1)."
    echo "     -c  - chunk number or jobindex (default: 1)."
    echo ""
    echo "General options:"
    echo "     -w  - working directory (default: current working directory)"
    echo "     -b  - Keep temporary files."
    echo "     -i  - p-value threshold (default: 1e-5) below which the temporary files will be saved."
    echo ""
    echo "Monster parameters:"
    echo "     -p  - phenotype (required, no default)"
    echo "     -P  - phenotype file (required, no default)"
    echo "     -K  - Kinship file (required, no default)"
    echo "     -V  - VCF file (required, no default, Use * character for chromosome name eg 'chr*.vcf.gz')"
    echo ""
    echo "Other options:"
    echo "     -h  - print help message and exit"
    echo ""
    echo "More information: ds26@sanger.ac.uk"
    echo ""
    echo ""
    
    exit 1
}

# --- If a gene fails at some point, generate a report, and save files in a failed folder.
function failed (){
    message="${1}"

    echo "[Error] Gene has failed: ${gene}"
    echo "[Error] Gene has failed: ${gene}" >&2 # Reporting to the error log as well.
    echo "[Error] ${message}"
    echo "[Error] folder moved to ${folder}/failed"
    echo ""

    # Adding NA-s to the output file:
    echo -e "${gene}\tNA\tNA" >> ${workingDir}/gene_set.${chunkNo}/results

    # Generating failed directory and compress gene dir:
    mkdir -p ${workingDir}/failed
    tar -czvf ${workingDir}/failed/${gene}.tar.gz -C ${workingDir}/gene_set.${chunkNo}/ ${gene}

    # Once the folder is compressed, we add to the pooled tar archive:
    if [[ ! -e ${workingDir}/failed/failed_genes.tar ]]; then
        tar -c --remove-files -f ${workingDir}/failed/failed_genes.tar -C ${workingDir}/failed  ${gene}.tar.gz # Create archive
    else
        tar -r --remove-files -f ${workingDir}/failed/failed_genes.tar -C ${workingDir}/failed  ${gene}.tar.gz  # Add to archive
        rm ${workingDir}/failed/${gene}.tar.gz
    fi

    # Removing directory:
    cd ${workingDir}
    rm -rf ${workingDir}/gene_set.${chunkNo}/${gene}
    return 0
}

# --- If a p-value is really low, the whole run backed up:
function savehit (){
    mkdir -p ${workingDir}/hits
    tar -czvf ${workingDir}/hits/${gene}.tar.gz -C ${workingDir}/gene_set.${chunkNo}/ ${gene}
    cd ${workingDir}
    return 0
}

# --- Capture command line options --------------------------------------------

# Help message is printed if no command line arguments has been passed:
if [ $# == 0 ]; then display_help; fi

# Looping through all command line options:
OPTIND=1
while getopts ":hL:c:d:p:P:K:V:bi:g:m:s:l:e:x:k:t:ofw:" optname; do
    case "$optname" in
      # Gene list related parameters:
        "L") geneListFile=${OPTARG} ;;
        "c") export chunkNo=${OPTARG} ;;
        "d") chunkCount=${OPTARG} ;;

      # MONSTER input files:
        "p" ) phenotype=${OPTARG} ;;
        "P" ) phenotypeFile=${OPTARG} ;;
        "K" ) kinshipFile=${OPTARG} ;;
        "V" ) vcfFile="${OPTARG}" ;;

      # Wrapper option:
        "b") keep_temp="yes";;
        "i") signifThreshold=${OPTARG};;

      # variant filter parameters:
        "g") gencode=${OPTARG} ;;
        "m") MAF=${OPTARG} ;;
        "s") score=${OPTARG} ;;
        "l") overlap=${OPTARG} ;;
        "e") gtex=${OPTARG} ;;
        "x") xtend=${OPTARG} ;;
        "k") cutoff=${OPTARG} ;;
        "t") scoreshift=${OPTARG} ;;
        "o") lof=1 ;;
        "f") loftee=1 ;;

      # Other parameters:
        "w") rootDir=${OPTARG} ;;
        "h") display_help ;;
        "?") display_help "[Error] Unknown option $OPTARG" ;;
        ":") display_help "[Error] No argument value for option $OPTARG";;
        *) display_help "[Error] Unknown error while processing options";;
    esac
done

#--- checking input files - if any of the tests fails, the script exits.---------

# Phenotype file (Is it set? Does it exists?):
if [[ -z "${phenotypeFile}" ]]; then
    display_help "[Error] Phenotype file has to be specified!";
elif [[ ! -e "${phenotypeFile}" ]]; then
    echo "[Error] Phenotype file could not be opened: $phenotypeFile";
    exit;
fi

# Kinship file (Is it set? Does it exists?):
if [[ -z "${kinshipFile}" ]]; then
    display_help "[Error] Kinship file has to be specified!";
elif [[ ! -e "${kinshipFile}" ]]; then
    echo "[Error] Kiship file could not be opened: $kinshipFile";
    exit;
fi

# vcf file (Is it set? Does it exists?):
if [[ -z "${vcfFile}" ]]; then
    display_help "[Error] VCF file has to be specified!";
elif [[ ! -e $( echo "${vcfFile}" | sed -e 's/\%/12/') ]]; then
    echo "[Error] VCF file could not be opened: $vcfFile";
    exit;
else
    commandOptions=" -vcfFile ${vcfFile} "
fi

# Gene list file (Is it set? Does it exists?):
if [[ ! -e "${geneListFile}" ]]; then
    echo "[Error] Gene list file could not be opened: $geneListFile";
    exit;
fi

# Working directory (Is it set? Does it exists?):
if [[ -z "${rootDir}" ]]; then
    rootDir=$(pwd)
elif [[ ! -d "${rootDir}" ]]; then
    echo "[Error] The directory does not exists: $rootDir";
    exit;
fi

# GENCODE -expecting a list of feature names separated by a comma.
if [[ ! -z "${gencode}" ]]; then
    folder="${gencode}"
    commandOptions="${commandOptions} -G ${gencode}"
fi

# GTEx - expecting a list of feature names separeted by comma.
if [[ ! -z "$gtex" ]]; then
    folder="${folder}.GTEX.${gtex}"
    commandOptions="${commandOptions} -E ${gtex}"
fi

# Overlap - expecting a list of features separeated by comma.
if [[ ! -z "${overlap}" ]]; then
    folder="${folder}.Overlap.${overlap}"
    commandOptions="${commandOptions} -L ${overlap}"
fi

# MAF - expecting a float between 0 and 0.5
if [[ ! -z "$MAF" ]]; then
    folder="${folder}.MAF.${MAF}"
    commandOptions="${commandOptions} --maf ${MAF}"
fi

# IF loftee is set, only those variants are included in the test that were replrted to be LC or HC lof variants by loftee:
if [[ ! -z "$loftee" ]]; then
    folder="${folder}.loftee"
    commandOptions="${commandOptions} --loftee "
fi

# If lof is set, only variants with severe consequences will be selected.
if [[ ! -z "$lof" ]]; then
    folder="${folder}.lof"
    commandOptions="${commandOptions} --lof "
fi

# Score - If score is not given we apply no score. Otherwise we test the submitted value:
# Accepted scores:
if [[ ! -z "${score}" ]]; then
    score="${score^^}"
    case "${score}" in
        EIGEN )        score="Eigen";;
        EIGENPC )      score="EigenPC";;
        EIGENPHRED )   score="EigenPhred";;
        EIGENPCPHRED ) score="EigenPCPhred";;
        CADD )         score="CADD";;
        LINSIGHT )     score="Linsight";;
        MIXED )        score="Mixed";;
        * )            score="noweight";;
    esac
else
    echo "[Warning] Submitted score is not recognized! Accepted scores: CADD, Eigen, EigenPC, EigenPhred, EigenPCPhred, Linsight or Mixed."
    echo "[Warning] No scores are being applied."
    score="noweight"
fi

# Only adding score to command line if score is requested:
if [[ "${score}" != "noweight" ]]; then
    commandOptions="${commandOptions} -s ${score}";
fi

# If Eigen score is applied, we shift the scores by 1, if no other value is specified:
if [[ ("${score}" == "Eigen") && (-z "${scoreshift}" ) ]]; then scoreshift=1; fi

# Adjust dir name, noweight if no score is specified:
folder="${folder}.${score}"

# Exons might be etended with a given number of residues:
if [[ ! -z "${xtend}" ]]; then
    folder="${folder}.Xtend.${xtend}"
    commandOptions="${commandOptions} --extend ${xtend}"
fi

# Setting score cutoff, below which the variant will be removed from the test:
if [[ ! -z "${cutoff}" ]]; then
    folder="${folder}.cutoff.${cutoff}"
    commandOptions="${commandOptions} --cutoff ${cutoff}"
fi

# Setting score shift, a number that will be added to every scores (MONSTER does not accept scores <= 0!!):
if [[ ! -z "${scoreshift}" ]]; then
    folder="${folder}.shift.${scoreshift}"
    commandOptions="${commandOptions} --shift ${scoreshift}"
fi

# Checking if phenotype is provided and if that phenotype file:
if [[ -z "${phenotype}" ]]; then
    echo "[Error] Phenotype was not set! Exiting."
    exit 1
fi

# Based on the size of the list and the chunk count we determine the size of each chunk:
chunkSize=$(wc -l "${geneListFile}" | perl -lane 'print int($F[0]/$ENV{"chunkCount"}+0.5)')

# Once everything is checked, print out the parameters, and ready to go:


# Updating working dir, and create folder:
folder=$( echo $folder | perl -lane '$_ =~ s/^\.//; print $_')
workingDir="${rootDir}/${folder}/Pheno.${phenotype}"
mkdir -p "${rootDir}"

# --- Reporting parameters ------------------------------------------------------
echo "##"
echo "## Genome-wide Monster wrapper version ${version}"
echo "## Date: ${today}"
echo "##"
echo ""

echo "[Info] General options:"
echo -e "\tVariant selector: ${regionSelector}"
echo -e "\tScript dir: ${scriptDir}"
echo -e "\tSignificance threshold: ${signifThreshold}"
echo -e "\tWorking directory: ${workingDir}/gene_set.${chunkNo}"
echo -e "\tKeeping temporary files: ${keep_temp}"
echo ""

echo "[Info] Gene list options:"
echo -e "\tGene list file: ${geneListFile}"
echo -e "\tNumber of chunks the gene list is split into: ${chunkCount}"
echo -e "\tCurrent chunk: ${chunkNo}"
echo -e "\tNumber of genes in one chunk: ${chunkSize}"
echo ""

echo "[Info] Variant filtering options:"
echo -e "\tvcf file: ${vcfFile}"
echo -e "\tGENCODE feaures: ${gencode:--}"
echo -e "\tGTEx feaures: ${gtex:--}"
echo -e "\tOverlapping reg.feaures: ${overlap:-NA}"
echo -e "\tFeatures are extended by ${xtend:-0}bp"
echo -e "\tUpper minor allele frequency: ${MAF:-1}"
echo ""

echo "[Info] Weighting options:"
echo -e "\tWeighting: ${score}"
echo -e "\tScore cutoff: ${cutoff:-0}"
echo -e "\tScores shifted by: ${scoreshift:-0}"
echo -e "\tOutput folder: ${workingDir}"
echo ""

echo "[Info] command line options for burden get region: ${commandOptions}"
echo ""

echo "[Info] MONSTER options:"
echo -e "\tMONSTER executable: ${MONSTER}"
echo -e "\tMissingness: ${missing_cutoff}"
echo -e "\tImputation method: ${imputation_method:-BLUP}"
echo -e "\tKinship matrix: ${kinshipFile}"
echo -e "\tPhenotype file: ${phenotypeFile}"
echo -e "\tPhenotype: ${phenotype}"
echo ""

# --- Main loop executed for all genes --------------------------------------------

# Creating results file with header:
mkdir -p ${workingDir}/gene_set.${chunkNo}
if [[ ! -e ${workingDir}/gene_set.${chunkNo}/results ]]; then 
    echo -e "GeneName\tp-value\tSNP_count" >> ${workingDir}/gene_set.${chunkNo}/results
fi

# Looping through all genes in the gene set:
awk -v cn="${chunkNo}" -v cs="${chunkSize}" 'NR > (cn-1)*cs && NR <= cn*cs' ${geneListFile} | while read gene ; do

    # Removing colons if genomic coordinates are submitted:
    gene=${gene/:/_}

    echo ""
    echo "[Info] Processing $gene"
    (>&2 echo "[Info] Processing $gene" ) # Adding mark into the error file, so we'll know in which gene the error comes up.
    echo "[Info] Step 1.: extracting variants given the list of parameters."

    # Creating folders for the gene:
    mkdir -p ${workingDir}/gene_set.${chunkNo}/${gene}
    cd ${workingDir}/gene_set.${chunkNo}/${gene}

    # Reporting call:
    echo "${scriptDir}/${regionSelector}  -i $gene -o ${gene}_output ${commandOptions} -v > ${gene}_log"
    ${scriptDir}/${regionSelector}  -i $gene -o ${gene}_output ${commandOptions} -v > ${gene}_log

    # We have to test if the output files are OK, then we go to the next gene:
    if [[ ! -e ${gene}_output_genotype ]]; then
        failed "${gene} has failed, no genotype output generated. Gene skipped."
        continue
    fi
    genoLines=($(wc -l ${gene}_output_genotype))

    # Checking if output files are created and has enough lines to run MONSTER:
    if [[ ! ${genoLines[0]} -ge 0 ]]; then
        failed "${gene} has failed, genotype file is empty. Gene skipped."
        continue
    fi
    if [[ ! -e ${gene}_output_variants ]]; then
        failed "${gene} has failed, phenotype file was not generated. Gene skipped."
        continue
    fi

    # The variant IDs have to be adjusted: the non-alphanumeric characters are also removed:
    cat ${gene}_output_variants | perl -lane '$_ =~ s/[^0-9a-z\t\.]//gi; print $_'  > snpfile.mod.txt
    
    # The sample IDs also have to be mapped to numbers, so a mapping file has to be created when the first gene of the chunk is created:
    if [[ ! -e ${workingDir}/gene_set.${chunkNo}/sample_mapping.sed ]]; then
        cat <(head -n1 ${gene}_output_genotype | tr "\t" "\n") \
            <(tail -n+2  ${phenotypeFile} | cut -f1) \
            | sort -u | awk '{printf "s/%s/%s/g\n", $1, NR}' > ${workingDir}/gene_set.${chunkNo}/sample_mapping.sed
    fi
    
    # The IDs have to be adjusted, and the non-alphanumeric characters are also removed:
    sed -f ${workingDir}/gene_set.${chunkNo}/sample_mapping.sed ${gene}_output_genotype | perl -lane '$F[0] =~ s/[^0-9a-z]//gi; print join "\t", @F'  > genotype.mod.txt

    # Preparing phenotype file
    sed -f ${workingDir}/gene_set.${chunkNo}/sample_mapping.sed  ${phenotypeFile} | perl -lane 'next if $. == 1; next unless $_ =~ /^\d+/;
        printf "1\t%s\t0\t0\t%s\t%s\t%s\n", $F[0], $F[1], $F[3]' > pheno.mod.txt # Adding sex as coveriate

    # Adjusting the order of the phenotype file according to the samples in the genotype file:
    head -n1 genotype.mod.txt | tr "\t" "\n" | tail -n+2 | perl -lane 'BEGIN {open $pf, "< pheno.mod.txt";
            while ($l = <$pf>){chomp $l;@a = split(/\s/, $l);$h{$a[1]} = $l;}}{print $h{$F[0]} if exists $h{$F[0]}}' > pheno.mod.ordered.txt

    # Now we have to go back to the genotype file and remove those samples that do not have
    # phenotype information:
    export samples=$(grep -w -vf <( cat pheno.mod.ordered.txt | cut -f2 ) \
                   <(head -n1 genotype.mod.txt | tr "\t" "\n" | tail -n+2))

    # cuting out the samples from genotype file:
    head -n1 genotype.mod.txt | tr "\t" "\n" | perl -lane 'BEGIN {foreach $s ( split /\s/, $ENV{"samples"}){$h{$s} = 1;}}{
        push @a, $. unless exists $h{$F[0]}} END{$s = sprintf("cut -f%s genotype.mod.txt > genotype.mod.filtered.txt",join(",", @a));`$s`}'

    export samples=$( cat pheno.mod.ordered.txt | cut -f2 )
    LC_ALL=C cat ${kinshipFile} | perl -lane 'BEGIN{foreach $s ( split /\s/, $ENV{"samples"}){$h{$s} = 1;}
        }{print $_ if exists $h{$F[1]} and exists $h{$F[2]};}' > kinship.mod.filtered.txt

    # Calling MONSTER:
    echo "[Info] MONSTER call: ${MONSTER} \
        -k kinship.mod.filtered.txt \
        -p pheno.mod.ordered.txt \
        -m ${missing_cutoff} \
        -g genotype.mod.filtered.txt \
        -s snpfile.mod.txt \
        ${imputation_method}"

    ${MONSTER} \
        -k kinship.mod.filtered.txt \
        -p pheno.mod.ordered.txt \
        -m ${missing_cutoff} \
        -g genotype.mod.filtered.txt \
        -s snpfile.mod.txt \
        ${imputation_method}
    
    # Checking if output was generated:
    if [[ ! -e MONSTER.out ]]; then
        failed "MONSTER has failed for ${gene}, no output has been generated. Temporary files are kept for troubleshooting."
        continue;
    fi

    # Once MONSTER run is complete, we have to rename output files:
    mv MONSTER.out MONSTER_out_${gene}.out
    mv MONSTER.param MONSTER_out_${gene}.param

    # Checking if the p value is also in the file:
    if [[ $(tail -1 MONSTER_out_${gene}.out | awk '{print NF}') -ne 5 ]]; then
        failed "MONSTER has failed for ${gene}, p-value is missing!. Temporary files are kept for troubleshooting."
        continue;
    fi

    # Pool results together:
    varcnt=$(cut -f3- *output_variants | head -n1 | tr "\t" "\n" | wc -l )
    pval=$(tail -1 MONSTER_out_${gene}.out | cut -f5)

    # Adding p-values to file:
    echo -e "${gene}\t${pval:-NA}\t${varcnt:-NA}" >> ${workingDir}/gene_set.${chunkNo}/results

    # If p-value is really low save into a folder:
    if [[ $(echo $pval | perl -lane 'print $F[0] < $ENV{"signifThreshold"} ? 1 : 0') == 1 ]]; then
        echo "[Info] Found a hit! ${gene} p-value: ${pval}";
        savehit;
    fi

    cd ${workingDir};
    
    # Cleaning up: removing all files or just the mod files if requested.
    if [[ ${keep_temp} == "not" ]]; then
        rm -rf ${workingDir}/gene_set.${chunkNo}/${gene};
    else
        tar -czvf ${workingDir}/gene_set.${chunkNo}/${gene}.tar.gz -C ${workingDir}/gene_set.${chunkNo}/ ${gene}
        rm -rf ${workingDir}/gene_set.${chunkNo}/${gene}
    fi

done

exit 0;
