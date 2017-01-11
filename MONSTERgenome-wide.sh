#!/usr/local/bin/bash

# To be fixed in a later version: 2017.01.05
    # Constant values are sourced from a config file.
    #

# A wrapper script to automate genome-wide burden testing using MONSTER.
# For more information on the applied method see: http://www.stat.uchicago.edu/~mcpeek/software/MONSTER/

version="v9.2 Last modified: 2017.01.05"
today=$(date "+%Y.%m.%d") # Get the date

# The variant selector script, that generates input for MONSTER:
regionSelector=burden_get_regions.pl

# Folder with all the scripts:
export scriptDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Folder with all the phenotypes:
phenotypeDir=/lustre/scratch115/projects/t144_helic_15x/analysis/HA/phenotypes

# File with all the gene names (only on autosomes):
geneListFile=${scriptDir}/gene_list.lst

# MONSTER link that has to be adjusted:

# Kinship matrix file:
kinshipMatrix=/nfs/team144/ds26/burden_testing/kinship/2016.10.20_fix_diagonal/kinship.fixdiag.txt

# By default, the first chunk is read, or it can be submitted using -d or form jobindex:
# By default we process all genes at once.
export chunkNo=1 # Which chunk we are processing?
export chunkCount=1 # How many chunks is the gene list split into?

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
    echo "     -s  - turn weights on. Arguments: CADD, Eigen or GWAVA."
    echo "     -x  - extend genomic regions (bp)."
    echo "     -o  - exclude all non loss-of-function variants from the test."
    echo ""
    echo "Gene filter options:"
    echo "     -L  - list file with the gene names."
    echo "     -d  - chunk count: how many chunks the gene list should be split into."
    echo "     -c  - chunk number or jobindex"
    echo ""
    echo "General options:"
    echo "     -w  - working directory or actual directory"
    echo "     -b  - Keep temporary files."
    echo ""
    echo "Monster parameters:"
    echo "     -p  - phenotype"
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
    ( >&2 echo "[Error] Gene has failed: ${gene}" ) # Reporting to the error log as well.
    echo "[Error] ${message}"
    echo "[Error] folder moved to ${folder}/failed"
    echo ""

    # Generating failed directory and compress gene dir:
    mkdir -p ${workingDir}/failed
    tar -czvf ${workingDir}/failed/${gene}.tar.gz -C ${workingDir}/gene_set.${chunkNo}/ ${gene}

    # Removing directory:
    rm -rf ${workingDir}/gene_set.${chunkNo}/${gene}
}

# --- If a p-value is really low, the whole run backed up:
function savehit (){
    mkdir -p ${workingDir}/hits
    tar -czvf ${workingDir}/hits/${gene}.tar.gz -C ${workingDir}/gene_set.${chunkNo}/ ${gene}
}


# --- Capture command line options --------------------------------------------

# If help message is needed:
if [ $# == 0 ]; then display_help; fi

# Looping through all command line options:
OPTIND=1
while getopts ":hg:m:w:c:l:e:p:s:x:d:k:t:oL:b" optname; do
    case "$optname" in
      # gene filter parameters:
        "L") geneListFile=${OPTARG} ;;
        "c") export chunkNo=${OPTARG} ;;
        "d") chunkCount=${OPTARG} ;;

      # MONSTER parameters:
        "p") phenotype=${OPTARG} ;;

      # Wrapper option:
        "b") keep_temp=1;;

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

      # Other parameters:
        "w") workingDir=${OPTARG} ;;
        "h") display_help ;;
        "?") display_help "[Error] Unknown option $OPTARG" ;;
        ":") display_help "[Error] No argument value for option $OPTARG";;
        *) display_help "[Error] Unknown error while processing options";;
    esac
done

# --- Process and test parameters --------------------------------------------

# Reporting parameters:
echo "##"
echo "## Genome-wide Monster wrapper version ${version}"
echo "## Run date: ${today}"
echo "##"
echo ""
echo "[Info] Files used by the script:"
echo -e "\tVariant selector: ${regionSelector}"
echo -e "\tGene list: ${geneListFile}"
echo -e "\tPhenotype folder: ${phenotypeDir}"
echo -e "\tScript dir: ${scriptDir}"
echo -e "\tKinship matrix: ${kinshipMatrix}"
echo ""


# Initializing folder names and command line parameters:
folder=""
commandOptions=""

# Exiting if the gene list could not opened:
if [[ ! -e $geneListFile ]]; then
    echo "[Error] $geneListFile could not opened. Exiting."
    exit;
fi

echo "[Info] The following parameters were submitted:"
# 1. GENCODE -expecting a list of feature names separated by a comma.
if [ ! -z $gencode ]; then
    folder=${gencode}
    commandOptions=" -G ${gencode}"
    echo -e "\tGENCODE feaures: ${gencode}"
fi

# If there is no specified output folder, use current working dir:
if [ -z $workingDir ]; then
    workingDir=$(pwd)
    echo -e "\tWorking directory: $workingDir"
fi

# Reporting if temporary files are saved:
if [[ ${keep_temp} -eq 1 ]]; then echo -e "\tTemporary files will not be deleted!"; fi

# 2. GTEx - expecting a list of feature names separeted by comma.
if [ ! -z $gtex ]; then
    folder=${folder}".GTEX."${gtex}
    commandOptions=${commandOptions}" -E ${gtex}"
    echo -e "\tGTEx feaures: ${gtex}"
fi

# 3. Overlap - expecting a list of features separeated by comma.
if [ ! -z $overlap ]; then
    folder=${folder}".Overlap."${overlap}
    commandOptions=${commandOptions}" -L ${overlap}"
    echo -e "\tOverlapping reg.feaures: ${overlap}"
fi

# 4. MAF - expecting a float between 0 and 0.5
if [ ! -z $MAF ]; then
    folder=${folder}".MAF."${MAF}
    commandOptions=${commandOptions}" --maf ${MAF}"
    echo -e "\tUpper minor allele frequency: ${MAF}"
fi

# 5. If lof is set, we only need loss of function variants:
if [ ! -z $lof ]; then
    folder=${folder}".loftee"
    commandOptions=${commandOptions}" --loftee "
    echo -e "\tOnly loss-of-function variants will be considered (loftee HC and LC)."
fi

# 5. Score - If score is not given we apply no score. Otherwise we test the submitted value:
if [[ ! -z ${score} ]]; then
    score=${score^^}
    case ${score} in
        EIGEN) score="Eigen";;
        CADD) score="CADD";;
        GWAVA) score="GWAVA";;
        *) score="noweight";;
    esac
else
    score="noweights"
fi
echo -e "\tWeighting: ${score}"

# Only adding score to command line if score is requested:
if [ ${score} != "noweights" ]; then commandOptions=${commandOptions}" -s ${score}"; fi
folder=${folder}"."${score}

# 6. Extend regions:
if [ ! -z $xtend ]; then
    folder=${folder}".Xtend."${xtend}
    commandOptions=${commandOptions}" --extend ${xtend}"
    echo -e "\tFeatures are extended by ${xtend}bp"
fi

# 7. Set score cutoff:
if [ ! -z ${cutoff} ]; then
    folder=${folder}".cutoff."${cutoff}
    echo -e "\tScore cutoff: ${cutoff}"
    commandOptions=${commandOptions}" --cutoff ${cutoff}"
fi
if [ ! -z ${scoreshift} ]; then
    folder=${folder}".shift."${scoreshift}
    echo -e "\tScores shifted by: ${scoreshift}"
    commandOptions=${commandOptions}" --shift ${scoreshift}"
fi

# 8. If phenotype is not given the script dies:
if [ -z $phenotype ]; then
    echo "[Error] Phenotype was not set! Exiting"
    exit 1
elif [ ! -e "${phenotypeDir}/MANOLIS.${phenotype}.txt" ]; then
    echo "[Error] Phenotype file was not found! ${phenotypeDir}/MANOLIS.${phenotype}.txt"
    exit 1
else
    folder=${folder}".Pheno."${phenotype}
fi
echo -e "\tPhenotype: ${phenotype}"
echo -e "\tPhenotype file: ${phenotypeDir}/MANOLIS.${phenotype}.txt"

# Chunk number can be explicitly submitted, or read from the LSB_JOBINDEX
if [[ ! -z ${LSB_JOBINDEX} ]]; then chunkNo=${LSB_JOBINDEX}; fi;
chunkSize=$(wc -l ${geneListFile} | perl -lane 'print int($F[0]/$ENV{"chunkCount"}+0.5)')

# Reporting chunk size and count
echo -e "\n[Info] Gene slice report:"
echo -e "\tNumber of chunks the gene list is split into: ${chunkCount}"
echo -e "\tCurrent chunk: ${chunkNo}"
echo -e "\tNumber of genes in one chunk: ${chunkSize}\n"

# printing out report:
echo -e "[Info] output folder: ${workingDir}/${folder}"
echo "[Info] command line options for burden get region: ${commandOptions}"

# Updating working dir, and create folder:
workingDir=${workingDir}/${folder}
mkdir -p ${workingDir}

# Testing kinship matrix:
if [[ ! -e ${kinshipMatrix} ]]; then
    echo "[Error] Kinship could not be opened! ${kinshipMatrix}";
    exit 1;
fi

# --- Main loop executed for all genes --------------------------------------------

# Looping through all genes in the gene set:
awk -v cn="${chunkNo}" -v cs="${chunkSize}" 'NR > (cn-1)*cs && NR <= cn*cs' ${geneListFile} | while read gene ; do

    # Testing if a run was already completed:
    #flag=$( if [[ -e ${workingDir}/gene_set.${chunkNo}/results ]]; then grep -w ${gene} ${workingDir}/gene_set.${chunkNo}/results; fi)
    #if [[ ! -z "${flag}" ]]; then echo "$gene is done! Next."; continue; fi

    echo ""
    echo "[Info] Processing $gene"
    (>&2 echo "[Info] Processing $gene" ) # Adding mark into the error file, so we'll know in which gene the error comes up.
    echo "[Info] Step 1.: extracting variants given the list of parameters."

    # Creating folders for the parameters:
    mkdir -p ${workingDir}/gene_set.${chunkNo}/${gene}
    cd ${workingDir}/gene_set.${chunkNo}/${gene}

    # Reporting call:
    echo "${scriptDir}/${regionSelector}  -i $gene -o ${gene}_output ${commandOptions} -v > ${gene}_log"
    perl ${scriptDir}/${regionSelector}  -i $gene -o ${gene}_output ${commandOptions} -v > ${gene}_log

    # We have to test if the output files are OK, then we go to the next gene:
    if [[ ! -e ${gene}_output_genotype ]]; then
        failed "${gene} has failed, no genotype output generated. Gene skipped."
        cd .. && rm -rf ${gene}
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

    # Assigning default cutoff values for monster (will be moved up somewhere..):
    if [[ ! -n ${missing_cutoff} ]]; then missing_cutoff=0; fi

    # Echoing MONSTER parameters: (Will be located above somewhere...)
    echo "Optional MONSTER parameters: missingness: ${missing_cutoff}, imputation method: ${imputation_method}"


    # The IDs have to be adjusted: the non-alphanumeric characters are also removed:
    cat ${gene}_output_variants | perl -lane '$_ =~ s/[^0-9a-z\t\.]//gi; print $_'  > snpfile.mod.txt

    # The IDs have to be adjusted, and the non-alphanumeric characters are also removed:
    sed -f ${scriptDir}/HELIC.to.Num_genotype.sed ${gene}_output_genotype | perl -lane '$F[0] =~ s/[^0-9a-z]//gi; print join "\t", @F'  > genotype.mod.txt

    ## Preparing phenotype file
    export PhenoFile=${phenotypeDir}/MANOLIS.${phenotype}.txt
    sed -f ${scriptDir}/HELIC.to.Num.sed  ${PhenoFile}| perl -lane 'next if $. == 1; next unless $_ =~ /^\d+/;
        printf "1\t%s\t0\t0\t%s\t%s\t%s\n", $F[0], $F[1], $F[3], $F[1]' > pheno.mod.txt # Adding sex as coveriate

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
    LC_ALL=C cat ${kinshipMatrix} | perl -lane 'BEGIN{foreach $s ( split /\s/, $ENV{"samples"}){$h{$s} = 1;}
        }{print $_ if exists $h{$F[1]} and exists $h{$F[2]};}' > kinship.mod.filtered.txt

    # Calling MONSTER:
    echo "[Info] MONSTER call: /nfs/team144/software/MONSTER/MONSTER -k kinship.mod.filtered.txt \
-p pheno.mod.ordered.txt \
-m ${missing_cutoff} \
-g genotype.mod.filtered.txt \
-s snpfile.mod.txt \
${imputation_method}"

    /nfs/team144/software/MONSTER/MONSTER -k kinship.mod.filtered.txt \
        -p pheno.mod.ordered.txt \
        -m ${missing_cutoff} \
        -g genotype.mod.filtered.txt \
        -s snpfile.mod.txt \
        ${imputation_method}

    # Once MONSTER run is complete, we have to rename output files:
    mv MONSTER.out MONSTER_out_${gene}.out
    mv MONSTER.param MONSTER_out_${gene}.param

    # Checking if output was generated:
    if [[ ! -e MONSTER_out_${gene}.out ]]; then
        failed "MONSTER has failed for ${gene}, no output has been generated. Temporary files are kept for troubleshooting."
        continue;
    fi

    # Checking if al the p value is also in the file:
    if [[ $(tail -1 MONSTER_out_${gene}.out | awk '{print NF}') -ne 5 ]]; then
        failed "MONSTER has failed for ${gene}, p-value calculated!. Temporary files are kept for troubleshooting."
        continue;
    fi

    # Pool results together:
    varcnt=$(cut -f3- *output_variants | head -n1 | tr "\t" "\n" | wc -l )
    pval=$(tail -1 MONSTER_out_${gene}.out | cut -f5)

    # If p-value is really low save into a folder:
    if [[ $(echo $pval | perl -lane 'print int(-log(abs($F[0]))/log(10))') -ge 6 ]]; then echo "[Info] Found a hit! ${gene} p-value: ${pval}"; savehit; fi

    # Saving file:
    echo -e "${gene}\t${pval}\t${varcnt}" >> ../results

    # Cleaning up: removing all files or just the mod files if requested.
    if [[ ${keep_temp} -eq 1 ]]; then rm *mod*; else cd .. && rm -rf ${gene}; fi

done

exit 0