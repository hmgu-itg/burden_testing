#!/usr/local/bin/bash

# A wrapper script to automate genome-wide burden testing using MONSTER.
# For more information on the applied method see: http://www.stat.uchicago.edu/~mcpeek/software/MONSTER/

# In this version when the script founds a significant gene, it tests if the association is
## driven by a single variant or more by repeating the test with removing a variant
## each time.

version="v10.3 Last modified: 2017.03.15"

# New in this version:
    # The test if a single variant drives the association is now optional
    # New scores are added: Linsight and mixed.
    # Directory structure is more stringent
    # Within a run, traits are separated in subfolders.

##
## todo list for the next release:
##  # 1. We have to modify the test function.
    #### Only one variant is tested the one with the lowest single point p-value.
    #### The genotype of variant will be added as covariate to the phenotype file.
    #### Re-run MONSTER once and report p-value.
    ####

today=$(date "+%Y.%m.%d") # Get the date

# The variant selector script, that generates input for MONSTER:
regionSelector=burden_get_regions.pl

# Folder with all the scripts:
export scriptDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Folder with all the phenotypes:
#phenotypeDir=/lustre/scratch115/projects/t144_helic_15x/analysis/HA/phenotypes
phenotypeDir=/lustre/scratch115/projects/t144_helic_15x/analysis/HA/phenotypes/correct_names.andmissing

# File with all the gene names (only on autosomes):
geneListFile=${scriptDir}/gene_list.lst

# MONSTER link that has to be adjusted:

# Kinship matrix file:
kinshipMatrix=/nfs/team144/ds26/burden_testing/kinship/2016.10.20_fix_diagonal/kinship.fixdiag.txt

# Single point p-values are read from here:
singlePointDir=/lustre/scratch115/projects/t144_helic_15x/analysis/HA/single_point/output

# Path to MONSTER executable:
MONSTER=/nfs/team144/software/MONSTER_v1.3/MONSTER

# By default we are not keeping temporary files:
keep_temp=0

# By default, the first chunk is read, or it can be submitted using -d or form jobindex:
# By default we process all genes at once.
export chunkCount=1 # How many chunks is the gene list split into?

# Which chunk we are processing?
export chunkNo=1 # Default value
if [[ $LSB_JOBINDEX > 0 ]]; then chunkNo=$LSB_JOBINDEX; fi # If jobindex is available, we set that as chunk number.


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
    echo "     -s  - turn weights on. Arguments: GWAVA, CADD, Eigen, EigenPC, EigenPhred, EigenPCPhred, Linsight, Mixed"
    echo "     -x  - extend genomic regions (bp) (de)."
    echo "     -o  - exclude all non loss-of-function variants from the test."
    echo ""
    echo "Gene filter options:"
    echo "     -L  - list file with the gene names (default: ${geneListFile})."
    echo "     -d  - chunk count: how many chunks the gene list should be split into (default: 1)."
    echo "     -c  - chunk number or jobindex (default: 1)."
    echo ""
    echo "General options:"
    echo "     -w  - working directory (default: current working directory)"
    echo "     -b  - Keep temporary files."
    echo "     -a  - Skip testing variants."
    echo ""
    echo "Monster parameters:"
    echo "     -p  - phenotype (required, no default)"
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

# --- if a gene has low p-value, we test if the association is driven by a single variant or not:
function testhit (){

    # Report the test:
    echo "[Info] Testing which variants were driving the association in gene ${gene} for trait ${phenotype}"

    # Printing header:
    echo -e "trait\tSNP\tsingle_pval\tburden_pval\trsID\tallele_string\tconsequence\tMAF" > ${gene}.pvalues.txt

    # Indexing single point output file if not indexed yet:
    if [[ ! -e ${singlePointDir}/${phenotype}/MANOLIS.${phenotype}.assoc.txt.bgz.tbi ]]; then
        echo "[Warning] Single point output file was not indexed. Indexing."
        zcat ${singlePointDir}/${phenotype}/MANOLIS.${phenotype}.assoc.txt.gz | bgzip > ${singlePointDir}/${phenotype}/MANOLIS.${phenotype}.assoc.txt.bgz
        tabix -f -s 1 -b 3 -e 3 -S 1 ${singlePointDir}/${phenotype}/MANOLIS.${phenotype}.assoc.txt.bgz
    fi

    # Now we have to get a list of all variants:
    for pos in $(seq 3 $(head -n1 ${gene}_output_variants | tr "\t" "\n" | wc -l )); do

        # Grepping p value from the corresponding single point associations:
        variant=$(cut -f $pos ${gene}_output_variants | head -n1)
        echo "[Info] ${gene} is a significant gene for ${phenotype}. Now re-testing association without ${variant}"

        tabixQuery=$( echo ${variant} | perl -F"_" -lane 'printf "%s:%s-%s", substr($F[0],3), $F[1], $F[1]' )

        singlePoint=$( tabix ${singlePointDir}/${phenotype}/MANOLIS.${phenotype}.assoc.txt.bgz ${tabixQuery} | cut -f14 )

        # Get variant info:
        variantString=$(grep $(echo ${variant} | perl -F"_" -lane 'printf $F[1]') ${gene}_output_SNP_info | cut -f3-6)

        # Call Monster if more than 2 variants are in the original set:
        if [[ $( cat ${gene}_output_genotype | wc -l ) -gt 3 ]]; then

            # Which is the variant we are excluding:
            variant=$(cut -f${pos} ${gene}_output_variants | head -n1 )

            start=$(( $pos - 1 ))
            end=$(( $pos + 1 ))

            # Exclude variant from the variant dataset:
            cut -f-${start},${end}- snpfile.mod.txt > snpfile.mod.${pos}.txt

            # Exclude variant from the genotype dataset:
            grep -w -v $variant genotype.mod.filtered.txt > genotype.mod.filtered.${pos}.txt

            # Calling Monster:
            #${MONSTER} -k kinship.mod.filtered.txt -p pheno.mod.ordered.txt -m ${missing_cutoff} \
            #    -g genotype.mod.filtered.${pos}.txt -s snpfile.mod.${pos}.txt

            # If everything went fine, we have to extract the p-value:
            pval_hit=$(cut -f5  MONSTER.out  | tail -1)
        fi;

        # Now print out result:
        echo -e "${phenotype}\t${variant}\t${singlePoint:--}\t${pval_hit:--}\t${variantString}" >> ${gene}.pvalues.txt
    done

    return 0
}

# --- Capture command line options --------------------------------------------

# If help message is needed:
if [ $# == 0 ]; then display_help; fi

# Looping through all command line options:
OPTIND=1
while getopts ":hg:m:w:c:l:e:p:s:x:d:ak:t:oL:b" optname; do
    case "$optname" in
      # gene filter parameters:
        "L") geneListFile=${OPTARG} ;;
        "c") export chunkNo=${OPTARG} ;; # If -c is specified it will overwrite the default values.
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
        "a") skiptest=1;;
        "o") lof=1 ;;

      # Other parameters:
        "w") rootDir=${OPTARG} ;;
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
echo -e "\tSingle point association results: ${singlePointDir}"
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
if [ -z $rootDir ]; then
    rootDir=$(pwd)
    echo -e "\tWorking directory: $rootDir"
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
    ## lines for loftee:
    # folder=${folder}".loftee"
    #commandOptions=${commandOptions}" --loftee "
    #echo -e "\tOnly loss-of-function variants will be considered (loftee HC and LC)."

    ## lines for severe variants:
    folder=${folder}".lof"
    commandOptions=${commandOptions}" --lof "
    echo -e "\tOnly variants with severe consequences will be considered (down to missense)."
fi

# 5. Score - If score is not given we apply no score. Otherwise we test the submitted value:
# Accepted scores:
if [[ ! -z ${score} ]]; then
    score=${score^^}
    case ${score} in
        EIGEN )        score="Eigen";;
        EIGENPC )      score="EigenPC";;
        EIGENPHRED )   score="EigenPhred";;
        EIGENPCPHRED ) score="EigenPCPhred";;
        CADD )         score="CADD";;
        GWAVA )        score="GWAVA";;
        LINSIGHT )     score="Linsight";;
        MIXED )        score="Mixed";;
        * )            score="noweight";;
    esac
else
    echo "[Warning] Submitted score is not recognized! Accepted scores: GWAVA, CADD, Eigen, EigenPC, EigenPhred, EigenPCPhred."
    echo "[Warning] No scores are being applied."
    score="noweights"
fi
echo -e "\tWeighting: ${score}"

# Only adding score to command line if score is requested:
if [ ${score} != "noweights" ]; then commandOptions=${commandOptions}" -s ${score}"; fi

# If Eigen score is applied, we automatically shift scores by 1:
if [ ${score} == "Eigen" ]; then commandOptions=${commandOptions}" --shift 1 "; fi

# Adjust dir name:
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
fi
echo -e "\tPhenotype: ${phenotype}"
echo -e "\tPhenotype file: ${phenotypeDir}/MANOLIS.${phenotype}.txt"

# Chunk number can be explicitly submitted, or read from the LSB_JOBINDEX
#if [[ ${LSB_JOBINDEX} -ne 0 ]]; then chunkNo=${LSB_JOBINDEX}; fi;
chunkSize=$(wc -l ${geneListFile} | perl -lane 'print int($F[0]/$ENV{"chunkCount"}+0.5)')

# Reporting chunk size and count
echo -e "\n[Info] Gene slice report:"
echo -e "\tNumber of chunks the gene list is split into: ${chunkCount}"
echo -e "\tCurrent chunk: ${chunkNo}"
echo -e "\tNumber of genes in one chunk: ${chunkSize}\n"

# Updating working dir, and create folder:
folder=$( echo $folder | perl -lane '$_ =~ s/^\.//; print $_')
workingDir=${rootDir}/${folder}/Pheno.${phenotype}
mkdir -p ${rootDir}

# printing out report:
echo -e "[Info] output folder: ${workingDir}"
echo "[Info] command line options for burden get region: ${commandOptions}"

# Testing kinship matrix:
if [[ ! -e ${kinshipMatrix} ]]; then
    echo "[Error] Kinship could not be opened! ${kinshipMatrix}";
    exit 1;
fi

# Do we want to perform test for all variations:
if [ -z $skiptest ]; then
    echo -e "[Info] Corresponding variants will be tested to see if a single variant drives the associations."
else
    echo -e "[Info] Variants will not be tested."
fi

# --- Main loop executed for all genes --------------------------------------------

# Creating results file with header:
mkdir -p ${workingDir}/gene_set.${chunkNo}
if [[ ! -e ${workingDir}/gene_set.${chunkNo}/results ]]; then 
    echo -e "GeneName\tp-value\tSNP_count" >> ${workingDir}/gene_set.${chunkNo}/results
fi

# Looping through all genes in the gene set:
awk -v cn="${chunkNo}" -v cs="${chunkSize}" 'NR > (cn-1)*cs && NR <= cn*cs' ${geneListFile} | while read gene ; do

    # Removing colons from genomic coordinates:
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

    # Assigning default cutoff values for monster (will be moved up somewhere..):
    if [[ ! -n ${missing_cutoff} ]]; then missing_cutoff=1; fi

    # Echoing MONSTER parameters: (Will be located above somewhere...)
    echo "Optional MONSTER parameters: missingness: ${missing_cutoff}, imputation method: ${imputation_method}"


    # The IDs have to be adjusted: the non-alphanumeric characters are also removed:
    cat ${gene}_output_variants | perl -lane '$_ =~ s/[^0-9a-z\t\.]//gi; print $_'  > snpfile.mod.txt

    # The IDs have to be adjusted, and the non-alphanumeric characters are also removed:
    sed -f ${scriptDir}/HELIC.to.Num_genotype.sed ${gene}_output_genotype | perl -lane '$F[0] =~ s/[^0-9a-z]//gi; print join "\t", @F'  > genotype.mod.txt

    # Preparing phenotype file
    export PhenoFile=${phenotypeDir}/MANOLIS.${phenotype}.txt
    sed -f ${scriptDir}/HELIC.to.Num.sed  ${PhenoFile}| perl -lane 'next if $. == 1; next unless $_ =~ /^\d+/;
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
    LC_ALL=C cat ${kinshipMatrix} | perl -lane 'BEGIN{foreach $s ( split /\s/, $ENV{"samples"}){$h{$s} = 1;}
        }{print $_ if exists $h{$F[1]} and exists $h{$F[2]};}' > kinship.mod.filtered.txt

    # Calling MONSTER:
    echo "[Info] MONSTER call: ${MONSTER} -k kinship.mod.filtered.txt \
-p pheno.mod.ordered.txt \
-m ${missing_cutoff} \
-g genotype.mod.filtered.txt \
-s snpfile.mod.txt \
${imputation_method}"

    ${MONSTER} -k kinship.mod.filtered.txt \
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

    # Adding p-values to file:
    echo -e "${gene}\t${pval:--}\t${varcnt:--}" >> ${workingDir}/gene_set.${chunkNo}/results

    # If p-value is really low save into a folder:
    if [[ $(echo $pval | perl -lane 'print int(-log(abs($F[0]))/log(10))') -ge 4 ]]; then
        echo "[Info] Found a hit! ${gene} p-value: ${pval}";
        if [ -z $skiptest ]; then
            echo "[Info] Testing if a single variant dirves the associations..."
            testhit; # Run test with removing variants.
        fi
        savehit; # Saving all intermediate files.
    fi

    # Cleaning up: removing all files or just the mod files if requested.
    if [[ ${keep_temp} -eq 1 ]]; then
        echo "cicaful"
        # rm ${workingDir}/gene_set.${chunkNo}/${gene}/*mod*;
    else
        cd ${workingDir};
        rm -rf ${workingDir}/gene_set.${chunkNo}/${gene};
    fi

done

exit 0
