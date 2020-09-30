#!/bin/bash

# if pattern corresponds to a filename (doesn't contain %), check if the file exists
# otherwise check if files for each chromosome (1-22) exist

function testVCFs {
    pattern=$1
    c=0
    if [[ $pattern =~ % ]];then
	for i in $(seq 1 22);do
	    fname=$(echo $pattern|sed "s/\%/$i/")
	    ifname=${fname}".tbi"
	    if [[ ! -e ${fname} ]];then
		echo `date "+%Y.%b.%d_%H:%M"` "[Warning] No VCF for chromosome $i"
	    else
		if [[ ! -e ${ifname} ]];then
		    echo `date "+%Y.%b.%d_%H:%M"` "[Error] No index file for $fname"
		    return 1
		fi
		c=$((c+1))
	    fi
	done
	if [[ $c -gt 0 ]];then # some VCFs exist; trying to work with them
	    return 0
	else
	    return 1
	fi
    else
	if [[ -e $pattern ]];then
	    return 0
	else
	    echo `date "+%Y.%b.%d_%H:%M"` "[Error] $pattern does not exist"
	    return 1
	fi
    fi
}


version="v12 Last modified: 2020.Feb.14"
today=$(date "+%Y.%b.%d")

# Folder with the variant selector script:
scriptDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
regionSelector=${scriptDir}/"Burden_testing.pl"

imputation_method='-A' # The default imputation method is BLUP, slowest, but the most accurate. For other options, see MONSTER documentation.
configFile=""

chunksTotal=1
chunkNo=""
MAF=0.05 # By default this is the upper minor allele frequency.

# --- print out help message and exit:
function display_help() {
    echo "$1"
    echo ""
    echo "Variant selector"
    echo "version: ${version}"
    echo ""
    echo "Variant filter options:"
    echo "     -C  - config file (reguired, no default)"
    echo "     -V  - VCF file(s) (required, no default; use % character for chromosome name eg 'chr%.vcf.gz')"
    echo "     -g  - comma separated list of GENCODE features (gene, exon, transcript, CDS or UTR)"
    echo "     -e  - comma separated list of GTEx features (promoter, CTCF, enhancer, promoterFlank, openChrom, TF_bind or allreg)"
    echo "     -l  - comma separated list of overlap features (promoter, CTCF, enhancer, promoterFlank, openChrom, TF_bind or allreg)"
    echo "     -m  - upper MAF thresholds (default: 0.05)"
    echo "     -x  - extend genomic regions by this amount (bp) (default: 0)"
    echo "     -o  - include variants with severe consequences only (more severe than missense)"
    echo "     -f  - include only HC and LC loftee variants"
    echo "     -j  - include only HC loftee variants"
    echo ""
    echo "Parameters to set up scores for variants:"
    echo "     -s  - apply weighting; accepted values: CADD, EigenPhred"
    echo "     -t  - the value with which the scores will be shifted (default value: if Eigen score weighting specified: 1, otherwise: 0)"
    echo "     -k  - below the specified cutoff value, the variants will be excluded (default: 0)"
    echo ""
    echo "Gene list and chunking:"
    echo "     -L  - file with gene IDs (if not specified all genes will be analyzed)"
    echo "     -d  - total number of chunks (default: 1)"
    echo "     -c  - chunk number (default: 1); takes precedence over SLURM_ARRAY_TASK_ID etc."
    echo ""
    echo "General options:"
    echo "     -w  - output directory where the chunk subdirectories will be created (required, no default)"
    echo ""
    echo "Other options:"
    echo "     -h  - print this message and exit"
    echo ""
    echo ""

    exit 1
}

# --- Capture command line options --------------------------------------------

if [ $# == 0 ]; then display_help; fi

OPTIND=1
score=""
geneListFile=""

while getopts ":hL:c:d:bg:m:s:l:e:x:k:t:ofw:jC:V:" optname; do
    case "$optname" in
      # Gene list related parameters:
        "L") geneListFile=${OPTARG} ;;
        "c") chunkNo=${OPTARG} ;;
        "d") chunksTotal=${OPTARG} ;;

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
        "j") lofteeHC=1 ;;
        "C") configFile=${OPTARG} ;;
        "V") vcfFile=${OPTARG} ;;

      # Other parameters:
        "w") outputDir=${OPTARG} ;;
        "h") display_help ;;
        "?") display_help "[Error] Unknown option $OPTARG" ;;
        ":") display_help "[Error] No argument value for option $OPTARG" ;;
        *) display_help "[Error] Unknown error while processing options" ;;
    esac
done

if [[ -z "${outputDir}" ]]; then
    echo `date "+%Y.%b.%d_%H:%M"` "[Error] Output directory not specified";
    exit 1
fi

#----------------------------------------- CHUNKS --------------------------------------------------------

chunk_warning=""
if [[ ! -z ${SLURM_ARRAY_TASK_ID} ]];then
    if [[ ! -z ${chunkNo} ]];then
	chunk_warning="WARNING: both SLURM_ARRAY_TASK_ID and chunkNo ( -c ) are defined; using chunkNo"
    else
	chunkNo=${SLURM_ARRAY_TASK_ID}
    fi
elif [[ ! -z ${LSB_JOBINDEX} ]];then
    if [[ ! -z ${chunkNo} ]];then
	chunk_warning="WARNING: both LSB_JOBINDEX and chunkNo ( -c ) are defined; using chunkNo"
    else
	chunkNo=${LSB_JOBINDEX}
    fi
else
    if [[ -z ${chunkNo} ]];then
	chunkNo=1 # default
    fi
fi

if [[ ${chunkNo} -gt ${chunksTotal} ]];then
    echo "ERROR: current chunk number ($chunkNo) is greater than the total number of chunks ($chunksTotal). EXIT"
    exit 1
fi

#---------------------------------------------------------------------------------------------------------

if [[ -z "${configFile}" ]]; then
    echo `date "+%Y.%b.%d_%H:%M"` "[Error] Config file was not specified"
    exit 1
fi

if [[ ! -e "${configFile}" ]]; then
    echo `date "+%Y.%b.%d_%H:%M"` "[Error] Config file does not exist: $configFile"
    exit 1
fi

commandOptions=" --config ${configFile} "

#------------------------------------ INPUT GENE LIST ---------------------------------------

no_list_warning=""
if [[ -z ${geneListFile} ]];then
    gencode_file=$(grep "^gencode_file" ${configFile} | cut -f 2 -d '=')
    zcat ${gencode_file} | cut -f 4 > temp_gene_list.txt
    geneListFile="temp_gene_list.txt"
    no_list_warning="[Warning] No gene list specified; using all genes from $gencode_file"
fi

if [[ ! -e "${geneListFile}" ]]; then
    echo `date "+%Y.%b.%d_%H:%M"` "[Error] Gene list file could not be opened: $geneListFile"
    exit 1
fi

#-------------------------------------- VCFs -------------------------------------------------

if [[ -z "${vcfFile}" ]]; then
    echo "[Error] No VCF file(s) specified!"
    exit 1
else
    testVCFs ${vcfFile}
    if [[ $? -eq 1 ]];then
	exit 1
    fi

    commandOptions="${commandOptions} --vcf ${vcfFile} "
fi

#-------------------------------- GENCODE features -------------------------------------------

if [[ ! -z "${gencode}" ]]; then
    folder="${gencode}"
    commandOptions="${commandOptions} --GENCODE ${gencode}"
fi

#--------------------------------- GTEx features ---------------------------------------------

if [[ ! -z "$gtex" ]]; then
    folder="${folder}.GTEX.${gtex}"
    commandOptions="${commandOptions} --GTEx ${gtex}"
fi

#------------------------------- Overlap features --------------------------------------------

if [[ ! -z "${overlap}" ]]; then
    folder="${folder}.Overlap.${overlap}"
    commandOptions="${commandOptions} --overlap ${overlap}"
fi

#--------------------------------- MAF threshold ---------------------------------------------

if [[ ! -z "$MAF" ]]; then
    folder="${folder}.MAF.${MAF}"
    commandOptions="${commandOptions} --maf ${MAF}"
fi

#---------------------------------- Cpnsequence options -------------------------------------

if [[ ! -z "$loftee" ]]; then
    folder="${folder}.loftee"
    commandOptions="${commandOptions} --loftee "
fi

if [[ ! -z "$lofteeHC" ]]; then
    folder="${folder}.lofteeHC"
    commandOptions="${commandOptions} --lofteeHC "
fi

if [[ ! -z "$lof" ]]; then
    folder="${folder}.severe"
    commandOptions="${commandOptions} --lof "
fi

#------------------------------------------ SCORES -------------------------------------------

warning1=""
warning2=""
score_tmp=""
# Score - If score is not given we apply no score. Otherwise we test the submitted value:
# Accepted scores:
if [[ ! -z "${score}" ]]; then
    score="${score^^}"
    case "${score}" in
#        EIGEN )        score="Eigen";;
#        EIGENPC )      score="EigenPC";;
        EIGENPHRED )   score="EigenPhred";;
#        EIGENPCPHRED ) score="EigenPCPhred";;
        CADD )         score="CADD";;
#       MIXED )        score="Mixed";;
        * )            score_tmp="noweight";;
    esac
else
    score="noweight"
fi

if [[ ! -z ${score_tmp} ]];then
    warning1="[Warning] Submitted score name ($score) is not recognized! Accepted scores: CADD, EigenPhred"
    warning2="[Warning] No scoring will be applied"
    score="noweight"
fi

# Only adding score to command line if score is requested:
if [[ "${score}" != "noweight" ]]; then
    commandOptions="${commandOptions} --score ${score}";
fi

# If Eigen score is applied, we shift the scores by 1, if no other value is specified:
if [[ ("${score}" == "Eigen") && (-z "${scoreshift}" ) ]]; then scoreshift=1; fi

# Adjust dir name, noweight if no score is specified:
folder="${folder}.${score}"

#----------------------------------------------------------------------------------------------

# GENCODE features might be extended with a given number of bps:
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

# Updating working dir, and creating folder:
folder=$( echo $folder | perl -lane '$_ =~ s/^\.//;$_ =~ s/,/_/g; print $_;')
outputDir="${outputDir}/${folder}"
outputDir=${outputDir}/gene_set.${chunkNo}
mkdir -p ${outputDir}
if [[ ! -d ${outputDir} ]]; then
    echo `date "+%Y.%b.%d_%H:%M"` "[Error] Chunk directory (${outputDir}) could not be created. Exiting."
    exit 1
fi

LOGFILE=${outputDir}/"variant_selector.log"

if [[ ! -z ${warning1} ]];then
    echo `date "+%Y.%b.%d_%H:%M"` ${warning1} >> ${LOGFILE}
fi

if [[ ! -z ${warning2} ]];then
    echo `date "+%Y.%b.%d_%H:%M"` ${warning2} >> ${LOGFILE}
fi

if [[ ! -z ${chunk_warning} ]];then
    echo `date "+%Y.%b.%d_%H:%M"` ${chunk_warning} >> ${LOGFILE}
fi

if [[ ! -z ${no_list_warning} ]];then
    echo `date "+%Y.%b.%d_%H:%M"` ${no_list_warning} >> ${LOGFILE}
fi

# --------------------------------------------- GENE LIST FOR THE GIVEN CHUNK ------------------------------

totalGenes=$(cat ${geneListFile} | wc -l)
if [ $totalGenes -lt $chunksTotal ];then
    echo `date "+%Y.%b.%d_%H:%M"` "[Warning] Number of chunks ($chunksTotal) is larger than number of genes in the gene list ($totalGenes) "  >> ${LOGFILE}
    echo `date "+%Y.%b.%d_%H:%M"` "[Warning] Analyzing all genes in one chunk"  >> ${LOGFILE}
    chunkNo=1
    cat ${geneListFile} > ${outputDir}/input_gene.list    
else
    rem=$(( totalGenes % chunksTotal ))
    chunkSize=$(( totalGenes / chunksTotal ))
    lastChunkSize=$(( chunkSize + rem ))
    if [[ ${chunkNo} -eq ${chunksTotal} ]];then
	tail -n ${lastChunkSize} ${geneListFile} > ${outputDir}/input_gene.list
    else
	awk -v cn="${chunkNo}" -v cs="${chunkSize}" 'NR > (cn-1)*cs && NR <= cn*cs' ${geneListFile} > ${outputDir}/input_gene.list
    fi

    n=$( cat ${outputDir}/input_gene.list | wc -l)
    if [[ $n -eq 0 ]];then
	echo "Chunk ${chunkNo} is empty; EXIT" >> ${LOGFILE}
	echo `date "+%Y.%b.%d_%H:%M"` "[Info] VARIANT SELECTION DONE" >> ${LOGFILE}
	exit 0
    fi
fi

# ------------------------------------------------ Reporting parameters ------------------------------------

echo `date "+%Y.%b.%d_%H:%M"` "##"  >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"` "## Variant selector version ${version}" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"` "## Date: ${today}" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"` "##" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"` "" >> ${LOGFILE}

echo `date "+%Y.%b.%d_%H:%M"` "[Info] General options:" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"`  "Variant selector: ${regionSelector}" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"`  "Script dir: ${scriptDir}" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"`  "Output directory: ${outputDir}" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"` "" >> ${LOGFILE}

echo `date "+%Y.%b.%d_%H:%M"` "[Info] Gene list options:" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"`  "Gene list file: ${geneListFile}" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"`  "Number of chunks the gene list is split into: ${chunksTotal}" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"`  "Current chunk: ${chunkNo}" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"`  "Number of genes in one chunk: ${chunkSize}" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"` "" >> ${LOGFILE}

echo `date "+%Y.%b.%d_%H:%M"` "[Info] Variant filtering options:" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"`  "vcf file: ${vcfFile}" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"`  "GENCODE feaures: ${gencode:--}" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"`  "GTEx feaures: ${gtex:--}" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"`  "Overlapping reg.features: ${overlap:-NA}" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"`  "Features are extended by ${xtend:-0}bp" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"`  "Upper minor allele frequency: ${MAF:-1}" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"` "" >> ${LOGFILE}

echo `date "+%Y.%b.%d_%H:%M"` "[Info] Weighting options:" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"`  "Weighting: ${score}" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"`  "Score cutoff: ${cutoff:-0}" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"`  "Scores shifted by: ${scoreshift:-0}" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"` "" >> ${LOGFILE}

echo `date "+%Y.%b.%d_%H:%M"` "[Info] command line options for the variant selecor: ${commandOptions}" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"` "" >> ${LOGFILE}

#------------------------------------------------- CALLING VARIANT SELECTOR ------------------------------------------------

selectorLog=${outputDir}/chunk_${chunkNo}.output.log
echo `date "+%Y.%b.%d_%H:%M"` "Calling ${regionSelector}  --input ${outputDir}/input_gene.list --output gene_set_output --output-dir ${outputDir} ${commandOptions} --verbose > ${selectorLog} 2 > ${selectorLog}"  >> ${LOGFILE}
${regionSelector} --input ${outputDir}/input_gene.list --output gene_set_output --output-dir ${outputDir} ${commandOptions} --verbose > ${selectorLog} 2 > ${selectorLog}

#------------------------------------------------------- CHECKING OUTPUT ---------------------------------------------------

# We expect to get 2 files: gene_set_output_genotype_file.txt and gene_set_output_SNPinfo_file.txt

echo `date "+%Y.%b.%d_%H:%M"` "[Info] Checking output..." >> ${LOGFILE}

cd ${outputDir}

# We have to check if both files are generated AND they have enough lines.
gene_notenough=$(cat ${selectorLog} | grep -c NOT_ENOUGH_VAR)
gene_toomany=$(cat ${selectorLog} | grep -c TOO_MANY_VAR)
gene_noremain=$(cat ${selectorLog} | grep -c NO_VAR_REMAIN)
gene_absent=$(cat ${selectorLog} | grep -c NO_GENE)
region_absent=$(cat ${selectorLog} | grep -c NO_REGION)

echo `date "+%Y.%b.%d_%H:%M"` -e "[Info] WARNINGS/ERRORS FROM VARIANT SELECTOR" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"` -e "[Info] =====================================" >> ${LOGFILE}

if [[ "$gene_notenough" -ne 0 ]]; then
        echo `date "+%Y.%b.%d_%H:%M"` -e "[Warning] Not enough variants [NOT_ENOUGH_VAR]: $(cat ${selectorLog} | grep NOT_ENOUGH_VAR | sed 's/.*Gene.//;s/ .*//' | tr '\n' ' ')" >> ${LOGFILE}
fi
if [[ "$gene_toomany" -ne 0 ]]; then
        echo `date "+%Y.%b.%d_%H:%M"` -e "[Warning] Too many variants [TOO_MANY_VAR]: $(cat ${selectorLog} | grep TOO_MANY_VAR | sed 's/.*Gene.//;s/ .*//'| tr '\n' ' ')" >> ${LOGFILE}
fi
if [[ "$gene_noremain" -ne 0 ]]; then
        echo `date "+%Y.%b.%d_%H:%M"` -e "[Warning] [Warning] No variants after scoring [NO_VAR_REMAIN] for genes:: $(cat ${selectorLog} | grep NO_VAR_REMAIN | sed 's/.*Gene.//;s/ .*//'| tr '\n' ' ')" >> ${LOGFILE}
fi
if [[ "$gene_absent" -ne 0 ]]; then
        echo `date "+%Y.%b.%d_%H:%M"` -e "[Warning] Gene name unknown [NO_GENE]: $(cat ${selectorLog} | grep NO_GENE | sed 's/.*Gene.//;s/ .*//'| tr '\n' ' ')" >> ${LOGFILE}
fi
if [[ "$region_absent" -ne 0 ]]; then
        echo `date "+%Y.%b.%d_%H:%M"` -e "[Warning] No region in gene [NO_REGION]: $(cat ${selectorLog} | grep NO_REGION | sed 's/.*Gene.//;s/ .*//'| tr '\n' ' ')" >> ${LOGFILE}
fi

if [[ ! -e gene_set_output_genotype_file.txt ]]; then
    echo `date "+%Y.%b.%d_%H:%M"` "[Error] Gene set ${chunkNo} has failed. No genotype file has been generated." >> ${LOGFILE}
#    exit 1
elif [[ $(cat gene_set_output_genotype_file.txt | wc -l ) -lt 2 ]]; then
    echo `date "+%Y.%b.%d_%H:%M"` "[Error] Gene set ${chunkNo} has failed, genotype file is empty." >> ${LOGFILE}
#    exit 1
elif [[ ! -e gene_set_output_variant_file.txt ]]; then
    echo `date "+%Y.%b.%d_%H:%M"` "[Error] Gene set ${chunkNo} has failed, SNP file was not generated." >> ${LOGFILE}
#    exit 1
elif [[ $( cat gene_set_output_variant_file.txt | wc -l ) -lt 1 ]]; then
    echo `date "+%Y.%b.%d_%H:%M"` "[Error] Gene set ${chunkNo} has failed, SNP file is empty." >> ${LOGFILE}
#    exit 1
fi

echo `date "+%Y.%b.%d_%H:%M"` "[Info] VARIANT SELECTION DONE" >> ${LOGFILE}
