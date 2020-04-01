#!/bin/bash


version="v12 Last modified: 2020.Mar.30"
today=$(date "+%Y.%b.%d")

# The variant selector script:
regionSelector="Burden_testing.pl"

# Folder with the variant selector script:
#scriptDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
#missing_cutoff=1 # Missingness threshold, individuals having missingness higher than this threshold will be excluded.
configFile=""

chunksTotal=1
chunkNo=""
#MAF=1.0

# --- print out help message and exit:
function display_help() {
    echo "$1"
    echo ""
    echo "Variant selector"
    echo "version: ${version}"
    echo ""
    echo "Usage: $0 <parameters>"
    echo ""
    echo "Variant filter options:"
    echo "     -g  - list of gencode features."
    echo "     -e  - list of linked GTEx featuress"
    echo "     -l  - list of linked overlapping features."
#    echo "     -m  - upper maf thresholds"
    echo "     -x  - extend genomic regions (bp)."
    echo "     -o  - include variants with severe consequences only (more severe than missense)."
#    echo "     -f  - include only HC and LC loftee variants."
#    echo "     -j  - include only HC loftee variants."
    echo "     -C  - config file"
    echo "     -i  - input variant list, tab-separated, bgzipped and indexed (required, no default, each line has to have 5 fields: chr,pos,ID,ref,alt)"
    echo ""
    echo "Parameters to set up scores for variants:"
    echo "     -s  - turn weights on. Arguments: CADD, EigenPhred"
    echo "     -t  - the value with which the scores will be shifted (default value: if Eigen score weighting specified: 1, otherwise: 0)"
    echo "     -k  - below the specified cutoff value, the variants will be excluded (default: 0)"
    echo ""
    echo "Gene list and chunking:"
    echo "     -L  - file with gene IDs (if not specified all genes will be analyzed)."
    echo "     -d  - total number of chunks (default: 1)."
    echo "     -c  - chunk number (default: 1). Takes precedence over SLURM_ARRAY_TASK_ID etc."
    echo ""
    echo "General options:"
    echo "     -w  - output directory (required, no default)"
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

#while getopts ":hL:c:d:bg:m:s:l:e:x:k:t:ofw:jC:V:" optname; do
while getopts ":hL:c:d:bg:s:l:e:x:k:t:ow:C:i:" optname; do
    case "$optname" in
      # Gene list related parameters:
        "L") geneListFile=${OPTARG} ;;
        "c") chunkNo=${OPTARG} ;;
        "d") chunksTotal=${OPTARG} ;;

      # variant filter parameters:
        "g") gencode=${OPTARG} ;;
#        "m") MAF=${OPTARG} ;;
        "s") score=${OPTARG} ;;
        "l") overlap=${OPTARG} ;;
        "e") gtex=${OPTARG} ;;
        "x") xtend=${OPTARG} ;;
        "k") cutoff=${OPTARG} ;;
        "t") scoreshift=${OPTARG} ;;
        "o") lof=1 ;;
#        "f") loftee=1 ;;
#        "j") lofteeHC=1 ;;
        "C") configFile=${OPTARG} ;;
        "i") inputFile=${OPTARG} ;;

      # Other parameters:
        "w") outputDir=${OPTARG} ;;
        "h") display_help ;;
        "?") display_help "[Error] Unknown option $OPTARG" ;;
        ":") display_help "[Error] No argument value for option $OPTARG" ;;
        *) display_help "[Error] Unknown error while processing options" ;;
    esac
done


if [[ -z "${inputFile}" ]]; then
    "[Error] input file not specified!"
    exit 1
fi

if [[ -z "${outputDir}" ]]; then
    echo `date "+%Y.%b.%d_%H:%M"` "[Error] Output directory not specified";
    exit 1
fi

# remove trailing slash
outputDir=${outputDir%/}

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

#--- checking input files - if any of the tests fails, the script exits.---------

if [[ -z "${configFile}" ]]; then
    echo `date "+%Y.%b.%d_%H:%M"` "[Error] Config file was not specified";
    exit;
fi

if [[ ! -e "${configFile}" ]]; then
    echo `date "+%Y.%b.%d_%H:%M"` "[Error] Config file does not exist: $configFile";
    exit;
fi

commandOptions=" --config ${configFile} "
# -----------------------------------------------------------------------------------------------------------------------------

outprefix="group_file"
# GENCODE -expecting a list of feature names separated by a comma.
if [[ ! -z "${gencode}" ]]; then
    commandOptions="${commandOptions} --GENCODE ${gencode}"
    str=$( echo "${gencode}" | perl -lane '$_ =~ s/^\.//;$_ =~ s/,/_/g; print $_;')
    outprefix=${outprefix}"_GENCODE_"${str}
fi

# GTEx - expecting a list of feature names separeted by comma.
if [[ ! -z "$gtex" ]]; then
    commandOptions="${commandOptions} --GTEx ${gtex}"
    str=$( echo "${gtex}" | perl -lane '$_ =~ s/^\.//;$_ =~ s/,/_/g; print $_;')
    outprefix=${outprefix}"_GTEx_"${str}
fi

# Overlap - expecting a list of features separeated by comma.
if [[ ! -z "${overlap}" ]]; then
    commandOptions="${commandOptions} --overlap ${overlap}"
    str=$( echo "${overlap}" | perl -lane '$_ =~ s/^\.//;$_ =~ s/,/_/g; print $_;')
    outprefix=${outprefix}"_overlap_"${str}
fi

# If lof is set, only variants with severe consequences will be selected.
if [[ ! -z "$lof" ]]; then
    commandOptions="${commandOptions} --lof "
    outprefix=${outprefix}"_severe"
fi

warning1=""
warning2=""
score_tmp=""
# Score - If score is not given we apply no score. Otherwise we test the submitted value:
# Accepted scores:
if [[ ! -z "${score}" ]]; then
    score="${score^^}"
    case "${score}" in
        EIGENPHRED )   score="EigenPhred";;
        CADD )         score="CADD";;
        * )            score_tmp="noweight";;
    esac
else
    score="noweight"
fi

if [[ ! -z ${score_tmp} ]];then
    warning1="[Warning] Submitted score name ($score) is not recognized! Accepted scores: CADD, EigenPhred."
    warning2="[Warning] No scoring will be applied."
    score="noweight"
fi

# Only adding score to command line if score is requested:
if [[ "${score}" != "noweight" ]]; then
    commandOptions="${commandOptions} --score ${score}";
fi
outprefix=${outprefix}"_score_"${score}

# If Eigen score is applied, we shift the scores by 1, if no other value is specified:
#if [[ ("${score}" == "Eigen") && (-z "${scoreshift}" ) ]]; then scoreshift=1; fi

# Exons might be extended with a given number of bps:
if [[ ! -z "${xtend}" ]]; then
    commandOptions="${commandOptions} --extend ${xtend}"
    outprefix=${outprefix}"_extend_"${xtend}
fi

# Setting score cutoff, below which the variant will be removed from the test:
if [[ ! -z "${cutoff}" ]]; then
    commandOptions="${commandOptions} --cutoff ${cutoff}"
fi

# Setting score shift, a number that will be added to every scores (MONSTER does not accept scores <= 0!!):
if [[ ! -z "${scoreshift}" ]]; then
    commandOptions="${commandOptions} --shift ${scoreshift}"
fi

outFile="group_file_gene_set."${chunkNo}
outputDir2=${outputDir}"/"${outprefix}
commandOptions="${commandOptions} --smmat ${inputFile} --output-dir ${outputDir2} --output ${outFile}"

LOGFILE=${outputDir2}/"make_group_file_gene_set.${chunkNo}.log"

#--------------------------------------------------------------------------------------------
no_list_warning=""
if [[ -z ${geneListFile} ]];then
    gencode_file=$(grep "^gencode_file" ${configFile} | cut -f 2 -d '=')
    zcat ${gencode_file} | cut -f 4 > "${outputDir2}/temp_gene_list.txt"
    geneListFile="${outputDir2}/temp_gene_list.txt"
    no_list_warning="[Warning] No gene list specified; using all genes from $gencode_file"
fi
#--------------------------------------------------------------------------------------------

if [[ ! -e "${geneListFile}" ]]; then
    echo `date "+%Y.%b.%d_%H:%M"` "[Error] Gene list file could not be opened: $geneListFile"
    exit 1
fi

# -----------------------------------------------------------------------------------------------------------------------------
# creating gene list

totalGenes=$(cat ${geneListFile} | wc -l)
if [ $totalGenes -lt $chunksTotal ];then
    echo `date "+%Y.%b.%d_%H:%M"` "[Warning] Number of chunks ($chunksTotal) is larger than number of genes in the gene list ($totalGenes) "  >> ${LOGFILE}
    echo `date "+%Y.%b.%d_%H:%M"` "[Warning] Analyzing all genes in one chunk"  >> ${LOGFILE}
    chunkNo=1
    cat ${geneListFile} > ${outputDir2}/input_gene.list    
else
    rem=$(( totalGenes % chunksTotal ))
    chunkSize=$(( totalGenes / chunksTotal ))
    lastChunkSize=$(( chunkSize + rem ))
    if [[ ${chunkNo} -eq ${chunksTotal} ]];then
	tail -n ${lastChunkSize} ${geneListFile} > ${outputDir2}/input_gene.list
    else
	awk -v cn="${chunkNo}" -v cs="${chunkSize}" 'NR > (cn-1)*cs && NR <= cn*cs' ${geneListFile} > ${outputDir2}/input_gene.list
    fi

    n=$( cat ${outputDir2}/input_gene.list | wc -l)
    if [[ $n -eq 0 ]];then
	echo "Chunk ${chunkNo} is empty; EXIT" >> ${LOGFILE}
	exit 0
    fi
fi

# -----------------------------------------------------------------------------------------------------------------------------
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

# --- Reporting parameters ------------------------------------------------------
echo `date "+%Y.%b.%d_%H:%M"` "##"  >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"` "## Version ${version}" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"` "## Date: ${today}" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"` "##" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"` "" >> ${LOGFILE}

echo `date "+%Y.%b.%d_%H:%M"` "[Info] General options:" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"`  "Variant selector: ${regionSelector}" >> ${LOGFILE}
#echo `date "+%Y.%b.%d_%H:%M"`  "Script dir: ${scriptDir}" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"`  "Output directory: ${outputDir}" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"` "" >> ${LOGFILE}

echo `date "+%Y.%b.%d_%H:%M"` "[Info] Gene list options:" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"`  "Gene list file: ${geneListFile}" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"`  "Number of chunks the gene list is split into: ${chunksTotal}" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"`  "Current chunk: ${chunkNo}" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"`  "Number of genes in one chunk: ${chunkSize}" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"` "" >> ${LOGFILE}

echo `date "+%Y.%b.%d_%H:%M"` "[Info] Variant filtering options:" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"`  "input file: ${inputFile}" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"`  "GENCODE feaures: ${gencode:--}" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"`  "GTEx feaures: ${gtex:--}" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"`  "Overlapping reg.features: ${overlap:-NA}" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"`  "Features are extended by ${xtend:-0}bp" >> ${LOGFILE}
#echo `date "+%Y.%b.%d_%H:%M"`  "Upper minor allele frequency: ${MAF:-1}" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"` "" >> ${LOGFILE}

echo `date "+%Y.%b.%d_%H:%M"` "[Info] Weighting options:" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"`  "Weighting: ${score}" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"`  "Score cutoff: ${cutoff:-0}" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"`  "Scores shifted by: ${scoreshift:-0}" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"` "" >> ${LOGFILE}

echo `date "+%Y.%b.%d_%H:%M"` "[Info] command line options for the selector script: ${commandOptions}" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"` "" >> ${LOGFILE}

selectorLog=${outputDir2}/"selector_chunk_"${chunkNo}.log

echo `date "+%Y.%b.%d_%H:%M"` "Calling ${regionSelector} --input ${outputDir2}/input_gene.list ${commandOptions} --verbose > ${selectorLog} 2 > ${selectorLog}"  >> ${LOGFILE}
${regionSelector} --input ${outputDir2}/input_gene.list ${commandOptions} --verbose > ${selectorLog} 2 > ${selectorLog}

echo `date "+%Y.%b.%d_%H:%M"` "[Info] Checking output..." >> ${LOGFILE}

cd ${outputDir2}

# We have to check if both files are generated AND they have enough lines.
gene_notenough=$(cat ${selectorLog} | grep -c NOT_ENOUGH_VAR)
gene_toomany=$(cat ${selectorLog} | grep -c TOO_MANY_VAR)
gene_noremain=$(cat ${selectorLog} | grep -c NO_VAR_REMAIN)
gene_absent=$(cat ${selectorLog} | grep -c NO_GENE)
region_absent=$(cat ${selectorLog} | grep -c NO_REGION)

echo `date "+%Y.%b.%d_%H:%M"` -e "[Warning] ERROR REPORTING FROM VARIANT SELECTOR" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"` -e "[Warning] =====================================" >> ${LOGFILE}

if [[ "$gene_notenough" -ne 0 ]]; then
        echo `date "+%Y.%b.%d_%H:%M"` -e "[Warning] Not enough variants [NOT_ENOUGH_VAR]: $(cat ${selectorLog} | grep NOT_ENOUGH_VAR | sed 's/.*Gene.//;s/ .*//' | tr '\n' ' ')" >> ${LOGFILE}
fi
if [[ "$gene_toomany" -ne 0 ]]; then
        echo `date "+%Y.%b.%d_%H:%M"` -e "[Warning] Too many variants [TOO_MANY_VAR]: $(cat ${selectorLog} | grep TOO_MANY_VAR | sed 's/.*Gene.//;s/ .*//'| tr '\n' ' ')" >> ${LOGFILE}
fi
if [[ "$gene_noremain" -ne 0 ]]; then
        echo `date "+%Y.%b.%d_%H:%M"` -e "[Warning] All scoring failed [NO_VAR_REMAIN]: $(cat ${selectorLog} | grep NO_VAR_REMAIN | sed 's/.*Gene.//;s/ .*//'| tr '\n' ' ')" >> ${LOGFILE}
fi
if [[ "$gene_absent" -ne 0 ]]; then
        echo `date "+%Y.%b.%d_%H:%M"` -e "[Warning] Gene name unknown [NO_GENE]: $(cat ${selectorLog} | grep NO_GENE | sed 's/.*Gene.//;s/ .*//'| tr '\n' ' ')" >> ${LOGFILE}
fi
if [[ "$region_absent" -ne 0 ]]; then
        echo `date "+%Y.%b.%d_%H:%M"` -e "[Warning] No region in gene [NO_REGION]: $(cat ${selectorLog} | grep NO_REGION | sed 's/.*Gene.//;s/ .*//'| tr '\n' ' ')" >> ${LOGFILE}
fi

if [[ ! -e ${outFile} ]]; then
    echo `date "+%Y.%b.%d_%H:%M"` "[Error] Gene set ${chunkNo} has failed. No group file has been generated. Exiting." >> ${LOGFILE}
    exit 1
elif [[ $(cat ${outFile} | wc -l ) -lt 1 ]]; then
    echo `date "+%Y.%b.%d_%H:%M"` "[Error] Gene set ${chunkNo} has failed, group file is empty. Exiting." >> ${LOGFILE}
    exit 1
fi
