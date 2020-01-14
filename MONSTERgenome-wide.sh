#!/bin/bash

# if pattern corresponds to a filename (doesn't contain %), check if the file exists
# otherwise check if files for each chromosome (1-22) exist
function testVCFs {
    pattern=$1
    c=0
    if [[ $pattern =~ % ]];then
	for i in $(seq 1 22);do
	    if [[ ! -e $(echo $pattern|sed "s/\%/$i/") ]];then
		echo `date "+%Y.%b.%d_%H:%M"` "[Warning] No VCF for chromosome $i"
	    else
		c=$((c+1))
	    fi
	done
	if [[ $c -gt 0 ]];then
	    return 0
	else
	    return 1
	fi
    else
	if [[ -e $pattern ]];then
	    return 0
	else
	    return 1
	fi
    fi
}

version="v12 Last modified: 2019.Dec.21"
today=$(date "+%Y.%b.%d-%H_%M")

# The variant selector script, that generates snp and genotype input for MONSTER:
regionSelector="Burden_testing.pl"

# Folder with the variant selector script:
scriptDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
MONSTER=$(which MONSTER)
missing_cutoff=1 # Missingness threshold, individuals having missingness higher than this threshold will be excluded.
imputation_method='-A' # The default imputation method is BLUP, slowest, but the most accurate. For other options, see MONSTER documentation.
configFile=""

chunksTotal=1
chunkNo=""
MAF=0.05 # By default this is the upper minor allele frequency.

# --- print out help message and exit:
function display_help() {
    echo "$1"
    echo ""
    echo "Genome-wide Monster wrapper"
    echo "version: ${version}"
    echo ""
    echo "This script was written to run MONSTER genome wide. This script takes a series of arguments
          based on which it calls downstream helper scripts, and generates specific directory for the output files.
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
    echo "     -o  - include variants with severe consequences only (more severe than missense)."
    echo "     -f  - include only HC and LC loftee variants."
    echo "     -j  - include only HC loftee variants."
    echo "     -C  - config file for Burden_testing.pl."
    echo ""
    echo "Parameters to set up scores for variants:"
    echo "     -s  - turn weights on. Arguments: CADD, Eigen, EigenPC, EigenPhred, EigenPCPhred, Mixed"
    echo "     -t  - the value with which the scores will be shifted (default value: if Eigen score weighting specified: 1, otherwise: 0)"
    echo "     -k  - below the specified cutoff value, the variants will be excluded (default: 0)"
    echo ""
    echo "Gene list and chunking:"
    echo "     -L  - file with gene IDs (required, no default)."
    echo "     -d  - total number of chunks (default: 1)."
    echo "     -c  - chunk number (default: 1)."
    echo ""
    echo "General options:"
    echo "     -w  - working directory where the output files of \"prepare-regions\" are (required, no default)"
    echo ""
    echo "Monster parameters:"
    echo "     -p  - phenotype name (required, no default)"
    echo "     -P  - phenotype file (required, no default)"
    echo "     -K  - kinship matrix (required, no default)"
    echo "     -V  - VCF file(s) (required, no default; use % character for chromosome name eg 'chr%.vcf.gz')"
    echo ""
    echo "Other options:"
    echo "     -h  - print this message and exit"
    echo ""
    echo ""

    exit 1
}

# --- Capture command line options --------------------------------------------

if [ $# == 0 ]; then display_help; fi

zipout="yes" # by default gzip output results
OPTIND=1
while getopts ":hL:c:d:p:P:K:V:bg:m:s:l:e:x:k:t:ofw:jC:z" optname; do
    case "$optname" in
      # Gene list related parameters:
        "L") geneListFile=${OPTARG} ;;
        "c") chunkNo=${OPTARG} ;;
        "d") chunksTotal=${OPTARG} ;;

      # MONSTER input files:
        "p" ) phenotype=${OPTARG} ;;
        "P" ) phenotypeFile=${OPTARG} ;;
        "K" ) kinshipFile=${OPTARG} ;;
        "V" ) vcfFile=${OPTARG} ;;

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
        "z") zipout="no" ;;
        "f") loftee=1 ;;
        "j") lofteeHC=1 ;;
        "C") configFile=${OPTARG} ;;

      # Other parameters:
        "w") rootDir=${OPTARG} ;;
        "h") display_help ;;
        "?") display_help "[Error] Unknown option $OPTARG" ;;
        ":") display_help "[Error] No argument value for option $OPTARG" ;;
        *) display_help "[Error] Unknown error while processing options" ;;
    esac
done

if [[ -z "${rootDir}" ]]; then
    echo `date "+%Y.%b.%d_%H:%M"` "[Error] Working directory not specified";
    exit;
fi

if [[ ! -d "${rootDir}" ]]; then
    echo `date "+%Y.%b.%d_%H:%M"` "[Error] Working directory does not exist: $rootDir";
    exit;
fi

chunk_warning=""
if [[ ! -z ${SLURM_ARRAY_TASK_ID} ]];then
    if [[ ! -z ${chunkNo} ]];then
	chunk_warning="WARNING: both SLURM_ARRAY_TASK_ID and chunkNo ( -c ) are defined; using chunkNo"
    else
	chunkNo=${SLURM_ARRAY_TASK_ID}
    fi
else
    if [[ -z ${chunkNo} ]];then
	chunkNo=1 # default
    fi
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

if [[ -z "${phenotypeFile}" ]]; then
    echo `date "+%Y.%b.%d_%H:%M"` "[Error] Phenotype file has to be specified!"
    exit 1
elif [[ ! -e "${phenotypeFile}" ]]; then
    echo `date "+%Y.%b.%d_%H:%M"` "[Error] Phenotype file could not be opened: $phenotypeFile"
    exit 1
fi

if [[ -z "${kinshipFile}" ]]; then
    echo `date "+%Y.%b.%d_%H:%M"` "[Error] Kinship file has to be specified!"
    exit 1
elif [[ ! -e "${kinshipFile}" ]]; then
    echo `date "+%Y.%b.%d_%H:%M"` "[Error] Kiship file could not be opened: $kinshipFile"
    exit;
fi

if [[ -z "${vcfFile}" ]]; then
    "[Error] VCF file(s) not specified!"
    exit 1
else
    testVCFs ${vcfFile}
    if [[ $? -eq 1 ]];then
	echo `date "+%Y.%b.%d_%H:%M"` "[Error] VCF file could not be opened: $vcfFile"
	exit 1
    fi

    commandOptions="${commandOptions} --vcf ${vcfFile} "
fi

if [[ ! -e "${geneListFile}" ]]; then
    echo `date "+%Y.%b.%d_%H:%M"` "[Error] Gene list file could not be opened: $geneListFile"
    exit 1
fi

commandOptions="${commandOptions} --working-dir ${rootDir} "

# GENCODE -expecting a list of feature names separated by a comma.
if [[ ! -z "${gencode}" ]]; then
    folder="${gencode}"
    commandOptions="${commandOptions} --GENCODE ${gencode}"
fi

# GTEx - expecting a list of feature names separeted by comma.
if [[ ! -z "$gtex" ]]; then
    folder="${folder}.GTEX.${gtex}"
    commandOptions="${commandOptions} --GTEx ${gtex}"
fi

# Overlap - expecting a list of features separeated by comma.
if [[ ! -z "${overlap}" ]]; then
    folder="${folder}.Overlap.${overlap}"
    commandOptions="${commandOptions} --overlap ${overlap}"
fi

# MAF - expecting a float between 0 and 0.5
if [[ ! -z "$MAF" ]]; then
    folder="${folder}.MAF.${MAF}"
    commandOptions="${commandOptions} --maf ${MAF}"
fi

# IF loftee is set, only those variants are included in the test that were predicted to be LC or HC lof variants by loftee:
if [[ ! -z "$loftee" ]]; then
    folder="${folder}.loftee"
    commandOptions="${commandOptions} --loftee "
fi

# IF lofteeHC is set, only those variants are included in the test that were predicted to be HC lof variants by loftee:
if [[ ! -z "$lofteeHC" ]]; then
    folder="${folder}.lofteeHC"
    commandOptions="${commandOptions} --lofteeHC "
fi

# If lof is set, only variants with severe consequences will be selected.
if [[ ! -z "$lof" ]]; then
    folder="${folder}.severe"
    commandOptions="${commandOptions} --lof "
fi

warning1=""
warning2=""
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
        MIXED )        score="Mixed";;
        * )            score="noweight";;
    esac
else
    warning1="[Warning] Submitted score name is not recognized! Accepted scores: CADD, Eigen, EigenPC, EigenPhred, EigenPCPhred or Mixed."
    warning2="[Warning] No scoring will be applied."
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

# Exons might be extended with a given number of bps:
if [[ ! -z "${xtend}" ]]; then
    folder="${folder}.Xtend.${xtend}"
    commandOptions="${commandOptions} --extend ${xtend}"
fi

# Setting score cutoff, below which the variant will be removed from the test:
if [[ ! -z "${cutoff}" ]]; then
    folder="${folder}.cutoff.${cutoff}"
    commandOptions="${commandOptions} --cutoff ${cutoff}"
fi

# TODO: scoreshift, general vs Eigen ?
# Setting score shift, a number that will be added to every scores (MONSTER does not accept scores <= 0!!):
if [[ ! -z "${scoreshift}" ]]; then
    folder="${folder}.shift.${scoreshift}"
    commandOptions="${commandOptions} --shift ${scoreshift}"
fi

# Checking if phenotype is provided and if that phenotype file:
if [[ -z "${phenotype}" ]]; then
    echo `date "+%Y.%b.%d_%H:%M"` "[Error] Phenotype was not set! Exiting."
    exit 1
fi

# Updating working dir, and creating folder:
folder=$( echo $folder | perl -lane '$_ =~ s/^\.//;$_ =~ s/,/_/g; print $_;')
workingDir="${rootDir}/${folder}/Pheno.${phenotype}"
outDir=${workingDir}/gene_set.${chunkNo}
mkdir -p ${outDir}
if [[ ! -d ${outDir} ]]; then
    echo `date "+%Y.%b.%d_%H:%M"` "[Error] Chunk directory (${outDir}) could not be created. Exiting."
fi

LOGFILE=${outDir}/"MONSTER-"${today}."chunk_$chunkNo".log

if [[ ! -z ${warning1} ]];then
    echo ${warning1} >> ${LOGFILE}
fi

if [[ ! -z ${warning2} ]];then
    echo ${warning2} >> ${LOGFILE}
fi

if [[ ! -z ${chunk_warning} ]];then
    echo ${chunk_warning} >> ${LOGFILE}
fi

# Creating gene set:
totalGenes=$(cat ${geneListFile} | wc -l)
rem=$(( totalGenes % chunksTotal ))
chunkSize=$(( totalGenes / chunksTotal ))
if [[ $rem -ne 0 ]];then
    chunkSize=$(( chunkSize + 1 ))
fi

# --- Reporting parameters ------------------------------------------------------
echo `date "+%Y.%b.%d_%H:%M"` "##"  >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"` "## Genome-wide Monster wrapper version ${version}" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"` "## Date: ${today}" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"` "##" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"` "" >> ${LOGFILE}

echo `date "+%Y.%b.%d_%H:%M"` "[Info] General options:" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"`  "Variant selector: ${regionSelector}" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"`  "Script dir: ${scriptDir}" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"`  "Working directory: ${workingDir}/gene_set.${chunkNo}" >> ${LOGFILE}
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
echo `date "+%Y.%b.%d_%H:%M"`  "Output folder: ${workingDir}" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"` "" >> ${LOGFILE}

echo `date "+%Y.%b.%d_%H:%M"` "[Info] command line options for burden get region: ${commandOptions}" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"` "" >> ${LOGFILE}

echo `date "+%Y.%b.%d_%H:%M"` "[Info] MONSTER options:" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"`  "MONSTER executable: ${MONSTER}" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"`  "Missingness: ${missing_cutoff}" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"`  "Imputation method: ${imputation_method:-BLUP}" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"`  "Kinship matrix: ${kinshipFile}" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"`  "Phenotype file: ${phenotypeFile}" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"`  "Phenotype: ${phenotype}" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"` "" >> ${LOGFILE}

# --- Main loop executed for all genes --------------------------------------------

awk -v cn="${chunkNo}" -v cs="${chunkSize}" 'NR > (cn-1)*cs && NR <= cn*cs' ${geneListFile} > ${outDir}/input_gene.list
n=$( cat ${outDir}/input_gene.list | wc -l)
if [[ $n -eq 0 ]];then
    echo "Chunk ${chunkNo} is empty; EXIT" >> ${LOGFILE}
    exit 0
fi

selectorLog=${outDir}/chunk_${chunkNo}.output.log
echo `date "+%Y.%b.%d_%H:%M"` "Calling ${scriptDir}/${regionSelector}  --input ${outDir}/input_gene.list --output gene_set_output --output-dir ${outDir} ${commandOptions} --verbose > ${selectorLog} 2 > ${selectorLog}"  >> ${LOGFILE}
${scriptDir}/${regionSelector} --input ${outDir}/input_gene.list --output gene_set_output --output-dir ${outDir} ${commandOptions} --verbose > ${selectorLog} 2 > ${selectorLog}

# We are expecting to get 2 files: gene_set_output_genotype_file.txt & gene_set_output_SNPinfo_file.txt
echo `date "+%Y.%b.%d_%H:%M"` "[Info] Checking output..." >> ${LOGFILE}

cd ${outDir}

# We have to check if both files are generated AND they have enough lines.
gene_notenough=$(cat ${selectorLog} | grep -c NOT_ENOUGH_VAR)
gene_toomany=$(cat ${selectorLog} | grep -c TOO_MANY_VAR)
gene_noremain=$(cat ${selectorLog} | grep -c NO_VAR_REMAIN)
gene_absent=$(cat ${selectorLog} | grep -c NO_GENE)
region_absent=$(cat ${selectorLog} | grep -c NO_REGION)

echo `date "+%Y.%b.%d_%H:%M"` -e "[Warning] ERROR REPORTING FROM VARIANT SELECTOR" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"` -e "[Warning] =====================================" >> ${LOGFILE}

if [[ "$gene_notenough" -ne 0 ]]; then
        echo `date "+%Y.%b.%d_%H:%M"` -e "[Warning] Not enough variants [NOT_ENOUGH_VAR]:\t $(cat ${selectorLog} | grep NOT_ENOUGH_VAR | sed 's/.*Gene.//;s/ .*//' | tr '\n' ' ')" >> ${LOGFILE}
fi
if [[ "$gene_toomany" -ne 0 ]]; then
        echo `date "+%Y.%b.%d_%H:%M"` -e "[Warning] Too many variants [TOO_MANY_VAR]:\t $(cat ${selectorLog} | grep TOO_MANY_VAR | sed 's/.*Gene.//;s/ .*//'| tr '\n' ' ')" >> ${LOGFILE}
fi
if [[ "$gene_noremain" -ne 0 ]]; then
        echo `date "+%Y.%b.%d_%H:%M"` -e "[Warning] All scoring failed [NO_VAR_REMAIN]:\t $(cat ${selectorLog} | grep NO_VAR_REMAIN | sed 's/.*Gene.//;s/ .*//'| tr '\n' ' ')" >> ${LOGFILE}
fi
if [[ "$gene_absent" -ne 0 ]]; then
        echo `date "+%Y.%b.%d_%H:%M"` -e "[Warning] Gene name unknown [NO_GENE]:\t $(cat ${selectorLog} | grep NO_GENE | sed 's/.*Gene.//;s/ .*//'| tr '\n' ' ')" >> ${LOGFILE}
fi
if [[ "$region_absent" -ne 0 ]]; then
        echo `date "+%Y.%b.%d_%H:%M"` -e "[Warning] No region in gene [NO_REGION]:\t $(cat ${selectorLog} | grep NO_REGION | sed 's/.*Gene.//;s/ .*//'| tr '\n' ' ')" >> ${LOGFILE}
fi

if [[ ! -e gene_set_output_genotype_file.txt ]]; then
    echo `date "+%Y.%b.%d_%H:%M"` "[Error] Gene set ${chunkNo} has failed. No genotype file has been generated. Exiting." >> ${LOGFILE}
    exit 1
elif [[ $(cat gene_set_output_genotype_file.txt | wc -l ) -lt 2 ]]; then
    echo `date "+%Y.%b.%d_%H:%M"` "[Error] Gene set ${chunkNo} has failed, genotype file is empty. Exiting." >> ${LOGFILE}
    exit 1
elif [[ ! -e gene_set_output_variant_file.txt ]]; then
    echo `date "+%Y.%b.%d_%H:%M"` "[Error] Gene set ${chunkNo} has failed, SNP file was not generated. Exiting." >> ${LOGFILE}
    exit 1
elif [[ $( cat gene_set_output_variant_file.txt | wc -l ) -lt 1 ]]; then
    echo `date "+%Y.%b.%d_%H:%M"` "[Error] Gene set ${chunkNo} has failed, SNP file is empty. Exiting." >> ${LOGFILE}
    exit 1
fi

# At this point the we have to process the above created files to syncronize wi/nfs/team144/ds26/scripts/burden_testing/MONSTERgenome-wide_updated.sh -g exon -x 50 -e promoter,enhancer,TF_bind -l promoter,enhancer,TF_bind -s EigenPCPhred -c 3 -P /lustre/scratch115/projects/t144_helic_15x/analysis/HA/phenotypes/correct_names.andmissing/MANOLIS.HDL.txt -w /lustre/scratch115/realdata/mdt0/projects/t144_helic_15x/analysis/HA/burdentesting/arthur_rerun -V /lustre/scratch115/realdata/mdt0/projects/t144_helic_15x/analysis/HA/release/postrelease_missingnessfilter/chr%.missingfiltered-0.01_consequences.lof.HWE.vcf.gz -K /lustre/scratch115/realdata/mdt0/projects/t144_helic_15x/analysis/HA/relmat/final/burden/matrix.monster.txt -p HDL -d 50th the phenotype file and the kinship matrix.

# ASSUMING PHENOTYPE FILE HAS 2 COLUMNS: ID PHENOTYPE, TAB DELIMITED, NO HEADER
echo `date "+%Y.%b.%d_%H:%M"` "[Info] Extracting phenotype." >> ${LOGFILE}
cat ${phenotypeFile} | awk 'BEGIN{FS="\t";OFS="\t";}{if ($2!="NA") {printf "1\t%s\t0\t0\t0\t%s\n", $1, $2} }' > pheno.txt # subset of samples in pheno.txt (without NA)

# Order the sample IDs in the phenotype file:
echo `date "+%Y.%b.%d_%H:%M"` "[Info] Re-ordering samples in the phenotype file." >> ${LOGFILE}
head -n 1 gene_set_output_genotype_file.txt | tr "\t" "\n" | tail -n +2 |sort| perl -lane 'BEGIN {open $pf, "< pheno.txt";while ($l = <$pf>){chomp $l;@a = split(/\s/, $l);$h{$a[1]} = $l;}}{print $h{$F[0]} if exists $h{$F[0]}}' > pheno.ordered.txt

# Get the list of samples that are not in the pheno file:
export samples=$(grep -v -w -f <(cut -f 2 pheno.ordered.txt) <(head -n 1 gene_set_output_genotype_file.txt | cut -f 2- | tr "\t" "\n"))

# From the genotype file, extract only those samples that are present in the pheno file:
echo `date "+%Y.%b.%d_%H:%M"` "[info] Extracting un-used samples from the genotype file." >> ${LOGFILE}
head -n 1 gene_set_output_genotype_file.txt | tr "\t" "\n" | perl -lane 'BEGIN {foreach $s ( split /\s/, $ENV{"samples"}){$h{$s} = 1;}}{push @a, $. unless exists $h{$F[0]}} END{$s = sprintf("cut -f%s gene_set_output_genotype_file.txt > genotype.filtered.txt", join(",", @a));`$s`}'

# Generate a mapping file that helps to convert IDs to numbers:
echo `date "+%Y.%b.%d_%H:%M"` "[Info] Generate sample mapping file." >> ${LOGFILE}
cut -f 2 pheno.ordered.txt | awk '{printf "s/%s/%s/g\n", $1, NR+2 }' > sample.map.sed

# Generate an inclusion list with the samples to be kept:
cut -f 2 pheno.ordered.txt > samples.to.keep.txt

# Get the kinship matrix:
echo `date "+%Y.%b.%d_%H:%M"` "[Info] Processing kinship file." >> ${LOGFILE}
R --slave -e 'library(data.table); mlong=fread("'$kinshipFile'"); tokeep=fread("samples.to.keep.txt", header=F)$V1; direct=mlong[(mlong$V2 %in% tokeep) & (mlong$V3 %in% tokeep),];  mapping = fread("sample.map.sed", sep="/", header=FALSE);direct$V3 = mapping[match(direct$V3, mapping$V2),]$V3; direct$V2 = mapping[match(direct$V2, mapping$V2),]$V3;write.table(direct, file="filtered_kinship.txt", quote=FALSE, sep=" ", col.names = FALSE, row.names=FALSE)'

# Remap IDs and remove special characters from the snp, phenotype and genotype files:
echo `date "+%Y.%b.%d_%H:%M"` "[Info] Changing IDs and variant names." >> ${LOGFILE}
sed -i -f sample.map.sed pheno.ordered.txt
sed -i -f sample.map.sed genotype.filtered.txt
cat genotype.filtered.txt | perl -lane '$_ =~ s/[^0-9A-Za-z\-\t\._]//gi; print $_'  > genotype.filtered.mod.txt
cat gene_set_output_variant_file.txt | perl -lane '$_ =~ s/[^0-9A-Za-z\-\t\._]//gi; $_ =~ s/Inf/0.0001/g; ;print $_'  > snpfile.mod.txt

# Filter out genes which have only monomorphic variants, as it might cause a crash:
echo `date "+%Y.%b.%d_%H:%M"` "[Info] Looking for monomorphic variants..." >> ${LOGFILE}
#tail -n+2 genotype.filtered.mod.txt | while read snp genotype ; do if [[ -z $( echo $genotype | awk '$0 ~ 1' ) ]]; then echo $snp; fi; done | awk '{printf "s/%s//g\n", $1}' > mono_remove.sed
tail -n+2 genotype.filtered.mod.txt | perl -lne '@f=split(/\s+/);$\="\n";$s=shift(@f);foreach (@f){$H{$_}=1;}if (scalar(keys(%H))==1){print $s;}' | awk '{printf "s/%s//g\n", $1}' > mono_remove.sed

echo `date "+%Y.%b.%d_%H:%M"` "[Info] Removing monomorphic variants from the SNPs file." >> ${LOGFILE}
sed -f mono_remove.sed snpfile.mod.txt > snpfile.nomono

echo `date "+%Y.%b.%d_%H:%M"` "[Info] Get genes where only monomorphics remain. Exclude them." >> ${LOGFILE}
grep -v -w -f <(cat snpfile.nomono | awk 'NF  == 2 {print $1; a+=1}END{if(a == 0){print "noting to remove"}}') snpfile.mod.txt > snpfile.mod.nomono.txt

# Calling MONSTER
echo `date "+%Y.%b.%d_%H:%M"` "[Info] MONSTER call: MONSTER -k filtered_kinship.txt -p pheno.ordered.txt -m 1 -g genotype.filtered.mod.txt  -s snpfile.mod.nomono.txt ${imputation_method}" >> ${LOGFILE}
while true; do
    MONSTER -k filtered_kinship.txt -p pheno.ordered.txt -m 1 -g genotype.filtered.mod.txt  -s snpfile.mod.nomono.txt ${imputation_method}

    # We break the loop if the run was successful.
    if [[ $? -eq 0 ]]; then break; fi
    
    # Test if we've analyzed all genes
    if [[ $(awk 'NF == 5' MONSTER.out | cut -f1 | tail -n+2 | wc -l ) == $(cut -f1 snpfile.mod.txt | sort| uniq -u | wc -l) ]]; then break; fi

    if [[ ! -e MONSTER.out ]]; then
        echo `date "+%Y.%b.%d_%H:%M"` "[Error] MONSTER failed before creating the output file. Cannot be resolved. Exiting" >> ${LOGFILE}
        break;
    elif [[ $( cat MONSTER.out | wc -l) -eq 1 ]]; then
        firstGene=$(cut -f1 snpfile.mod.nomono.txt | head -n1)
        echo `date "+%Y.%b.%d_%H:%M"` "[Warning] It seems that the first gene (${firstGene}) has failed. Re-running MONSTER." >> ${LOGFILE}
        grep -vw $(firstGene) snpfile.mod.nomono.txt | sponge snpfile.mod.nomono.txt
    else
        lastGene=$(awk 'NF == 5' MONSTER.out | tail -n1 | cut -f1 )
        failedGene=$(grep -A1 -w ${lastGene} snpfile.mod.nomono.txt | cut -f1 | tail -n1 )
        echo `date "+%Y.%b.%d_%H:%M"` "[Warning] Monster has failed after ${lastGene}, next gene (${failedGene}) is removed and re-run." >> ${LOGFILE}
        grep -vw ${failedGene} snpfile.mod.nomono.txt | sponge snpfile.mod.nomono.txt
    fi
done

# Once MONSTER is finished, we remove the un-used temporary files:
rm genotype.filtered.txt

# Moving MONSTER.out to the root directory:
if [[ -e MONSTER.out ]]; then
    cp MONSTER.out ../MONSTER.${phenotype}.${chunkNo}.out
else
    echo `date "+%Y.%b.%d_%H:%M"` "[Error] MONSTER.out file was not found. Something went wrong."  >> ${LOGFILE}
fi

cp ${selectorLog} ..

if [[ ${zipout} = "yes" ]];then
    # Compress folder:
    echo `date "+%Y.%b.%d_%H:%M"` "[Info] Compressing and removing files." >> ${LOGFILE}
    tar -zcvf gene_set.${chunkNo}.tar.gz *
    mv gene_set.${chunkNo}.tar.gz ..
    cd .. && rm -rf gene_set.${chunkNo}
fi
echo `date "+%Y.%b.%d_%H:%M"` "DONE"  >> ${LOGFILE}
