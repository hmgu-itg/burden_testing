#!/usr/local/bin/bash

## This script was written to automate conditional testing of potential associations 
## found by MONSTER using MONSTER_wrapper.sh

## This script takes variants (rsID or SNPID (GRCh38), or a file with variants), and a hit gene folder and
## extracts that variant from the 15x vcf file adds as a covariate to the phenotype
## file and re-run the MONSTER annotation.

## There's an option to provide variants by specify only a gene around which the variants
## with phenotype annotations will be extracted from Ensmbl.
## The window around the gene name can also be set by specifying with -w switch.


## v.1.0 Last modified: 2017.04.24

# Step -1: Preparation, setting up environment:

# Directory of the vcf files:
vcfDir=/lustre/scratch115/realdata/mdt0/projects/t144_helic_15x/analysis/HA/release/postrelease_missingnessfilter/

# Single point dir: 
singlePointDir=/lustre/scratch115/projects/t144_helic_15x/analysis/HA/single_point/output/missing_chunks/

# genome-wide plink file:
plinkFile=/lustre/scratch115/projects/t144_helic_15x/analysis/HA/single_point/input/whole_genome/autosomal

# Setting up working dir:
workingDir=$(pwd)

# Distance around the gene we fetch variants with phenotypes:
export window=500000

# Path to the burden p-values:


# Path to MONSTER executable:
MONSTER=/nfs/team144/software/MONSTER_v1.3/MONSTER

# Folder with all the scripts:
export scriptDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# --- print out help message and exit
function display_help() {
    echo ""
    echo "${1}"
    echo ""
    echo "Usage: $0 <parameters>"
    echo ""
    echo "     -g  - Folder with gene hit, or the compressed tar.gz file. (Obligatory)"
    echo "     -s  - A file with a list of snps or rsIDs to condition on"
    echo "     -x  - Gene name."
    echo "     -d  - check LD: the tested variants will be grouped based on LD (only if -o is defined.)"
    echo "     -o  - output file name (stdout by default)"
    echo "     -w  - Window length to fetch variants with phenotypes (bp)."
    echo "     (gene name or snp list has to be specified!)"
    echo ""
    echo "     -h  - Print out this help message."
    echo ""
    echo "More information: ds26@sanger.ac.uk"
    echo ""
    echo ""

    exit 1
}

function test_variant(){
    local chr="${1}"
    local pos="${2}"
    local variant="${3}"
    local ID MAF consequences

    # Extract variant:
    START_TIME=$SECONDS
    read ID MAF consequences <<< $(tabix ${vcfDir}/chr${chr}.missingfiltered-0.01_consequences.vcf.gz chr${chr}:${pos}-${pos} | cut -f-8 | perl -lane '$_ =~ /AF=(.+?);.+consequence=(\S+)/; print "$F[2] $1 $2"')
    ELAPSED_TIME=$(($SECONDS - $START_TIME))
    
    if [[ ${ELAPSED_TIME} -gt 1 ]]; then echo "[Warning] A single tabix took ${ELAPSED_TIME} sec." >&2; fi # report if the tabix query takes unusually long... 
    
    # Checking if extraction was successful:
    if [[ -z $ID ]]; then
        echo "[Warning] The queried variant ${variant} was not found in the HA15X dataset." >&2;
        return 1;
    else
        echo "[Info] Testing ${variant} (chr${chr}:${pos})..." >&2
    fi
    
    # If it was successful, we have to extract the genotypes:
    paste <( zgrep -m1 "#CHROM" ${vcfDir}/chr11.missingfiltered-0.01_consequences.vcf.gz | cut -f10- | tr "\t" "\n" )\
        <( tabix ${vcfDir}/chr${chr}.missingfiltered-0.01_consequences.vcf.gz chr${chr}:${pos}-${pos} \
            | cut -f10- | tr "\t" "\n" | cut -f1 -d":" | perl -F"/" -lane '$F[0] eq "." ? print "-9" : print $F[0] + $F[1]') \
            | sed -f ${scriptDir}/HELIC.to.Num.sed | tr " " "\t" | cut -f1,3 | perl -lane '
        $h{$F[0]} = $F[1];
        END{
            open($ph, "< pheno.mod.ordered.txt");
            while ( $l = <$ph> ){
                chomp $l;
                @a = split("\t", $l);
                next unless exists $h{$a[1]};
                next if $h{$a[1]} == -9;
                print "$l\t$h{$a[1]}"
            }
        }' > pheno.mod.ordered_conditioned_chr${chr}.${pos}.txt
    
    # Now we have to run MONSTER
    ${MONSTER} -k kinship.mod.filtered.txt \
        -p pheno.mod.ordered_conditioned_chr${chr}.${pos}.txt \
        -m 1 \
        -g genotype.mod.filtered.txt \
        -s snpfile.mod.txt \
        ${imputation_method} > /dev/null
    
    # Checking if the output is exists or not:
    if [[ ! -e MONSTER.out ]]; then
        echo "[Warning] MONSTER has failed for ${variant}. No output." >&2;
        return 1;
    fi        
    
    # re-name output:
    mv MONSTER.out MONSTER_cond_chr${chr}.${pos}.out
    
    # Testing MONSTER output:
    if [[ $( cat MONSTER_cond_chr${chr}.${pos}.out | wc -l  ) -eq 1 ]]; then
        echo "[Warning] ${variant} has failed. No MONSTER output!" >&2;
        return 1;
    fi
    
    # Extract p-values:
    tail -n1 MONSTER_cond_chr${chr}.${pos}.out | awk -v snpID="chr${chr}:${pos}" '{print snpID, $0}'
    
    return 0;
}

function rsID2chrPos(){
    local rsID="${1}"
    wget -q --header='Content-type:application/json' "http://rest.ensembl.org/variation/human/${rsID}?phenotypes=1"  -O - | perl -MJSON \
        -lane 'BEGIN{ ($chr, $pos) = ("NA", "NA")}{
            $h = decode_json($_); $rsID = $h->{"name"}; $var_type = $h->{"var_class"};
            foreach $m (@{$h->{"mappings"}}){
                if (length($m->{"seq_region_name"}) < 3 ){
                    $chr = $m->{"seq_region_name"};
                    $pos = $h->{"var_class"} eq "SNP" ? $m->{"start"} : $m->{"start"} - 1;
                }
            }
            print "$chr $pos $rsID $var_type";
            }'
    return 0
}

function getTrait(){

    if [[ ! -z $( ls *.pvalues.txt ) ]]; then    
        trait=$(tail -n1  *pvalues.txt | cut -f1)
        case "$trait" in
            "adiponectin") trait=ADIPONECTIN;;
            "BGP") trait=OSTEOCALCIN;;
            "Bilirubin") trait=BILIRUBIN;;
            "Fe_iron") trait=FEIRON;;
            "Ferritin") trait=FERRITIN;;
            "gammaGT") trait=GAMMAGT;;
            "leptin") trait=LEPTIN;;
        esac
        echo ${trait}
    else
        echo "[Warning] \${gene}.pvalues.txt does not exits\!" >&2
        echo "NA"
    fi
    return 0
}

function getSinglePoint(){
    trait="${1}"
    snp="${2}"
    
    # Use tabix:
    pvals=$( tabix ${singlePointDir}/${trait}/MANOLIS.${trait}.assoc.missingfilled.txt.bgz "${snp}" | head -n1 | cut -f14 )
    
    if [[ -z "${pvals}" ]]; then
        pvals="NA"
    fi
    
    echo "${pvals}";
    return 0;
}

function parseSNPID(){
    local snpid="${1}"
    local chr pos ID type
    read chr pos <<< $(echo "${snpid}" | perl -lane '$_ =~ /chr(.+)\:(\d+)/i; print "$1 $2"')
    echo ${chr:-"NA"} ${pos:-"NA"} ${snpid} "SNPID"
}

function get_coordinates(){
    local variant="${1}"
    local chr pos ID type
    
    
    # checking if an rsID was submitted:
    if [[ $( echo "${variant}" | grep -i rs ) ]]; then
        read chr pos ID type <<< $( rsID2chrPos "${variant}" )
    elif [[ $( echo "${variant}" | grep -i chr ) ]]; then
        read chr pos ID type <<< $( parseSNPID "${variant}" )
    else
        echo "[Warning] The submitted variant (${variant}) is not recognized. Submit rsID or SNPID." >&2
        read chr pos ID type <<< $( echo "NA NA NA NA" )
    fi

    # Returning variant data:
    echo ${chr:-"NA"} ${pos:-"NA"} ${ID:-"NA"} ${type:-"NA"}
}

function get_phenotypes_for_gene(){
    
    geneName="${1}"
    wget -q --header='Content-type:application/json' "http://rest.ensembl.org/lookup/symbol/homo_sapiens/${geneName}?" -O - \
        | perl -MJSON -lane '$h = decode_json($_); printf "%s:%s-%s\n", $h->{"seq_region_name"}, $h->{"start"}-$ENV{"window"}, $h->{"end"}+$ENV{"window"}' \
        | while read coord ; do  wget -q --header='Content-type:application/json' "http://rest.ensembl.org/phenotype/region/homo_sapiens/${coord}?feature_type=Variation" -O - ; done \
        | perl -MJSON -lane 'foreach $a (@{decode_json($_)}){ next unless $a->{"id"} =~ /^rs/; print $a->{"id"}}' \
        > ${geneName}_phenotypes.rsIDs.lst
    varNo=$( cat  ${geneName}_phenotypes.rsIDs.lst | wc -l )
    echo "[Info] For ${geneName} we have found ${varNo} variants with phenotype annotations within 1Mbp window." >&2
    echo  "$(cd "$(dirname "${geneName}_phenotypes.rsIDs.lst")"; pwd)/$(basename "${geneName}_phenotypes.rsIDs.lst")"
}

# Once the run is ready, we extract all the variants that have a high impact on the burden p-value.
function getLDblocks(){
    outfile="$1"
    export gene="$2"
    
    # Get a list of snps:
    awk '$6 != "NA" && NR != 1 {print $4 }' "${outfile}"| tr ":" " " | sort -k1,1 -k2,2n | tr " " ":" | awk '{print $0 "[b38]"}'  > ${gene}_variants;
    
    # Get variants in LD:
    echo "[Info] Running plink to find correlated variants."
    
    plink --memory --bfile ${plinkFile} --extract ${gene}_variants -r2 square gz yes-really --out ${gene}_R2  > /dev/null
    
    # Checking if plink output file was generated. If not, the log file will be displayed.
    if [[ -e ${gene}_R2.ld.gz ]];then
        echo "[Info] plink finished successfully."
    else
        echo "[Warning] plink has failed. Dumping log file:"
        cat  ${gene}_R2.log
    fi
    
    # Parsing plink file:
    echo "[Info] plink finished, parsing output."
    zcat ${gene}_R2.ld.gz | perl -MData::Dumper -lane 'BEGIN {
        open ($snp, "sort $ENV{gene}_variants |");
        $i = 0;
        while ( $l = <$snp>){
            chomp $l; $h{$i} = $l; $i++;
        }}{
            $seen = "pocok";
            $current = $h{$. - 1};
            $i = $. - 1;
            for ($j = 0; $j < scalar(@F); $j ++){
                $pair = $h{$j};
                if(($F[$j] >= 0.2) && ($i <= $j)){
                    if ( $seen == "pocok"){
                        $seen = $j;
                    }
                    $k{$seen}{$pair}++;
                }
                if(($F[$j] >= 0.2) && ($i > $j)){
                    if ( $seen == "pocok"){
                        $seen = $j;
                    }
                    $k{$seen}{$current}++;
                }        
            }
        }; END {
            $LD_count = 1;
            foreach $b (values %k){
                @vars = keys %{$b};
                if (scalar @vars > 1){
                    foreach $v (@vars){
                        print "$v\t$LD_count";
                    }
                    $LD_count ++;
                }
                else {
                    print  "$vars[0]\t0";
                }
            }
        
        }' > ${gene}_LD_gropus.tsv
    
    # Adding extra column to the tsv file:
    cat "${outfile}" | perl -MData::Dumper -lane 'BEGIN {
            open($ld, " < $ENV{gene}_LD_gropus.tsv");
            while ($l = <$ld>){
                chomp $l;
                ($v, $g) = $l =~ /(.+?)\[b38\]\t(.+)/;
                $h{$v} = $g unless exists $h{$v} or $h{$v} > $g;
            }
        }{
            if ($. == 1){ print "$_\tLD"; next; }
            $LD = exists $h{$F[3]} ? $h{$F[3]} : "NA";
            print "$_\t$LD"
        }' > ${outfile}_ext
}

# If help message is needed:
if [ $# == 0 ]; then display_help; fi

# Default output file:
outputFile=/dev/stdout

# Step 0: Processing and checking input parameters:
OPTIND=1
while getopts ":hg:s:dw:x:o:" optname; do
    case "$optname" in
        "g") hitFile=${OPTARG} ;;
        "s") SNPFile=${OPTARG} ;;
        "x") geneName=${OPTARG} ;;
        "h") display_help ;;
        "o") outputFile=${OPTARG} ;;
        "w") window=$OPTARG ;;
        "d") checkLD=1;;
        "?") display_help "[Error] Unknown option $OPTARG" ;;
        ":") display_help "[Error] No argument value for option $OPTARG";;
        *) display_help "[Error] Unknown error while processing options";;
    esac
done

# If output file is specified we have to update with the absolute path:
if [[ $outputFile != "/dev/stdout" ]]; then outputFile=${workingDir}/${outputFile}; fi

# Testing files:
if [[ -z "${hitFile}" ]]; then display_help "[Error] Specify MONSTER output with -g."; fi
if [[ ( -z "${SNPFile}" ) && ( -z "${geneName}" ) ]]; then display_help "[Error] Specify SNP list (-s) or gene name (-x)."; fi

# Testing hit file:
if [[ $hitFile =~ .*.tar.gz ]]; then
    echo "[Info] Input file is compressed. Decompressing." >&2
    tar -xvzf $hitFile > /dev/null
    
    
    # checking if extraction went well:
    if [ $? -ne 0 ]; then echo "[Error] Extracting ${hitFile} failed. Exiting" >&2; exit 1; fi

    # Testing and entering working directory:
    workingDir=${workingDir}/$(basename "${hitFile}" | sed -e 's/\.tar\.gz//')
    if [ ! -d ${workingDir} ]; then
        echo "[Error] Entering directory has failed. Exiting $workingDir" >&2; exit 1;
    fi
    
else # If directory is submitted:
    if [ -d ${hitFile} ]; then workingDir=${hitFile};
    else echo "[Error] Input hit file is not a compressed file or a directory. Exiting" >&2; exit 1; fi
fi

# If the provided SNP input is a file, then generate absolute path:
if [[ -e "${SNPFile}" ]]; then SNPFile=$(echo "$(cd "$(dirname "${SNPFile}")"; pwd)/$(basename "${SNPFile}")"); fi

# If we are giving a gene name, we have to extract all variants with phenotypes:
if [[ ! -z ${geneName} ]]; then
    SNPFile=$( get_phenotypes_for_gene "$geneName" );
fi

# Entering directory:
cd $workingDir

echo -e "Gene\tTrait\tBurden_pVal\tVariant\tSNPID\tSpPval\tSampleCnt\tSNPCnt\tCond_pVal" > "${outputFile}"

# Testing if all the required files are located in the working directory:
for file in "pheno.mod.ordered.txt" "genotype.mod.filtered.txt" "snpfile.mod.txt" "kinship.mod.filtered.txt"; do
    if [[ ! -e "${file}" ]]; then "[Error] $file is required but missing! Exiting." >&2; exit 1; fi;
done

# Get basic information:
trait=$( getTrait )
gene=$(ls *bed | head -n1 | cut -f1 -d "_")
basePval=$( tail -n1  *_${gene}.out  | cut -f5 )

# Step 1: Processing variants:
if [[ -e "${SNPFile}" ]]; then # A valid file has been submitted.
    echo "[Info] A file with variants has been submitted." >&2
    
    cat "${SNPFile}" | while read variant; do
        read chr pos ID type <<< $( get_coordinates "${variant}" )
        # Expected output: APOC3	1445	5	0.1	6.3049e-16
        singlePval=$( getSinglePoint ${trait} ${chr}:${pos}-${pos} )
        read SNPID genex sampleCnt snpCNT x pval <<< $( test_variant ${chr} ${pos} ${variant} )
        echo -e "${gene:-NA}\t${trait:-NA}\t${basePval:-NA}\t${variant:-NA}\t${SNPID:-NA}\t${singlePval:-NA}\t${sampleCnt:-NA}\t${snpCNT:-NA}\t${pval:-NA}" >> "${outputFile}" 
    done
else
    read chr pos ID type <<< $( get_coordinates "${SNPFile}" )
    singlePval=$( getSinglePoint ${trait} ${chr}:${pos}-${pos} )
    read SNPID genex sampleCnt snpCNT x pval <<< $( test_variant ${chr} ${pos} ${variant} )
    echo -e "${gene:-NA}\t${trait:-NA}\t${basePval:-NA}\t${variant:-NA}\t${SNPID:-NA}\t${singlePval:-NA}\t${sampleCnt:-NA}\t${snpCNT:-NA}\t${pval:-NA}" >> "${outputFile}"
fi

# If output filename h0as been specified and LD info is requested:
if [[ ("${outputFile}" != "/dev/stdout") &&  ($checkLD -eq 1 ) ]] ; then
    echo "[Info] Adding LD information to the tested variants..."
    getLDblocks "${outputFile}" "${gene}"
fi

