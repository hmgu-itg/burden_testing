#!/usr/bin/python3

import argparse
import logging

from burden import config
from burden import utils

# ------------------------------------------------------------------------------------------------------------------------

# default
verbosity=logging.INFO

# ------------------------------------------------------------------------------------------------------------------------

parser=argparse.ArgumentParser()
input_options=parser.add_argument_group('Input options')
mutex_options=parser.add_mutually_exclusive_group(required=True)
mutex_options_loftee=parser.add_mutually_exclusive_group()

# Required options
input_options.add_argument("--input", "-i",action="store",type=str,help="Required: input gene list",required=True)
input_options.add_argument('--output-dir','-output-dir',action="store",type=str,help="Required: output directory where output and temporary files will be created",required=True)
input_options.add_argument('--output-prefix','-output-prefix',action="store",type=str,help="Required: output filename prefix",required=True)

input_options.add_argument('--gencode','-gencode',action="store",type=str,help="Optional: comma separated list of GENCODE features (gene, exon, transcript, CDS or UTR)",required=False)
input_options.add_argument('--gtex','-gtex',action="store",type=str,help="Optional: comma separated list of GTEx features (promoter, CTCF, enhancer, promoterFlank, openChrom, TF_bind or allreg)",required=False)
input_options.add_argument('--overlap','-overlap',action="store",type=str,help="Optional: comma separated list of overlap features (promoter, CTCF, enhancer, promoterFlank, openChrom, TF_bind or allreg)",required=False)
input_options.add_argument('--extend','-extend',action="store",type=int,help="Optional: by how many basepairs the GENCODE features should be extended",default=config.DEFAULT_EXTENSION,required=False)
input_options.add_argument('--skipminor','-skipminor',action="store_true",help="Optional: skip minor APPRIS transcripts",required=False)
input_options.add_argument("--verbose", "-v", help="Optional: verbosity level",required=False,choices=("debug","info","warning","error"),default="info")
input_options.add_argument('--score','-score',action="store",type=str,help="Optional: which score to use to weight variants",required=False,choices=("CADD","EigenPhred"))
input_options.add_argument('--lof','-lof',action="store_true",help="Optional: only select high impact variants",required=False)

# vcf subcommand
mutex_options.add_argument('--vcf','-vcf',action="store",type=str,help="Input VCF(s)")
mutex_options_loftee.add_argument('--loftee','-loftee',action="store_true",help="Optional: only select LOFTEE variants",required=False)
mutex_options_loftee.add_argument('--lofteeHC','-lofteeHC',action="store_true",help="Optional: only select high confidence LOFTEE variants",required=False)
input_options.add_argument('--maxvars','-maxvars',action="store",type=int,help="Optional: max number of variants in a gene",default=config.DEFAULT_MAXVARS,required=False)
input_options.add_argument('--maf','-maf',action="store",type=float,help="Optional: MAF upper threshold",default=config.DEFAULT_MAF,required=False)
input_options.add_argument('--mac','-mac',action="store",type=float,help="Optional: MAC lower threshold",default=config.DEFAULT_MAC,required=False)
input_options.add_argument('--missingness','-missingness',action="store",type=float,help="Optional: upper missingness threshold",default=config.DEFAULT_MISSINGNESS,required=False)
input_options.add_argument('--chr-prefix','-chr-prefix',action="store",type=str,help="Optional: chromosome prefix in VCF files",default=config.CHR_PREFIX,required=False)

# smmat subcommamd
mutex_options.add_argument('--smmat','-smmat',action="store",type=str,help="5 column tab-delimited input list of variants")


args=parser.parse_args()

# ------------------------------- READ COMMANDLINE ARGUMENTS


# ------------------------------- CREATE FILTER DICTS





# ------------------------------- BEDTOOLS



# ------------------------------- FILTER



# ------------------------------- SELECT VARIANTS


# ------------------------------- FILTER VARIANTS



# ------------------------------- ADD SCORES














