#!/usr/bin/python3

import argparse
import logging
import os

from burden import config
from burden import utils
from burden import filters
from burden import io

# ------------------------------------------------------------------------------------------------------------------------

parser=argparse.ArgumentParser()

# common options
parser.add_argument("--input",action="store",type=str,help="Required: input gene list",required=True)
parser.add_argument("--config",action="store",type=str,help="Required: config file",required=True)
parser.add_argument("--output-dir",action="store",type=str,help="Required: output directory",required=True)
parser.add_argument("--output-prefix",action="store",type=str,help="Required: output filename prefix",required=True)
parser_smmat.add_argument("--smmat",action="store",type=str,help="Required: 5 column tab-delimited input list of variants, bgzipped and tabixed",required=True)

parser.add_argument("--gencode",action="store",type=str,help="Optional: comma separated list of GENCODE features (%s)" %(",".join(config.GENCODE_FEATURES+"[all]")),required=False)
parser.add_argument("--gtex",action="store",type=str,help="Optional: comma separated list of regulatory features (%s)" %(",".join(config.REG_FEATURES+"[all]")),required=False)
parser.add_argument("--overlap",action="store",type=str,help="Optional: comma separated list of regulatory features (%s)" %(",".join(config.REG_FEATURES+"[all]")),required=False)
parser.add_argument("--extend",action="store",type=int,help="Optional: by how many basepairs the GENCODE features should be extended; default: %d" %(config.DEFAULT_EXTENSION),default=config.DEFAULT_EXTENSION,required=False)
parser.add_argument("--skipminor",action="store_true",help="Optional: skip minor APPRIS transcripts; default: False",required=False)
parser.add_argument("--verbose",help="Optional: verbosity level; default: info",required=False,choices=("debug","info","warning","error"),default="info")
parser.add_argument("--score",action="store",type=str,help="Optional: which score to use to weight variants; default: none",required=False,default="none",choices=("CADD","EigenPhred","none"))
parser.add_argument("--lof",action="store_true",help="Optional: only select high impact variants; default: False",required=False)

args=parser.parse_args()

# ------------------------------------------------------------------------------------------------------------------------

record_filters=[]
gencode_opts=None
if not args.gencode is None:
    if "all" in args.gencode.split(","):
        gencode_opts=config.GENCODE_FEATURES
    else:
        gencode_opts=args.gencode.split(",")
    record_filters["GENCODE"]=gencode_opts
gtex_opts=None
if not args.gtex is None:
    if "all" in args.gtex.split(","):
        gtex_opts=config.REG_FEATURES
    else:
        gtex_opts=args.gtex.split(",")
    record_filters["GTEx"]=gtex_opts
overlap_opts=None
if not args.overlap is None:
    if "all" in args.overlap.split(","):
        overlap_opts=config.REG_FEATURES
    else:
        overlap_opts=args.overlap.split(",")
    record_filters["overlap"]=overlap_opts

variant_filters=[]
if args.lof:
    variant_filters.append(filters.createLofFilter())
        
# -------------------------------------------------------- LOGGING -------------------------------------------------------

verbosity=config.LOGGING[args.verbose]
outdir=args.output_dir
if outdir.endswith("/"):
    outdir=outdir[:-1]
if not os.path.exists(outdir):
    os.makedirs(outdir)
logfile=outdir+"/variant_selector.log"

LOGGER=logging.getLogger("Variant Selector")
LOGGER.setLevel(verbosity)
#ch=logging.StreamHandler()
ch=logging.FileHandler(logfile,'w')
ch.setLevel(verbosity)
formatter=logging.Formatter('%(levelname)s - %(name)s - %(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
ch.setFormatter(formatter)
LOGGER.addHandler(ch)
LOGGER.addHandler(logging.StreamHandler(sys.stdout))

for m in ["utils","filters","io"]:
    logging.getLogger("burden."+m).addHandler(ch)
    logging.getLogger("burden."+m).addHandler(logging.StreamHandler(sys.stdout))
    logging.getLogger("burden."+m).setLevel(verbosity)

# ------------------------------------------------------------------------------------------------------------------------

for arg in vars(args):
    LOGGER.info("INPUT OPTIONS: %s : %s" % (arg, getattr(args, arg)))
LOGGER.info("")

# ------------------------------------------------------------------------------------------------------------------------
        
io.readConfig(args.config)
GENCODE=io.readGencode()

with open(args.input) as F:
    for line in F:
        current_gene=line.rstrip()
        if not current_gene in GENCODE:
            if current_gene in GENCODE["duplicates"]:
                LOGGER.warning("Gene %s occurs multiple times in GENCODE; skipping" %(current_gene))
            else:
                LOGGER.warning("Gene %s not found in GENCODE; skipping" %(current_gene))
            continue
        rec=GENCODE[current_gene]


        
# ------------------------------- BEDTOOLS



# ------------------------------- FILTER



# ------------------------------- SELECT VARIANTS


# ------------------------------- FILTER VARIANTS



# ------------------------------- ADD SCORES














