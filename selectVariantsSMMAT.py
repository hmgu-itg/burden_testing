#!/usr/bin/python3

import argparse
import logging

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
parser.add_argument("--gencode",action="store",type=str,help="Optional: comma separated list of GENCODE features (%s)" %(",".join(config.GENCODE_FEATURES+"[all]")),required=False)
parser.add_argument("--gtex",action="store",type=str,help="Optional: comma separated list of regulatory features (%s)" %(",".join(config.REG_FEATURES+"[all]")),required=False)
parser.add_argument("--overlap",action="store",type=str,help="Optional: comma separated list of regulatory features (%s)" %(",".join(config.REG_FEATURES+"[all]")),required=False)
parser.add_argument("--extend",action="store",type=int,help="Optional: by how many basepairs the GENCODE features should be extended; default: %d" %(config.DEFAULT_EXTENSION),default=config.DEFAULT_EXTENSION,required=False)
parser.add_argument("--skipminor",action="store_true",help="Optional: skip minor APPRIS transcripts; default: False",required=False)
parser.add_argument("--verbose",help="Optional: verbosity level; default: info",required=False,choices=("debug","info","warning","error"),default="info")
parser.add_argument("--score",action="store",type=str,help="Optional: which score to use to weight variants; default: none",required=False,default="none",choices=("CADD","EigenPhred","none"))
parser.add_argument("--lof",action="store_true",help="Optional: only select high impact variants; default: False",required=False)

# sub-commands
subparsers=parser.add_subparsers(dest="subcommand",help="sub-command help")
parser_monster=subparsers.add_parser("monster", help="MONSTER help")
parser_smmat=subparsers.add_parser("smmat", help="SMMAT help")

# MONSTER options
parser_monster.add_argument("--vcf",action="store",type=str,help="Input VCF(s)",required=True)
parser_monster.add_argument("--maxvars",action="store",type=int,help="Optional: max number of variants in a gene",default=config.DEFAULT_MAXVARS,required=False)
parser_monster.add_argument("--maf",action="store",type=float,help="Optional: MAF upper threshold",default=config.DEFAULT_MAF,required=False)
parser_monster.add_argument("--mac",action="store",type=float,help="Optional: MAC lower threshold",default=config.DEFAULT_MAC,required=False)
parser_monster.add_argument("--missingness",action="store",type=float,help="Optional: upper missingness threshold",default=config.DEFAULT_MISSINGNESS,required=False)
parser_monster.add_argument("--chr-prefix",action="store",type=str,help="Optional: chromosome prefix in VCF files",default=config.CHR_PREFIX,required=False)
mutex_options_loftee=parser_monster.add_mutually_exclusive_group()
mutex_options_loftee.add_argument("--loftee",action="store_true",help="Optional: only select LOFTEE variants",required=False)
mutex_options_loftee.add_argument("--lofteeHC",action="store_true",help="Optional: only select high confidence LOFTEE variants",required=False)

# SMMAT options
parser_smmat.add_argument("--smmat",action="store",type=str,help="5 column tab-delimited input list of variants")

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
if args.subcommand=="monster":
    variant_filters.append(filters.createMAFFilter(args.maf))
    variant_filters.append(filters.createMACFilter(args.mac))
    variant_filters.append(filters.createMISSFilter(args.missingness))
    if args.loftee:
        variant_filters.append(filters.createLofteeFilter())
    if args.lofteeHC:
        variant_filters.append(filters.createLofteeHCFilter())
        
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


# ------------------------------- BEDTOOLS



# ------------------------------- FILTER



# ------------------------------- SELECT VARIANTS


# ------------------------------- FILTER VARIANTS



# ------------------------------- ADD SCORES














