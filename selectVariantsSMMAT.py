#!/usr/bin/python3

import argparse
import logging
import os
import sys

from burden import config
from burden import utils
from burden import filters
from burden import io

# ------------------------------------------------------------------------------------------------------------------------

parser=argparse.ArgumentParser()

# required
parser.add_argument("--input",action="store",type=str,help="Required: input gene list",required=True)
parser.add_argument("--config",action="store",type=str,help="Required: config file",required=True)
parser.add_argument("--output-dir",action="store",type=str,help="Required: output directory",required=True)
parser.add_argument("--smmat",action="store",type=str,help="Required: 5 column tab-delimited input list of variants, bgzipped and tabixed",required=True)

# optional
parser.add_argument("--gencode",action="store",type=str,help="Optional: comma separated list of GENCODE features (%s)" %(",".join(config.GENCODE_FEATURES+["all"])),required=False)
parser.add_argument("--gtex",action="store",type=str,help="Optional: comma separated list of regulatory features (%s)" %(",".join(list(config.REG_FEATURES.keys())+["all"])),required=False)
parser.add_argument("--overlap",action="store",type=str,help="Optional: comma separated list of regulatory features (%s)" %(",".join(list(config.REG_FEATURES.keys())+["all"])),required=False)
parser.add_argument("--extend",action="store",type=int,help="Optional: by how many basepairs the GENCODE features should be extended; default: %d" %(config.DEFAULT_EXTENSION),default=config.DEFAULT_EXTENSION,required=False)
parser.add_argument("--skipminor",action="store_true",help="Optional: skip minor APPRIS transcripts; default: False",required=False)
parser.add_argument("--verbose",help="Optional: verbosity level; default: info",required=False,choices=("debug","info","warning","error"),default="info")
parser.add_argument("--score",action="store",type=str,help="Optional: which score to use to weight variants; default: none",required=False,default=None,choices=("CADD","EigenPhred"))
parser.add_argument("--lof",action="store_true",help="Optional: only select high impact variants; default: False",required=False)
parser.add_argument("--log",action="store",help="Optional: log file; default: \"variant_selector.log\" in the output directory",required=False)

args=parser.parse_args()

# -------------------------------------------------------- LOGGING -------------------------------------------------------

verbosity=config.LOGGING[args.verbose]
outdir=args.output_dir
if outdir.endswith("/"):
    outdir=outdir[:-1]
if not os.path.exists(outdir):
    os.makedirs(outdir)
if args.log is None:
    logfile=outdir+"/variant_selector.log"
else:
    logfile=args.log
outfile=outdir+"/group_file.txt"

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

# ------------------------------------------------- RECORD AND VARIANT FILTERS -------------------------------------------

record_filters=dict()
gencode_opts=list()
if not args.gencode is None:
    if "all" in args.gencode.split(","):
        for f in args.gencode.split(","):
            if f=="all":
                continue
            if not f in config.GENCODE_FEATURES:
                LOGGING.warning("Provided GENCODE feature \"%s\" is not valid" %(f))
        LOGGING.info("Using all GENCODE features")
        gencode_opts=config.GENCODE_FEATURES
    else:
        for f in args.gencode.split(","):
            if f in config.GENCODE_FEATURES:
                if not f in gencode_opts:
                    gencode_opts.append(f)
            else:
                LOGGING.warning("Provided GENCODE feature \"%s\" is not valid; skipping" %(f))
    if len(gencode_opts)!=0:
        record_filters["GENCODE"]=gencode_opts
gtex_opts=list()
if not args.gtex is None:
    if "all" in args.gtex.split(","):
        for f in args.gtex.split(","):
            if f=="all":
                continue
            if not f in config.REG_FEATURES:
                LOGGING.warning("Provided GTEx feature \"%s\" is not valid" %(f))
        LOGGING.info("Using all regulatory features for \"GTEx\" source")
        gtex_opts=list(config.REG_FEATURES.keys())
    else:
        for f in args.gtex.split(","):
            if f in config.REG_FEATURES:
                if not config.REG_FEATURES[f] in gtex_opts:
                    gtex_opts.append(config.REG_FEATURES[f])
            else:
                LOGGING.warning("Provided GTEx feature \"%s\" is not valid; skipping" %(f))
    if len(gtex_opts)!=0:
        record_filters["GTEx"]=gtex_opts
overlap_opts=list()
if not args.overlap is None:
    if "all" in args.overlap.split(","):
        for f in args.overlap.split(","):
            if f=="all":
                continue
            if not f in config.REG_FEATURES:
                LOGGING.warning("Provided overlap feature \"%s\" is not valid" %(f))
        LOGGING.info("Using all regulatory features for \"overlap\" source")
        overlap_opts=list(config.REG_FEATURES.keys())
    else:
        for f in args.overlap.split(","):
            if f in config.REG_FEATURES:
                if not config.REG_FEATURES[f] in overlap_opts:
                    overlap_opts.append(config.REG_FEATURES[f])
            else:
                LOGGING.warning("Provided GTEx feature \"%s\" is not valid; skipping" %(f))
    if len(overlap_opts)!=0:
        record_filters["overlap"]=overlap_opts

score=args.score
variant_filters=list()
if args.lof:
    variant_filters.append(filters.createLofFilter())
if not score is None:
    variant_filters.append(filters.createIndelFilter())
        
# ------------------------------------------------------------------------------------------------------------------------

bp_extension=args.extend
in_list=args.smmat
io.readConfig(args.config)
for x in config.CONFIG:
    LOGGER.info("CONFIG: %s=%s" % (x,config.CONFIG[x]))
LOGGER.info("")

GENCODE=io.readGENCODE()

if os.path.isfile(outfile):
    os.remove(outfile)

# ------------------------------------------------------ MAIN LOOP -------------------------------------------------------

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
        LOGGER.info("CURRENT GENE: %s %s" %(current_gene,str(rec)))
        regions=utils.queryLinkedFeatures(rec["chr"],str(int(rec["start"])-1),rec["end"],rec["ID"])
        for r in regions:
            LOGGER.debug("REGIONS: %s:%s-%s" %(r["chr"],r["start"],r["end"]))
        filtered_regions=utils.filterRecords(regions,record_filters,"source","class")
        for r in filtered_regions:
            LOGGER.debug("FILTERED REGIONS: %s:%s-%s" %(r["chr"],r["start"],r["end"]))
        merged_regions=utils.mergeRecords(filtered_regions,bp_extension)
        for r in merged_regions:
            LOGGER.debug("MERGED REGIONS: %s:%s-%s" %(r["chr"],r["start"],r["end"]))
        variants=utils.selectVariants(merged_regions,in_list)
        for v in variants:
            LOGGER.debug("VARIANTS: %s\t%s\t%s\t%s" %(v["chr"],v["pos"],v["ref"],v["alt"]))        
        if len(variants)>0:
            if args.lof:
                utils.addConsequences(variants,rec["ID"])
                for v in variants:
                    LOGGER.debug("CONSEQUENCES: %s\t%s\t%s\t%s\t%s" %(v["chr"],v["pos"],v["ref"],v["alt"],v["consequence"]))        
            if score is None:
                utils.addScore(variants)
            else:
                utils.addScore(variants,config.SCORE_SPECS[score],config.CONFIG[config.SCORE_FILES[score]])
            for v in variants:
                LOGGER.debug("SCORES: %s\t%s\t%s\t%s\t%s\t%s" %(v["chr"],v["pos"],v["ref"],v["alt"],v["consequence"],v["score"]))        
            filtered_variants=utils.filterVariants(variants,variant_filters)
            for v in filtered_variants:
                LOGGER.debug("FILTERED VARIANTS: %s\t%s\t%s\t%s\t%s\t%s" %(v["chr"],v["pos"],v["ref"],v["alt"],v["consequence"],v["score"]))        
            if len(filtered_variants)>0:
                io.writeOutput(filtered_variants,current_gene,outfile)
        LOGGER.info("")












