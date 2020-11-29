import logging
import re
import sys
import gzip

import config

LOGGER=logging.getLogger(__name__)

# ==============================================================================================================================

def readGENCODE():
    D=dict()
    D["duplicates"]=[]
    par_re=re.compile("^ENSG\d+.*_PAR_Y$")
    ID_re=re.compile("^(ENSG\d+)\.")
    with gzip.open(config.CONFIG["gencode_file"]) as F:
        for line in F:
            (chrom,start,end,name,ID)=line.rstrip().split("\t")
            if par_re.match(ID):
                continue
            if name in D:
                LOGGER.warning("duplicate gene name: %s" %(name))
                del D[name]
                D["duplicates"].append(name)
            else:
                rec={"chr":chrom,"start":start,"end":end,"name":name,"ID":ID}
                D[name]=rec
                D[ID]=rec
                if ID_re.match(ID):
                    D[ID_re.match(ID).group(1)]=rec
    return D

# ==============================================================================================================================

def readConfig(fname):
    D=dict()
    with open(fname) as F:
        for line in F:
            (key,val)=line.rstrip().split("=")
            D[key]=val
        config.CONFIG=D
    if not all(x in config.CONFIG for x in config.CONFIG_KEYS):
        LOGGER.error("some keys in config (%s) are missing" %(fname))
        LOGGER.error("keys present: %s" %(",".join(list(config.CONFIG.keys()))))
        LOGGER.error("required keys: %s" %(",".join(list(config.CONFIG_KEYS))))
        sys.exit(1)

# ==============================================================================================================================

def writeSMMAT(variants,gene,fname):
    with open(fname,"w+") as F:
        for v in variants:
            if not v["score"] is None:
                F.write("\t".join([gene,v["chr"],v["pos"],v["ref"],v["alt"],v["score"]])+"\n")
