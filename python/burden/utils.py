import logging
import subprocess
import json
import re
import tempfile as tf
import os

import config

LOGGER=logging.getLogger(__name__)

# ==============================================================================================================================

def runLiftOver(input_data,build="38"):
    if build!="38" and build!="37":
        LOGGER.error("provided build: %s; build should be either 37 or 38" % build)
        return None
    
    L=list()

    chain="/usr/local/bin/hg38ToHg19.over.chain.gz"
    if build=="37":
        chain="/usr/local/bin/hg19ToHg38.over.chain.gz"

    in_bed=tf.NamedTemporaryFile(delete=False,mode="w",prefix="annotator_")
    LOGGER.debug("Input bed file: %s" % (in_bed.name))
    out_fname=tf.mktemp(prefix="annotator_")
    LOGGER.debug("Output bed file: %s" % (out_fname))
    unmapped_fname=tf.mktemp(prefix="annotator_")
    LOGGER.debug("Unmapped file: %s" % unmapped_fname)

    for x in input_data:
        in_bed.write("%s\t%d\t%d\t%s\n" %(x["chr"],int(x["start"]),int(x["end"]),x["id"]))
    in_bed.close()

    LOGGER.debug("Input: %d record(s)" % len(input_data))
    
    if which("liftOver") is None:
        LOGGER.error("liftOver was not found in PATH")
        return None
        
    LOGGER.debug("Calling: liftOver %s %s %s %s" %(in_bed.name,chain,out_fname,unmapped_fname))
    cmdline="liftOver %s %s %s %s" %(in_bed.name,chain,out_fname,unmapped_fname)
    subprocess.Popen(cmdline,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT).communicate()

    if not os.path.isfile(out_fname):
        LOGGER.error("liftOver failed to create output file %s" % out_fname)
        return None

    if os.path.getsize(out_fname)==0:
        LOGGER.warning("liftOver produced empty output file %s" % out_fname)
        return None

    count=0
    with open(out_fname) as F:
        for line in F:
            (chrom,start,end,ID)=line.rstrip().split("\t")
            L.append({"chr":chrom,"start":start,"end":end,"id":ID})
            count+=1

    LOGGER.debug("Remapped: %d record(s)" % count)
    unmapped=sum(1 for line in open(unmapped_fname) if not line.startswith("#"))
    LOGGER.debug("Unmapped: %d record(s)" % unmapped)
    
    LOGGER.debug("Removing temporary files")
    if os.path.isfile(in_bed.name):
        os.remove(in_bed.name)
    if os.path.isfile(out_fname):
        os.remove(out_fname)
    if os.path.isfile(unmapped_fname):
        os.remove(unmapped_fname)
        
    return L

# ==============================================================================================================================

def selectLines(cmd):
    L=list()
    LOGGER.debug("Command line: %s" %(cmd))
    output = subprocess.Popen(query,shell=True,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    for line in output.stdout.readlines():
        L.append(line)
    return L

# ==============================================================================================================================

def queryLinkedFeatures(chrom,start,end,gene_ID=None):
    L=list()
    LOGGER.debug("chrom: %s, start: %s, end: %s, filename: %s" %(chrom,start,end,config.CONFIG["Linked_features"]))
    cmd="intersectBed -wb -a <(echo -e \"%s\\t%s\\t%s\") -b %s -sorted 2>/dev/null | cut -f9-" %(chrom,start,end,config.CONFIG["Linked_features"])
    for line in selectLines(cmd):
        if gene_ID is None:
            d=json.decode(line)
            L.append(d)
        else:
            if gene_ID==d["gene_ID"]:
                L.append(d)
    return L

# ==============================================================================================================================

# key1: "source"
# key2: "class"
# filters: "GENCODE" --> ["exon", ... ], "overlap" --> ["promoter", ... ]
def filterRecords(records,filters,key1,key2):
    L=list()
    for r in records:
        if r[key1] in filters:
            if r[key2] in filters[r[key1]]:
                L.append(r)
    return L

# ==============================================================================================================================

def mergeRecords(records,extension):
    L=list()
    tmpfile=tf.NamedTemporaryFile(delete=False,mode="w",prefix="burden_")
    for x in sorted(records,key=lambda x:int(x["start"])):
        tmpfile.write("%s\t%d\t%d\n" %(x["chr"],int(x["start"])-int(extension),int(x["end"])+int(extension)))
    tmpfile.close()
    for line in selectLines("mergeBed -i %s" % (tmpfile.name)):
        (chrom,start,end)=line.split("\t")
        L.append({"chr":chrom,"start":start,"end":end})
    if os.path.isfile(tmpfile.name):
        os.remove(tmpfile.name)
    return L

# ==============================================================================================================================

# TODO: use temp file for tabix
def selectVariants(records,fname,prefix,vcf=False):
    L=list()
    if vcf:
        p_ac=re.compile("^AC=(\d+)")
        p_an=re.compile("^AN=(\d+)")
        p_loftee=re.compile("^LoF_conf=([^;]+)")
    cmd="tabix %s" %(fname)
    for r in records:
        cmd+=" %s%s:%s-%s" %(prefix,r["chr"],r["start"],r["end"])
    for line in selectLines(cmd):
        fields = line.split("\t")
        r={"chr":fields[0],"pos":fields[1],"ref":fields[3],"alt":fields[4],"id":"_".join([fields[0],fields[1],fields[3],fields[4]])}
        if vcf:
            f=fields[7].split(";")
            r["loftee"]=next((p_loftee.match(s).group(1) for s in f if p_loftee.match(s)),None)
            r["AC"]=next((p_ac.match(s).group(1) for s in f if p_ac.match(s)),None)
            r["AN"]=next((p_an.match(s).group(1) for s in f if p_an.match(s)),None)
            r["NS"]=len(fields)-9
        else:
            r["loftee"]=None
            r["AC"]=None
            r["AN"]=None
            r["NS"]=None
        L.append(r)
    return L

# ==============================================================================================================================

def filterVariants(variants,filters):
    return [s for s in filter( lambda x: all(f(x) for f in filters), variants )]

# ==============================================================================================================================

def getVepRecord(v):
    l0=len(v["ref"])
    l=len(v["alt"])
    if l0==1 and l==1:
        return "\t".join([v["chr"],v["pos"],v["pos"],v["ref"]+"/"+v["alt"],"+",v["id"]])
    if l0>l:
        return "\t".join([v["chr"],str(int(v["pos"])+l),str(int(v["pos"])+l0-1),v["ref"][l:]+"/-","+",v["id"]])
    if l0<l:
        return "\t".join([v["chr"],str(int(v["pos"])+l0),str(int(v["pos"])+l0-1),"-/"+v["alt"][l0:],"+",v["id"]])
    return None

# ==============================================================================================================================

def getMostSevereEffect(effects,severity_scores):
    max_score=0
    max_effect=None
    for e in effects.split(","):
        if e in severity_scores:
            if severity_scores[e]>max_score:
                max_score=severity_scores[e]
                max_effect=e
        else:
            if max_effect is None:
                max_effect=e
    return max_effect
                
# ==============================================================================================================================

def addConsequences(variants,geneID,severity_scores):
    ID=geneID
    m=re.match("(ENSG\d+)\.",geneID)
    if m:
        ID=m.group(1)
    vepdir=config.VEPDIR
    vepexec=config.VEPEXEC
    count=0
    vep_input=tf.NamedTemporaryFile(delete=False,mode="w",prefix="burden_vep_")
    for v in variants:
        s=getVepRecord(v)
        if s in None:
            LOGGER.warning("Wrong variant type for %s; skipping" %(v["id"]))
            continue
        else:
            vep_input.write("%s\n" %(s))
            count+=1
    vep_input.close()
    D=dict()
    if count==0:
        LOGGER.warning("No variants left")
    else:
        cmd=" ".join(["PERL5LIB=\$PERL5LIB:%s" %(vepdir),vepexec,"-i",vep_input.name,"--dir",vepdir,"--dir-cache",vepdir,"-o STDOUT --offline --no_stats | grep -v \"^#\" | awk -v g=%s" %(ID),"'BEGIN{FS=\"\\t\";}\$4==g{print \$0;}' | cut -f 1,7"])
        for line in selectLines(cmd):
            fields = line.split("\t")
            vID=fields[0]
            effs=fields[1]
            D[vID]=getMostSevereEffect(effs,severity_scores)
        for f in variants:
            if v["id"] in D:
                v["consequence"]=D[v["id"]]
            else:
                v["consequence"]=None
    if os.path.isfile(vep_input.name):
        os.remove(vep_input.name)

# ==============================================================================================================================

def liftOver(variants,build="37"):
    L=list()
    lifted=runLiftOver([{"chr":x["chr"],"start":str(int(x["pos"])-1),"end":x["pos"],"":x["id"]} for x in variants],build) # 0-based
    for v in variants:
        new_chrpos=next(({"chr":x["chr"],"pos":x["end"]} for x in lifted if x["id"]==v["id"]),None)
        if not new_chrpos is None:
            x=v.copy()
            x["chr"]=new_chrpos["chr"]
            x["pos"]=new_chrpos["pos"]
            L.append(x)
    return L

# ==============================================================================================================================

def addScore(variants,score_file_specs=None,score_file=None):
    if score_file_specs is None or score_file is None:
        for v in variants:
            v["score"]=1.0
    else:
        tmpfile=tf.NamedTemporaryFile(delete=False,mode="w",prefix="burden_score_")
        for v in variants:
            tmpfile.write("%s\t%s\n" %(v["chr"],v["pos"])) # 1-based
        tmpfile.close()
        L=list()
        for line in selectLines("tabix -R %s %s | cut -f %s,%s,%s,%s" % (tmpfile.name,score_file,score_specs["CHR"],score_specs["POS"],score_specs["ALT"],score_specs["SCORE"])):
            (chrom,pos,alt,score)=line.split("\t")
            L.append({"chr":chrom,"pos":pos,"alt":alt,"score":score})
        for v in variants:
            v["score"]=next((x["score"] for x in L if x["chr"]==v["chr"] and x["pos"]==v["pos"] and x["alt"]==v["alt"]),None)
        if os.path.isfile(tmpfile.name):
            os.remove(tmpfile.name)
