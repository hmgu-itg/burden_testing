import logging
import re

LOGGER=logging.getLogger(__name__)

# ==============================================================================================================================

def createMAFFilter(f):
    def filter(v):
        maf=v["MAF"]
        if maf is None:
            return False
        return maf<f
    return filter

# ==============================================================================================================================

def createMACFilter(c):
    def filter(v):
        ac=int(v["AC"])
        if ac is None:
            return False        
        an=int(v["AN"])
        if an is None:
            return False
        maf=ac/an
        if maf>0.5:
            ac=an-ac
        return ac<c
    return filter

# ==============================================================================================================================

def createMISSFilter(miss):
    def filter(v):
        ms=v["MISS"]
        if ms is None:
            return False
        return ms<miss
    return filter

# ==============================================================================================================================

def createLofteeFilter:
    def filter(v):
        lt=v["loftee"]
        if lt is None:
            return False
        return lt!="-"
    return filter

# ==============================================================================================================================

def createLofteeHCFilter:
    def filter(v):
        lt=v["loftee"]
        if lt is None:
            return False
        return lt=="HC"
    return filter
