import logging

LOGGING={
    "debug":logging.DEBUG,
    "info":logging.INFO,
    "warning":logging.WARNING,
    "error":logging.ERROR
}

LOF_SEVERITY={
    "transcript_ablation"      : 10,
    "splice_acceptor_variant"  : 9,
    "splice_donor_variant"     : 8,
    "stop_gained"              : 7,
    "frameshift_variant"       : 6,
    "stop_lost"                : 5,
    "start_lost"               : 4,
    "transcript_amplification" : 3,
    "inframe_insertion"        : 2,
    "inframe_deletion"         : 1
}

EIGEN_SPECS={
    "CHR"=1,
    "POS"=2,
    "REF"=3,
    "ALT"=4,
    "SCORE"=7
}

CADD_SPECS={
    "CHR"=1,
    "POS"=2,
    "REF"=3,
    "ALT"=4,
    "SCORE"=6
}

GENCODE_FEATURES=["gene","exon","transcript","CDS","UTR"]
REG_FEATURES=["promoter","CTCF","enhancer","promoterFlank","openChrom","TF_bind"]

COFIG_KEYS=["VEPdir","VEPexec","Linked_features","gencode_file","EigenPath","caddPath"]

DEFAULT_EXTENSION=0
