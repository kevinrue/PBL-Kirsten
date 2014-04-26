########
# test for most efficient  annotation
########    

# List of all informative probesets in the microarray data: 5082
test.probe.list = PBL.toptable$ID
length(test.probe.list)

# Obtain the ensembl ids for the informative probesets
test.probe.ensembl.table = toTable(bovineENSEMBL[test.probe.list])
# Number of unique probesets mapped : 3360
length(unique(test.probe.ensembl.table$probe_id))
# Number of unique ensembl ids identified from the list of probesets: 3043
length(unique(test.probe.ensembl.table$ensembl_id))
# Number of final ensembl ids compared to number of initial probe sets: 2039
5082-3043
# Number of probe sets lost in the process: 1722 informative probesets do not map to any ensembl id in bovine.db
5082-3360


# Obtain the gene symbols for the informative probesets
test.probe.symbol.table = toTable([test.probe.list])
nrow(test.probe.symbol.table)
# Number of unique probesets mapped : 3607
length(unique(test.probe.symbol.table$probe_id))
# Number of unique ensembl ids identified from the list of probesets: 3223
length(unique(test.probe.symbol.table$symbol))
# Number of final gene symbol compared to number of initial probe sets: 1859
5082-3223
# Number of probe sets lost in the process: 1475 informative probesets do not map to any gene symbol in bovine.db
5082-3607