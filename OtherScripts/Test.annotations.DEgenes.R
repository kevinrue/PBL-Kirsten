########
# test for most efficient  annotation for DE genes
########    

# List of all informative probesets in the microarray data: 5082
DE.probe.list = PBL.toptable[PBL.toptable$adj.P.Val < 0.05,]$ID
length(DE.probe.list)

library(bovine.db)

table = toTable(bovineENSEMBL[DE.probe.list])
length(unique(table$probe_id)) # 1834 unique probe ids are found in the bovine.db annotation
length(unique(table$ensembl_id)) # 1705 unique ensembl ids of the bovine.db package map to the 1834 probe ids in our dataset
2808 - 1834 # 974 DE probe sets are not annotated in the bovine.db package 
1834/2808 # 65% of our the probesets in our dataset have an annotation in the bovine.db package

# How many annotated probesets in the affymetrix most recent annotation?
affy_annot = read.table(file="ProbeAnnotation/Bovine.na33.annot_filtered.csv", header=T, sep=",", stringsAsFactors=F)

affy_annot$Probe.Set.ID
affy_annot$Ensembl

# How many of the probesets
sum(affy_annot[affy_annot$Probe.Set.ID %in% DE.probe.list,]$Ensembl != "---")