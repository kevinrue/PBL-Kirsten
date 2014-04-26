
# Understand CellMix Abbas cell signature ---------------------------------

library(CellMix)

# ExpressionSet containing the markers human probesets for each pure cell type
packageData("Abbas")
Abbas
pData(Abbas)

# Extract the list of cell type markers from the Abbas expression set
markerList = MarkerList(Abbas)
summary(markerList)

# We find the 359 features again in the lists of markers for the 17 cell types
# Note: 2 cell types have 0 markers, what are they doing here?
geneIds(markerList)
sapply(X=markerList, FUN=length)
sum(sapply(X=markerList, FUN=length))


# Create a table of two columns: Cell type, Human probeset ----------------

cellType = c()
humanProbeset = c()
# For each cell type
for(ct in names(geneIds(markerList))){
  # For each probeset marker of that cell type
  for(hps in geneIds(markerList)[[ct]]){
    # add an entry in the mapping table
    # (prepared using one list per column of the mapping table)
    cellType = c(cellType, ct)
    humanProbeset = c(humanProbeset, hps)
  }
  rm(ct, hps)
}
mapping = data.frame(CellType=cellType, humanProbeset=humanProbeset, stringsAsFactors=F)
rm(cellType, humanProbeset)

# Check that each probeset is a marker for at most one cell type
# Counts number of cell types per probeset, and summarises the range: 1-1
range(tapply(X=mapping$CellType, INDEX=mapping$humanProbeset, FUN=function(x){length(x)}))


# Create a table of two columns: human probeset, ensembl gene id ----------

library(hgu133a.db)
library(hgu133b.db)
# Use the id converter from CellMix to obtain the ensembl gene id for each probeset
hgu133a_ensembl = convertIDs(object=mapping$humanProbeset, to="ENSEMBL", from="hgu133a.db", unlist=F)
hgu133b_ensembl = convertIDs(object=mapping$humanProbeset, to="ENSEMBL",from="hgu133b.db")
# Show that each probeset is uniquely annotated by either hgu133a or hgu133b
sum(!is.na(hgu133a_ensembl) & !is.na(hgu133b_ensembl))
# How many don't have ensembl gene id in either array? 41
sum(is.na(hgu133a_ensembl) & is.na(hgu133b_ensembl))
# NOTE: We will have to remove those 41 probesets
# Merge lists of ensembl, leave NA where no annotation
mapping2 = mapping
mapping2$humanENS = NA
mapping2$humanENS[which(!is.na(hgu133a_ensembl))] = hgu133a_ensembl[which(!is.na(hgu133a_ensembl))]
mapping2$humanENS[which(!is.na(hgu133b_ensembl))] = hgu133b_ensembl[which(!is.na(hgu133b_ensembl))]
rm(hgu133a_ensembl, hgu133b_ensembl)
# The column in mapping2 is still a list at the moment
mapping2$humanENS = unlist(mapping2$humanENS, use.names=F)
# NOTE: the table is 359 rows, meaning that there has been no duplications, let's check: 0-1
# NOTE: Indeed, either one annotation, or NA (zero)
range(tapply(X=mapping2$humanENS, INDEX=mapping2$humanProbeset, FUN=function(x){sum(!is.na(x))}))
# NOTE: this means that each probeset got either 1 or 0 ensembl gene id
# In the mapping2 table, there are now 41 human probesets without annotation
sum(is.na(mapping2$humanENS))
# NOTE: this mean 359-41= 318 probesets have an ensembl gene id
# How many unique ensembl gene ids do we have in fact? 262
sum(!is.na(unique(mapping2$humanENS)))

# How many ensembl gene id does each tissue have?
tapply(X=mapping2$humanENS, INDEX=mapping2$CellType, FUN=function(x){sum(!is.na(x))})
# Can compare to the initial number of probeset per cell type in Abbas
sapply(X=geneIds(markerList), FUN=length)

# How many probesets lead to the same ensembl gene id? 1-2
# NOTE: no need for unique(x) because we know there has been no row duplication so far
range(tapply(X=mapping2$humanProbeset, INDEX=mapping2$humanENS, FUN=function(x){sum(!is.na(x))}))
# How many cell types is each human ensemble gene id a marker for? 1-2
# NOTE: we know already that different probesets lead to the same human ensembl gene id
# if those probesets as markers of the same tissue, the tissue will be counted multiple times
# to prevent that, use unique()
range(tapply(X=mapping2$CellType, INDEX=mapping2$humanENS, FUN=function(x){sum(!is.na(unique(x)))}))


# Create a table of two columns: human ensembl id, bovine ensembl  --------

library(biomaRt)
# Check which Biomart web services are available
listMarts()
# Loads the set of ensembl mart
mart = useMart(biomart="ensembl")
# Look at the datasets in the mart
listDatasets(mart)
# Search for mention of "Homo sapiens" in description
listDatasets(mart)[grep(pattern="Homo sapiens", x=listDatasets(mart)[,"description"]),]
# redefine mart to the human information
mart = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
# We look for filters allowing to select the human_ensembl_id  we have 
listFilters(mart)
names(listFilters(mart))
listFilters(mart)[grep(pattern="ensembl", x=listFilters(mart)[,"name"]),] # ensembl_gene_id 

# We look for attribute fields containing the bovine orthologous ensembl_id
listAttributes(mart)
# Looks for mention of "Cow" in the "description" column
listAttributes(mart)[grep(pattern="Cow", x=listAttributes(mart)[,"description"]),] # btaurus_homolog_ensembl_gene
# Looks for mention of "ensembl" in the "name" column
listAttributes(mart)[grep(pattern="ensembl", x=listAttributes(mart)[,"name"]),] # ensembl_gene_id 

# Get the human and the bovine ensembl id for all ensembl gene ids in our 
ensemblMap = getBM(attributes=c("ensembl_gene_id", "btaurus_homolog_ensembl_gene"), filters="ensembl_gene_id", values=mapping2$humanENS, mart=mart, bmHeader=TRUE)
# Check that we have again our 262 human ensembl gene ids
length(unique(ensemblMap$"Ensembl Gene ID"))

# To avoid a table with empty characters as the bovine, we could filter biomart for all the ones with orthologs in cow
# NOTE: we will filter for the human ensembl gene ids in our data when we merge the tables
listFilters(mart)[grep(pattern="Cow", x=listFilters(mart)[,"description"]),] # with_homolog_btau 
# Get the probe id, ensembl id, and gene symbol for the top 10 DE genes (adj. p-value)
ensemblMap = getBM(attributes=c("ensembl_gene_id", "btaurus_homolog_ensembl_gene"), filters="with_homolog_btau", values=TRUE, mart=mart, bmHeader=TRUE)
rm(mart)
# Merges the bovine homolog ensembl gene ids with the existing mapping2 table
mapping3 = merge(x=mapping2, y=ensemblMap, by.x="humanENS", by.y="Ensembl Gene ID", all.x=T, all.y=F, sort=F)
rm(ensemblMap)
# NOTE: The number of rows has increased from 359 to 380
# NOTE: this means that some human ensembl gene ids have multi-mapped to bovine ensembl gene ids
# NOTE: this means that the same human probeset (each of them were initially unique to one cell type)
# Rename the cow ensembl gene id to a name without spaces
colnames(mapping3)[colnames(mapping3) == "Cow Ensembl Gene ID"] = "bovineENS"

# How many ensembl gene id does each tissue have left?
tapply(X=mapping3$bovineENS, INDEX=mapping3$CellType, FUN=function(x){sum(!is.na(unique(x)))})
# Can compare to the initial number of markers per cell type in Abbas
sapply(X=geneIds(markerList), FUN=length)

# How many bovine ensembl do each human ensembl annotate to? 0-7
range(tapply(X=mapping3$bovineENS, INDEX=mapping3$humanENS, FUN=function(x){sum(!is.na(unique(x)))}))
# For instance, a DC marker (single human probeset,annotated by single human ensembl) converts to 7 bovine ensembl
mapping3[mapping3$humanENS == "ENSG00000158485",]
# How many human ensembl lead to a same bovine ensembl gene id? 1
# NOTE: No two human ensembl gene ids lead to the same bovine ensembl gene id
range(tapply(X=mapping3$humanENS, INDEX=mapping3$bovineENS, FUN=function(x){sum(!is.na(unique(x)))}))

# How many cell types is each bovine ensembl gene id associated with? 1-2
range(tapply(X=mapping3$CellType, INDEX=mapping3$bovineENS, FUN=function(x){sum(!is.na(unique(x)))}))


# Create a table of two columns: bovine ensembl gene id, bovine probeset --------
library(bovine.db)
bovineMAP = toTable(bovineENSEMBL2PROBE)
rm(bovineENS)
# Merges the conversion table with the existing mapping3 table
mapping4 = merge(x=mapping3, y=bovineMAP, by.x="bovineENS", by.y="ensembl_id", all.x=T, all.y=F, sort=F)
rm(bovineMAP)
# Rename the cow probeset column to a unambiguous name
colnames(mapping4)[colnames(mapping4) == "probe_id"] = "bovineProbeset"

# How many cow probesets does each tissue have left?
tapply(X=mapping4$bovineProbeset, INDEX=mapping4$CellType, FUN=function(x){sum(!is.na(unique(x)))})
# Can compare to the initial number of markers per cell type in Abbas
sapply(X=geneIds(markerList), FUN=length)

# How many probesets do each bovine ensembl gene id annotate to? 0-6
range(tapply(X=mapping4$bovineProbeset, INDEX=mapping4$bovineENS, FUN=function(x){sum(!is.na(unique(x)))}))
# How many bovine ensembl lead to a same probeset? 1-2
range(tapply(X=mapping4$bovineENS, INDEX=mapping4$bovineProbeset, FUN=function(x){sum(!is.na(unique(x)))}))

# NOTE: following all 3 successive conversions, lack of mapping and multimapping
# have led to some cell types having (1) less probesets in cow than human, and others
# having (2) more probesets in cow compared to human.

# Exports this critical table into a file
setwd("PhD/Human2BovineProbesets")
save(mapping4, file="mapping4.RData")


# Analysis of the probeset conversion table -------------------------------

# Keep only the rows which reached a bovine probeset
# The other rows don't have anything to offer anyway
noNA.mapping = mapping4[!is.na(mapping4$bovineProbeset),]

# How many cell types is each bovine probeset a marker for? 1-2
range(tapply(X=noNA.mapping$CellType, INDEX=noNA.mapping$bovineProbeset, FUN=function(x){sum(!is.na(unique(x)))}))

# Check that the human probesets have remained cell type specific (otherwise, I messed up somewhere)
# NOTE: each human probeset is still associated with at most 1 cell type
# therefore, the problem is due to the conversion, where two probesets from different cell types
# map to the same human ensembl gene id,bovine ensembl gene id, or bovine probeset
range(tapply(X=noNA.mapping$CellType, INDEX=noNA.mapping$humanProbeset, FUN=function(x){sum(!is.na(unique(x)))}))

# How many cow probesets are unique to a cell type? 218
length(which(tapply(X=noNA.mapping$CellType, INDEX=noNA.mapping$bovineProbeset, FUN=function(x){sum(!is.na(unique(x)))}) == 1))


# Keep only those 218 bovine probesets unique to a cell type: leaves 266 rows (multimapping)
mappingBovineUniqueCelltype = noNA.mapping[noNA.mapping$bovineProbeset %in% names(which(tapply(X=noNA.mapping$CellType, INDEX=noNA.mapping$bovineProbeset, FUN=function(x){sum(!is.na(unique(x)))}) == 1)),]
# Check that now, all remaining bovine probesets uniquely mark on cell type: 1-1
range(tapply(X=mappingBovineUniqueCelltype$CellType, INDEX=mappingBovineUniqueCelltype$bovineProbeset, FUN=function(x){sum(!is.na(unique(x)))}))
# Which  cell types are still represented in that mapping table, and how many bovine probesets do they have?
tapply(X=mappingBovineUniqueCelltype$bovineProbeset, INDEX=mappingBovineUniqueCelltype$CellType, FUN=function(x){sum(!is.na(unique(x)))})
# Compare to the initial list of cell types and human markers
sapply(X=geneIds(markerList), FUN=length)
# In this dataset how many cell types do the different oher identifiers mark?
# 1 bovine ensembl gene id uniquely marks one cell type
range(tapply(X=mappingBovineUniqueCelltype$CellType, INDEX=mappingBovineUniqueCelltype$bovineENS, FUN=function(x){sum(!is.na(unique(x)))}))
# 1 human ensembl gene id uniquely marks one cell type
range(tapply(X=mappingBovineUniqueCelltype$CellType, INDEX=mappingBovineUniqueCelltype$humanENS, FUN=function(x){sum(!is.na(unique(x)))}))
# 1 human probeset uniquely marks one cell type (as it always did)
range(tapply(X=mappingBovineUniqueCelltype$CellType, INDEX=mappingBovineUniqueCelltype$humanProbeset, FUN=function(x){sum(!is.na(unique(x)))}))

# How many bovine probesets can each human probeset be replaced by? 1-6
range(tapply(X=mappingBovineUniqueCelltype$bovineProbeset, INDEX=mappingBovineUniqueCelltype$humanProbeset, FUN=function(x){sum(!is.na(unique(x)))}))
# How many human probesets can each bovine probeset replace? 1-4
range(tapply(X=mappingBovineUniqueCelltype$humanProbeset, INDEX=mappingBovineUniqueCelltype$bovineProbeset, FUN=function(x){sum(!is.na(unique(x)))}))

# How many unique human probesets do we have? 177
length(unique(mappingBovineUniqueCelltype$humanProbeset))
# How many unique bovine probesets do we have? 218
length(unique(mappingBovineUniqueCelltype$bovineProbeset))
# NOTE: We have more bovine probeset than human probesets
# We will be able to test different replacements of human probesets by bovine probesets
# NOTE: my filters have ensured that one bovine probeset can only mark one cell type
# therefore, it may replace different human probesets, as long as they all mark that same cell type.
# How many mapping relationshipsdo we have between them?

# Knowing this, we can randomly replace human probesets by one the possible corresponding bovine probesets
# We know that each human probeset in the table has one or more bovine probeset replacements in the table
# we know that each bovine probeset can only be used to replace human probesets from a single cell type



# Save this critical mapping table
save(mappingBovineUniqueCelltype, file="mappingBovineUniqueCelltype.RData")


# Gene Symbol annotations readily available -------------------------------

# Obtains mapping between human probeset and gene symbol from Abbas
symbolMapping = data.frame(humanProbeset=featureNames(Abbas), symbol=featureData(Abbas)$SYMBOL)
# Merges the information with the existing table
mappingBovineUniqueCelltypeSymbol = merge(x=mappingBovineUniqueCelltype, y=symbolMapping, by.x="humanProbeset", by.y="humanProbeset", all.x=T, all.y=F, sort=F)
# Saves this nice table
save(mappingBovineUniqueCelltypeSymbol, file="mappingBovineUniqueCelltypeSymbol.RData")


# Modification of the Abbas signature -------------------------------------

load(file="mappingBovineUniqueCelltypeSymbol.RData")
map = mappingBovineUniqueCelltypeSymbol
# Re-initialise the Abbas variable
packageData("Abbas")
# Edit 
annotation(Abbas) = "bovine.db"
# We need first to remove from Abbas the human probesets without match in our table
# (or keep the ones present in our table)
Abbas = Abbas[featureNames(Abbas) %in% unique(mappingBovineUniqueCelltypeSymbol$humanProbeset),]
# Then for each human probese, we randomly replace it by one of its bovine counterpart
# NOTE: in case 2 or mroe human probesets conflict for the same bovine probeset
# the first one gets annotated, and the following one is remove from the ExpressionSet
# But to avoid being always the same ones removed from the eSet, the replacement will be performed
# in a random fashion
# For a random feature
rowToRemove = c()
for(i in sample(x=1:nrow(Abbas), size=nrow(Abbas))){
  bovine = sample(x=map[map$humanProbeset == featureNames(Abbas)[i],]$bovineProbeset, size=1)
  # checking that this counterpart hasn't been used to replace another human probeset
  if(bovine %in% featureNames(Abbas)){
    if(all(map[map$humanProbeset == featureNames(Abbas)[i],]$bovineProbeset %in% featureNames(Abbas))){
      rowToRemove = c(rowToRemove, i)
    }
    else {
      repeat{
        bovine = sample(x=map[map$humanProbeset == featureNames(Abbas)[i],]$bovineProbeset, size=1)
        if(!bovine %in% featureNames(Abbas)){
          featureNames(Abbas)[i] = bovine
          break
        }
      }
    }
  }
  else{
    featureNames(Abbas)[i] = bovine
  }
  rm(bovine)
}
if(length(rowToRemove) > 0){
  Abbas = Abbas[-rowToRemove,]
}
cat(length(rowToRemove), "rows removed. No more bovine probeset available to replace those human probesets.")
rm(map, rowToRemove)


# Deconvolution! ----------------------------------------------------------

# Custom command using the modified signature
res = gedProportions(object=farms, x=Abbas, method="lsfit", CLsubset=sampleNames(Abbas), normalize=F, log=F, verbose=2)

# result object
res
# proportions are stored in the coefficient matrix
dim(coef(res))
coef(res)[1:3, 1:4]

# cell type names
basisnames(res)
# basis signatures (with converted IDs)
basis(res)[1:5, 1:3]

# plot against actual CBC (ask Kate for sample-wise data)
##profplot(acr, cbc)
# plot cell proportion differences between groups
boxplotBy(res, pData(farms)$Treatment, main = "Cell proportions vs Transplant status")

# aggregate into CBC
cbc <- asCBC(res)
dim(cbc)

coef(cbc)
# plot cell proportion differences between groups ()
boxplotBy(cbc, pData(farms)$Treatment, main = "Cell proportions vs Transplant status")

basisnames(cbc)

# A few t-tests
t.test(x=coef(cbc)["Lymphocytes",which(pData(farms)$Treatment == "infected")], y=coef(cbc)["Lymphocytes", which(pData(farms)$Treatment == "control")], alternative="two.sided")
t.test(x=coef(cbc)["Neutrophils",which(pData(farms)$Treatment == "infected")], y=coef(cbc)["Neutrophils", which(pData(farms)$Treatment == "control")], alternative="two.sided")
t.test(x=coef(cbc)["Monocytes",which(pData(farms)$Treatment == "infected")], y=coef(cbc)["Monocytes", which(pData(farms)$Treatment == "control")], alternative="two.sided")

# I export one cell-type signature that worked relatively well
#save(Abbas, file="Abbas_bovine1.RData")

