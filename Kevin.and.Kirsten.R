# Moves to the directory with the data
setwd("C:/Users/Kevin RUE/Documents/PBL-Kirsten/")

# Saves the path to the current directory
workDir = getwd()

# lists files in the current folder
dir()
grep(pattern=".CEL$", x=dir())

# Saves the CEL file names
celFiles = dir()[grep(pattern=".CEL$", x=dir())]
celFiles

write.csv(celFiles, file="celFiles.csv", row.names=F)

# Methods for Affymetrix Oligonucleotide Arrays 
library(affy)

###Read annotated data file "TB_aug_cattle_adf.txt" into a tragets file
targets = read.AnnotatedDataFrame ("celFiles.txt", header = TRUE, as.is = TRUE)
pData(targets)

# get sample names from targets file
SampleNames <- sampleNames(targets)
SampleNames


###Read all CEL files into an affybatch object called rawdata###
rawdata = ReadAffy(filenames = celFiles, phenoData = targets)
rawdata

# Update the sample names attached to the raw data with the information saved earlier
sampleNames(rawdata) = SampleNames

# Look at phenodata attached to the raw data. Make sure the information matches the right samples
pData(rawdata)


###Create boxplot of rawdata###
boxplot(rawdata, col = "magenta", las=2)

#
library(limma)

# MDS of samples base on raw data
plotMDS(x=assayData(rawdata)$exprs, dim.plot=c(1,2))
plotMDS(x=assayData(rawdata)$exprs, dim.plot=c(1,3))
plotMDS(x=assayData(rawdata)$exprs, dim.plot=c(2,3))


# Methods for fitting probe-level models 
library(affyPLM)
# Very simple high level analysis of Affymetrix data
library(simpleaffy)

###Change parameters of image plots to show 4 x 4 arrays###
par(mfrow=c(4,4))
par(mar=c(2,2,2,2))

###Generate image plots of rawdata###
image(rawdata)

par(mfrow=c(1,1))
par(mar=c(0,0,0,0))

pData(rawdata)
nrow(pData(rawdata))
###Generate histogram of rawdata###
hist(rawdata, col = 1:nrow(pData(rawdata)), lty=1)
legend(x="topright", legend=rownames(pData(rawdata)), col=1:nrow(pData(rawdata)), lty=1)

# Factor Analysis for Robust Microarray Summarization 
library(farms)

###Normalize the object rawdata using the FARMS normalization method. Call the normalized object farms###
farms <- qFarms(rawdata)

###Write out farms normalised data###
write.table(assayData(farms)$exprs, sep="\t", file="farms.txt")

###Generate boxplot of normalised data###
boxplot(assayData(farms)$exprs, col = "magenta", las=2)

###Generate histogram of normalised data###
hist(assayData(farms)$exprs, col = "magenta")

###Find informative genes from the farms normalised genes using INIcalls###
INIs <- INIcalls(farms)
str(INIs)

farms_informative <-getI_Eset(INIs)

save(farms_informative, file = "farms_informative.RData")

write.table(assayData(farms_informative)$exprs, file = "PBL_farms_informative.txt")

###Take the infection status of each animal and call it group###
unique(targets$Treatment)
groups <- unique(targets$Treatment)
groups

###Let f be the infection status of each array###
f <- factor(targets$Treatment, levels = groups)
f

###Design the matrix for the 16 arrays, with the columns "Control" and "Infected"###
design <- model.matrix(~0 + f)
design
colnames(design) <- c("Control", "Infected")
design

###Make contrast matrix showing infected vs control###
contrasts <- makeContrasts(IvsC=Infected-Control, levels=design)
contrasts

###Fit a linear model to the farms normalized data, call it fit###
fit = lmFit(object=farms, design=design)
fit

###Apply contrasts to the linear model fit. Call this cotrast_fit###
contrast_fit = contrasts.fit(fit=fit, contrasts=contrasts)
contrast_fit

###Filter contrast_fit using the informative genes found earlier (farms_informative).
# Call the resulting gene list PBL_filtered ###
filter <- rownames(contrast_fit) %in% rownames(exprs(farms_informative))
PBL_filtered <- contrast_fit[filter,]
# Empirical Bayes Statistics for Differential Expression
ebayes_filtered <- eBayes(PBL_filtered)

###Write Aug08_DEgenes to an outfile###
write.table(ebayes_filtered, file = "PBL_DEgenes.txt", quote = FALSE, sep = "\t", eol = "\n", na = "NA", row.names = TRUE, col.names = TRUE) 

# 
PBL.toptable <- topTable(ebayes_filtered, coef="IvsC", number=nrow(ebayes_filtered), adjust.method="BH", sort.by="P")
write.table(PBL.toptable, file = "PBL_DEgenes_topTable.txt", quote = FALSE, sep = "\t", eol = "\n", na = "NA", row.names = TRUE, col.names = TRUE) 

# MA plot
plot(logFC~AveExpr, data=PBL.toptable, col=(PBL.toptable$adj.P.Val < 0.05)+1, pch=20, cex=0.70)
# How many DE genes
sum(PBL.toptable$adj.P.Val < 0.05)

# Function to annotate a dataset where column ID contains the probe set ID
AnnotateProbeTable = function(df, annotPkg="bovine.db", probeCol="ID",
         gene.symbol=TRUE, ENSEMBL=TRUE, ENTREZ=TRUE)
{
  # Annotates a data frame containing probe IDs with Official gene symbol,
  # ENTREZ gene ID, and ENSEMBL gene ID.
  #
  # Args:
  #   df: Data frame containing the probe ID.
  #   annotPkg: Annotation package to use. Default is "bovine.db"
  #   probeCol: identifier of the column containing the probe IDs. Default is "ID".
  #   symbol: If TRUE, annotates with official gene symbol.
  #   ENSEMBL: If TRUE, annotates with ENSEMBL gene ID.
  #   ENTREZ: If TRUE, annotates with ENTREZ gene ID.
  #
  # Returns:
  #   The input data.frame with annotation columns appended.
  #
  n = nrow(df)
  df$rankFC <- 1:n # Keeps track of original row order
  # Packages required to annotate
  library(annotate, quietly=TRUE)
  # Fancy way to pass the library to load as an argument
  library(annotPkg, character.only=TRUE, quietly=TRUE)
  # Error handling
  if (!probeCol %in% colnames(df)){
    stop("Column does not exist in data.frame: ",
         probeCol, ".")
  }
  # Obtains the gene symbol annotations
  if (gene.symbol){
    df = cbind(df, gene.symbol=getSYMBOL(df[,probeCol], annotPkg))
  }
  # Obtains the ENSEMBL annotations
  if (ENSEMBL){
    df = merge(df, toTable(bovineENSEMBL[df[,probeCol]]), by.x=probeCol, by.y="probe_id",
               all.x=TRUE)
  }
  # Obtains the ENTREZ annotations
  if (ENTREZ){
    df = merge(df, toTable(bovineENTREZID[df[,probeCol]]), by.x=probeCol, by.y="probe_id",
               all.x=TRUE)
  }
  # Restores the original row order and rownames (messed up by "merge")
  df = df[order(df$rankFC),]
  # Discards the ranking column
  df$rankFC = NULL
  return(df)
}

PBL.toptable$ID = rownames(PBL.toptable)

# Annotates the microarray top table with the known gene symbol for each probe set
annotated.PBL.topTable = AnnotateProbeTable(df=PBL.toptable, annotPkg="bovine.db",probeCol="ID", gene.symbol=T, ENSEMBL=F, ENTREZ=F)

head(annotated.PBL.topTable)
sum(!is.na(annotated.PBL.topTable$gene.symbol))/nrow(annotated.PBL.topTable)


# Loads the RNA differential expression result file
RNA_PBL_table = read.table(file="RNA_allGenes_PBL_final.txt", header=T, sep="\t", stringsAsFactors=F, na.strings="")

# 
names(RNA_PBL_table)

# We want the DE gene symbols in RNA seq where FDR < 0.05
RNA.DE.symbols = RNA_PBL_table$external_gene_id[RNA_PBL_table$FDR < 0.05]
sum(is.na(RNA.DE.symbols))

# Remove the NA gene symbols
RNA.DE.noNA = RNA.DE.symbols[which(!is.na(RNA.DE.symbols))]
RNA.DE.noNA


# Look at the annotated microarray results
head(annotated.PBL.topTable)
# Filter the microarray results where adjusted p value is below 0.05
micro.DE.symbols = as.character(annotated.PBL.topTable$gene.symbol[annotated.PBL.topTable$adj.P.Val < 0.05])
# Get rid of the NAs
micro.DE.noNA = micro.DE.symbols[which(!is.na(micro.DE.symbols))]

plot(Venn(Sets=list(micro.DE.noNA, RNA.DE.noNA), SetNames=c("Microarray", "RNA-Seq")))


# Repeat the Venn analysis using the Ensembl id
annotated.PBL.topTable = AnnotateProbeTable(df=PBL.toptable, gene.symbol=T, ENSEMBL=T, ENTREZ=F)

# Filter the microarray results where adjusted p value is below 0.05
micro.DE.ensembl = annotated.PBL.topTable$ensembl_id[annotated.PBL.topTable$adj.P.Val < 0.05]
# Get rid of the NAs
sum(is.na(micro.DE.ensembl)) / length(micro.DE.ensembl)
micro.DE.noNA = micro.DE.ensembl[which(!is.na(micro.DE.ensembl))]

# We want the DE gene symbols in RNA seq where FDR < 0.05
RNA.DE.ensembl = RNA_PBL_table$ensembl_gene_id[RNA_PBL_table$FDR < 0.05]

# There are no NA in ensembl ids for RNA seq because counts are based on ensembl ids

# Venn diagram
length(unique(RNA.DE.ensembl))
length(unique(micro.DE.noNA))
head(RNA.DE.ensembl)
head(micro.DE.noNA)


plot(Venn(Sets=list(micro.DE.noNA, RNA.DE.ensembl), SetNames=c("Microarray", "RNA-Seq")))



##########
# Overlapping set of ups and down
##########
# Microarray
# We want the list of significant up and down-regulated ensembl ids
# Filter the microarray results where adjusted p value is below 0.05
micro.DE.ensembl.up = annotated.PBL.topTable$ensembl_id[annotated.PBL.topTable$adj.P.Val < 0.05 & annotated.PBL.topTable$logFC > 0]
micro.DE.ensembl.down = annotated.PBL.topTable$ensembl_id[annotated.PBL.topTable$adj.P.Val < 0.05 & annotated.PBL.topTable$logFC < 0]
# We need to remove NAs for micrarrays
micro.DE.up.noNA = micro.DE.ensembl.up[which(!is.na(micro.DE.ensembl.up))]
micro.DE.down.noNA = micro.DE.ensembl.down[which(!is.na(micro.DE.ensembl.down))]

# We want the DE gene symbols in RNA seq where FDR < 0.05 and
RNA.DE.ensembl.up = RNA_PBL_table$ensembl_gene_id[RNA_PBL_table$FDR < 0.05 & RNA_PBL_table$logFC > 0]
RNA.DE.ensembl.down = RNA_PBL_table$ensembl_gene_id[RNA_PBL_table$FDR < 0.05 & RNA_PBL_table$logFC < 0]

Venn(Sets=list(micro.DE.up.noNA, micro.DE.down.noNA, RNA.DE.ensembl.up, RNA.DE.ensembl.down), SetNames=c("Up microarray", "Down microarray", "Up RNA-seq", "Down RNA-seq"))
plot(Venn(Sets=list(micro.DE.up.noNA, micro.DE.down.noNA, RNA.DE.ensembl.up, RNA.DE.ensembl.down), SetNames=c("Up microarray", "Down microarray", "Up RNA-seq", "Down RNA-seq")))

# What are the ensembl ids that are up and down at the same time in microarray?
ensembl.conflicts = intersect(micro.DE.up.noNA, micro.DE.down.noNA)
tmp = annotated.PBL.topTable[annotated.PBL.topTable$ensembl_id %in% ensembl.conflicts,]
tmp[order(tmp$ensembl_id),]

# Remove the inconsistent ensembl ids from the lists for the venn diagram
RNA.DE.up.filtered = RNA.DE.ensembl.up[!RNA.DE.ensembl.up %in% ensembl.conflicts]
RNA.DE.down.filtered = RNA.DE.ensembl.down[!RNA.DE.ensembl.down %in% ensembl.conflicts]
micro.DE.noNA.up.filtered = micro.DE.up.noNA[!micro.DE.up.noNA %in% ensembl.conflicts]
micro.DE.noNA.down.filtered = micro.DE.down.noNA[!micro.DE.down.noNA %in% ensembl.conflicts]

# Venn of the filtered data (removed inconsistent ensembl ids in microarray and corresponding rna seq rows)
Venn(Sets=list(micro.DE.noNA.up.filtered, micro.DE.noNA.down.filtered, RNA.DE.up.filtered, RNA.DE.down.filtered), SetNames=c("Up microarray", "Down microarray", "Up RNA-seq", "Down RNA-seq"))
plot(Venn(Sets=list(micro.DE.noNA.up.filtered, micro.DE.noNA.down.filtered, RNA.DE.up.filtered, RNA.DE.down.filtered), SetNames=c("Up microarray", "Down microarray", "Up RNA-seq", "Down RNA-seq")))


##########
# Correlation of fold change between techologies
###########
# RNA-seq has 12k ensembl ids post-filtering, and each of them is unique in the result file
# microarray has 5k informative probe sets, which we have converted into ensembl ids
head(annotated.PBL.topTable)
# We remove the NAs
micro.topTable.noNA = annotated.PBL.topTable[!is.na(annotated.PBL.topTable$ensembl_id),]
# then we summarise the duplicates in the microarray by taking the median
summarised.micro.topTable = aggregate(x=micro.topTable.noNA$logFC, by=list(ensembl_id=micro.topTable.noNA$ensembl_id), FUN=median)
# We rename the x column with its real content: median logFC
colnames(summarised.micro.topTable)[2] = "logFC"
head(summarised.micro.topTable)
# The datasetis 3043 rows and 2 columns (ensembl id and median logFC)
dim(summarised.micro.topTable)
# We replace the numbered rownames by the actual ensembl id
rownames(summarised.micro.topTable) = summarised.micro.topTable$ensembl_id


# then we look for the ensembl ids present in both microarray and rna seq datasets
intersect.logFC = intersect(summarised.micro.topTable$ensembl_id, RNA_PBL_table$ensembl_gene_id)
length(intersect.logFC) # 2965 ensembl ids overlap between the filtered RNA-seq data and informative microarray data
head(RNA_PBL_table)
rownames(RNA_PBL_table) = RNA_PBL_table$ensembl_gene_id
# and we plot those to compare the fold change in microarray and rna seq
plot(x=summarised.micro.topTable[intersect.logFC,]$logFC, y=RNA_PBL_table[intersect.logFC,]$logFC)

nrow(summarised.micro.topTable)
length(unique(summarised.micro.topTable$ensembl_id))
## We can color code the DE genes called by both technologies, but we need to identify those
# We take the list of ensembl ids used in the plot
head(intersect.logFC)
# We identify which of those are DE in microarray
intersect.logFC.DE.micro = intersect.logFC %in% micro.DE.noNA
sum(intersect.logFC.DE.micro)
# We identify which of those are DE in RNA seq
intersect.logFC.DE.RNA = intersect.logFC %in% RNA.DE.ensembl
sum(intersect.logFC.DE.RNA)
# The color code will be 1 for not DE, 2 for microarray specific, 3 for RNA seq specific, and 4 for DE if both
color.logFC = 1+1*intersect.logFC.DE.micro+2*intersect.logFC.DE.RNA
# Let's plot again with the colors
plot(x=summarised.micro.topTable[intersect.logFC,]$logFC, y=RNA_PBL_table[intersect.logFC,]$logFC, col=color.logFC, pch=20, cex=0.8)
legend(x="topleft", legend=c("None", "Microarray", "RNA-seq", "Both"), col=1:4, pch=20)
# In ths plot, we have more microarray specific DE genes 
# this is due to the double filter applied to the data, the low expression for RNA seq and the informative for microarray
# the informative filter was more stringent, leaving only 5,000 genes for microarray
# some of those being DE in RNA seq
# On the other hand, low expression filtering left 12,000 genes for RNA seq
# most of those being DE in microarray (which is poorly sensitive to low expressed genes anyway)
# Therefore, microarray technology ends up favored in the above plot, with more microarra-DE genes left to plot
# and a lot of RNA seq DE genes discarded before the plot

# THe problem is that we only have DE data for 5,000 genes for the microarray platform
# and therefore we cannot consider more than those genes for the plot
# Otherwise, we would need to apply the eBayes correction on the unfiltered 24,000 genes microarray dataset
# which would alter the statistics and particularly the adjusted p-value.

# Have a look at the correlation (spearman for ranked correlation)
cor.test(x=summarised.micro.topTable[intersect.logFC,]$logFC, y=RNA_PBL_table[intersect.logFC,]$logFC, col=color.logFC, method="spearman")
# Have a look at the pearson correlation (fold change value correlation)
cor.test(x=summarised.micro.topTable[intersect.logFC,]$logFC, y=RNA_PBL_table[intersect.logFC,]$logFC, col=color.logFC, method="pearson")
# ranked correlation is a bit better
# In both cases, the correlation is really significant, meaning that the fold change for the genes
# that passed filtering is correlated between the two platforms.

##########
# Biomart query to convert affymetrix ids into ensembl ids
#########
# loads the library
library("biomaRt")
# Check which Biomart web services are available
listMarts()

# Loads the set of ensembl mart
mart = useMart(biomart="ensembl")
# Lookat the datasets in the mart
listDatasets(mart)
# redefine mart to the bovine information
mart = useMart(biomart="ensembl", dataset="btaurus_gene_ensembl")

# We look for filters allowing to search affymetrix ids
listFilters(mart)
# Looks for mention of Affy in the "description" column
grep(pattern="Affy", x=listFilters(mart)[,2]) # 12, 79
# What's on rows 12, and 79?
listFilters(mart)[c(12,79),]
# The affymetrix probe id is in field "affy_bovine"

# We look for attributes ensembl_id
listAttributes(mart) 
# Where is the ensembl gene id attribute?
grep(pattern="ensembl_gene_id", x=listAttributes(mart)[,1])
# 
listAttributes(mart)[grep(pattern="ensembl_gene_id", x=listAttributes(mart)[,1]),]
# The ensembl gene id is in field "ensembl_gene_id"

# where is the gene symbol field?
grep(pattern="symbol", x=listAttributes(mart)[,2])
# what's on row 36?
listAttributes(mart)[36,]

# Get the ensembl and affy probe id for all the informative microarray probe sets
getBM(attributes=c("ensembl_gene_id", "affy_bovine"), filters="affy_bovine", values=head(annotated.PBL.topTable,n=10)$ID, mart=mart, bmHeader=TRUE)

# Get the probe id, ensembl id, and gene symbol for the top 10 DE genes (adj. p-value)
getBM(attributes=c("affy_bovine", "ensembl_gene_id", "hgnc_symbol"), filters="affy_bovine", values=head(annotated.PBL.topTable,n=10)$ID, mart=mart, bmHeader=TRUE)





##########
# Correlation of expression between techologies
###########
length(intersect.logFC) # 2965 ensembl ids overlap between the filtered RNA-seq data and informative microarray data
length(unique(annotated.PBL.topTable$ensembl_id)) # 3,044 ensembl were in the microarray topTable (passed informative)
nrow(RNA_PBL_table) # 12,294 


##########
# Cell type deconvolution
#########
# Package for expression 
library(CellMix)
# Information about the markers
cellMarkersInfo()
# The Abbas marke looks good
m <- MarkerList("Abbas")
# Summary of the marker
summary(m)
barplot(m)
MarkerList(m)
details(MarkerList(m))
str(m)
# Gene ids of each cell type in the Abbas signature
geneIds(m)
# NetAffx file for best match between human probeset and bovine probeset
bestMatch = read.table(file="ProbeAnnotation/nfs/netaffx/arraycompare/data_dumps/U133PlusVsBovine_BestMatch.txt", header=T, sep="\t")
head(bestMatch)
head(bestMatch[,c("A.Probe.Set.Name","B.Probe.Set.Name")])
# How many of the probesets inthe human signature are present in the best match file?
geneIds(m) %in% bestMatch$A.Probe.Set.Name
lapply(X=geneIds(m), FUN=function(x){x %in% bestMatch$A.Probe.Set.Name })

# The NetAffx file has very few markers in it
# Let's use the annotation packages to try and convert the human probesets to 
# an cross reference identifiers, and then back to the bovine identifer
library(hgu133a.db, hgu133b.db)


lapply(X=geneIds(m), FUN=function(x){hgu133aUNIGENE[x]})
lapply(X=geneIds(m), FUN=function(x){hgu133bUNIGENE[x]})

hgu133bUNIGENE[geneIds(m)$neutro]




# Creates an ExpressionMix object from farms normalised microarray data
cellmix.eSet = ExpressionMix(farms)
cellmix.eSet
pData(cellmix.eSet) # pData works even on the ExpressionMix object, we have preserved group labels
featureNames(cellmix.eSet)




gedBlood(object=cellmix.eSet)

??.getMappingData()

