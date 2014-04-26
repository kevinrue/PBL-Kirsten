# Understand CellMix Abbas cell signature
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

# Create a table of two columns: Cell type, Human probeset
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

# Create a table of two columns: human probeset, gene symbol
library(hgu133a.db)
library(hgu133b.db)
probesetsA = toTable(hgu133aSYMBOL)$probe_id
probesetsB = toTable(hgu133bSYMBOL)$probe_id
humanProbeset = c()
geneSymbol = c()
# For each human probeset
for(hps in mapping$humanProbeset){
  print(hps)
  # if the probeset is associated with a gene symbol in the hgu133a array
  if(hps %in% probesetsA){
    print("A")w
    for(symbol in toTable(hgu133aSYMBOL[hps])$symbol){
      humanProbeset = c(humanProbeset, hps)
      geneSymbol = c(geneSymbol, symbol)
    }
  }
  # if the probeset is associated with a gene symbol in the hgu133b array
  else if(hps %in% probesetsB){
    print("B")
    for(symbol in toTable(hgu133bSYMBOL)$symbol){
      humanProbeset = c(humanProbeset, hps)
      geneSymbol = c(geneSymbol, symbol)
    }
  }
  else{
    print("Nope")
    symbol=NULL
    # the probeset is not associated with any gene symbol in any of the human arrays
  }
  # Clears the temporary variable
  rm(hps, symbol)
}

rm(symbolsA, symbolsB)