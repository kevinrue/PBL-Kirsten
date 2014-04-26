source("http://bioconductor.org/biocLite.R")
biocLite("GEOquery")


# load series and platform data from GEO
library(GEOquery)
# Command to directly download the PBL dataset from the GEO website
#gset <- getGEO("GSE33359", GSEMatrix =TRUE)

gset <- getGEO(filename="C:/Users/Kevin RUE/Downloads/GSE33359_RAW.tar", GSEMatrix =TRUE)