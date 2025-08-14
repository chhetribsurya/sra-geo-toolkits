#!/usr/bin/env Rscript


############################################################
# Install and load GEOquery
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("GEOquery")
############################################################


library(GEOquery)

# Download GEO dataset GSE188486
gse <- getGEO("GSE188486", GSEMatrix = TRUE)

# If the dataset has more than one series, select the first one
if (length(gse) > 1) gse <- gse[[1]]

# Extract metadata
metadata <- pData(phenoData(gse))
print(metadata)

# Extract expression data
exprData <- exprs(gse[[1]])
head(exprData)

# Download all supplementary files associated with GSE188486
sapply(gse@supplementary_file, function(x) {
  download.file(x, destfile = basename(x))
})
