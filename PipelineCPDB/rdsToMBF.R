#!/usr/bin/env Rscript

###############################################################################
# Author: Andre de la Rambelje                                                #
# Version: 1                                                                  #
# Script for extracting sparse matrix from RDS with the barcodes and features #
# for CellPhoneDB                                                             #
###############################################################################

library(Matrix)
library("Seurat")
# Creating input for script
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

print(paste0("Working in DIR: ", args[1]))
print(paste0("RDS file selected: ", args[2]))
# Extracting sparse matrix from RDS
subset.to.tear <- readRDS(paste0(args[2]))
counts.to.save <- GetAssayData(subset.to.tear,slot = "counts", assay = "RNA")
writeMM(counts.to.save, paste0(args[1],"/","matrix.mtx"))

# Extracting barcodes & Features from RDS
barcodes.to.save <- colnames(counts.to.save)
colnames(barcodes.to.save) <- NULL
write.table(barcodes.to.save, paste0(args[1],"/","barcodes.tsv"), sep = "\t", row.names = FALSE, col.names = F)
features.to.save <- rownames(counts.to.save)
colnames(features.to.save) <- NULL
write.table(features.to.save, paste0(args[1],"/","features.tsv"), sep = "\t", row.names = FALSE, col.names = F)

# Saving metadata
meta.data.to.save <- subset.to.tear@meta.data
meta.data.to.save$Barcodes <- rownames(meta.data.to.save)
meta.B.C <- meta.data.to.save[,c("Barcodes","Cluster")]
colnames(meta.B.C) <- NULL
write.csv(meta.B.C, paste0(args[1],"/","metadata.csv"), row.names = F)

