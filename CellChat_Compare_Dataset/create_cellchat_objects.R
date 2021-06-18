# Create CellChat objects from given subsets:
# Nick Veltmaat
# 16-5-2021

# Convert Seurat to Cellchat Objects
cellchat1 <- createCellChat(object = subset1, group.by = "Cluster") # E.g. Healthy subset
cellchat2 <- createCellChat(object = subset2, group.by = "Cluster") # E.g. Inflamed subset

# Subset 1 (Healthy subset):
cellchat1 <- setIdent(cellchat1, ident.use = "Cluster") # set "labels" as default cell identity
levels(cellchat1@idents) # show factor levels of the cell labels
groupSize_subset1 <- as.numeric(table(cellchat1@idents)) # number of cells in each cell group
# Subset 2 (Inflamed subset):
cellchat2 <- setIdent(cellchat2, ident.use = "Cluster") # set "labels" as default cell identity
levels(cellchat2@idents) # show factor levels of the cell labels
groupSize_subset2 <- as.numeric(table(cellchat2@idents)) # number of cells in each cell group

# RUN DATABASE.R to load either the original CellChatDB or the ICELLNET based DB
#source('database.R')

# Compute Interactions:
#### (Healthy) subset 1:
cellchat1 <- subsetData(cellchat1) # subset the expression data of signaling genes for saving computation cost
future::plan("multiprocess", workers = 8) # do parallel
cellchat1 <- identifyOverExpressedGenes(cellchat1)
cellchat1 <- identifyOverExpressedInteractions(cellchat1)
cellchat1 <- projectData(cellchat1, PPI.human)
cellchat1 <- computeCommunProb(cellchat1, raw.use = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat1 <- filterCommunication(cellchat1, min.cells = 10)
cellchat1 <- computeCommunProbPathway(cellchat1)
cellchat1 <- aggregateNet(cellchat1)
cellchat1 <- netAnalysis_computeCentrality(cellchat1)
saveRDS(cellchat1, file = "cellchat1.rds")
# Compute Interactions:
#### (Inflamed) subset 2:
cellchat2 <- subsetData(cellchat2) # subset the expression data of signaling genes for saving computation cost
future::plan("multiprocess", workers = 8) # do parallel
cellchat2 <- identifyOverExpressedGenes(cellchat2)
cellchat2 <- identifyOverExpressedInteractions(cellchat2)
cellchat2 <- projectData(cellchat2, PPI.human)
cellchat2 <- computeCommunProb(cellchat2, raw.use = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat2 <- filterCommunication(cellchat2, min.cells = 10)
cellchat2 <- computeCommunProbPathway(cellchat2)
cellchat2 <- aggregateNet(cellchat2)
cellchat2 <- netAnalysis_computeCentrality(cellchat2)
saveRDS(cellchat2, file = "cellchat2.rds")
