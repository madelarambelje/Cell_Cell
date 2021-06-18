# Analyzing 1 seurat sub-dataset at a time: 
# Nick Veltmaat
# 21 - 5 - 2021
library(CellChat)
library(NMF)
library(ggalluvial)
library(Seurat)
setwd('~/Desktop/Master Files/Cell_Cell_Project/Cellchat')

## SETTINGS: 
# Set tutorial data you want to load: 
tutorial_data <- "CellPhoneDB" # CellChat, CellPhoneDB or ICELLNET

# Set Permutated Data to TRUE or FALSE:
permutated_dataset <- FALSE


# Read RDS seurat object
if (permutated_dataset == TRUE) {
  print(paste("Loading the Permutated dataset of the ", tutorial_data, " tutorial", sep = ""))
  seurat.object <- readRDS(paste("Tutorial_Datasets/data/Permutation",tutorial_data,"/PermutatedCombined", tutorial_data,".rds", sep=""))
} else if (permutated_dataset == FALSE){
  print(paste("Loading Normal (non-permutated) dataset of the ", tutorial_data, " tutorial", sep = ""))
  seurat.object <- readRDS(paste("Tutorial_Datasets/data/", tutorial_data, "_TutorialData_SeuratObject.rds", sep=""))
}

# Normalize ICELLNET data and rename clusters
if (tutorial_data == "ICELLNET") {
  seurat.object <- NormalizeData(seurat.object)
  seurat.object <- ScaleData(seurat.object)
  # Rename Cell Cluster Names
  celltypes_list <- list("CB0" = "Act B-Cells", "CB1" = "Plasma B-Cells", "CB2a" = "Naive B-Cells", "CB2b" = "pDCs", 
                         "CB3" = "ISG-high B-Cells", "CD0" = "Dividing Cells", "CE0" = "Epithelial Cells", 
                         "CM0" = "Inflam. CD16+ Macrophage", "CM1" = "Phagocytic CD16+ Macrophage", 
                         "CM2" = "Tissue-resident macrophage", "CM3" = "cDCs", "CM4" = "M2-like CD16+ Macrophage", 
                         "CT0a" = "Eff Mem CD4+ T-Cell", "CT0b" = "Central Mem CD4+ T-Cells", "CT1" = "CD56dim CD16+ NK", 
                         "CT2" = "CTLs","CT3a" = "Tregs", "CT3b" = "TFG-Like Cells", "CT4" = "GZMK+ CD8+ Cells", 
                         "CT5a" = "Res Mem CD8+ T-Cells", "CT5b" = "CD56bright CD16- NK", "CT6" = "ISGhigh CD4+ T-Cells")
  
  for (i in 1:length(seurat.object@meta.data$author_annotation)){
    oldval <- seurat.object@meta.data$author_annotation[i]
    newval <- celltypes_list[[oldval]]
    seurat.object@meta.data$author_annotation[i] <- newval
  }
  
}


# Load CellChat full DataBase
CellChatDB <- CellChatDB.human
#showDatabaseCategory(CellChatDB)

# Create CellChat object from Seurat
if (tutorial_data == "CellChat") {
  cellchat <- createCellChat(object = seurat.object, group.by = "Cluster")
  CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
} else if (tutorial_data == "ICELLNET") {
  cellchat <- createCellChat(object = seurat.object, group.by = "author_annotation")
  CellChatDB.use <- CellChatDB # simply use the default CellChatDB
} else {
  cellchat <- createCellChat(object = seurat.object, group.by = "Cluster")
  CellChatDB.use <- CellChatDB # simply use the default CellChatDB
}
cellchat@DB <- CellChatDB.use

# Preprocessing the expression data for cell-cell communication analysis
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 16) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# project gene expression data onto PPI network (optional)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 1)
cellchat <- computeCommunProbPathway(cellchat, )
cellchat <- aggregateNet(cellchat)






###### 
#VISUALIZATIONS for CellChat tutorial Data: (For CellChatDB or ICELLNET run other scripts)
#

#pdf(paste(tutorial_data, "_TutorialData_CellChat.pdf", sep = ""), width = 14, height = 8) 
pdf(paste(tutorial_data, "_Permutated_Combined_TutorialData_CellChat.pdf", sep = ""), width = 14, height = 8) 


## PART I
# visualize the aggregated cell-cell communication network
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

netVisual_heatmap(object = cellchat, color.heatmap = c('white', 'red'))


cellchat@netP$pathways
# Part III: Visualization of cell-cell communication network
pathways.show <- c("CXCL") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)

# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
#> Note: The first link end is drawn out of sector 'Inflam. FIB'.
# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
#> Do heatmap based on a single object


# Chord diagram
group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4)) # grouping cell clusters into fibroblast, DC and TC cells
names(group.cellType) <- levels(cellchat@idents)
netVisual_chord_cell(cellchat, signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network"))

netAnalysis_contribution(cellchat, signaling = pathways.show)

pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)

# Circle plot
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
# Chord diagram
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")

# Visualize cell-cell communication mediated by multiple ligand-receptors or signaling pathways
# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), remove.isolate = FALSE)

# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), signaling = c("CCL","CXCL"), remove.isolate = FALSE)

# show all the significant interactions (L-R pairs) based on user's input (defined by `pairLR.use`)
pairLR.use <- extractEnrichedLR(cellchat, signaling = c("CCL","CXCL","FGF"))
netVisual_bubble(cellchat, sources.use = c(3,4), targets.use = c(5:8), pairLR.use = pairLR.use, remove.isolate = TRUE)

# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
# show all the interactions sending from Inflam.FIB
netVisual_chord_gene(cellchat, sources.use = 4, targets.use = c(5:11), lab.cex = 0.5,legend.pos.y = 30)
# show all the interactions received by Inflam.DC
netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = 8, legend.pos.x = 15)
# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = c(5:11), signaling = c("CCL","CXCL"),legend.pos.x = 8)
# show all the significant signaling pathways from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = c(5:11), slot.name = "netP", legend.pos.x = 10)


# Plot the signaling gene expression distribution using violin/dot plot
plotGeneExpression(cellchat, signaling = "CXCL")
plotGeneExpression(cellchat, signaling = "CXCL", enriched.only = FALSE)


# Part IV: Systems analysis of cell-cell communication network
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

# Visualize the dominant senders (sources) and receivers (targets) in a 2D space
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("CXCL", "CCL"))
#> Signaling role analysis on the cell-cell communication network from user's input
gg1 + gg2

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2

# Signaling role analysis on the cell-cell communication networks of interest
ht <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("CXCL", "CCL"))
ht


# Identify and visualize outgoing communication pattern of secreting cells
# Here we run selectK to infer the number of patterns.
selectK(cellchat, pattern = "outgoing")
nPatterns = 3
plot.new()
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
# river plot
netAnalysis_river(cellchat, pattern = "outgoing")
# dot plot
netAnalysis_dot(cellchat, pattern = "outgoing")

# Identify and visualize incoming communication pattern of target cells
selectK(cellchat, pattern = "incoming")
nPatterns = 4
plot.new()
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
# river plot
netAnalysis_river(cellchat, pattern = "incoming")
# dot plot
netAnalysis_dot(cellchat, pattern = "incoming")

# Identify signaling groups based on their functional similarity
cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "functional")
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)
netVisual_embeddingZoomIn(cellchat, type = "functional", nCol = 2)

# Identify signaling groups based on structure similarity
cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "structural")
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "structural", label.size = 3.5)
netVisual_embeddingZoomIn(cellchat, type = "structural", nCol = 2)

dev.off()
######

