#Comparing the smillie subsets one by one
#Nick Veltmaat
#25-4-2021
pacman::p_load(ComplexHeatmap, circlize, NMF, CellChat, Matrix, Seurat)
setwd('~/Desktop/Master Files/Cell_Cell_Project/Cellchat')
#Load seurat objects / subsets
t_cell <- readRDS('smilliedata/t_cell_subset.rds')
natural_composition <- readRDS('smilliedata/natural_composition_subset.rds')
uniform <- readRDS('smilliedata/uniform_subset.rds')

#Convert seurat to Cellchat Objects
cellchat_t_cell <- createCellChat(object = t_cell, group.by = "ident")
cellchat_natural_comp <- createCellChat(object = natural_composition, group.by = "ident")
cellchat_uniform <- createCellChat(object = uniform, group.by = "ident")
#cellchat.list <- list(
  #cellchat_t_cell, 
  #cellchat_natural_comp, cellchat_uniform)
#cellchat.list

#Natural composition subset:
cellchat_natural_comp <- setIdent(cellchat_natural_comp, ident.use = "Cluster") # set "labels" as default cell identity
levels(cellchat_natural_comp@idents) # show factor levels of the cell labels
groupSize_nat <- as.numeric(table(cellchat_natural_comp@idents)) # number of cells in each cell group

#Uniform subset:
cellchat_uniform <- setIdent(cellchat_uniform, ident.use = "Cluster") # set "labels" as default cell identity
levels(cellchat_uniform@idents) # show factor levels of the cell labels
groupSize_uni <- as.numeric(table(cellchat_uniform@idents)) # number of cells in each cell group

#Set database
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
#showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
cellchat_natural_comp@DB <- CellChatDB.use
cellchat_uniform@DB <- CellChatDB.use

#Natural composition:
cellchat_natural_comp <- subsetData(cellchat_natural_comp) # subset the expression data of signaling genes for saving computation cost
future::plan("multiprocess", workers = 4) # do parallel
cellchat_natural_comp <- identifyOverExpressedGenes(cellchat_natural_comp)
cellchat_natural_comp <- identifyOverExpressedInteractions(cellchat_natural_comp)
cellchat_natural_comp <- projectData(cellchat_natural_comp, PPI.human)
cellchat_natural_comp <- computeCommunProb(cellchat_natural_comp, raw.use = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_natural_comp <- filterCommunication(cellchat_natural_comp, min.cells = 10)
cellchat_natural_comp <- computeCommunProbPathway(cellchat_natural_comp)
cellchat_natural_comp <- aggregateNet(cellchat_natural_comp)

#Uniform Distribution: 
cellchat_uniform <- subsetData(cellchat_uniform) # subset the expression data of signaling genes for saving computation cost
future::plan("multiprocess", workers = 4) # do parallel
cellchat_uniform <- identifyOverExpressedGenes(cellchat_uniform)
cellchat_uniform <- identifyOverExpressedInteractions(cellchat_uniform)
cellchat_uniform <- projectData(cellchat_uniform, PPI.human)
cellchat_uniform <- computeCommunProb(cellchat_uniform, raw.use = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_uniform <- filterCommunication(cellchat_uniform, min.cells = 10)
cellchat_uniform <- computeCommunProbPathway(cellchat_uniform)
cellchat_uniform <- aggregateNet(cellchat_uniform)


###########################################################
#### IMPORTANT ###
#Switch out 'cellchat_natural_comp' or 'cellchat_uniform' for the rest of the script
cellchat <- cellchat_natural_comp
###########################################################

#Visualizing interactions natural comp
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

#Check cell interactions group by group
mat <- cellchat@net$weight #Change subset here
par(mfrow = c(2,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

#PATHWAY VISUALISATIONS
unique(cellchat@netP$pathways) #Select a pathway from this list

pathways.show <- c("MK") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,3) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
#Circle plot:
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
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("MK", "CD99"))
#> Signaling role analysis on the cell-cell communication network from user's input
gg1 + gg2


# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2

# Identify signaling groups based on their functional similarity
cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)


##########
##########
##########
##########
##########
##########

#Inflamed vs Healthy
# 2-5-2021
# Nick Veltmaat


#Load seurat objects / subsets
healthy_seurat <- readRDS('smilliedata/healthy_natural_composition.rds')
inflamed_seurat <- readRDS('smilliedata/inflamed_natural_composition.rds')

#Creating Cellchat objects
healthy_cellchat <- createCellChat(object = healthy_seurat, group.by = "ident")
inflamed_cellchat <- createCellChat(object = inflamed_seurat, group.by = "ident")

#Healthy = Natural comp
healthy_cellchat <- setIdent(healthy_cellchat, ident.use = "Cluster") # set "labels" as default cell identity
levels(healthy_cellchat@idents) # show factor levels of the cell labels
groupSize_healthy <- as.numeric(table(healthy_cellchat@idents)) # number of cells in each cell group

#Inflamed = Uniform comp:
inflamed_cellchat <- setIdent(inflamed_cellchat, ident.use = "Cluster") # set "labels" as default cell identity
levels(inflamed_cellchat@idents) # show factor levels of the cell labels
groupSize_inflamed <- as.numeric(table(inflamed_cellchat@idents)) # number of cells in each cell group



#Set database
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
#showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
healthy_cellchat@DB <- CellChatDB.use
inflamed_cellchat@DB <- CellChatDB.use

#Natural composition:
healthy_cellchat <- subsetData(healthy_cellchat) # subset the expression data of signaling genes for saving computation cost
future::plan("multiprocess", workers = 4) # do parallel
healthy_cellchat <- identifyOverExpressedGenes(healthy_cellchat)
healthy_cellchat <- identifyOverExpressedInteractions(healthy_cellchat)
healthy_cellchat <- projectData(healthy_cellchat, PPI.human)
healthy_cellchat <- computeCommunProb(healthy_cellchat, raw.use = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
healthy_cellchat <- filterCommunication(healthy_cellchat, min.cells = 10)
healthy_cellchat <- computeCommunProbPathway(healthy_cellchat)
healthy_cellchat <- aggregateNet(healthy_cellchat)


#Uniform Distribution: 
inflamed_cellchat <- subsetData(inflamed_cellchat) # subset the expression data of signaling genes for saving computation cost
future::plan("multiprocess", workers = 4) # do parallel
inflamed_cellchat <- identifyOverExpressedGenes(inflamed_cellchat)
inflamed_cellchat <- identifyOverExpressedInteractions(inflamed_cellchat)
inflamed_cellchat <- projectData(inflamed_cellchat, PPI.human)
inflamed_cellchat <- computeCommunProb(inflamed_cellchat, raw.use = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
inflamed_cellchat <- filterCommunication(inflamed_cellchat, min.cells = 10)
inflamed_cellchat <- computeCommunProbPathway(inflamed_cellchat)
inflamed_cellchat <- aggregateNet(inflamed_cellchat)



#### IMPORTANT ###
#Switch out 'cellchat_natural_comp' or 'cellchat_uniform' for the rest of the script
cellchat <- inflamed_cellchat
###########################################################

#Visualizing interactions natural comp
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

#Check cell interactions group by group
mat <- cellchat@net$weight #Change subset here
par(mfrow = c(2,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

#PATHWAY VISUALISATIONS
unique(cellchat@netP$pathways) #Select a pathway from this list

pathways.show <- c("MK") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,3) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
#Circle plot:
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
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("MK", "CD99"))
#> Signaling role analysis on the cell-cell communication network from user's input
gg1 + gg2


# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2

# Identify signaling groups based on their functional similarity
cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)

