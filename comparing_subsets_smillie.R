#Comparing the smillie subsets
#Nick Veltmaat
#3-5-2021
pacman::p_load(ComplexHeatmap, circlize, NMF, CellChat, Matrix, Seurat)
#BiocManager::install('UMAP')
setwd('~/Desktop/Master Files/Cell_Cell_Project/Cellchat')

#OPTIONAL: load files and compute interactions manually to create cellchat objects: 
###------
#Load seurat objects / subsets
inflamed_natural <- readRDS('smilliedata/inflamed_natural_composition.rds')
healthy_natural <- readRDS('smilliedata/natural_composition_subset.rds')

#Convert seurat to Cellchat Objects
cellchat_inflamed <- createCellChat(object = inflamed_natural, group.by = "ident")
cellchat_healthy <- createCellChat(object = healthy_natural, group.by = "ident")

#Healthy subset:
cellchat_healthy <- setIdent(cellchat_healthy, ident.use = "Cluster") # set "labels" as default cell identity
levels(cellchat_healthy@idents) # show factor levels of the cell labels
groupSize_healthy <- as.numeric(table(cellchat_healthy@idents)) # number of cells in each cell group
#Inflamed subset:
cellchat_inflamed <- setIdent(cellchat_inflamed, ident.use = "Cluster") # set "labels" as default cell identity
levels(cellchat_inflamed@idents) # show factor levels of the cell labels
groupSize_inflamed <- as.numeric(table(cellchat_inflamed@idents)) # number of cells in each cell group

#Set database
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
#showDatabaseCategory(CellChatDB)
#dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
cellchat_healthy@DB <- CellChatDB.use    #Set DB for Healthy Subset
cellchat_inflamed@DB <- CellChatDB.use   #Set DB for Inflamed Subset

#Compute Interactions:
#### Healthy subset:
cellchat_healthy <- subsetData(cellchat_healthy) # subset the expression data of signaling genes for saving computation cost
future::plan("multiprocess", workers = 4) # do parallel
cellchat_healthy <- identifyOverExpressedGenes(cellchat_healthy)
cellchat_healthy <- identifyOverExpressedInteractions(cellchat_healthy)
cellchat_healthy <- projectData(cellchat_healthy, PPI.human)
cellchat_healthy <- computeCommunProb(cellchat_healthy, raw.use = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_healthy <- filterCommunication(cellchat_healthy, min.cells = 10)
cellchat_healthy <- computeCommunProbPathway(cellchat_healthy)
cellchat_healthy <- aggregateNet(cellchat_healthy)
cellchat_healthy <- netAnalysis_computeCentrality(cellchat_healthy)
saveRDS(cellchat_healthy, file = "cellchat_healthy.rds")
#Compute Interactions:
#### Healthy subset:
cellchat_inflamed <- subsetData(cellchat_inflamed) # subset the expression data of signaling genes for saving computation cost
future::plan("multiprocess", workers = 4) # do parallel
cellchat_inflamed <- identifyOverExpressedGenes(cellchat_inflamed)
cellchat_inflamed <- identifyOverExpressedInteractions(cellchat_inflamed)
cellchat_inflamed <- projectData(cellchat_inflamed, PPI.human)
cellchat_inflamed <- computeCommunProb(cellchat_inflamed, raw.use = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_inflamed <- filterCommunication(cellchat_inflamed, min.cells = 10)
cellchat_inflamed <- computeCommunProbPathway(cellchat_inflamed)
cellchat_inflamed <- aggregateNet(cellchat_inflamed)
cellchat_inflamed <- netAnalysis_computeCentrality(cellchat_inflamed)
saveRDS(cellchat_inflamed, file = "cellchat_inflamed.rds")
###--------


##OR LOAD CELLCHAT OBJECTS FROM RDS:
cellchat_healthy <- readRDS("cellchat_healthy.rds")
cellchat_inflamed <- readRDS("cellchat_inflamed.rds")

### Merge Cellchat objects
object.list <- list(healthy = cellchat_healthy, 
                    inflamed = cellchat_inflamed)
#  object.list
cellchat <- mergeCellChat(object.list, 
                          add.names = names(object.list))
cellchat


###############
#VISUALIZATIONS

#Compare the total number of interactions and interaction strength
gg1 <- compareInteractions(cellchat, 
                           show.legend = F, 
                           group = c(1,2)) 
gg2 <- compareInteractions(cellchat, 
                           show.legend = F, 
                           group = c(1,2), 
                           measure = "weight")
gg1 + gg2

#Differentially expressed interactions circle plot
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

#Heatmap: differential number of interactions or interaction strength
gg3 <- netVisual_heatmap(cellchat)
gg4 <- netVisual_heatmap(cellchat, measure = "weight")
gg3 + gg4

# Identify signaling networks with larger (or less) difference as well as 
# signaling groups based on their functional similarity
#Reset cellchat object
#object.list <- list(healthy = cellchat_healthy, inflamed = cellchat_inflamed)
#cellchat <- mergeCellChat(object.list, add.names = names(object.list))
#
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
#> 2D visualization of signaling networks from datasets 1 2
netVisual_embeddingPairwiseZoomIn(cellchat, type = "functional", nCol = 2)
#> 2D visualization of signaling networks from datasets 1 2

# Identify signaling groups based on structure similarity
cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "structural")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "structural")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)
#> 2D visualization of signaling networks from datasets 1 2
netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2)
#> 2D visualization of signaling networks from datasets 1 2

rankSimilarity(cellchat, type = "functional")
rankSimilarity(cellchat, type = "structural")
#> Compute the distance of signaling networks between datasets 1 2

#Identify and visualize the conserved and context-specific signaling pathways
#Compare the overall information flow of each signaling pathway
gg5 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg6 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg5 + gg6

i = 1
# combining all the identified signaling pathways from different datasets 
#Outgoing signals
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
#Incoming Signals
ht3 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "GnBu")
ht4 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "GnBu")
draw(ht3 + ht4, ht_gap = unit(0.5, "cm"))
#Overall signaling patterns: 
ht5 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "OrRd")
ht6 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "OrRd")
draw(ht5 + ht6, ht_gap = unit(0.5, "cm"))

#Part III: Identify the upgulated and down-regulated signaling ligand-receptor pairs
#Identify dysfunctional signaling by comparing the communication probabities
#Reset cellchat object
#object.list <- list(healthy = cellchat_healthy, inflamed = cellchat_inflamed)
#cellchat <- mergeCellChat(object.list, add.names = names(object.list))
netVisual_bubble(cellchat, sources.use =2, targets.use = c(0:1000),  
                 comparison = c(1, 2), angle.x = 45)
#> Comparing communications on a merged object

gg7 <- netVisual_bubble(cellchat, sources.use = 5, targets.use = c(0:100),  
                        comparison = c(1, 2), max.dataset = 2, 
                        title.name = "Increased signaling in inflamed", 
                        angle.x = 45, remove.isolate = T)
gg8 <- netVisual_bubble(cellchat, sources.use = 5, targets.use = c(0:100),  
                        comparison = c(1, 2), max.dataset = 1, 
                        title.name = "Decreased signaling in inflamed", 
                        angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg7 + gg8


###
#Part IV: Visually compare cell-cell communication using Hierarchy plot, Circle plot or Chord diagram
#PATHWAYS
###
cellchat@netP[["healthy"]][["pathways"]]
cellchat@netP[["inflamed"]][["pathways"]]
pathways.show <- c("LAMININ") 

#Circleplot
weight.max <- getMaxWeight(object.list, 
                           slot.name = c("netP"), 
                           attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, 
                      layout = "circle", edge.weight.max = weight.max[1], 
                      edge.width.max = 10, signaling.name = paste(pathways.show, 
                                                                  names(object.list)[i]))
}
#Heatmap
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
# Chord diagram
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}


#Part V: Compare the signaling gene expression distribution between different datasets
#plot the gene expression distribution of signaling genes related to L-R pairs or signaling pathway using a Seurat wrapper function
####
cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("healthy", "inflamed")) # set factor level
plotGeneExpression(cellchat, signaling = pathways.show, split.by = "datasets", colors.ggplot = T)




