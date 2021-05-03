#Comparing the smillie subsets
#Nick Veltmaat
#3-5-2021
pacman::p_load(ComplexHeatmap, circlize, NMF, CellChat, Matrix, Seurat, UMAP)
BiocManager::install('UMAP')
nsetwd('~/Desktop/Master Files/Cell_Cell_Project/Cellchat')

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

### Merge Cellchat objects
object.list <- list(healthy = cellchat_healthy, inflamed = cellchat_inflamed)
#  object.list
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat

###############
#VISUALIZATIONS

#Compare the total number of interactions and interaction strength
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2)) #or with cellchat_merged
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

#Differentially expressed interactions circle plot
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

#Heatmap: differential number of interactions or interaction strength
gg3 <- netVisual_heatmap(cellchat)
gg4 <- netVisual_heatmap(cellchat, measure = "weight")
gg3 + gg4



cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "functional")
#----------
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
#> 2D visualization of signaling networks from datasets 1 2
#-----




