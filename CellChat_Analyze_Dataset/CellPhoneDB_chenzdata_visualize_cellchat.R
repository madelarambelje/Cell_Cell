
###### 
#VISUALIZATIONS
#
pdf(paste(tutorial_data, "_TutorialData_CellChat.pdf", sep = ""), width = 14, height = 8) 
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

gg1 <- netVisual_heatmap(cellchat, color.heatmap = 'Reds')
gg1


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
group.cellType <- c(rep("Macrophages", 5), rep("T/NK's", 10), rep("B / Ep Cells", 7)) # grouping cell clusters into fibroblast, DC and TC cells
names(group.cellType) <- levels(cellchat@idents)
netVisual_chord_cell(cellchat, signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network"))

netAnalysis_contribution(cellchat, signaling = pathways.show)

pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[5,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
# Circle plot
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
# Chord diagram
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")

# Visualize cell-cell communication mediated by multiple ligand-receptors or signaling pathways
# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')

# show all the significant interactions (L-R pairs) associated with certain signaling pathways
bb1 <- netVisual_bubble(cellchat, sources.use = 7, targets.use = c(1:8), signaling = c("CCL","CXCL"), remove.isolate = FALSE)
bb2 <- netVisual_bubble(cellchat, sources.use = c(1:8), targets.use = 7, signaling = c("CCL","CXCL"), remove.isolate = FALSE)
bb1+bb2


bb3 <- netVisual_bubble(cellchat, sources.use = 7, targets.use = c(1:8), signaling = c("VEGF","FGF", "IGF","PDGF", "HGF"), remove.isolate = FALSE)
bb4 <- netVisual_bubble(cellchat, sources.use = c(1:8), targets.use = 7, signaling = c("VEGF","FGF", "IGF","PDGF", "HGF"), remove.isolate = FALSE)
bb3+bb4
# show all the significant interactions (L-R pairs) based on user's input (defined by `pairLR.use`)
#pairLR.use <- extractEnrichedLR(cellchat, signaling = c("CCL","CXCL","TNF", "VEGF","FGF", "IGF","PDGF", "HGF"))
#netVisual_bubble(cellchat, sources.use = 7, targets.use = c(1:8), pairLR.use = pairLR.use, remove.isolate = TRUE)

# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_chord_gene(cellchat, sources.use = 7, targets.use = c(1:8), signaling = c("CCL","CXCL"),
                     lab.cex = 0.5,legend.pos.y = 30)
netVisual_chord_gene(cellchat, sources.use = 7, targets.use = c(1:8), signaling = c("VEGF","FGF", "IGF","PDGF", "HGF"),
                     lab.cex = 0.5,legend.pos.y = 30)

# show all the interactions received by Inflam.DC


# Plot the signaling gene expression distribution using violin/dot plot
plotGeneExpression(cellchat, signaling = c("CCL","CXCL"))
plotGeneExpression(cellchat, signaling = c("CCL","CXCL"), enriched.only = FALSE)


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
cellchat <- netClustering(cellchat, type = "functional", k = 4)
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)
netVisual_embeddingZoomIn(cellchat, type = "functional", nCol = 2)

# Identify signaling groups based on structure similarity
cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "structural", k = 7)
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "structural", label.size = 3.5)
netVisual_embeddingZoomIn(cellchat, type = "structural", nCol = 2)


dev.off()
######



#netVisual_embedding(cellchat, type = "structural", label.size = 3.5)
