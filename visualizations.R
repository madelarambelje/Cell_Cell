#Visualizations: 

# Part I: Predict general principles of cell-cell communication
# Compare the total number of interactions and interaction strength
PLOT.compare.interactions <- function(x = cellchat) {
  gg1 <- compareInteractions(cellchat, 
                             show.legend = F, 
                             group = c(1,2)) 
  gg2 <- compareInteractions(cellchat, 
                             show.legend = F, 
                             group = c(1,2), 
                             measure = "weight")
  return (gg1 + gg2)
}

# Differentially expressed interactions circle plot
PLOT.compare.interactions.circle <- function(x = cellchat) {
  par(mfrow = c(1,2), xpd=T)
  netVisual_diffInteraction(x, weight.scale = T)
  netVisual_diffInteraction(x, weight.scale = T, measure = "weight")
}

# Heatmap: differential number of interactions or interaction strength
PLOT.compare.interactions.HM <- function(x=cellchat) {
  gg3 <- netVisual_heatmap(x)
  gg4 <- netVisual_heatmap(x, measure = "weight")
  return (gg3 + gg4)
}


# Part II: Identify the conserved and context-specific signaling pathways
#> Compute the distance of signaling networks between datasets 1 2
PLOT.rank.similarity <- function(x = cellchat) {
  par(mfrow = c(1,2), xpd=F)
  g <- rankSimilarity(x, type = "functional")
  g2 <- rankSimilarity(x, type = "structural")
  return (g+g2)
}

# Identify and visualize the conserved and context-specific signaling pathways
# Compare the overall information flow of each signaling pathway
PLOT.comparison.info.flow <- function(x = cellchat) {
  par(mfrow = c(1,2), xpd=F)
  gg5 <- rankNet(x, mode = "comparison", stacked = T, do.stat = TRUE)
  gg6 <- rankNet(x, mode = "comparison", stacked = F, do.stat = TRUE)
  return (gg5 + gg6)
}

# combining all the identified signaling pathways from different datasets 
# Outgoing signals
PLOT.outgoing.signals.HM<- function(i=1) {
  pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
  ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", 
                                          signaling = pathway.union, 
                                          title = names(object.list)[i], 
                                          width = 7, height = 10)
  ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", 
                                          signaling = pathway.union, 
                                          title = names(object.list)[i+1], 
                                          width = 7, height = 10)
  return (draw(ht1 + ht2, ht_gap = unit(0.5, "cm")))
}

# Incoming Signals
PLOT.incoming.signals.HM<- function(i=1) {
  pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
  ht3 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", 
                                          signaling = pathway.union, 
                                          title = names(object.list)[i], width = 7, 
                                          height = 10, color.heatmap = "GnBu")
  ht4 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", 
                                          signaling = pathway.union, 
                                          title = names(object.list)[i+1], width = 7, 
                                          height = 10, color.heatmap = "GnBu")
  return (draw(ht3 + ht4, ht_gap = unit(0.5, "cm")))
}

# Overall signaling patterns: 
PLOT.overall.signals.HM<- function(i=1) {
  pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
  ht5 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", 
                                          signaling = pathway.union, 
                                          title = names(object.list)[i], width = 7, 
                                          height = 10, color.heatmap = "OrRd")
  ht6 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", 
                                          signaling = pathway.union, 
                                          title = names(object.list)[i+1], width = 7, 
                                          height = 10, color.heatmap = "OrRd")
  return (draw(ht5 + ht6, ht_gap = unit(0.5, "cm")))
}


# Part III: Identify the upgulated and down-regulated signaling ligand-receptor pairs
#> Comparing communications on a merged object
PLOT.comm.prob.BP<- function(x = cellchat, sources.use = 5, targets.use = c(0:100),  
                                   comparison = c(1, 2)) {
  gg7 <- netVisual_bubble(x, sources.use = sources.use, targets.use,  
                          comparison = comparison, max.dataset = 2, 
                          title.name = "Increased signaling in inflamed", 
                          angle.x = 45, remove.isolate = T)
  gg8 <- netVisual_bubble(x, sources.use = sources.use, targets.use,  
                          comparison = comparison, max.dataset = 1,  
                          title.name = "Decreased signaling in inflamed", 
                          angle.x = 45, remove.isolate = T)
  #> Comparing communications on a merged object
  return (gg7 + gg8)
}


###
# Part IV: Visually compare cell-cell communication using 
# Hierarchy plot, Circle plot or Chord diagram
# PATHWAYS
###
# Circleplot
PLOT.pathway.circle<- function(i = 1, pathways.show = pathways.show) {
  weight.max <- getMaxWeight(object.list, 
                             slot.name = c("netP"), 
                             attribute = pathways.show) # control the edge weights across different datasets
  par(mfrow = c(1,2), xpd=TRUE)
  #> Comparing communications on a merged object
  return (for (i in 1:length(object.list)) {
    netVisual_aggregate(object.list[[i]], signaling = pathways.show, 
                        layout = "circle", edge.weight.max = weight.max[1], 
                        edge.width.max = 10, signaling.name = paste(pathways.show, 
                                                                    names(object.list)[i]))
  })
}

# Heatmap
PLOT.pathway.HM <- function(i = 1, pathways.show = pathways.show) {
  par(mfrow = c(1,2), xpd=TRUE)
  ht <- list()
  for (i in 1:length(object.list)) {
    ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, 
                                 color.heatmap = "Reds",
                                 title.name = paste(pathways.show, "signaling ",
                                                    names(object.list)[i]))
  }
  return (ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm")))
}

# Chord diagram
PLOT.pathway.chord <- function(i = 1, pathways.show = pathways.show) {
  par(mfrow = c(1,2), xpd=TRUE)
  return (for (i in 1:length(object.list)) {
    netVisual_aggregate(object.list[[i]], signaling = pathways.show, 
                        layout = "chord", signaling.name = paste(pathways.show, 
                                                                 names(object.list)[i]))
  })
}


# Part V: Compare the signaling gene expression distribution between different datasets
# plot the gene expression distribution of signaling genes related to L-R pairs or 
# signaling pathway using a Seurat wrapper function
####
PLOT.expression.distribution <- function(x = cellchat, pathways.show = pathways.show) {
  cellchat@meta$datasets = factor(x@meta$datasets, 
                                  levels = c(names(object.list)[1], 
                                             names(object.list)[2])) # set factor level
  return (plotGeneExpression(x, signaling = pathways.show, 
                             split.by = "datasets", colors.ggplot = T))
}

# End of visualizations
