#cellchat 1
pacman::p_load(ComplexHeatmap, circlize, NMF, CellChat, Matrix)
setwd('~/Desktop/Master Files/Cell_Cell_Project/Cellchat')

# Here we load a scRNA-seq data matrix and its associated cell meta data
#--test
#load(url("https://ndownloader.figshare.com/files/25950872")) # This is a combined data from two biological conditions: normal and diseases
#testdata.input = data_humanSkin$data # normalized data matrix

#Meta
cell_subsets = read.table('smilliedata/cell_subsets.txt', sep='\t', header=F, stringsAsFactors=F)
meta = read.table('smilliedata/all.meta2.txt', sep='\t', header=T, row.names=1, stringsAsFactors=F)
meta = na.omit(meta[meta$Location == 'Epi',])

#Matrix
epi.counts = as(readMM('smilliedata/gene_sorted-Epi.matrix.mtx'), "dgCMatrix")
rownames(epi.counts) = readLines('smilliedata/Epi.genes.tsv')
colnames(epi.counts) = readLines('smilliedata/Epi.barcodes2.tsv')

#--test
#testmeta = data_humanSkin$meta # a dataframe with rownames containing cell meta data
#testcell.use = rownames(testmeta)[testmeta$condition == "LS"] # extract the cell names from disease data
#testcell.use

#subset inflamed
cell.use = rownames(meta)[meta$Health == "Inflamed"] # extract the cell names from disease data
cell.use

# Prepare input data for CelChat analysis
data.input = epi.counts[, cell.use]
meta = meta[cell.use, ]
#testmeta = data.frame(labels = meta$Cluster[cell.use], row.names = colnames(data.input)) # manually create a dataframe consisting of the cell labels
meta = data.frame(labels = meta$Cluster, row.names = colnames(data.input))
unique(meta$labels) # check the cell labels


#Create cellChat object
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
future::plan("multiprocess", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

###########################################################
###########################################################
###########################################################
#Visualizing interactions
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
#
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

unique(cellchat@netP$pathways) #Select a pathway from this list
pathways.show <- c("MK") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
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


#Comparing multiple pathways
unique(cellchat@netP$pathways) #Select a pathway from this list

pathways.show <- c("MK") 
netAnalysis_contribution(cellchat, signaling = pathways.show)
pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)

# Circle plot
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")


# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), remove.isolate = FALSE)
#> Comparing communications on a single object

# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), signaling = c("EGF","MK"), remove.isolate = FALSE)
#> Comparing communications on a single object

# show all the significant interactions (L-R pairs) based on user's input (defined by `pairLR.use`)
pairLR.use <- extractEnrichedLR(cellchat, signaling = c("EGF","MK","CDH"))
netVisual_bubble(cellchat, sources.use = c(3,4), targets.use = c(5:8), pairLR.use = pairLR.use, remove.isolate = TRUE)
#> Comparing communications on a single object
# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
# show all the interactions sending from Inflam.FIB
netVisual_chord_gene(cellchat, sources.use = 4, targets.use = c(5:11), lab.cex = 0.5,legend.pos.y = 30)
#> Note: The first link end is drawn out of sector 'MIF'.

# show all the interactions received by Inflam.DC
netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = 8, legend.pos.x = 15)

plotGeneExpression(cellchat, signaling = "MK")

