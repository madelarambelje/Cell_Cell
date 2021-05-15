# Run this script:
# Nick Veltmaat
# 15-5-2021
# In this example, an Inflamed and Healthy subset from Smillie's data (https://github.com/cssmillie/ulcerative_colitis)
# with natural cell compositions will be compared using CellChat (https://github.com/sqjin/CellChat)
setwd('~/Desktop/Master Files/Cell_Cell_Project/Cellchat')

# Firstly, install/load the required packages using the 'requisites.r' script
source('Vignette/requisites.R')



# Give the script the Seurat objects that you want to compare / analyze
subset1 <- readRDS('Vignette/data/healthy_subset_natural.rds') # E.g. Healthy subset
subset2 <- readRDS('Vignette/data/inflamed_subset_natural.rds') # E.g. Inflamed subset



# Create Cellchat Objects (takes time)
source('Vignette/create_cellchat_objects.R')
## OR LOAD CELLCHAT OBJECTS FROM RDS to skip line above:
# subset1 <- readRDS("cellchat1.rds")
# subset2 <- readRDS("cellchat2.rds")

### Merge Cellchat objects, set condition for subsets: (Replace Healthy / Inflamed as preferred)
object.list <- list(Healthy = cellchat1, 
                    Inflamed = cellchat2)
cellchat <- mergeCellChat(object.list, 
                          add.names = names(object.list))
cellchat



# Identify the conserved and context-specific signaling pathways
source('Vignette/similarity.R')

# Save merged cellchatfile as RDS:
# saveRDS(cellchat, file = "merged_cellchat.rds")

# Skip all previous code by loading merged cellchat object RDS file:
# cellchat <- readRDS("merged_cellchat.rds")



# Create Comparison visualization functions:
source('Vignette/visualizations.R')



## PART I
PLOT.compare.interactions(x = cellchat)
PLOT.compare.interactions.circle(x = cellchat)
PLOT.compare.interactions.HM(x = cellchat)


## PART II
# Identify signaling networks with larger (or less) difference as well as 
# signaling groups based on their functional similarity
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5) # Functional similarity
netVisual_embeddingPairwiseZoomIn(cellchat, type = "functional", nCol = 2) # Zoom in on Clusters
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5) # Structural similarity
netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2) # Zoom in on Clusters
PLOT.rank.similarity(x = cellchat)
PLOT.comparison.info.flow(x = cellchat)
PLOT.outgoing.signals.HM(i = 1)
PLOT.incoming.signals.HM(i = 1)
PLOT.overall.signals.HM(i = 1)


## PART III
# Identify dysfunctional signaling by comparing the communication probabities
netVisual_bubble(cellchat, sources.use =2, targets.use = c(0:1000),  
                 comparison = c(1, 2), angle.x = 45)
PLOT.comm.prob.BP(x = cellchat, sources.use = 5, targets.use = c(0:1000),
                  comparison = c(1,2))


## PART IV
#Select a pathway from these lists:
cellchat@netP[[names(object.list)[1]]][["pathways"]]
cellchat@netP[[names(object.list)[2]]][["pathways"]]
pathways.show <- c("IL1") #Set pathway to analyze here

PLOT.pathway.circle(i = 1, pathways.show = pathways.show)
PLOT.pathway.HM(i = 1, pathways.show = pathways.show)
PLOT.pathway.chord(i = 1, pathways.show = pathways.show)


## PART V
PLOT.expression.distribution(x = cellchat, pathways.show = pathways.show)
