# Author: Andre de la Rambelje

library("Seurat")

setwd("/students/2020-2021/master/CellPhoneDB/")

load(file = "article_cell.rds")

rownames(log2cpm) <- featuredata$Associated.Gene.Name


# Numbers connected to celltypes retrieved from the READme 
tsne.data$dbCluster[tsne.data$dbCluster == 1] <- "Myeloid"
tsne.data$dbCluster[tsne.data$dbCluster == 2] <- "Epithelial"
tsne.data$dbCluster[tsne.data$dbCluster == 3] <- "mCAF"
tsne.data$dbCluster[tsne.data$dbCluster == 4] <- "iCAF"
tsne.data$dbCluster[tsne.data$dbCluster == 5] <- "Endothelial"
tsne.data$dbCluster[tsne.data$dbCluster == 6] <- "T cell"
tsne.data$dbCluster[tsne.data$dbCluster == 7] <- "B cell"
tsne.data$dbCluster[tsne.data$dbCluster == 8] <- "Mast cell"
meta.data <- tsne.data$dbCluster
seurat.object <- CreateSeuratObject(log2cpm)
seurat.object@meta.data$Cluster <- meta.data





saveRDS(seurat.object, "Article_CellPhoneDB.rds")
