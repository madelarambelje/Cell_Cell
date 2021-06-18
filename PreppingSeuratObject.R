# Author: Andre de la Rambelje

Lupus <- readRDS("seurat_object.rds")
#Locate the column where cell labels are stored. Change this column name to "Cluster"
# In this case the column name was found at index four in the meta data.
colnames(seurat.object@meta.data)[4] <- "Cluster"
# Save this adjusted seurat.object
saveRDS(seurat.object, "seurat_object.rds")
