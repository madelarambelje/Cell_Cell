# Creating Seurat object (as RDS) from example CellChat tutorial data
# https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html
# Nick Veltmaat
# 21-5-2021

setwd('~/Desktop/Master Files/Cell_Cell_Project/Cellchat/Vignette/data')

# Load the raw data and metadata as a list
load("data_humanSkin_CellChat.rda")

# Extract data and metadata from list
data.input = data_humanSkin$data # normalized data matrix
meta = data_humanSkin$meta # a dataframe with rownames containing cell meta data

cell.use = rownames(meta)[meta$condition == "LS"] # extract the cell names from disease data


# Prepare input data for converting to seurat object
data.input = data.input[, cell.use]
meta = meta[cell.use, ]

# Create Seurat object
seurat.object <- CreateSeuratObject(data.input, assay = "RNA", meta.data = meta)

# Save Seurat object as RDS file
saveRDS(seurat.object, file = "CellChat_TutorialData_SeuratObject.rds")
