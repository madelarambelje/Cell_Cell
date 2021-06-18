# Author: Andre de la Rambelje

library("Seurat")
seurat.for.genes <- readRDS("CellChat_TutorialData_SeuratObject.rds")
seurat.for.cells <- readRDS("CellChat_TutorialData_SeuratObject.rds")

colnames(seurat.for.genes@meta.data)[6] <- "Cluster"
colnames(seurat.for.cells@meta.data)[6] <- "Cluster"
# Extracting counts and metadata
counts <- GetAssayData(seurat.for.genes,slot = "counts")
meta.data.for.genes <- seurat.for.genes@meta.data 
# Shuffle the Gene names
for (n in 1:100) {
  print(n)
  rownames(counts) <- sample(rownames(counts))
  }

permutated.genes.to.save <- CreateSeuratObject(counts = counts, meta.data= meta.data ) 
saveRDS(permutated.genes.to.save, "permutated_GeneLabels.rds")

# Permutate the cell labels
for (n in 1:100) {
  print(n)
  seurat.for.cells@meta.data$Cluster <- sample(seurat.for.cells@meta.data$Cluster)
  }

saveRDS(seurat.for.cells, "permutated_CellLabels.rds")

# Combined Permutation
seurat.for.combined <- CreateSeuratObject(counts = counts, meta.data = seurat.for.cells@meta.data)
saveRDS(seurat.for.combined,"permutated_Combined.rds")

