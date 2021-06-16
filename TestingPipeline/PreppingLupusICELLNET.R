library("Seurat")

getwd()
Lupus <- readRDS("TestingPipeline/Lupus_Seurat_SingleCell_Landscape.Rds")
colnames(Lupus@meta.data)[4] <- "Cluster"

saveRDS(Lupus, "TestingPipeline/Lupus_Seurat_SingleCell_Landscape.Rds")


