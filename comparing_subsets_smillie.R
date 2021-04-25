pacman::p_load(ComplexHeatmap, circlize, NMF, CellChat, Matrix, Seurat)
setwd('~/Desktop/Master Files/Cell_Cell_Project/Cellchat')

#Load seurat objects / subsets
t_cell <- readRDS('smilliedata/t_cell_subset.rds')
natural_composition <- readRDS('smilliedata/natural_composition_subset.rds')
uniform <- readRDS('smilliedata/uniform_subset.rds')

#Merge Seurat objects
#Test
#all.counts <- merge(t_cell,
#                   y= c(natural_composition, uniform),
#                    add.cell.ids = c("t_cells", "natural_composition", "uniform_distribution"))

#Convert seurat to Cellchat Objects
cellchat_t_cell <- createCellChat(object = t_cell, group.by = "ident")
cellchat_natural_comp <- createCellChat(object = natural_composition, group.by = "ident")
cellchat_uniform <- createCellChat(object = uniform, group.by = "ident")
#cellchat_merged <- createCellChat(object = all.counts, group.by = "ident")


#Merge Cellchat objects
object.list <- list(t_cell = cellchat_t_cell ,
                    nat_comp = cellchat_natural_comp
                    ,uniform = cellchat_uniform
                    )

cellchat <- mergeCellChat(object.list, add.names = names(object.list),cell.prefix = TRUE)
cellchat


#Compare the total number of interactions and interaction strength
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2)) #or with cellchat_merged
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2
