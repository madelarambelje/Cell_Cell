setwd("~/Desktop/ulcerative_colitus")
library("Seurat")
library("Matrix")
library(ggplot2)
# Loading smillie preprocessed
epi.seur <- readRDS('train.Epi.seur.rds')
fib.seur <- readRDS('train.Fib.seur.rds')
imm.seur <- readRDS('train.Imm.seur.rds')


# Merging all celltypes
all.counts <- merge(epi.seur, 
                    y= c(fib.seur, imm.seur),
                    add.cell.ids = c("epi", "fib", "imm"))
# Get number of celltypes 
length(unique(all.counts@meta.data$Cluster))
# 51

#Setting Ident on Health
Idents(all.counts) <- "Health"
#Extract healthy counts 
healthy.counts <- subset(all.counts, idents = "Healthy")
# Get number of  celltypes in healthy
length(unique(healthy.counts@meta.data$Cluster))
# 51

#Extract inflamed counts
inflamed.counts <- subset(all.counts, idents = "Inflamed")

# Get total number of cells in inflamed and healthy
sum.cells.inflamed <- sum(table(inflamed.counts$Cluster))
sum.cells.health <- sum(table(healthy.counts$Cluster))
# Getting percentage of celltypes present for how https://satijalab.org/howmanycells/
# Inflamed
inflamed.cell.freq <- table(inflamed.counts$Cluster)/sum.cells.inflamed
inflamed.cell.freq[order(inflamed.cell.freq)]
# Healthy
healthy.cell.freq <- table(healthy.counts$Cluster)/sum.cells.health
healthy.cell.freq[order(healthy.cell.freq)]

# Actual subsetting samplesize retrieved from https://satijalab.org/howmanycells/, 
# minimum cell type fraction extracted from fractions healthy
subset_healthy <- healthy.counts[,sample(colnames(healthy.counts), size = 37546)]
freq.subset.healthy <- table(subset_healthy$Cluster) / sum(table(subset_healthy$Cluster))
freq.subset.healthy[order(freq.subset.healthy)]

saveRDS(subset_healthy, "healthy_subset_natural.rds")
# Subsetting not required for inflamed because of low cell count
saveRDS(inflamed.counts, "inflamed_subset_natural.rds")
