library(Seurat)
library(tidyverse)
library(Matrix)
library(RCurl)
library(scales)
library(cowplot)
library(SingleCellExperiment)
library(AnnotationHub)

reading_mtx <- function(mtx, genes, cells) {
  epi.counts = readMM(mtx)
  rownames(epi.counts) = readLines(genes)
  colnames(epi.counts) = readLines(cells)
  return(epi.counts)
}

# Reading in MTX
epi.counts <- reading_mtx('gene_sorted-Epi.matrix.mtx','Epi.genes.tsv','Epi.barcodes2.tsv')


imm.counts <- reading_mtx('gene_sorted-Imm.matrix.mtx', 'Imm.genes.tsv', 'Imm.barcodes2.tsv')
fib.counts <- reading_mtx('gene_sorted-Fib.matrix.mtx', 'Fib.genes.tsv', 'Fib.barcodes2.tsv')

# Creating Seurat Object
epi.counts <- CreateSeuratObject(counts = epi.counts, project="epithelium", min.features = 100)
imm.counts <- CreateSeuratObject(counts = imm.counts, project="Immuno", min.features = 100)
fib.counts <- CreateSeuratObject(counts = fib.counts, project="Fibroblasts", min.features = 100)
head(imm.counts$nCount_RNA)
head(fib.counts$orig.ident, n=1000)
head(imm.counts$orig.ident, n=1000)
head(epi.counts$orig.ident, n=1000)
# Get Log10 ratio genes per UMI
epi.counts$log10GenesPerUMI <- log10(epi.counts$nFeature_RNA) / log10(epi.counts$nCount_RNA)
# Get Mito ratio
epi.counts$mitoRatio <- PercentageFeatureSet(object = epi.counts, pattern = "^MT-")
epi.counts$mitoRatio <- epi.counts@meta.data$mitoRatio / 100

metadata <- epi.counts@meta.data
metadata$cells <- rownames(metadata)
# Renaming columns to more convenient names
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)
# Extracting type of sample
metadata$sample <- NA
metadata$sample <- substr(gsub("^.+?\\.(.+?)\\..*$", "\\1", metadata$cells), start = 1, stop = 3)

epi.counts@meta.data <- metadata
#  Sum of Cells per Sample
metadata %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

# Frequency of transcripts per cell
metadata %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Log10 Cell density") + 
  geom_vline(xintercept = 500) + 
  ggtitle("Frequency of Transcripts Detected in Cells")


# Frequency of genes detected per cell 
metadata %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  ylab("Cell Density") +
  geom_vline(xintercept = 300) +
  ggtitle("Frequency of Unique Genes Detected in Cells")

# Genes per Cell  
metadata %>% 
  ggplot(aes(x=sample, y=log10(nGene), fill=sample)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")

# UMI(Transcripts) in relation to number of unique genes
metadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)

# Mitochondrial ratio frequency 
metadata %>% 
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)

#Filtering:
#Log10()
filtered_seurat <- subset(x = epi.counts, 
                          subset= (nUMI >= 500) & 
                            (nGene >= 250) & 
                            (log10GenesPerUMI > 0.80) & 
                            (mitoRatio < 0.20))

counts <- GetAssayData(object = filtered_seurat, slot = "counts")

nonzero <- counts > 0

keep_genes <- Matrix::rowSums(nonzero) >= 10

filtered_counts <- counts[keep_genes, ]

filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)

metadata_clean <- filtered_seurat@meta.data

#  Sum of Cells per Sample
metadata_clean %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

# Frequency of transcripts per cell
metadata_clean %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Log10 Cell density") + 
  geom_vline(xintercept = 500) + 
  ggtitle("Frequency of Transcripts Detected in Cells")


# Frequency of genes detected per cell 
metadata_clean %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  ylab("Cell Density") +
  geom_vline(xintercept = 300) +
  ggtitle("Frequency of Unique Genes Detected in Cells")
# Genes per Cell  
metadata_clean %>% 
  ggplot(aes(x=sample, y=log10(nGene), fill=sample)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")

# UMI(Transcripts) in relation to number of unique genes
metadata_clean %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)

# Mitochondrial ratio frequency 
metadata_clean %>% 
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)

normalized_epi <- NormalizeData(filtered_seurat)

s.genes <- cc.genes$s.genes
g2m <- cc.genes$g2m.genes
normalized_epi_phase <- CellCycleScoring(normalized_epi, 
                                         g2m.features =  g2m.genes, 
                                         s.features =  s.genes)

View(normalized_epi_phase@meta.data)

normalized_epi_phase <- FindVariableFeatures(normalized_epi_phase,
                                             selection.method = "vst",
                                             nfeatures = 2000,
                                             verbose = FALSE)
normalized_epi_phase <- ScaleData(normalized_epi_phase)

normalized_epi_phase <- RunPCA(normalized_epi_phase)
DimPlot(normalized_epi_phase,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase")
options(future.globals.maxSize = 4000 * 1024^2)

split_seurat <- SplitObject(filtered_seurat, split.by = "sample")
head(split_seurat$)
split_seurat <- split_seurat[c("Epi","LPA","LPB")]
for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- NormalizeData(split_seurat[[i]], verbose = TRUE)
  split_seurat[[i]] <- CellCycleScoring(split_seurat[[i]], g2m.features=g2m.genes, s.features=s.genes)
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("mitoRatio"))
}


integ_features <- SelectIntegrationFeatures(object.list = split_seurat, 
                                            nfeatures = 3000) 

split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = integ_features)

integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)

seurat_integrated <- IntegrateData(anchorset = integ_anchors, 
                                   normalization.method = "SCT")


