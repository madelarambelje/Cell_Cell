# Install / Load required packages for cellchat analysis
# Nick Veltmaat
# 16-5-2021

install.packages('pacman')
install.packages('BiocManager')

pacman::p_load(ComplexHeatmap, circlize, NMF, Matrix, Seurat)
devtools::install_github("sqjin/CellChat")
pacman::p_load('CellChat')

system("pip3 install umap-learn")
