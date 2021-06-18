install.packages("pacman")
pacman::p_load(devtools, curl, jetset, readxl, psych, GGally, gplots, ggplot2, RColorBrewer, data.table, gridExtra, ggthemes, scales, rlist, Seurat,dplyr, NormalizeData) 

setwd('~/Desktop/Masters/Cell_Cell_Project/ICELLNET')

# Installs Bioconductor dependancies 
if (!requireNamespace("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager")
BiocManager::install(c("BiocGenerics", "org.Hs.eg.db", "hgu133plus2.db", "annotate"))

library(devtools)
library(BiocGenerics)
library("org.Hs.eg.db")
library("hgu133plus2.db")
library(jetset)
library(ggplot2)
library(dplyr)
library(icellnet)
library(gridExtra)
library(Seurat)

#installs ICELLNET package from github
install_github("soumelis-lab/ICELLNET",ref="master", subdir="icellnet")


# creating a dataframe of ICELLNET database creating couples based on interaction pairs present
db=as.data.frame(read.csv(curl::curl(url="https://raw.githubusercontent.com/soumelis-lab/ICELLNET/master/data/ICELLNETdb.tsv"), sep="\t",header = T, check.names=FALSE, stringsAsFactors = FALSE, na.strings = ""))
#Here a database couple is created based on the interaction pairs present in ICELLNET DB based on Family of molecules 
db.name.couple=name.lr.couple(db, type="Family")
