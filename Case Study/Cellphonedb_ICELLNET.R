
# This script loads the dataset for CellPhoneDB and compares the interaction between ligand receptor pairs
# of cytokines and interaction between growth factors between iCAF and other cell types.

# Author: Shreya Dey
# University: Hanzehogeschool Groningen,Netherlands
# Date: June 15, 2021


setwd('~/Desktop/Masters/Cell_Cell_Project/ICELLNET')

#Run prerequisite.r to load the necessary libraries and functions
source("prerequisite.r")


cpdb <- readRDS(file = "Article_CellPhoneDB.rds")

#only for UMAP visualization, not for ICELLNET purpose
cpdb <- FindVariableFeatures(cpdb, selection.method = "vst", nfeatures = 2000)
cpdb <- RunPCA(cpdb)
cpdb <- RunUMAP(cpdb, dims = 1:50)
DimPlot(cpdb, reduction = 'umap', group.by = 'Cluster', label = T)


#Average.clean function does not work for the datasets hence the manual calculatation 
#mentioned by author is used to get calculations

#calculation of average gene expression
data <- as.data.frame(GetAssayData(cpdb, slot = "data")) #or other matrix for which features expression values # are scaled by the total expression in each cell


##---------------------------------- Important Note ---------------------------------------------------------##
## Replace the name of column in which Cell label is defined as per your seurat object in calculations       ##
## You can do this by writing name of your seurat object followed by two square brackets as shown belowthis  ##
## cpdb[[]]                                                                                                  ##
##-----------------------------------------------------------------------------------------------------------##

target.cpdb <- cpdb@meta.data
target.cpdb$Class=target.cpdb$Cluster
target.cpdb$Cell=rownames(target.cpdb)

average.manual.cpdb=matrix(ncol=length(unique(target.cpdb$Cluster)), nrow=length(rownames(data)))
colnames(average.manual.cpdb)=unique(target.cpdb$Cluster)
rownames(average.manual.cpdb)=rownames(data)
dim(average.manual.cpdb)
for (cell in unique(target.cpdb$Cluster)){
  cells.clust=target.cpdb$Cell[which(target.cpdb$Cluster==cell)]
  average.manual.cpdb[,cell]=apply(data[,which(colnames(data)%in%cells.clust)], 1, mean)
}

average.clean.cpdb=average.manual.cpdb #average gene expression

data.icell.cpdb=as.data.frame(gene.scaling(as.data.frame(average.clean.cpdb), n=1, db=db))

# label and color label if you are working families of molecules already present in the database

my.family=c("Growth factor","Chemokine","Checkpoint","Cytokine","Notch family","Antigen binding")
family.col = c( "Growth factor"= "#AECBE3", "Chemokine"= "#66ABDF", "Checkpoint"= "#1D1D18"  ,
                "Cytokine"="#156399", "Notch family" ="#676766", "Antigen binding" = "#12A039",  "other" = "#908F90",  "NA"="#908F90")


recpex = colnames(data.icell.cpdb[, !names(data.icell.cpdb) %in% c("Symbol")]) #receptors expressed

PC.data.cpdb=as.data.frame(data.icell.cpdb[], row.names = rownames(data.icell.cpdb))

#my.colname = c(recpex)
PC.target.cpdb=data.frame("Class"=c(recpex), "ID"= c(recpex), "Cell_type"=c(recpex))
rownames(PC.target.cpdb)=c(recpex)
my.selection.cpdb=c(recpex)


outdirection = "out"
indirection = "in"

#iCAF selected as central cell for the comparison
centralcell = "iCAF"
#selecting cytokines interactions for icaf cell avaialble in ICELLNET
interactions_icaf = c("CXCL12 / CXCR4","CX3CL1 / CX3CR1","CXCL12 / ACKR3","CXCL8 / ACKR1","CXCL16 / CXCR6","CXCL12 / CXCR4","CXCL12 / ACKR3","CXCL1 / ACKR1","CCL5 / ACKR1","CCL5 / CCR5","CCL5 / CCR1","CCL21 / CCR7","CCL2 / ACKR1","CCL17 / ACKR1")
#seelcting icaf growth factor interactions available in ICELLNET
interactions_icaf_gf = c("HGF / MET","PDGFC / PDGFRA","VEGFA / FLT1","VEGFA / KDR","VEGFA / NRP1","VEGFB / FLT1")


score.computation.cpdbin= icellnet.score(direction=indirection, PC.data=PC.data.cpdb, 
                                          CC.data= as.data.frame(data.icell.cpdb[,c(centralcell)], row.names = rownames(data.icell.cpdb)),  
                                          PC.target = PC.target.cpdb, PC=my.selection.cpdb, CC.type = "RNAseq", 
                                          PC.type = "RNAseq",  db = db)
score1.cpdbin=as.data.frame(score.computation.cpdbin[[1]])
lr1.cpdbin=score.computation.cpdbin[[2]]

score.computation.cpdbout= icellnet.score(direction=outdirection, PC.data=PC.data.cpdb, 
                                           CC.data= as.data.frame(data.icell.cpdb[,c(centralcell)], row.names = rownames(data.icell.cpdb)),  
                                           PC.target = PC.target.cpdb, PC=my.selection.cpdb, CC.type = "RNAseq", 
                                           PC.type = "RNAseq",  db = db)
score1.cpdbout=as.data.frame(score.computation.cpdbout[[1]])
lr1.cpdbout=score.computation.cpdbout[[2]]

Scoresio=cbind(score1.cpdbin,score1.cpdbout)
colnames(Scoresio)=c(paste(paste(centralcell),paste("(R)-in")),paste(paste(centralcell),paste("(L)-out")))
Scoresio #gives in and out scores in one go

Scores.norm=(Scoresio-min(Scoresio))/(max(Scoresio)-min(Scoresio))*9+1
#Color scheme for cells available in the cluster
PC.col = c("Myeloid"="#88b04b","Epithelial"="#0ee6ec","mCAF"="#e861f0","iCAF"="#0cf628",
            "Endothelial"="#ff962c","T cell"="#e9f60c", "B cell"="#0cc1f6","Mast cell"="#f60c76")

network.plot1 = network.create(icn.score = Scores.norm[1], scale = c(round(min(Scores.norm)),round(max(Scores.norm))), direction = "in", PC.col)
network.plot2 = network.create(icn.score =Scores.norm[2], scale = c(round(min(Scores.norm)),round(max(Scores.norm))), direction = "out",PC.col)
gridExtra::grid.arrange(network.plot1, network.plot2, ncol=2, nrow=1)

#ymax for barplot
ymax=round(max(Scoresio))+1


#code for box plot
#Boxplot shows global scores at level of molecules of class, there is no need to select scores for specific interactions
contrib.family.1= LR.family.score(lr=lr1.cpdbin, my.family=my.family, db.couple=db.name.couple, plot= T, family.col=family.col, title=paste(paste(centralcell),paste("(R)-in")))
contrib.family.2= LR.family.score(lr=lr1.cpdbout, my.family=my.family, db.couple=db.name.couple, plot= T, family.col=family.col, title=paste(paste(centralcell),paste("(L)-out")))
grid.arrange(contrib.family.1, contrib.family.2, ncol=2, nrow=1)

#LR tables
LR.family.score(lr=lr1.cpdbin, my.family=my.family, db.couple=db.name.couple, plot=F)
LR.family.score(lr=lr1.cpdbout, my.family=my.family, db.couple=db.name.couple, plot=F)

#summing ligand-receptor in and out score 
lr_icafsum = lr1.cpdbin + lr1.cpdbout

#balloon plot
LR.balloon.plot(lr = lr1.cpdbin[interactions_icaf,], thresh = 10 , topn=30 , sort.by="sum",  db.name.couple=db.name.couple, family.col=family.col, title="Most contributing interactions of iCAF by sum/in")
LR.balloon.plot(lr = lr1.cpdbout[interactions_icaf,], thresh = 10 , topn=30 , sort.by="sum",  db.name.couple=db.name.couple, family.col=family.col, title="Most contributing interactions of iCAF by sum/out")
LR.balloon.plot(lr = lr_icafsum[interactions_icaf,], thresh = 10 , topn=11 , sort.by="sum",  db.name.couple=db.name.couple, family.col=family.col, title="Most contributing interactions of iCAF by sum")


LR.balloon.plot(lr = lr1.cpdbin[interactions_icaf_gf,], thresh = 10 , topn=30 , sort.by="sum",  db.name.couple=db.name.couple, family.col=family.col, title="Most contributing interactions of iCAF by sum/in")
LR.balloon.plot(lr = lr1.cpdbout[interactions_icaf_gf,], thresh = 10 , topn=30 , sort.by="sum",  db.name.couple=db.name.couple, family.col=family.col, title="Most contributing interactions of iCAF by sum/out")
LR.balloon.plot(lr = lr_icafsum[interactions_icaf_gf,], thresh = 10 , topn=11 , sort.by="sum",  db.name.couple=db.name.couple, family.col=family.col, title="Most contributing interactions of iCAF by sum")


#Heatmap

 # ---------Important Note---------------------------------------------------------------------------------------------------------#
 # Comment below 2 lines before running the code, this function is not yet published by the author, here *xarguments represents    #
 # arguments which were passed to create the heatmaps                                                                              #
 #---------------------------------------------------------------------------------------------------------------------------------#
  
 
LR.heatmap(*xarguments, title = paste(paste(centralcell),paste("-IN")))  
LR.heatmap(*xarguments, title = paste(paste(centralcell),paste("-OUT")))  
