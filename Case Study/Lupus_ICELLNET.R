
# This file contains the code for comparing ICELLNET on tutorial article to compare the interactions between
# specific cell of interest (CE0 and CM4)

# Author: Shreya Dey
# University: Hanzehogeschool, Groningen
# Date: June 18, 2021


setwd('~/Desktop/Masters/Cell_Cell_Project/ICELLNET')

#Run prerequisite.r to load the necessary libraries and functions
source("prerequisite.r")

#Load data
#Change the file names here
lup <- readRDS(file = "Lupus_Seurat_SingleCell_Landscape.Rds") 
lup <- NormalizeData(lup)
lup <- ScaleData(lup)

#only for UMAP visualization, not for ICELLNET purpose
lup <- FindVariableFeatures(lup, selection.method = "vst", nfeatures = 2000)
lup <- RunPCA(lup)
lup <- RunUMAP(lup, dims = 1:50)
DimPlot(lup, reduction = 'umap', group.by = 'author_annotation', label = T)

# The lupus dataset is by default annoted by clusters as compared to what was used in the actual tutorial,
# if any error is encountered while running the sc.data.cleaning function, run the next line before running the rest of the code

Idents(lup)=lup$author_annotation

filter.perc=0
average.clean.lup= sc.data.cleaning(object = lup, db = db, filter.perc = filter.perc, save_file = T, path=NULL, force.file = F)

# If the above code for calculating average clean does not work, uncomment the below code and try to run

# target <- lup@meta.data
# target$Class=target$author_annotation
# target$Cell=rownames(target)
# 
# average.manual=matrix(ncol=length(unique(target$author_annotation)), nrow=length(rownames(data)))
# colnames(average.manual)=unique(target$author_annotation)
# rownames(average.manual)=rownames(data)
# dim(average.manual)
# for (cell in unique(target$author_annotation)){
#   cells.clust=target$Cell[which(target$author_annotation==cell)]
#   average.manual[,cell]=apply(data[,which(colnames(data)%in%cells.clust)], 1, mean)
# }
# 
# average.clean=average.manual

data.icell.lup=as.data.frame(gene.scaling(as.data.frame(average.clean.lup), n=1, db=db))

# label and color label if you are working families of molecules already present in the database
#66ABDF e11e6d
my.family=c("Growth factor","Chemokine","Checkpoint","Cytokine","Notch family","Antigen binding")
family.col = c( "Growth factor"= "#AECBE3", "Chemokine"= "#66ABDF", "Checkpoint"= "#1D1D18"  ,
                "Cytokine"="#156399", "Notch family" ="#676766", "Antigen binding" = "#12A039",  "other" = "#908F90",  "NA"="#908F90")


recpex = colnames(data.icell.lup[, !names(data.icell.lup) %in% c("Symbol")]) #receptors expressed

PC.data.lup=as.data.frame(data.icell.lup[], row.names = rownames(data.icell.lup))

#my.colname = c(recpex)
PC.target.lup=data.frame("Class"=c(recpex), "ID"= c(recpex), "Cell_type"=c(recpex))
rownames(PC.target.lup)=c(recpex)
my.selection.lup=c(recpex)


outdirection = "out"
indirection = "in"

#Replace the central cell with CM4 for comparing interactions for CM4 cell for Lupus dataset
centralcell = "CE0"
interactions = c("CXCL12 / CXCR4","CX3CL1 / CX3CR1")

score.computation.lupin= icellnet.score(direction=indirection, PC.data=PC.data.lup, 
                                          CC.data= as.data.frame(data.icell.lup[,c(centralcell)], row.names = rownames(data.icell.lup)),  
                                          PC.target = PC.target.lup, PC=my.selection.lup, CC.type = "RNAseq", 
                                          PC.type = "RNAseq",  db = db)
score1.lupin=as.data.frame(score.computation.lupin[[1]])
lr1.lupin=score.computation.lupin[[2]]

score.computation.lupout= icellnet.score(direction=outdirection, PC.data=PC.data.lup, 
                                           CC.data= as.data.frame(data.icell.lup[,c(centralcell)], row.names = rownames(data.icell.lup)),  
                                           PC.target = PC.target.lup, PC=my.selection.lup, CC.type = "RNAseq", 
                                           PC.type = "RNAseq",  db = db)
score1.lupout=as.data.frame(score.computation.lupout[[1]])
lr1.lupout=score.computation.lupout[[2]]

Scoresio=cbind(score1.lupin,score1.lupout)
colnames(Scoresio)=c(paste(paste(centralcell),paste("(R)-in")),paste(paste(centralcell),paste("(L)-out")))
Scoresio #combined input and output score

Scores.norm=(Scoresio-min(Scoresio))/(max(Scoresio)-min(Scoresio))*9+1

#Color scheme for cells of Lupus dataset
PC.col = c("CE0" ="#0cc1f6", "CD0" ="#ff962c", "CM0"="#0cf628","CM1"="#0cf628",  "CM2"="#0cf628" , "CM3"="#0cf628",
           "CM4"="#0cf628" ,"CB2b"="#0cf628","CT1" = "#e9f60c", "CT2"="#e9f60c", "CT5b"="#e9f60c", "CT4" ="#e9f60c",
           "CT5a"="#e9f60c","CT0a"="#e9f60c","CT0b" ="#e9f60c", "CT3a" = "#e9f60c","CT6" ="#e9f60c","CT3b" ="#e9f60c",
           "CB1" ="#c100b9", "CB3" = "#f60c76","CB0"="#f60c76","CB2a" ="#f60c76")

network.plot1 = network.create(icn.score = Scores.norm[1], scale = c(round(min(Scores.norm)),round(max(Scores.norm))), direction = "in", PC.col)
network.plot2 = network.create(icn.score =Scores.norm[2], scale = c(round(min(Scores.norm)),round(max(Scores.norm))), direction = "out",PC.col)
gridExtra::grid.arrange(network.plot1, network.plot2, ncol=2, nrow=1)


#ymax for barplot
ymax=round(max(Scoresio))+1


#code for box plot
#Boxplot shows global scores at level of molecules of class, there is no need to select scores for specific interactions
contrib.family.1= LR.family.score(lr=lr1.lupin, my.family=my.family, db.couple=db.name.couple, plot= T, family.col=family.col, title=paste(paste(centralcell),paste("(R)-in")))
contrib.family.2= LR.family.score(lr=lr1.lupout, my.family=my.family, db.couple=db.name.couple, plot= T, family.col=family.col, title=paste(paste(centralcell),paste("(L)-out")))
grid.arrange(contrib.family.1, contrib.family.2, ncol=2, nrow=1)

#LR tables
LR.family.score(lr=lr1.lupin, my.family=my.family, db.couple=db.name.couple, plot=F)
LR.family.score(lr=lr1.lupout, my.family=my.family, db.couple=db.name.couple, plot=F)

#For balloon plot and heatmap we are selecting specific interactions as we want to check if interaction is present or not
#balloon plot
LR.balloon.plot(lr = lr1.lupin[interactions,], thresh = 10 , topn=11 , sort.by="sum",  db.name.couple=db.name.couple, family.col=family.col, title=paste(paste("Most contributing interactions of"),paste(centralcell),paste("by sum/in")))
LR.balloon.plot(lr = lr1.lupout[interactions,], thresh = 0 , topn=30 , sort.by="sum",  db.name.couple=db.name.couple, family.col=family.col, title=paste(paste("Most contributing interactions of"),paste(centralcell),paste("by sum/out")))
lrsum = lr1.lupin + lr1.lupout
LR.balloon.plot(lr = lrsum[interactions,], thresh = 10 , topn=11 , sort.by="sum",  db.name.couple=db.name.couple, family.col=family.col, title=paste(paste("Most contributing interactions of"),paste(centralcell),paste("by sum(in+out)")))

#heatmap
LR.heatmap(lr = lr1.lupin[interactions,], thresh = 0, topn = NULL, sort.by = "sum",  db.name.couple = db.name.couple, title = paste(paste(centralcell),paste("-IN")), family.col = family.col, value_display=50)  
LR.heatmap(lr = lr1.lupout[interactions,], thresh = 0, topn = NULL, sort.by = "sum",  db.name.couple = db.name.couple, title = paste(paste(centralcell),paste("-OUT")), family.col = family.col, value_display=50)  


