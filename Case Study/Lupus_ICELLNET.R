
#This file contains the code for comparing ICELLNET on tutorial article to compare the interactions between
#specific cell of interest 
#Author: Shreya Dey


#Load data
#Change the file names here
per.cc <- readRDS(file = "Lupus_Seurat_SingleCell_Landscape.Rds") 
per.cc <- NormalizeData(per.cc)
per.cc <- ScaleData(per.cc)

#only for UMAP visualization, not for ICELLNET purpose
seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)
seurat <- RunPCA(seurat)
seurat <- RunUMAP(seurat, dims = 1:50)
DimPlot(seurat, reduction = 'umap', group.by = 'author_annotation', label = T)

filter.perc=0
average.clean.percc.ic= sc.data.cleaning(object = per.cc, db = db, filter.perc = filter.perc, save_file = T, path=NULL, force.file = F)

data.icell.percc=as.data.frame(gene.scaling(as.data.frame(average.clean.percc), n=1, db=db))

# label and color label if you are working families of molecules already present in the database
#66ABDF e11e6d
my.family=c("Growth factor","Chemokine","Checkpoint","Cytokine","Notch family","Antigen binding")
family.col = c( "Growth factor"= "#AECBE3", "Chemokine"= "#66ABDF", "Checkpoint"= "#1D1D18"  ,
                "Cytokine"="#156399", "Notch family" ="#676766", "Antigen binding" = "#12A039",  "other" = "#908F90",  "NA"="#908F90")


recpex = colnames(data.icell.percc[, !names(data.icell.percc) %in% c("Symbol")]) #receptors expressed

PC.data.percc=as.data.frame(data.icell.percc[], row.names = rownames(data.icell.percc))

#my.colname = c(recpex)
PC.target.percc=data.frame("Class"=c(recpex), "ID"= c(recpex), "Cell_type"=c(recpex))
rownames(PC.target.percc)=c(recpex)
my.selection.percc=c(recpex)


outdirection = "out"
indirection = "in"

#Replace the central cell with CM4 for comparing interactions for CM4 cell for Lupus dataset
centralcell = "CE0"
interactions = c("CXCL12 / CXCR4","CX3CL1 / CX3CR1")

score.computation.perccin= icellnet.score(direction=indirection, PC.data=PC.data.percc, 
                                          CC.data= as.data.frame(data.icell.percc[,c(centralcell)], row.names = rownames(data.icell.percc)),  
                                          PC.target = PC.target.percc, PC=my.selection.percc, CC.type = "RNAseq", 
                                          PC.type = "RNAseq",  db = db)
score1.perccin=as.data.frame(score.computation.perccin[[1]])
lr1.perccin=score.computation.perccin[[2]]

score.computation.perccout= icellnet.score(direction=outdirection, PC.data=PC.data.percc, 
                                           CC.data= as.data.frame(data.icell.percc[,c(centralcell)], row.names = rownames(data.icell.percc)),  
                                           PC.target = PC.target.percc, PC=my.selection.percc, CC.type = "RNAseq", 
                                           PC.type = "RNAseq",  db = db)
score1.perccout=as.data.frame(score.computation.perccout[[1]])
lr1.perccout=score.computation.perccout[[2]]

Scoresio=cbind(score1.perccin,score1.perccout)
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
contrib.family.1= LR.family.score(lr=lr1.perccin, my.family=my.family, db.couple=db.name.couple, plot= T, family.col=family.col, title=paste(paste(centralcell),paste("(R)-in")))
contrib.family.2= LR.family.score(lr=lr1.perccout, my.family=my.family, db.couple=db.name.couple, plot= T, family.col=family.col, title=paste(paste(centralcell),paste("(L)-out")))
grid.arrange(contrib.family.1, contrib.family.2, ncol=2, nrow=1)

#LR tables
LR.family.score(lr=lr1.perccin, my.family=my.family, db.couple=db.name.couple, plot=F)
LR.family.score(lr=lr1.perccout, my.family=my.family, db.couple=db.name.couple, plot=F)

#For balloon plot and heatmap we are selecting specific interactions as we want to check if interaction is present or not
#balloon plot
LR.balloon.plot(lr = lr1.perccin[interactions,], thresh = 10 , topn=11 , sort.by="sum",  db.name.couple=db.name.couple, family.col=family.col, title=paste(paste("Most contributing interactions of"),paste(centralcell),paste("by sum/in")))
LR.balloon.plot(lr = lr1.perccout[interactions,], thresh = 0 , topn=30 , sort.by="sum",  db.name.couple=db.name.couple, family.col=family.col, title=paste(paste("Most contributing interactions of"),paste(centralcell),paste("by sum/out")))
lrsum = lr1.perccin + lr1.perccout
LR.balloon.plot(lr = lrsum[interactions,], thresh = 10 , topn=11 , sort.by="sum",  db.name.couple=db.name.couple, family.col=family.col, title=paste(paste("Most contributing interactions of"),paste(centralcell),paste("by sum(in+out)")))

#heatmap
LR.heatmap(lr = lr1.perccin[interactions,], thresh = 0, topn = NULL, sort.by = "sum",  db.name.couple = db.name.couple, title = paste(paste(centralcell),paste("-IN")), family.col = family.col, value_display=50)  
LR.heatmap(lr = lr1.perccout[interactions,], thresh = 0, topn = NULL, sort.by = "sum",  db.name.couple = db.name.couple, title = paste(paste(centralcell),paste("-OUT")), family.col = family.col, value_display=50)  


