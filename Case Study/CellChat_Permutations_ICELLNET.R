# This script contains code for running ICELLNET pipeline through CellChat tutorial datasets,
# We used this script for permuted CellChat datasets as well as the actual cellchat turorial datasets.
# Filenames need to be changed for visualizing differences.


# Author: Shreya Dey
# University: Hanzehogeschool Groningen,Netherlands
# Date: June 3, 2021



#pdf("CellChatresults.pdf") # To save all plots in a pdf file

#Run prerequisite.r to load the necessary libraries and functions
source(prerequisite.r)


#Load data
#Change the file names here
per.cc <- readRDS(file = "PermutationCombinedCellChat.rds") 
#per.cc <- NormalizeData(per.cc) #data is already normalised
#per.cc <- ScaleData(per.cc)

#umap visualization not for icellnet purpose, to check the clusters and a precheck
per.cc <- FindVariableFeatures(per.cc, selection.method = "vst", nfeatures = 2000)
per.cc <- RunPCA(per.cc)
per.cc <- RunUMAP(per.cc, dims = 1:50)
DimPlot(per.cc, reduction = 'umap', group.by = 'Cluster', label = T)

#Average.clean function does not work for the datasets hence the manual calculatation 
#mentioned by author is used to get calculations

#calculation of average gene expression
data <- as.data.frame(GetAssayData(per.cc, slot = "data")) #or other matrix for which features expression values # are scaled by the total expression in each cell


##---------------------------------- Important Note ----------------------------------##
## Replace the name of Cluster column as per your seurat object in calculations below ##
##------------------------------------------------------------------------------------##

target.percc <- per.cc@meta.data
target.percc$Class=target.percc$Cluster
target.percc$Cell=rownames(target.percc)

average.manual.percc=matrix(ncol=length(unique(target.percc$Cluster)), nrow=length(rownames(data)))
colnames(average.manual.percc)=unique(target.percc$Cluster)
rownames(average.manual.percc)=rownames(data)
dim(average.manual.percc)
for (cell in unique(target.percc$Cluster)){
  cells.clust=target.percc$Cell[which(target.percc$Cluster==cell)]
  average.manual.percc[,cell]=apply(data[,which(colnames(data)%in%cells.clust)], 1, mean)
}

average.clean.percc=average.manual.percc #average gene expression

#Applying Icellnet pipeline on cluster of interest
data.icell.percc=as.data.frame(gene.scaling(as.data.frame(average.clean.percc), n=1, db=db))

# label and color label if you are working families of molecules already present in the database

my.family=c("Growth factor","Chemokine","Checkpoint","Cytokine","Notch family","Antigen binding")
family.col = c( "Growth factor"= "#AECBE3", "Chemokine"= "#66ABDF", "Checkpoint"= "#1D1D18"  ,
                "Cytokine"="#156399", "Notch family" ="#676766", "Antigen binding" = "#12A039",  "other" = "#908F90",  "NA"="#908F90")

recpex = colnames(data.icell.percc[, !names(data.icell.percc) %in% c("Symbol")]) #receptors expressed
PC.data.percc=as.data.frame(data.icell.percc[], row.names = rownames(data.icell.percc)) #Peripheral cells

#my.colname = c(recpex)
PC.target.percc=data.frame("Class"=c(recpex), "ID"= c(recpex), "Cell_type"=c(recpex))
rownames(PC.target.percc)=c(recpex)
my.selection.percc=c(recpex)



for (labelname in recpex){
 
  PC.data.percc=as.data.frame(data.icell.percc[], row.names = rownames(data.icell.percc))
  PC.target.percc=data.frame("Class"=c(recpex), "ID"= c(recpex), "Cell_type"=c(recpex))
  rownames(PC.target.percc)=c(recpex)
  
  my.selection.percc=c(recpex)
  outdirection = "out"
  indirection = "in"
  
  #Score is computed for each cells, one cell at a time
  
  score.computation.perccin= icellnet.score(direction=indirection, PC.data=PC.data.percc, 
                                            CC.data= as.data.frame(data.icell.percc[,c(labelname)], row.names = rownames(data.icell.percc)),  
                                            PC.target = PC.target.percc, PC=my.selection.percc, CC.type = "RNAseq", 
                                            PC.type = "RNAseq",  db = db)
  score1.perccin=as.data.frame(score.computation.perccin[[1]])
  lr1.perccin=score.computation.perccin[[2]]
  
  score.computation.perccout= icellnet.score(direction=outdirection, PC.data=PC.data.percc, 
                                             CC.data= as.data.frame(data.icell.percc[,c(labelname)], row.names = rownames(data.icell.percc)),  
                                             PC.target = PC.target.percc, PC=my.selection.percc, CC.type = "RNAseq", 
                                             PC.type = "RNAseq",  db = db)
  score1.perccout=as.data.frame(score.computation.perccout[[1]])
  lr1.perccout=score.computation.perccout[[2]]
  
  Scoresio=cbind(score1.perccin,score1.perccout)
  colnames(Scoresio)=c(paste(paste(labelname),paste("(R)-in")),paste(paste(labelname),paste("(L)-out")))
  
  #prints in and out score 
  print(Scoresio) 
  
  #network plot
  
  #using combined scores for network plot
  #Score scaling
  Scores.norm=(Scoresio-min(Scoresio))/(max(Scoresio)-min(Scoresio))*9+1
  # Display intercellular communication networks
  
  #color scheme for network plot
  PC.col = c("cDC2"="#88b04b","LC"="#88b04b","Inflam. DC"="#88b04b","cDC1"="#88b04b",
             "Inflam. FIB"="#ff962c", "FBN1+ FIB"="#ff962c","APOE+ FIB"="#ff962c", "COL11A1+ FIB"="#ff962c",
             "CD40LG+ TC"="#c100b9","Inflam. TC"="#c100b9","TC"="#c100b9","NKT"="#c100b9")
  
  network.plot1 = network.create(icn.score = Scores.norm[1], scale = c(round(min(Scores.norm)),round(max(Scores.norm))), direction = "in", PC.col)
  network.plot2 = network.create(icn.score =Scores.norm[2], scale = c(round(min(Scores.norm)),round(max(Scores.norm))), direction = "out",PC.col)
  gridExtra::grid.arrange(network.plot1, network.plot2, ncol=2, nrow=1)
  
  
  
  #to define the y axis range of the barplot
  ymax=round(max(Scoresio))+1 
  
  #Compute the contribution of each family of molecules to the global communication scores + barplot (plot=T)
  contrib.family.1= LR.family.score(lr=lr1.perccin, my.family=my.family, db.couple=db.name.couple, plot= T, family.col=family.col, title=paste(paste(labelname),paste("(R)-in")))
  contrib.family.2= LR.family.score(lr=lr1.perccout, my.family=my.family, db.couple=db.name.couple, plot= T, family.col=family.col, title=paste(paste(labelname),paste("(L)-out")))
  grid.arrange(contrib.family.1, contrib.family.2, ncol=2, nrow=1)
  
  tabin = LR.family.score(lr=lr1.perccin, my.family=my.family, db.couple=db.name.couple, plot=F)
  tabout = LR.family.score(lr=lr1.perccout, my.family=my.family, db.couple=db.name.couple, plot=F)
  print("--------------------------LR family score(in)-----------------------------")
  print(tabin)
  print("--------------------------LR family score(out)----------------------------")
  print(tabout)
  
  #-------
  #balloon plot
  #For our code we only plotted based on sum as we wanted to check the most contributing factors 
  BP1 = LR.balloon.plot(lr = lr1.perccin, thresh = 0 , topn=30 , sort.by="sum",  db.name.couple=db.name.couple, family.col=family.col, title=paste("Most contributing interactions of",paste(labelname),paste("(by sum/in)"))) 
  #BP2 = LR.balloon.plot(lr = lr1.perccin, thresh = 0 , topn=30 , sort.by="var",  db.name.couple=db.name.couple, family.col=family.col, title=paste("Most contributing interactions of",paste(labelname),paste("(by var/in)")))
  grid.arrange(BP1)
  #grid.arrange(BP2)
  
  BP3 = LR.balloon.plot(lr = lr1.perccout, thresh = 0 , topn=30 , sort.by="sum",  db.name.couple=db.name.couple, family.col=family.col, title=paste("Most contributing interactions of",paste(labelname),paste("(by sum/out)")))
  #BP4 =LR.balloon.plot(lr = lr1.perccout, thresh = 0 , topn=30 , sort.by="var",  db.name.couple=db.name.couple, family.col=family.col, title=paste("Most contributing interactions of",paste(labelname),paste("(by var/out)")))
  grid.arrange(BP3)
  #grid.arrange(BP4)
  
  #plotting the ballon plots with the sum of in and out direction scores
  lrsum = lr1.perccin + lr1.perccout
  BP5 = LR.balloon.plot(lr = lrsum, thresh = 0 , topn=30 , sort.by="sum",  db.name.couple=db.name.couple, family.col=family.col, title=paste("Most contributing interactions of",paste(labelname),paste("(by sum)")))
  grid.arrange(BP5)
  
  #function to plot heatmap with combined sum
  HP = LR.heatmap(lr = lrsum, thresh = 0, topn = 30, sort.by = "sum",  db.name.couple = db.name.couple, title = paste("Most contributing interactions of",paste(labelname),paste("(by sum)")), family.col = family.col, value_display=50)  
  grid.arrange(HP)
  
}


#dev.off()