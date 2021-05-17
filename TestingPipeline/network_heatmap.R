#!/usr/bin/env Rscript

library(igraph)
library(RColorBrewer)
library(matrixStats)
library(ggplot2)
library(scAI)

# Loading in Network
network = read.csv("count_network.txt", sep = "\t", header = TRUE)

# Converting to adjacency matrix
get.adjacency <- graph.data.frame(network)

adjaceny.m = get.adjacency(get.adjacency, sparse = FALSE, attr='count')
adjaceny.m <- adjaceny.m[c("Plasma","WNT2B+ Fos-hi", "CD4+ Memory", "Macrophages",   "Enterocytes")
                         ,c("Plasma","WNT2B+ Fos-hi", "CD4+ Memory",   "Macrophages",   "Enterocytes")]

print(adjaceny.m)
# From adjacency make graph
g<-graph_from_adjacency_matrix(adjaceny.m, mode = "undirected", weighted= T)

#Setting edge width, weight limits
edge.width.max <- 8
edge.weight.max <- max(E(g)$weight)
# Scaling edge width
E(g)$width<-0.3+E(g)$weight/edge.weight.max*edge.width.max
#Setting colors vertices
color.use <- scPalette(length(unique(network[,1])))
V(g)$color <- color.use
# Setting vertex size
deg <- degree(g)
print(deg)
#V(g)$size <- deg*5
deg <- rowSums2(adjaceny.m[V(g)$name,]) / 3
p <- layout.circle(g)

jpeg("network.jpg", 
     width=6.8, height=6.8, 
     units='in',res=300)
plot(g,main="Interaction Network",vertex.size = deg ,layout=p, 
     rescale=T, edge.curved=0.1, margin=0, 
     title.name="Lap", 
     vertex.label.dist=c(3,2.2,3,-4.5,-4),
     vertex.label.color="black" ,vertex.label.degree = c(-0.9,-2,-2,5,4), 
     vertex.label.cex=1.2, vertex.label.family="helvitica" )
dev.off()

ggplot(network, aes(x=SOURCE, y=TARGET, fill=log(count))) + 
  geom_tile() +scale_fill_gradient(low = "white", high = "steelblue") + 
  ggtitle("Number of Interactions Cell to Cell")+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
  axis.title.x=element_blank(),
  axis.text.x = element_text(angle =45,hjust=1),
  axis.title.y= element_blank())
ggsave("heatmap.pdf")


