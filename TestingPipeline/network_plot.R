library(igraph)
library(RColorBrewer)
library(matrixStats)
library(ggplot2)
library(scAI)

setwd("/students/2020-2021/master/CellPhoneDB/TestingPipeline/CellChat/out/out/")

# Loading in Network
network = read.csv("count_network.txt", sep = "\t", header = TRUE)

# Converting to adjacency matrix
get.adjacency <- graph.data.frame(network)
# Make adjacency
adjaceny.m = get.adjacency(get.adjacency, sparse = FALSE, attr='count')
# Make to graph
g <- graph_from_adjacency_matrix(adjaceny.m, mode = "undirected", weighted= T)
# Setting position loops
edge.start <- ends(g, es = E(g), names = F)
coords <- layout_(g, in_circle())
edge.start <- ends(g, es = E(g), names = F)


# This code section is pretty advanced, cannot really say whats happening here.
# This code adjusts the label position.
if(nrow(coords)!=1){
  coords_scale=scale(coords)
}else{
  coords_scale<-coords
}
loop.angle<-ifelse(coords_scale[igraph::V(g),1]>0,-atan(coords_scale[igraph::V(g),2]/coords_scale[igraph::V(g),1]),
                   pi-atan(coords_scale[igraph::V(g),2]/coords_scale[igraph::V(g),1]))
if(sum(edge.start[,2]==edge.start[,1])!=0){
  igraph::E(g)$loop.angle[which(edge.start[,2]==edge.start[,1])]<-loop.angle[edge.start[which(edge.start[,2]==edge.start[,1]),1]]
}
radian.rescale <- function(x, start=0, direction=1) {
  c.rotate <- function(x) (x + start) %% (2 * pi) * direction
  c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
}

vertex.size.max <- 15
vertex.weight <- 20
vertex.weight.max <- max(vertex.weight)
vertex.weight <- vertex.weight/vertex.weight.max*vertex.size.max+5
label.locs <- radian.rescale(x=1:length(igraph::V(g)), direction=-1, start=0)
label.dist <- vertex.weight/max(vertex.weight)+2
#Setting edge width, weight limits
edge.width.max <- 8
edge.weight.max <- max(E(g)$weight)
# Scaling edge width
E(g)$width<-0.3+(E(g)$weight/edge.weight.max*edge.width.max)
#Setting colors vertices
color.use <- scPalette(length(unique(network[,1])))
V(g)$color <- color.use
# Setting vertex size
deg <- degree(g)
print(deg)
# ADJUST THIS LINE FOR NODE SIZE
deg <- rowSums2(adjaceny.m[V(g)$name,])  /5


# Network
network_graph <- plot(g,main="Interaction Network",vertex.size = deg, 
                          layout=coords, 
                          rescale=T, edge.curved=0.1, margin=0, 
                          title.name="Lap", 
                          vertex.label.degree=label.locs, 
                          vertex.label.dist=label.dist,
                          vertex.label.family="helvitica")

ggsave("network_graph.png", plot(g,main="Interaction Network",vertex.size = deg, 
                                 layout=coords, 
                                 rescale=T, edge.curved=0.1, margin=0, 
                                 title.name="Lap", 
                                 #vertex.label.dist=c(3,2.2,3,-4.5,-4),
                                 #vertex.label.color="black" ,vertex.label.degree = c(-0.9,-2,-2,5,4), 
                                 vertex.label.degree=label.locs, 
                                 vertex.label.dist=label.dist,
                                 vertex.label.family="helvitica"))

# Heatmap log(count)
ggplot(network, aes(x=SOURCE, y=TARGET, fill=log(count))) + 
  geom_tile() +scale_fill_gradient(low = "white", high = "steelblue") + 
  ggtitle("Number of Interactions Cell to Cell")+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title.x=element_blank(),
        axis.text.x = element_text(angle =45,hjust=1),
        axis.title.y= element_blank())
ggsave("heatmap_interaction_log.png", width = 7)

# Heatmap count
ggplot(network, aes(x=SOURCE, y=TARGET, fill=count)) + 
  geom_tile() +scale_fill_gradient(low = "white", high = "steelblue") + 
  ggtitle("Number of Interactions Cell to Cell")+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title.x=element_blank(),
        axis.text.x = element_text(angle =45,hjust=1),
        axis.title.y= element_blank())
ggsave("heatmap_interaction_count.png", width = 7)


