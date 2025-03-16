
### SensorID dynamic interactomes Network analysis
## This R Scripts contains the analysis of the SensorID dynamic interactome as networks using the R package igraph
#####
#Install the required packages
install.packages(c("igraph", "bipartite", "assortnet", "ggplot2", "igraphdata", "RColorBrewer", "tidyverse"))

#Load Packages
require(igraph)
require(rnetcarto)
require(bipartite)
require(assortnet)
require(igraphdata)
require(RColorBrewer)
require(tidyverse)
require(ggplot2)
require(dplyr)
require(Matrix)
require(ggraph)

#Set the Working Dir
setwd("C:/Users/alfredog/Desktop/SensorID_Manuscript/GitHub_Scripts/SensorID_Network_Analysis")

#Load the required files
#Load the protein annotation data, including the FoldChange values per time point
Mat<-read.csv("./SensorID_Matrix.csv")
#Format a s factors the categorical variables needed for the analysis.
Mat$DNA_Repair<-as.factor(Mat$DNA_Repair)
Mat$Annotation<-as.factor(Mat$Annotation)
Mat$An.color<-as.factor(Mat$An.color)
str(Mat)
#Node files
Mat#This may be the best node attribute list that we may have

#Lets first start with the Full Nuclear Network (FNN), that contains all the nodes of the SensorID networks. 
#The FNN will allow to settle the possitions and layout for all the nodes, and use these positions and layouts to all the SensorID networks.
###SensorID FNN
#SensorID All nodes Network Edge Lists
SensorID_FNN<-read.delim("./SensorID_FNN.tsv",sep = "\t",col.names = )#Load the edge list of the SensorID FNN, previously downloaded from String database
str(SensorID_FNN)#Review the information contained in the edge list
head(SensorID_FNN)

#Generate the igraph network object
FNN<-graph_from_data_frame(d=SensorID_FNN,vertices = Mat, directed = F )#Generate the igraph object with the edge and node annotation attributes
str(vertex_attr(FNN))#Check node attributes
str(edge_attr(FNN))#Check for edge attributes
FNN<-igraph::simplify(FNN, remove.multiple=TRUE, remove.loops=TRUE) #Remove repetitive edges and self-connecting edges
FNN #Review the information contained in the FNN network
#Small pilot plot
plot(FNN)
#Generate an adjacency Matrix that will be used after
FNN_AdjMat<-as_adjacency_matrix(FNN,type="both",attr = "weight",names=TRUE,sparse = TRUE)
FNN_AdjMat # The adjacency matrix is another format of the networks, and can be used to easily visualize the edges, this matrix can be printed as a CSV

#Generate a vertex size based on degree and count for plotting of the networks
vsize<-igraph::degree(FNN)/vcount(FNN)
#Identify and remove disconnected nodes
isolated<-which(igraph::degree(FNN)==0)
isolated
length(isolated)# number of unconnected nodes
#remove isolated nodes
FNN2<-igraph::delete_vertices(FNN,isolated)

#Generate layouts for the network graphs, the FNN layouts should be kept for the rest of the networks, to maintain the same position to each node in all of them.
Lyt.Nice<-layout_nicely(FNN,dim = 3)
Lyt.GOpt<-layout_with_graphopt(FNN, niter=5000, mass = 30, charge = E(FNN)$weight)
Lyt.DrL<-layout_with_drl(FNN, weights = E(FNN)$weight)
Lyt.FR<-layout.fruchterman.reingold(FNN, niter=10000, weights=E(FNN)$weight)#This one looks like the best option
Lyt.KK<-layout_with_kk(FNN, maxiter=10000, weights=E(FNN)$weight) #very large separation of peripheral nodes
Lyt.GEM<-layout_with_gem(FNN, maxiter=10000) #Bad separation
Lyt.LGL<-layout.lgl(FNN, maxit=150,maxdelta=833,area=(length(E(FNN)))^2,coolexp=1.5,repulserad=((length(E(FNN)))^2)*(length(E(FNN))))#Do not work

#Star and Tree layouts, to place the DSB-sensor proteins in the center of the networks
V(FNN2)[name=="SIRT6"]
which(V(FNN2)[name=="SIRT6"])# 610 #Check what is the row index for your desired node
which(V(FNN2)$name=="SIRT6")# 610 #Check what is the row index for your desired node
which(V(FNN2)$name=="MRE11A")# 327
which(V(FNN2)$name=="XRCC5")# 722
Lyt.Star<-layout_as_star(FNN2,center=V(FNN2)[name=="SIRT6"]) #it has a special order argument where we can assign an order to the vertex 
Lyt.Tree<-layout_as_tree(FNN2,root=c(610))
Lyt.RT<-layout.reingold.tilford(FNN2,root=c(610))# This one looks good, but the indexing do not allow to highlight the correct nodes

#Check what layout looks fine
plot.igraph(FNN2,layout=Lyt.Tree[-isolated,], alpha=0.6,
            vertex.color=V(FNN)$An.color, vertex.frame.color=V(FNN)$An.color, vertex.shape="circle",
            vertex.size=5, vertex.label=V(FNN2)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(FNN)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)

#With ggraph, it uses ggplot grammar and graphics for network plotting

ggraph(FNN2,layout=Lyt.FR[-isolated,])+
  geom_edge_fan(color="gray85", aes(width=(E(FNN2)$weight)/3), alpha=0.4)+
  geom_node_point(aes(color=V(FNN2)$Annotation), size=10,alpha=0.75)+
  geom_node_text(aes(label=name),color="gray25",size=4)+
  scale_color_manual(limits=V(FNN2)$Annotation, values=V(FNN2)$An.color)+
  theme_void()

#Community Detection analysis
clp <- igraph::cluster_optimal(FNN2)
class(clp)
V(FNN2)$community <- clp$membership
colrs <- adjustcolor( c("gray50", "tomato", "gold", "yellowgreen"), alpha=.6)

#Plot the networks with the community analysis highlighting the associated clusters of nodes.
plot.igraph(clp,FNN2,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(FNN)$An.color, vertex.frame.color=V(FNN)$An.color, vertex.shape="circle",
            vertex.size=5, vertex.label=V(FNN)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(FNN2)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)

#Calculate the distances from all the nodes in the networks to the DSB-Sensors, a normalized number of links between each node and the DSB-sensor.
dist.from.Sirt6 <- distances(FNN2, v=V(FNN2)[name=="SIRT6"],
                             to=V(FNN2), weights=E(FNN2)$weight)
dist.from.Sirt6<-dist.from.Sirt6[,!is.infinite(colSums(dist.from.Sirt6))] # Remove those nodes whose distance result into an infinite value, disconnected nodes.
dist.from.Sirt6# Review the distance values, this results can be exported as a measure of the closseness of each node to the DSB-Sensors
range(dist.from.Sirt6) 

#Set colors to plot the distances:
oranges <- colorRampPalette(c("dark red", "gold"))
col <- oranges(max(dist.from.Sirt6))
col <- col[dist.from.Sirt6+1]
plot(FNN2, layout=Lyt.FR[-isolated,],vertex.color=col, vertex.label=dist.from.Sirt6,
     vertex.label.color="white",vertex.shape="circle",
     vertex.size=5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(FNN2)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)

#Shortest path from Sirt6 to Nucleolin, this can be extracted from any other pair of nodes in the networks. It will highlight what nodes and edges are implicated in connecting a pair of proteins.
S6Ncl.path <- shortest_paths(FNN2,
                             from = V(FNN2)[name=="SIRT6"],
                             to = V(FNN2)[name=="NCL"],
                             output = "both") # both nodes and edges in the path
#Generate edge color variable to plot the path:
ecol <- rep("gray80", ecount(FNN2))
ecol[unlist(S6Ncl.path$epath)] <- "orange"

# Generate edge width variable to plot the path:
ew <- rep(2, ecount(FNN2))
ew[unlist(S6Ncl.path$epath)] <- 8

# Generate node color variable to plot the path:
vcol <- rep("gray40", vcount(FNN2))
vcol[unlist(S6Ncl.path$vpath)] <- "gold"
plot(FNN2, vertex.color=vcol, edge.color=ecol,
     edge.width=ew,vertex.label.color="white",vertex.shape="circle",
     vertex.size=5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(FNN2)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)

#Plot incident edges and neighbor nodes to Sirt6. This can be done to any other node in the networks.
inc.edges <- incident(FNN2, V(FNN2)[name=="SIRT6"], mode="all")
# Set colors to plot the incident edges and neighboring nodes.
ecol <- rep("gray80", ecount(FNN2))
ecol[inc.edges] <- "orange"
vcol <- V(FNN2)$An.color
vcol[V(FNN2)$name=="SIRT6"] <- "gold"
ew <- rep(2, ecount(FNN2))
ew[inc.edges] <- 8

#Plot the networks highlighting 
plot(FNN2, vertex.color=vcol, edge.color=ecol, layout=Lyt.FR[-isolated,],
     edge.width=ew,vertex.label.color="white",vertex.shape="circle",
     vertex.size=vsize, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(FNN2)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)

#Plot Degree Distribution in the network
FNN.Deg<-degree_distribution(FNN, cumulative = T, mode="all")
plot( x=0:max(igraph::degree(FNN)), y=1-FNN.Deg, pch=19, cex=1.2, col="orange",
      xlab="Degree", ylab="Cumulative Frequency")

#Network descriptive metrics
#Number of nodes and edges
vcount(FNN2) # number of connected nodes
ecount(FNN2)# Number of interactions, edges.
#Density (Number of edges in relation to the total number of edges in a full connected graph)
edge_density(FNN,loops = F)
#Transitivity (Total number of connected triangles)
transitivity(FNN,type="global") 
#modularity, the tendency of the nodes with common traits (similar functional group, or Fuzzy C means cluster), to associate in a common module.
mod<-V(FNN)$An.number #lets check for functional annotation associated modularity.
mod<-mod+1 #The modularity algorithm do not process 0s, 0s needs to be transformed to 1, so we add 1 to every cluster to avoid errors.
modularity(FNN,membership = mod)#Calculate the modularity associated with the functional annotation.

#Assorativity, accounts for the tendency of nodes with similar properties to interact among them. Hence, we can try several characteristics of the nodes and measure the assortativity.
assortativity_nominal(FNN2,factor(V(FNN2)$Annotation),directed=F) # Measuring the assortativity in relationship to the functional annotations.

#Diameter (longest geodesic distance), shortest path between two most peripheral nodes
diameter(FNN2,directed = F, weights = E(FNN)$weights) #calculate the diameter using weighted edges.
diam<-get_diameter(FNN,directed = F,weights = E(FNN)$weights) # Retrieve the shortest path that represent this diameter. 
diam # Printing the diam file will show what nodes participate in this path 
#Plot Diameter
# highlight the nodes and edges in the diameter of the network
vcol <- rep("gray50", vcount(FNN2))
vcol[diam] <- "gold"
ecol <- rep("gray90", ecount(FNN2))
ecol[E(FNN2, path=diam)] <- "orange"
#Plot the networks with the diameter highlighted.
plot(FNN2, vertex.color=vcol, edge.color=ecol, vertex.label.color="black",edge.arrow.mode=0,layout=Lyt.FR)

# Calculate Node-based Network metrics.
degree<-igraph::centr_degree(FNN,mode = "all",normalized = T) #Calculate the degree, # Edges, for every node in the network.
degree<-degree$res# Recover the degree results. 
eigen_cent<-igraph::eigen_centrality(FNN) # Calculate the eigencentrality measurement per node
eigen_cent<-eigen_cent$vector # Recover the eigencentrality results
clos_cent <- igraph::closeness(FNN, normalized=T) # Calculate the closseness measurement per node
betw_cent <- igraph::betweenness(FNN, normalized=T) # Calculate the betweness measurement per node
hs <- hub_score(FNN, weights=E(FNN)$weight)$vector # Calculate the hub score  measurement per node# 
as <- authority_score(FNN, weights=E(FNN)$weight)$vector # Calculate the authority score measurement per node
Trans<-transitivity(FNN,type="local") # calculate the local transitivity of every node
Trans
# Construct a final table and print the results with all the node-level metrics. These metrics will be manually included into the SensorID interactome results tables.
FNN.Centrality<-data.frame(NAME=V(FNN)$name, UniprotID=V(FNN)$Uniprot, Degree=degree, Eigen= eigen_cent, Closeness= clos_cent, Betweeness= betw_cent, Hubs=hs, Authority=as,Transitivity=Trans)
FNN.Centrality
write.csv(FNN.Centrality,"./FNN.Centrality.csv")

#Cocitation, Number of shared interactions among nodes
Cocitation<-cocitation(FNN)
Cocitation#This matrix is an adjacency matrix showing how many shared interactors two proteins have. It can be saved as a csv file.

#Plot with degree as node size
ggraph(FNN2,layout=Lyt.FR[-isolated,])+
  geom_edge_fan(color="gray85", aes(width=(E(FNN2)$weight)/3), alpha=0.4)+
  geom_node_point(aes(color=V(FNN2)$Annotation, size=10+degree),alpha=0.75)+
  geom_node_text(aes(label=name),color="gray25",size=4)+
  scale_color_manual(limits=V(FNN2)$Annotation, values=V(FNN2)$An.color)+
  theme_void()

#Modularity Analysis: Edge-Betweeness and Louvain Modularity.
#These modularity analysis gave too many modules, and do not calculate participation and connectivity metrics. Hence the rNetcarto modularity was done in a separate scrtipt.

#Edge_Betweeness Modularity
FNN.ebmod=edge.betweenness.community(FNN2,weights=E(FNN2)$weight,directed = F) 
FNN.ebmod# detect 132 groups, just too many
# Plot the Edge-betweeness calculated modules
plot.igraph(FNN.ebmod,FNN2,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(FNN)$An.color, vertex.frame.color=V(FNN)$An.color, vertex.shape="circle",
            vertex.size=5, vertex.label=V(FNN)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(FNN)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)
#Louvain modularity
FNN.Louvain.mod<-cluster_louvain(FNN2,weights=E(FNN2)$weight)
FNN.Louvain.mod #The result shows a total of 17 clusters with 0.53 modularity value

#Plot the Louvain modules
plot(FNN.Louvain.mod, FNN2,layout=Lyt.FR[-isolated,], alpha=0.6, 
     vertex.color=V(FNN2)$An.color, vertex.frame.color=V(FNN2)$An.color, vertex.shape="circle", 
     vertex.size=5, vertex.label=V(FNN2)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black", 
     edge.width=E(FNN2)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout", 
     rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)

#Generate HeatMaps of the adjacency matrixes
#Let's use the neighbor list for generating the Adjacency Matrix
SensorID_FNN
v<-levels(as.factor(unlist(SensorID_FNN[,1:2]))) #get vertex names from the network list
n<-length(v)#get the number of vertex in the network
e<-matrix(match(as.character(unlist(SensorID_FNN[,1:2])),v),ncol = 2)#Edgelist with row index instead of Names
w<-SensorID_FNN$weight #Edge weights

# M is our co-expression network adjacency matrix.

# distringuish between i,j and j,i we want to make M symmetric (undirected) by considering an undirected network.
FNN.M <- matrix(0, n, n)                  # set up a co-expression matrix
FNN.M[e] <- w                             # fill it in with edge weights
FNN.M <- FNN.M + t(FNN.M)                         # make this symmetric
dimnames(FNN.M) <- list(v, v)             # label the vertices
# We let A_ij = 1 if there is an interaction. 
#A <- 1*(M > 0) #right now we keep the weights as they are
FNN.M #The adjacency Matrix not as an igraph object, easier to manipulate, and to subset for smaller Sensor specific networks
quantile(FNN.M,c(0,0.25,0.5,0.75,1))

#Creating the HeatMap
# make a heatmap image with colors: white below 0,
# grey from 0.4 to 1 with darker closer to 1
breaks <- c(0, seq(0.4, 1.98, length=11)) # breaks in the color bins
cols <- grey(1-c(0,seq(0.5,1,length=10))) # colors

# color the vertex for SIRT6, MRE11A and Ku80 blue, its neighbors red, and the others white
SensorIDBind <- which(v == "SIRT6"|v == "MRE11A"|v == "XRCC5") # index for DSB-Sensors
vcols <- rep("white",n)
vcols[SensorIDBind] <- "blue"
#Subsetting FNN.M by doing a loop asking for look into every row and column the index in SensorIDBind worked
for (i in SensorIDBind) {
  vcols[which(FNN.M[,i] > 0 | FNN.M[i,])] <- "red"
}

#Generate a dendogram for keeping the proteins in the sample place
FNN.dist<-dist(FNN.M,method = "euclidian")
FNN.Dend<-hclust(FNN.dist, method = "ward.D2" )
plot(FNN.Dend)
FNN.Dend#Save the FNN dendogram to fix the order of the columns and rows in all the posterior SensorID heatmaps

#Cut the dendogram if the branches are too many
cut.height <- 10 # try numbers from 4 to 7
FNN.hclusts<- cutree(FNN.Dend, cut.height)
plot(FNN.hclusts)

# now actually make the heat map
par(mar=rep(0,4)) # make plot margins 0
heatmap(FNN.M, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=vcols, RowSideColors=vcols,  col=0:1, frame=T)+
  legend("topleft", c("Neighbors(SIRT6, MRE11, KU80)","SIRT6, MRE11, XRCC5"), fill=c("red","blue"),  
         bty="n", inset=0, xpd=T,  border=F)
#
dev.off()

################################################Sirt6ID Networks
#Sirt6ID Networks

#Let's keep the Mat as the Node attribute table and load the Edge list as the networks
Mat# Node attribute matrix for all the networks

#Sirt6ID nonIR5 nodes Network Edge Lists
S6_nonIR5<-read.delim("./Sirt6_Files/S6_nonIR5.tsv",sep = "\t",col.names = )
str(S6_nonIR5)
head(S6_nonIR5)

#Generate the igraph network object
S6.nonIR5<-graph_from_data_frame(d=S6_nonIR5,vertices = Mat, directed = F )#Generate the igraph object with the edge and node annotation attributes
str(vertex_attr(S6.nonIR5))#Check node attributes
str(edge_attr(S6.nonIR5))#Check for edge attributes
S6.nonIR5<-igraph::simplify(S6.nonIR5, remove.multiple=TRUE, remove.loops=TRUE) #Remove repetitive edges and self-connecting edges
S6.nonIR5
#Small plot
plot(S6.nonIR5)

#Generate an adjacency Matrix that will be used after
S6.nonIR5_AdjMat<-as_adjacency_matrix(S6.nonIR5,type="both",attr = "weight",names=TRUE,sparse = TRUE)
S6.nonIR5_AdjMat<-as.matrix(S6.nonIR5_AdjMat)
write.csv(S6.nonIR5_AdjMat,"./S6.nonIR5_AdjMat.csv")

#Generate a vertex size based on degree and count
vsize<-igraph::degree(S6.nonIR5)/vcount(S6.nonIR5)
#Identify and remove disconnected nodes
isolated<-which(igraph::degree(S6.nonIR5)==0)
isolated
length(isolated)# Number of disconnected nodes and nodes not present in the network.
#remove isolated nodes
S6.nonIR5<-igraph::delete_vertices(S6.nonIR5,isolated)

#Generate layouts for the graph
#Do not run again, to not to change the layouts that we already have
Lyt.Nice<-layout_nicely(S6.nonIR5,dim = 3)
Lyt.GOpt<-layout_with_graphopt(S6.nonIR5, niter=5000, mass = 30, charge = E(S6.nonIR5)$weight)
Lyt.DrL<-layout_with_drl(S6.nonIR5, weights = E(S6.nonIR5)$weight)
Lyt.FR<-layout.fruchterman.reingold(S6.nonIR5, niter=10000, weights=E(S6.nonIR5)$weight)#This one looks like the best option
Lyt.KK<-layout_with_kk(S6.nonIR5, maxiter=10000, weights=E(S6.nonIR5)$weight) #very large separation of peripheral nodes
Lyt.GEM<-layout_with_gem(S6.nonIR5, maxiter=10000) #Bad separation
Lyt.LGL<-layout.lgl(S6.nonIR5, maxit=150,maxdelta=833,area=(length(E(S6.nonIR5)))^2,coolexp=1.5,repulserad=((length(E(S6.nonIR5)))^2)*(length(E(S6.nonIR5))))#Do not work

#Star and Tree layouts
V(S6.nonIR5)[name=="SIRT6"]
which(V(S6.nonIR5)[name=="SIRT6"])# 610 #Check what is the rowindex for your desired node
which(V(S6.nonIR5)$name=="SIRT6")# 610 #Check what is the rowindex for your desired node
which(V(S6.nonIR5)$name=="MRE11A")# 327
which(V(S6.nonIR5)$name=="XRCC5")# 722
Lyt.Star<-layout_as_star(S6.nonIR5,center=V(S6.nonIR5)[name=="SIRT6"]) #it has a special order argument where we can assign an order to the vertex, hence maybe the OPLS-DA:VIP value order may be used as ordering vector

Lyt.Tree<-layout_as_tree(S6.nonIR5,root=c(190))
Lyt.RT<-layout.reingold.tilford(S6.nonIR5,root=c(190))# This one looks good, but the indexing do not allow to highlight the correct nodes
Lyt.sugi<-layout_with_sugiyama(S6.nonIR5, maxiter=100, weights=E(S6.nonIR52)$weight)#This one looks good for undirected tree layouts

#Check what layout looks fine

plot.igraph(S6.nonIR5,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(S6.nonIR5)$An.color, vertex.frame.color=V(S6.nonIR5)$An.color, vertex.shape="circle",
            vertex.size=5, vertex.label=V(S6.nonIR5)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(S6.nonIR5)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)

#With ggraph
ggraph(S6.nonIR5,layout=Lyt.FR[-isolated,])+
  geom_edge_fan(color="gray65", width=1, alpha=0.4)+
  geom_node_point(aes(color=V(S6.nonIR5)$Annotation, size=V(S6.nonIR5)$S6.FC.nonIR5),alpha=0.75)+
  geom_node_text(aes(label=name),color="gray25",size=4)+
  scale_color_manual(limits=V(S6.nonIR5)$Annotation, values=V(S6.nonIR5)$An.color)+
  scale_size_continuous(range = c(5, 15)) +
  theme_void()+
  title("Sirt6 nonIR5 Network")+
  guides(color=guide_legend(override.aes = list(size=7)))+
  theme(legend.text=element_text(size=12,face = "bold"),
        legend.title = element_text(size=13,face = "bold"))

#Community Detection
clp <- igraph::cluster_optimal(S6.nonIR5)
class(clp)
V(S6.nonIR5)$community <- clp$membership
colrs <- adjustcolor( c("gray50", "tomato", "gold", "yellowgreen"), alpha=.6)

plot.igraph(clp,S6.nonIR5,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(S6.nonIR5)$An.color, vertex.frame.color=V(S6.nonIR5)$An.color, vertex.shape="circle",
            vertex.size==V(S6.nonIR5)$S6.FC.nonIR5, vertex.label=V(S6.nonIR5)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(S6.nonIR5)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)+
legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)

#Distances edges highlight
dist.from.Sirt6 <- distances(S6.nonIR5, v=V(S6.nonIR5)[name=="SIRT6"],
                             to=V(S6.nonIR5), weights=E(S6.nonIR5)$weight)
dist.from.Sirt6<-dist.from.Sirt6[,!is.infinite(colSums(dist.from.Sirt6))]
dist.from.Sirt6
range(dist.from.Sirt6)

#Set colors to plot the distances:
oranges <- colorRampPalette(c("dark red", "gold"))
col <- oranges(max(dist.from.Sirt6))
col <- col[dist.from.Sirt6+1]
#Plot the network with the color based on its distance to Sirt6
plot(S6.nonIR5, layout=Lyt.FR[-isolated,],vertex.color=col, vertex.label=dist.from.Sirt6,
     vertex.label.color="white",vertex.shape="circle",
     vertex.size=V(S6.nonIR5)$S6.FC.nonIR5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(S6.nonIR5)$weight, edge.lty=1, main="Sirt6 nonIR5 Full Nuclear Network", sub="FR Layout",
     rescale=T)

#Shortest path from Sirt6 to Nucleolin
S6Ncl.path <- shortest_paths(S6.nonIR5,
                             from = V(S6.nonIR5)[name=="SIRT6"],
                             to = V(S6.nonIR5)[name=="NCL"],
                             output = "both") # both path nodes and edges
# Generate edge color variable to plot the path:
ecol <- rep("gray80", ecount(S6.nonIR5))
ecol[unlist(S6Ncl.path$epath)] <- "orange"
# Generate edge width variable to plot the path:
ew <- rep(2, ecount(S6.nonIR5))
ew[unlist(S6Ncl.path$epath)] <- 4
# Generate node color variable to plot the path:
vcol <- rep("gray40", vcount(S6.nonIR5))
vcol[unlist(S6Ncl.path$vpath)] <- "gold"
plot(S6.nonIR5, layout=Lyt.FR[-isolated,], vertex.color=vcol, edge.color=ecol,
     edge.width=ew,vertex.label.color="white",vertex.shape="circle",
     vertex.size=5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(S6.nonIR5)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)

#Plot based on Sirt6 as incident
inc.edges <- incident(S6.nonIR5, V(S6.nonIR5)[name=="SIRT6"], mode="all")
# Set colors to plot the selected edges.
ecol <- rep("gray80", ecount(S6.nonIR5))
ecol[inc.edges] <- "orange"
vcol <- rep("gray40", vcount(S6.nonIR5))
vcol[V(S6.nonIR5)$name=="SIRT6"] <- "gold"

#Plot the network highlighting the incident edges and neighbor nodes
plot(S6.nonIR5, layout=Lyt.FR[-isolated,], vertex.color=vcol, edge.color=ecol, 
     edge.width=ew,vertex.label.color="white",vertex.shape="circle",
     vertex.size=V(S6.nonIR5)$S6.FC.nonIR5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(S6.nonIR5)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)

#Plot Degree Distribution
S6.nonIR5.Deg<-degree_distribution(S6.nonIR5, cumulative = T, mode="all")
plot( x=0:max(igraph::degree(S6.nonIR5)), y=1-S6.nonIR5.Deg, pch=19, cex=1.2, col="orange",
      xlab="Degree", ylab="Cumulative Frequency")

#Network Descriptives
#Number of nodes and edges
vcount(S6.nonIR5)
ecount(S6.nonIR5)
#Density (Number of edges in relation to the total number of edges in a full connected graph)
edge_density(S6.nonIR5,loops = F) 
#Transitivity (Total number of conected triangles)
transitivity(S6.nonIR5,type="global") 
#modularity
mod<-V(S6.nonIR5)$Annotation #lets check for fuzzy Cmeans modularity
mod<-mod+1 #The non clustered interactors that appears as 0, needs to be transformed to 1, so we add 1 to every cluster
modularity(S6.nonIR5,membership = mod) 

#Assorativity
assortativity_nominal(S6.nonIR5,factor(V(S6.nonIR5)$Annotation),directed=F)

#Diameter (longest geodesic distance), shortest path between two most peripheral nodes
diameter(S6.nonIR5,directed = F, weights = E(S6.nonIR5)$weights)#
diam<-get_diameter(S6.nonIR5,directed = F,weights = E(S6.nonIR5)$weights)
diam #
#Plot Diameter
vcol <- rep("gray50", vcount(S6.nonIR5))
vcol[diam] <- "gold"
ecol <- rep("gray90", ecount(S6.nonIR5))
ecol[E(S6.nonIR5, path=diam)] <- "orange"
# E(net, path=diam) finds edges along a path, here 'diam'
plot(S6.nonIR5, vertex.color=vcol, edge.color=ecol, vertex.label.color="black",edge.arrow.mode=0,layout=Lyt.FR[-isolated,],vertex.size=V(S6.nonIR5)$S6.FC.nonIR5)

#Node-level network metrics
degree<-igraph::centr_degree(S6.nonIR5,mode = "all",normalized = T)
degree<-degree$res
eigen_cent<-igraph::eigen_centrality(S6.nonIR5)
eigen_cent<-eigen_cent$vector
clos_cent <- igraph::closeness(S6.nonIR5, normalized=T)
betw_cent <- igraph::betweenness(S6.nonIR5, normalized=T)
hs <- hub_score(S6.nonIR5, weights=E(S6.nonIR5)$weight)$vector
as <- authority_score(S6.nonIR5, weights=E(S6.nonIR5)$weight)$vector
Trans<-transitivity(S6.nonIR5,type="local")
Trans
#write results
S6.nonIR5.Centrality<-data.frame(NAME=V(S6.nonIR5)$name, UniprotID=V(S6.nonIR5)$Uniprot, Degree=degree, Eigen= eigen_cent, Closeness= clos_cent, Betweeness= betw_cent, Hubs=hs, Authority=as,Transitivity=Trans)
S6.nonIR5.Centrality
write.csv(S6.nonIR5.Centrality,"./S6.nonIR5.Centrality.csv")

#Cocitation, Number of shared interactions among nodes
Cocitation<-cocitation(S6.nonIR5)
Cocitation# shared interactors among proteins
write.csv(Cocitation,"./S6.nonIR5.cocitation.csv")
#Plot with degree as node size
ggraph(S6.nonIR5,layout=Lyt.FR[-isolated,])+
  geom_edge_fan(color="gray85", width=1, alpha=0.4)+
  geom_node_point(aes(color=V(S6.nonIR5)$Annotation, size=degree),alpha=0.75)+
  geom_node_text(aes(label=name),color="gray25",size=4)+
  scale_color_manual(limits=V(S6.nonIR5)$Annotation, values=V(S6.nonIR5)$An.color)+
  scale_size_continuous(range = c(5, 15)) +
  theme_void()+
  title("Sirt6 nonIR5 Network")+
  guides(color=guide_legend(override.aes = list(size=7)))+
  theme(legend.text=element_text(size=12,face = "bold"),
        legend.title = element_text(size=13,face = "bold"))


#Modularity Analysis
#Edge_Betweeness Modularity
S6.nonIR5.ebmod=edge.betweenness.community(S6.nonIR5,weights=E(S6.nonIR5)$weight,directed = F) 
S6.nonIR5.ebmod #37 groups, good enough mod 0.47
#
plot(S6.nonIR5.ebmod, S6.nonIR5,layout=Lyt.FR[-isolated,], alpha=0.6, 
     vertex.color=V(S6.nonIR5)$An.color, vertex.frame.color=V(S6.nonIR5)$An.color, vertex.shape="circle", 
     vertex.size=5, vertex.label=V(S6.nonIR5)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black", 
     edge.width=E(S6.nonIR5)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout", 
     rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)
#Louvain modularity
S6.nonIR5.Louvain.mod<-cluster_louvain(S6.nonIR5,weights=E(S6.nonIR5)$weight)
S6.nonIR5.Louvain.mod #IGRAPH clustering multi level, groups: 7, mod: 0.58

#Plot the Louvain communities
plot(S6.nonIR5.Louvain.mod, S6.nonIR5,layout=Lyt.FR[-isolated,], alpha=0.6, 
     vertex.color=V(S6.nonIR5)$An.color, vertex.frame.color=V(S6.nonIR5)$An.color, vertex.shape="circle", 
     vertex.size=5, vertex.label=V(S6.nonIR5)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black", 
     edge.width=E(S6.nonIR5)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout", 
     rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)

#Generate HeatMaps of the adjacency matrixes
#Let's use the neighbor list for generating the Adjaency Matrix
SensorID_FNN #Let's use the edge list that contain the total amount of possible interactors for the three sensors
v<-levels(as.factor(unlist(SensorID_FNN[,1:2]))) #get vertex names from the network list
n<-length(v)#get the number of vertex in the network
e<-matrix(match(as.character(unlist(SensorID_FNN[,1:2])),v),ncol = 2)#Edgelist with row index instead of Names
w<-S6_nonIR5$weight #Edge weights, only from the S6_nonIR5 network to quit those edges that do not exist at S6_nonIR5

# distringuish between i,j and j,i we want to make M symmetric (undirected) by
# considering an undirected network.
S6.nonIR5.M <- matrix(0, n, n)                  # set up a co-expression matrix
S6.nonIR5.M[e] <- w                             # fill it in with edge weights
S6.nonIR5.M <- S6.nonIR5.M + t(S6.nonIR5.M)                         # make this symmetric
any(is.na(S6.nonIR5.M))#check no NAs in the matrix
dimnames(S6.nonIR5.M) <- list(v, v)             # label the vertices
# We let A_ij = 1 if there is an interaction. 
#A <- 1*(M > 0) #right now we keep the weights as they are
S6.nonIR5.M #The adjacency Matrix not as an igraph object, easier to manipulate, and to subset for smaller Sensor specific networks
quantile(S6.nonIR5.M,c(0,0.25,0.5,0.75,1))

#Creating the HeatMap
#Red Blue Colors
library(RColorBrewer)
breaks <- c(seq(-1, 1, length=10))
breaks
cols <- gray(1-c(0,seq(0,1,length=10))) # colors
cols2 <- hcl.colors(10 , palette = "Blues3" , alpha = 1 , rev = TRUE)
cols3 <- colorRamp("red" , space = "rgb" , interpolate = "linear")
HMCol<- brewer.pal(10,name = "RdBu")
HMCol<- rev(HMCol)
HMCol


# color the vertex for CCNB1 blue, its neighbors red, and the others white
SensorIDBind <- which(v == "SIRT6") # index for ccnb1
#SensorIDBind<-as.factor(SensorIDBind)
vcols <- rep("white",n)
vcols[SensorIDBind] <- "blue"
#Subsetting S6.nonIR5.M by doing a loop asking for look into every row and column the index in SensorIDBind worked
for (i in SensorIDBind) {
  vcols[which(S6.nonIR5.M[,i] > 0 | S6.nonIR5.M[i,])] <- "red"
}
##Color annotation for functional annotation
MatSet<-subset(Mat,Mat$NAME%in%v)# select only those proteins that are cointained in the node list of our matrixes
match(v,MatSet$NAME) # Check that the indexes of the Annotation Dataframe and the node list are in the same position 
match(rownames(S6.nonIR5.M),MatSet$NAME) #Check that the row names of the Adj Matrix and the Annotation dataframe are in the same order
match(colnames(S6.nonIR5.M),MatSet$NAME) #Check that the column names of the Adj Matrix and the Annotation dataframe are in the same order
AnCol<-MatSet$An.color#Vector of colors
AnCol<-as.character(AnCol)#transform it into a character vector, so the colors can be read

# now actually make the heat map
par(mar=rep(0,4)) # make plot margins 0
heatmap(S6.nonIR5.M, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=cols, breaks=breaks,frame=T, cexRow = 0.5 , cexCol = 0.5, main = "Sirt6 nonIR5 AdjMat", ylab="Total Protein Interactors",xlab="Total Protein Interactors")

#
legend("topleft", c("Cell Cycle","Chromatin Org","DNA Repair","Nuclear Pore","p53","Proteasome","Ribosome Biogenesis","RNA Metabolism","Transcription Reg","UPR"), 
       fill=c("coral1","darkviolet","firebrick1","deepskyblue3","tomato","grey","lightgreen","gold","slateblue1","maroon"),  
       bty="n", inset=0, xpd=T,  border=F, cex =2,y.intersp = 0.7,x.intersp = 0.2,pt.cex = 1)

############################################Sirt6 IR5 Networks
###
#Sirt6 IR5 Networks
#Let's keep the Mat as the Node attribute table and load the Neighbor list as the matrixes
Mat# Node attribute matrix for all the networks

#Sirt6ID IR5 nodes Network Edge Lists
S6_IR5<-read.delim("./Sirt6_Files/S6_IR5.tsv",sep = "\t",col.names = )
str(S6_IR5)
head(S6_IR5)

#Generate the igraph network object
S6.IR5<-graph_from_data_frame(d=S6_IR5,vertices = Mat, directed = F )#Generate the igraph object with the edge and node annotation attributes
str(vertex_attr(S6.IR5))#Check node attributes
str(edge_attr(S6.IR5))#Check for edge attributes
S6.IR5<-igraph::simplify(S6.IR5, remove.multiple=TRUE, remove.loops=TRUE) #Remove repetitive edges and self-connecting edges
S6.IR5
#Small plot
plot(S6.IR5)
#Generate an adjacency Matrix that will be used after
S6.IR5_AdjMat<-as_adjacency_matrix(S6.IR5,type="both",attr = "weight",names=TRUE,sparse = TRUE)
S6.IR5_AdjMat<-as.matrix(S6.IR5_AdjMat)
write.csv(S6.IR5_AdjMat,"./S6.IR5_AdjMat.csv")

#Generate a vertex size based on degree and count
vsize<-igraph::degree(S6.IR5)/vcount(S6.IR5)
#Identify and remove disconnected nodes
isolated<-which(igraph::degree(S6.IR5)==0)
isolated
length(isolated)#
#remove isolated nodes
S6.IR5<-igraph::delete_vertices(S6.IR5,isolated)

#Generate layouts for the graph
#Do not run again, to not to change the layouts that I already have
#Lyt.Nice<-layout_nicely(S6.IR5,dim = 3)
#Lyt.GOpt<-layout_with_graphopt(S6.IR5, niter=5000, mass = 30, charge = E(S6.IR5)$weight)
#Lyt.DrL<-layout_with_drl(S6.IR5, weights = E(S6.IR5)$weight)
#Lyt.FR<-layout.fruchterman.reingold(S6.IR5, niter=10000, weights=E(S6.IR5)$weight)#This one looks like the best option
#Lyt.KK<-layout_with_kk(S6.IR5, maxiter=10000, weights=E(S6.IR5)$weight) #very large separation of peripheral nodes
#Lyt.GEM<-layout_with_gem(S6.IR5, maxiter=10000) #Bad separation
#Lyt.DH<-layout_with_dh(S6.IR5, maxiter=150, fineiter=max(10, log2(833)))
#Lyt.LGL<-layout.lgl(S6.IR5, maxit=150,maxdelta=833,area=(length(E(S6.IR5)))^2,coolexp=1.5,repulserad=((length(E(S6.IR5)))^2)*(length(E(S6.IR5))))#Do not work


#Star and Tree layouts
#V(S6.IR5)[name=="SIRT6"]
which(V(S6.IR5)[name=="SIRT6"])
which(V(S6.IR5)$name=="SIRT6")
which(V(S6.IR5)$name=="MRE11A")
which(V(S6.IR5)$name=="XRCC5")
Lyt.Star<-layout_as_star(S6.IR5,center=V(S6.IR5)[name=="SIRT6"]) #it has a special order argument where we can assign an order to the vertex, hence maybe the OPLS-DA:VIP value order may be used as ordering vector
Lyt.Tree<-layout_as_tree(S6.IR5,root=c(191))
Lyt.RT<-layout.reingold.tilford(S6.IR5,root=c(191))# This one looks good, but the indexing do not allow to highlight the correct nodes
Lyt.sugi<-layout_with_sugiyama(S6.IR5, maxiter=100, weights=E(S6.IR52)$weight)#This one looks good for undirected tree layouts

#Check what layout looks fine

plot.igraph(S6.IR5,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(S6.IR5)$An.color, vertex.frame.color=V(S6.IR5)$An.color, vertex.shape="circle",
            vertex.size=V(S6.IR5)$S6.FC.IR5, vertex.label=V(S6.IR5)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(S6.IR5)$weight, edge.lty=1, main="Sirt6ID IR5 Network", sub="FR Tree Layout",
            rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)

#With ggraph

ggraph(S6.IR5,layout=Lyt.FR[-isolated,])+
  geom_edge_fan(color="gray65", width=1, alpha=0.4)+
  geom_node_point(aes(color=V(S6.IR5)$Annotation, size=V(S6.IR5)$S6.FC.IR5),alpha=0.75)+
  geom_node_text(aes(label=name),color="gray25",size=4)+
  scale_color_manual(limits=V(S6.IR5)$Annotation, values=V(S6.IR5)$An.color)+
  scale_size_continuous(range = c(5, 15)) +
  theme_void()+
  title("Sirt6 IR5 Network")+
  guides(color=guide_legend(override.aes = list(size=7)))+
  theme(legend.text=element_text(size=12,face = "bold"),
        legend.title = element_text(size=13,face = "bold"))


#Community Detection
clp <- igraph::cluster_optimal(S6.IR5)
class(clp)
V(S6.IR5)$community <- clp$membership
colrs <- adjustcolor( c("gray50", "tomato", "gold", "yellowgreen"), alpha=.6)
#Plot the communities detected
plot.igraph(clp,S6.IR5,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(S6.IR5)$An.color, vertex.frame.color=V(S6.IR5)$An.color, vertex.shape="circle",
            vertex.size==V(S6.IR5)$S6.FC.IR5, vertex.label=V(S6.IR5)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(S6.IR5)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)

#Distances from Sirt6, measure the distances between each node in the network to Sirt6
dist.from.Sirt6 <- distances(S6.IR5, v=V(S6.IR5)[name=="SIRT6"],
                             to=V(S6.IR5), weights=E(S6.IR5)$weight)
dist.from.Sirt6<-dist.from.Sirt6[,!is.infinite(colSums(dist.from.Sirt6))]
dist.from.Sirt6
range(dist.from.Sirt6)

#Set colors to plot the distances:
oranges <- colorRampPalette(c("dark red", "gold"))
col <- oranges(max(dist.from.Sirt6))
col <- col[dist.from.Sirt6+1]
plot(S6.IR5, layout=Lyt.FR[-isolated,],vertex.color=col, vertex.label=dist.from.Sirt6,
     vertex.label.color="white",vertex.shape="circle",
     vertex.size=V(S6.IR5)$S6.FC.IR5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(S6.IR5)$weight, edge.lty=1, main="Sirt6 IR5 Full Nuclear Network", sub="FR Layout",
     rescale=T)

#Shortest path from Sirt6 to Nucleolin
S6Ncl.path <- shortest_paths(S6.IR5,
                             from = V(S6.IR5)[name=="SIRT6"],
                             to = V(S6.IR5)[name=="NCL"],
                             output = "both") # both path nodes and edges
# Generate edge color variable to plot the path:
ecol <- rep("gray80", ecount(S6.IR5))
ecol[unlist(S6Ncl.path$epath)] <- "orange"
# Generate edge width variable to plot the path:
ew <- rep(2, ecount(S6.IR5))
ew[unlist(S6Ncl.path$epath)] <- 4
# Generate node color variable to plot the path:
vcol <- rep("gray40", vcount(S6.IR5))
vcol[unlist(S6Ncl.path$vpath)] <- "gold"
#Plot the path between Sirt6 and Nucleolin
plot(S6.IR5, layout=Lyt.FR[-isolated,], vertex.color=vcol, edge.color=ecol,
     edge.width=ew,vertex.label.color="white",vertex.shape="circle",
     vertex.size=5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(S6.IR5)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)

#Identify those edges incident to Sirt6 and its neighbor nodes
inc.edges <- incident(S6.IR5, V(S6.IR5)[name=="SIRT6"], mode="all")
# Set colors to plot the selected edges.
ecol <- rep("gray80", ecount(S6.IR5))
ecol[inc.edges] <- "orange"
vcol <- rep("gray40", vcount(S6.IR5))
vcol[V(S6.IR5)$name=="SIRT6"] <- "gold"
#Plot the incident edges and neighbor nodes to Sirt6
plot(S6.IR5, layout=Lyt.FR[-isolated,], vertex.color=vcol, edge.color=ecol, 
     edge.width=ew,vertex.label.color="white",vertex.shape="circle",
     vertex.size=V(S6.IR5)$S6.FC.IR5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(S6.IR5)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)

#Plot Degree Distribution
S6.IR5.Deg<-degree_distribution(S6.IR5, cumulative = T, mode="all")
plot( x=0:max(igraph::degree(S6.IR5)), y=1-S6.IR5.Deg, pch=19, cex=1.2, col="orange",
      xlab="Degree", ylab="Cumulative Frequency")

#Network based metrics
#Number of nodes and edges
vcount(S6.IR5) 
ecount(S6.IR5)
#Density (Number of edges in relation to the total number of edges in a full connected graph)
edge_density(S6.IR5,loops = F)
#Transitivity (Total number of conected triangles)
transitivity(S6.IR5,type="global") 
#modularity
mod<-V(S6.IR5)$Annotation #lets check for Annnotation based modularity
mod<-mod+1# 0s needs to be transformed to 1, so we add 1 to every cluster
modularity(S6.IR5,membership = mod)

#Assorativity
assortativity_nominal(S6.IR5,factor(V(S6.IR5)$Annotation),directed=F)

#Diameter (longest geodesic distance), shortest path between two most peripheral nodes
diameter(S6.IR5,directed = F, weights = E(S6.IR5)$weights)
diam<-get_diameter(S6.IR5,directed = F,weights = E(S6.IR5)$weights)
diam  
#Plot Diameter
vcol <- rep("gray50", vcount(S6.IR5))
vcol[diam] <- "gold"
ecol <- rep("gray90", ecount(S6.IR5))
ecol[E(S6.IR5, path=diam)] <- "orange"
# E(net, path=diam) finds edges along a path, here 'diam'
plot(S6.IR5, vertex.color=vcol, edge.color=ecol, vertex.label.color="black",edge.arrow.mode=0,layout=Lyt.FR[-isolated,],vertex.size=V(S6.IR5)$S6.FC.IR5)

#Node based metrics
degree<-igraph::centr_degree(S6.IR5,mode = "all",normalized = T)
degree<-degree$res
eigen_cent<-igraph::eigen_centrality(S6.IR5)
eigen_cent<-eigen_cent$vector
clos_cent <- igraph::closeness(S6.IR5, normalized=T)
betw_cent <- igraph::betweenness(S6.IR5, normalized=T)
hs <- hub_score(S6.IR5, weights=E(S6.IR5)$weight)$vector
as <- authority_score(S6.IR5, weights=E(S6.IR5)$weight)$vector
Trans<-transitivity(S6.IR5,type="local")
Trans
#write results
S6.IR5.Centrality<-data.frame(NAME=V(S6.IR5)$name, UniprotID=V(S6.IR5)$Uniprot, Degree=degree, Eigen= eigen_cent, Closeness= clos_cent, Betweeness= betw_cent, Hubs=hs, Authority=as,Transitivity=Trans)
S6.IR5.Centrality
write.csv(S6.IR5.Centrality,"./S6.IR5.Centrality.csv")

#Cocitation, Number of shared interactions among nodes
Cocitation<-cocitation(S6.IR5)
Cocitation#This matrix is similar to an adjacency matrix
write.csv(Cocitation,"./S6.IR5.cocitation.csv")

#Plot with degree as node size
ggraph(S6.IR5,layout=Lyt.FR[-isolated,])+
  geom_edge_fan(color="gray85", width=1, alpha=0.4)+
  geom_node_point(aes(color=V(S6.IR5)$Annotation, size=degree),alpha=0.75)+
  geom_node_text(aes(label=name),color="gray25",size=4)+
  scale_color_manual(limits=V(S6.IR5)$Annotation, values=V(S6.IR5)$An.color)+
  scale_size_continuous(range = c(5, 15)) +
  theme_void()+
  title("Sirt6 IR5 Network")+
  guides(color=guide_legend(override.aes = list(size=7)))+
  theme(legend.text=element_text(size=12,face = "bold"),
        legend.title = element_text(size=13,face = "bold"))


#Modularity Analysis
#Edge_Betweeness Modularity
S6.IR5.ebmod=edge.betweenness.community(S6.IR5,weights=E(S6.IR5)$weight,directed = F) 
S6.IR5.ebmod

# Plot the networks with the modules highlighted
plot(S6.IR5.ebmod, S6.IR5,layout=Lyt.FR[-isolated,], alpha=0.6, 
     vertex.color=V(S6.IR5)$An.color, vertex.frame.color=V(S6.IR5)$An.color, vertex.shape="circle", 
     vertex.size=5, vertex.label=V(S6.IR5)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black", 
     edge.width=E(S6.IR5)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout", 
     rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)

#Louvain modularity
S6.IR5.Louvain.mod<-cluster_louvain(S6.IR5,weights=E(S6.IR5)$weight)
S6.IR5.Louvain.mod

#Plot the Louvain communities
plot(S6.IR5.Louvain.mod, S6.IR5,layout=Lyt.FR[-isolated,], alpha=0.6, 
     vertex.color=V(S6.IR5)$An.color, vertex.frame.color=V(S6.IR5)$An.color, vertex.shape="circle", 
     vertex.size=5, vertex.label=V(S6.IR5)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black", 
     edge.width=E(S6.IR5)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout", 
     rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#Generate HeatMaps of the adjacency matrixes
#Let's use the neighbor list for generating the Adjaency Matrix
SensorID_FNN #Let's use the edge list that contain the total amount of possible interactors for the three sensors
v<-levels(as.factor(unlist(SensorID_FNN[,1:2]))) #get vertex names from the network list
n<-length(v)#get the number of vertex in the network
e<-matrix(match(as.character(unlist(SensorID_FNN[,1:2])),v),ncol = 2)#Edgelist with row index instead of Names
w<-S6_IR5$weight #Edge weights, only from the S6_IR5 network to quit those edges that do not exist at S6_IR5

# M is our co-expression network adjacency matrix.

# distringuish between i,j and j,i we want to make M symmetric (undirected) by
# considering an undirected network.
S6.IR5.M <- matrix(0, n, n)                  # set up a co-expression matrix
S6.IR5.M[e] <- w                             # fill it in with edge weights
S6.IR5.M <- S6.IR5.M + t(S6.IR5.M)                         # make this symmetric
any(is.na(S6.IR5.M))#check no NAs in the matrix
dimnames(S6.IR5.M) <- list(v, v)             # label the vertices
# We let A_ij = 1 if there is an interaction. 
#A <- 1*(M > 0) #right now we keep the weights as they are
S6.IR5.M #The adjacency Matrix not as an igraph object, easier to manipulate, and to subset for smaller Sensor specific networks
quantile(S6.IR5.M,c(0,0.25,0.5,0.75,1))

#Creating the HeatMap
# make a heatmap image with colors: white below 0,
# grey from 0.4 to 1 with darker closer to 1
breaks <- c(0, seq(0.4, 0.991, length=11)) # breaks in the color bins
cols <- grey(1-c(0,seq(0.5,1,length=10))) # colors

#Red Blue Colors
breaks <- c(seq(-1, 1, length=11))
breaks
cols <- gray(1-c(0,seq(0,1,length=10))) # colors
cols2 <- hcl.colors(10 , palette = "Blues3" , alpha = 1 , rev = TRUE)
cols3 <- colorRamp("red" , space = "rgb" , interpolate = "linear")
HMCol<- brewer.pal(10,name = "RdBu")
HMCol<- rev(HMCol)
HMCol


# color the vertex for CCNB1 blue, its neighbors red, and the others white
SensorIDBind <- which(v == "SIRT6") # index for ccnb1
vcols <- rep("white",n)
vcols[SensorIDBind] <- "blue"
#Subsetting S6.IR5.M by doing a loop asking for look into every row and column the index in SensorIDBind worked
for (i in SensorIDBind) {
  vcols[which(S6.IR5.M[,i] > 0 | S6.IR5.M[i,])] <- "red"
}
##Color annotation for functional annotation
MatSet<-subset(Mat,Mat$NAME%in%v)# select only those proteins that are cointained in the node list of our matrixes
match(v,MatSet$NAME) # Check that the indexes of the Annotation Dataframe and the node list are in the same position 
match(rownames(S6.IR5.M),MatSet$NAME) #Check that the row names of the Adj Matrix and the Annotation dataframe are in the same order
match(colnames(S6.IR5.M),MatSet$NAME) #Check that the column names of the Adj Matrix and the Annotation dataframe are in the same order
AnCol<-MatSet$An.color#Vector of colors
AnCol<-as.character(AnCol)#transform it into a character vector, so the colors can be read

# now actually make the heat map
par(mar=rep(0,4)) # make plot margins 0
heatmap(S6.IR5.M, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=cols, breaks=breaks,frame=T, cexRow = 0.5 , cexCol = 0.5, main = "Sirt6 IR5 AdjMat", ylab="Total Protein Interactors",xlab="Total Protein Interactors")
#
plot.new()
legend("topleft", c("Cell Cycle","Chromatin Org","DNA Repair","Nuclear Pore","p53","Proteasome","Ribosome Biogenesis","RNA Metabolism","Transcription Reg","UPR"), 
       fill=c("coral1","darkviolet","firebrick1","deepskyblue3","tomato","grey","lightgreen","gold","slateblue1","maroon"),  
       bty="n", inset=0, xpd=T,  border=F, cex =2,y.intersp = 0.7,x.intersp = 0.2,pt.cex = 1)


#Now we substract the S6nonIR5 AdjMatrix to the S6IR6 AdjMatrix
S6IR5.AdjMat.Change<-S6.IR5.M-S6.nonIR5.M
range(S6IR5.AdjMat.Change)#Let's check the magnitude of the substracted weights
#[1] -0.807  0.805
any(is.na(S6IR5.AdjMat.Change))
#Generate the breaks and color functions
breaks <- c(seq(-0.5, 0.5, length=10))
breaks
HMCol<- brewer.pal(9,name = "RdBu")
HMCol<- rev(HMCol)
HMCol
[1] "#2166AC" "#4393C3" "#92C5DE" "#D1E5F0" "#F7F7F7" "#FDDBC7" "#F4A582" "#D6604D" "#B2182B"

#Plot the HeatMap
par(mar=rep(0,4)) # make plot margins 0
heatmap(S6IR5.AdjMat.Change, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=HMCol, breaks=breaks,frame=T, cexRow = 0.5 , cexCol = 0.5, main = "Sirt6 IR5 AdjMat (minus nonIR5 edges)", ylab="Total Protein Interactors",xlab="Total Protein Interactors")
#
#check with the FNN
#now actually make the heat map
par(mar=rep(0,4)) # make plot margins 0
heatmap(FNN.M, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=0:1, frame=T, cexRow = 0.5 , cexCol = 0.5)
#
legend("topleft", c("Neighbors SIRT6","SIRT6"), fill=c("red","blue"),  
       bty="n", inset=0, xpd=T,  border=F)


#################### From this ponint onwards the scripts will follow the same pattern and order as the Sirt6ID nonIR5 and IR5 showed before, 
##########################if the analysis of a particular network is needed or is desired to consult, it will appear in the next sections of the script
##############################################s
##########################################Sirt6 IR30

#Sirt6 IR30 Networks


#Let's keep the Mat as the Node attribute table and load the Neighbor list as the matrixes
Mat# Node attribute matrix for all the networks

#Sirt6ID IR30 nodes Network Edge Lists
S6_IR30<-read.delim("./Sirt6_Files/S6_IR30.tsv",sep = "\t",col.names = )
str(S6_IR30)
head(S6_IR30)

#Generate the igraph network object
S6.IR30<-graph_from_data_frame(d=S6_IR30,vertices = Mat, directed = F )#Generate the igraph object with the edge and node annotation attributes
str(vertex_attr(S6.IR30))#Check node attributes
str(edge_attr(S6.IR30))#Check for edge attributes
S6.IR30<-igraph::simplify(S6.IR30, remove.multiple=TRUE, remove.loops=TRUE) #Remove repetitive edges and self-connecting edges
S6.IR30
#Small plot
plot(S6.IR30)

#Generate an adjacency Matrix that will be used after
S6.IR30_AdjMat<-as_adjacency_matrix(S6.IR30,type="both",attr = "weight",names=TRUE,sparse = TRUE)
S6.IR30_AdjMat<-as.matrix(S6.IR30_AdjMat)
write.csv(S6.IR30_AdjMat,"./S6.IR30_AdjMat.csv")

#Generate a vertex size based on degree and count
vsize<-igraph::degree(S6.IR30)/vcount(S6.IR30)
#Identify and remove disconnected nodes
isolated<-which(igraph::degree(S6.IR30)==0)
isolated
length(isolated)# [1] 605
#remove isolated nodes
S6.IR30<-igraph::delete_vertices(S6.IR30,isolated)


#Generate layouts for the graph
#Do not run again, to not to change the layouts that I already have
Lyt.Nice<-layout_nicely(S6.IR30,dim = 3)
Lyt.GOpt<-layout_with_graphopt(S6.IR30, niter=5000, mass = 30, charge = E(S6.IR30)$weight)
Lyt.DrL<-layout_with_drl(S6.IR30, weights = E(S6.IR30)$weight)
Lyt.FR<-layout.fruchterman.reingold(S6.IR30, niter=10000, weights=E(S6.IR30)$weight)#This one looks like the best option
Lyt.KK<-layout_with_kk(S6.IR30, maxiter=10000, weights=E(S6.IR30)$weight) #very large separation of peripheral nodes
Lyt.GEM<-layout_with_gem(S6.IR30, maxiter=10000) #Bad separation
#Lyt.DH<-layout_with_dh(S6.IR30, maxiter=150, fineiter=max(10, log2(833)))
Lyt.LGL<-layout.lgl(S6.IR30, maxit=150,maxdelta=833,area=(length(E(S6.IR30)))^2,coolexp=1.5,repulserad=((length(E(S6.IR30)))^2)*(length(E(S6.IR30))))#Do not work


#Star and Tree layouts

V(S6.IR30)[name=="SIRT6"]
which(V(S6.IR30)[name=="SIRT6"])# 228 #Check what is the rowindex for your desired node
which(V(S6.IR30)$name=="SIRT6")# 191 #Check what is the rowindex for your desired node
which(V(S6.IR30)$name=="MRE11A")# 327
which(V(S6.IR30)$name=="XRCC5")# 722
Lyt.Star<-layout_as_star(S6.IR30,center=V(S6.IR30)[name=="SIRT6"]) #it has a special order argument where we can assign an order to the vertex, hence maybe the OPLS-DA:VIP value order may be used as ordering vector

Lyt.Tree<-layout_as_tree(S6.IR30,root=c(187))
Lyt.RT<-layout.reingold.tilford(S6.IR30,root=c(187))# This one looks good, but the indexing do not allow to highlight the correct nodes
Lyt.sugi<-layout_with_sugiyama(S6.IR30, maxiter=100, weights=E(S6.IR302)$weight)#This one looks good for undirected tree layouts

#Check what layout looks fine

plot.igraph(S6.IR30,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(S6.IR30)$An.color, vertex.frame.color=V(S6.IR30)$An.color, vertex.shape="circle",
            vertex.size=V(S6.IR30)$S6.FC.IR30, vertex.label=V(S6.IR30)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(S6.IR30)$weight, edge.lty=1, main="Sirt6ID IR30 Network", sub="FR Tree Layout",
            rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#With ggraph

ggraph(S6.IR30,layout=Lyt.FR[-isolated,])+
  geom_edge_fan(color="gray65", width=1, alpha=0.4)+
  geom_node_point(aes(color=V(S6.IR30)$Annotation, size=V(S6.IR30)$S6.FC.IR30),alpha=0.75)+
  geom_node_text(aes(label=name),color="gray25",size=4)+
  scale_color_manual(limits=V(S6.IR30)$Annotation, values=V(S6.IR30)$An.color)+
  scale_size_continuous(range = c(5, 15)) +
  theme_void()+
  title("Sirt6 IR30 Network")+
  guides(color=guide_legend(override.aes = list(size=7)))+
  theme(legend.text=element_text(size=12,face = "bold"),
        legend.title = element_text(size=13,face = "bold"))



#Community Detection
#clp <- igraph::cluster_optimal(S6.IR30)
class(clp)
V(S6.IR30)$community <- clp$membership
colrs <- adjustcolor( c("gray50", "tomato", "gold", "yellowgreen"), alpha=.6)


plot.igraph(clp,S6.IR30,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(S6.IR30)$An.color, vertex.frame.color=V(S6.IR30)$An.color, vertex.shape="circle",
            vertex.size==V(S6.IR30)$S6.FC.IR30, vertex.label=V(S6.IR30)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(S6.IR30)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#Distances edges highlight
dist.from.Sirt6 <- distances(S6.IR30, v=V(S6.IR30)[name=="SIRT6"],
                             to=V(S6.IR30), weights=E(S6.IR30)$weight)
dist.from.Sirt6<-dist.from.Sirt6[,!is.infinite(colSums(dist.from.Sirt6))]

dist.from.Sirt6
range(dist.from.Sirt6)

#Set colors to plot the distances:
oranges <- colorRampPalette(c("dark red", "gold"))
col <- oranges(max(dist.from.Sirt6))
col <- col[dist.from.Sirt6+1]
plot(S6.IR30, layout=Lyt.FR[-isolated,],vertex.color=col, vertex.label=dist.from.Sirt6,
     vertex.label.color="white",vertex.shape="circle",
     vertex.size=V(S6.IR30)$S6.FC.IR30, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(S6.IR30)$weight, edge.lty=1, main="Sirt6 IR30 Full Nuclear Network", sub="FR Layout",
     rescale=T)


#Shortest path from Sirt6 to Nucleolin
S6Ncl.path <- shortest_paths(S6.IR30,
                             from = V(S6.IR30)[name=="SIRT6"],
                             to = V(S6.IR30)[name=="NCL"],
                             output = "both") # both path nodes and edges
# Generate edge color variable to plot the path:
ecol <- rep("gray80", ecount(S6.IR30))
ecol[unlist(S6Ncl.path$epath)] <- "orange"
# Generate edge width variable to plot the path:
ew <- rep(2, ecount(S6.IR30))
ew[unlist(S6Ncl.path$epath)] <- 4
# Generate node color variable to plot the path:
vcol <- rep("gray40", vcount(S6.IR30))
vcol[unlist(S6Ncl.path$vpath)] <- "gold"
plot(S6.IR30, layout=Lyt.FR[-isolated,], vertex.color=vcol, edge.color=ecol,
     edge.width=ew,vertex.label.color="white",vertex.shape="circle",
     vertex.size=5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(S6.IR30)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)



#Plot based on Sirt6 as incident
inc.edges <- incident(S6.IR30, V(S6.IR30)[name=="SIRT6"], mode="all")
# Set colors to plot the selected edges.
ecol <- rep("gray80", ecount(S6.IR30))
ecol[inc.edges] <- "orange"
vcol <- rep("gray40", vcount(S6.IR30))
vcol[V(S6.IR30)$name=="SIRT6"] <- "gold"

plot(S6.IR30, layout=Lyt.FR[-isolated,], vertex.color=vcol, edge.color=ecol, 
     edge.width=ew,vertex.label.color="white",vertex.shape="circle",
     vertex.size=V(S6.IR30)$S6.FC.IR30, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(S6.IR30)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)


#Plot Degree Distribution
S6.IR30.Deg<-degree_distribution(S6.IR30, cumulative = T, mode="all")
plot( x=0:max(igraph::degree(S6.IR30)), y=1-S6.IR30.Deg, pch=19, cex=1.2, col="orange",
      xlab="Degree", ylab="Cumulative Frequency")

#Network Descriptives
#Number of nodes and edges
vcount(S6.IR30) #[1] 220
ecount(S6.IR30)#[1]  1451
#Density (Number of edges in relation to the total number of edges in a full connected graph)
edge_density(S6.IR30,loops = F) #[1] 0.06023246
#Transitivity (Total number of conected triangles)
transitivity(S6.IR30,type="global") #[1] 0.5737691
#modularity
mod<-V(S6.IR30)$S6.FuzzyCMeans #lets check for fuzzy Cmeans modularity
mod<-mod+1#The non clustered interactors that appears as 0, needs to be transformed to 1, so we add 1 to every cluster
modularity(S6.IR30,membership = mod)
#S6.FuzzyCMeans modularity [1] 0.01282487
#Ku.FuzzyCMeans modularity [1] 0.02614348
#Mre.FuzzyCMeans modularity [1] 0.0196079

#Assorativity
assortativity_nominal(S6.IR30,factor(V(S6.IR30)$Mre.FuzzyCMeans),directed=F)
#Annotation [1] [1] 0.2118502
#S6.FuzzyCMeans [1] 0.0190546
#Ku.FuzzyCMeans [1] 0.04291517
#Mre.FuzzyCMeans [1] 0.02494721

#Diameter (longest geodesic distance), shortest path between two most peripheral nodes
diameter(S6.IR30,directed = F, weights = E(S6.IR30)$weights)#[1] 4.799
diam<-get_diameter(S6.IR30,directed = F,weights = E(S6.IR30)$weights)
diam #+ 8/220 vertices, named, from e95f181: [1] ACLY      CS        DRG1      RPS8      YBX1      HNRNPA2B1 PSMA3     PSMB3
#Plot Diameter
vcol <- rep("gray50", vcount(S6.IR30))
vcol[diam] <- "gold"
ecol <- rep("gray90", ecount(S6.IR30))
ecol[E(S6.IR30, path=diam)] <- "orange"
# E(net, path=diam) finds edges along a path, here 'diam'
plot(S6.IR30, vertex.color=vcol, edge.color=ecol, vertex.label.color="black",edge.arrow.mode=0,layout=Lyt.FR[-isolated,],vertex.size=V(S6.IR30)$S6.FC.IR30)

#Node topological descriptives
degree<-igraph::centr_degree(S6.IR30,mode = "all",normalized = T)
degree<-degree$res
eigen_cent<-igraph::eigen_centrality(S6.IR30)
eigen_cent<-eigen_cent$vector
clos_cent <- igraph::closeness(S6.IR30, normalized=T)
betw_cent <- igraph::betweenness(S6.IR30, normalized=T)
hs <- hub_score(S6.IR30, weights=E(S6.IR30)$weight)$vector
as <- authority_score(S6.IR30, weights=E(S6.IR30)$weight)$vector
Trans<-transitivity(S6.IR30,type="local")
Trans
#write results
S6.IR30.Centrality<-data.frame(NAME=V(S6.IR30)$name, UniprotID=V(S6.IR30)$Uniprot, Degree=degree, Eigen= eigen_cent, Closeness= clos_cent, Betweeness= betw_cent, Hubs=hs, Authority=as,Transitivity=Trans)
S6.IR30.Centrality
write.csv(S6.IR30.Centrality,"./S6.IR30.Centrality.csv")

#Cocitation, Number of shared interactions among nodes
Cocitation<-cocitation(S6.IR30)
Cocitation#This matrix is similar to an adjacency matrix
write.csv(Cocitation,"./S6.IR30.cocitation.csv")
#Plot with degree as node size
ggraph(S6.IR30,layout=Lyt.FR[-isolated,])+
  geom_edge_fan(color="gray85", width=1, alpha=0.4)+
  geom_node_point(aes(color=V(S6.IR30)$Annotation, size=degree),alpha=0.75)+
  geom_node_text(aes(label=name),color="gray25",size=4)+
  scale_color_manual(limits=V(S6.IR30)$Annotation, values=V(S6.IR30)$An.color)+
  scale_size_continuous(range = c(5, 15)) +
  theme_void()+
  title("Sirt6 IR30 Network")+
  guides(color=guide_legend(override.aes = list(size=7)))+
  theme(legend.text=element_text(size=12,face = "bold"),
        legend.title = element_text(size=13,face = "bold"))


#Modularity Analysis

#Edge_Betweeness Modularity
S6.IR30.ebmod=edge.betweenness.community(S6.IR30,weights=E(S6.IR30)$weight,directed = F) 
S6.IR30.ebmod# IGRAPH clustering edge betweenness, groups: 12, mod: 0.49

#
plot.igraph(S6.IR30.ebmod,S6.IR30,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(S6.IR30)$An.color, vertex.frame.color=V(S6.IR30)$An.color, vertex.shape="circle",
            vertex.size=5, vertex.label=V(S6.IR30)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(S6.IR30)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)

#
plot(S6.IR30.ebmod, S6.IR30,layout=Lyt.FR[-isolated,], alpha=0.6, 
     vertex.color=V(S6.IR30)$An.color, vertex.frame.color=V(S6.IR30)$An.color, vertex.shape="circle", 
     vertex.size=5, vertex.label=V(S6.IR30)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black", 
     edge.width=E(S6.IR30)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout", 
     rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)
#Dendogram of the Edge_Betweeness
dendPlot(S6.IR30.ebmod, mode="hclust")

#Louvain modularity
S6.IR30.Louvain.mod<-cluster_louvain(S6.IR30,weights=E(S6.IR30)$weight)
S6.IR30.Louvain.mod #IGRAPH clustering multi level, groups: 9, mod: 0.54

#Plot the Louvain communities
plot(S6.IR30.Louvain.mod, S6.IR30,layout=Lyt.FR[-isolated,], alpha=0.6, 
     vertex.color=V(S6.IR30)$An.color, vertex.frame.color=V(S6.IR30)$An.color, vertex.shape="circle", 
     vertex.size=5, vertex.label=V(S6.IR30)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black", 
     edge.width=E(S6.IR30)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout", 
     rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#Generate HeatMaps of the adjacency matrixes
#Let's use the neighbor list for generating the Adjaency Matrix
SensorID_FNN #Let's use the edge list that contain the total amount of possible interactors for the three sensors
v<-levels(as.factor(unlist(SensorID_FNN[,1:2]))) #get vertex names from the network list
n<-length(v)#get the number of vertex in the network
e<-matrix(match(as.character(unlist(SensorID_FNN[,1:2])),v),ncol = 2)#Edgelist with row index instead of Names
w<-S6_IR30$weight #Edge weights, only from the S6_IR30 network to quit those edges that do not exist at S6_IR30

# M is our co-expression network adjacency matrix.

# distringuish between i,j and j,i we want to make M symmetric (undirected) by
# considering an undirected network.
S6.IR30.M <- matrix(0, n, n)                  # set up a co-expression matrix
S6.IR30.M[e] <- w                             # fill it in with edge weights
S6.IR30.M <- S6.IR30.M + t(S6.IR30.M)                         # make this symmetric
any(is.na(S6.IR30.M))#check no NAs in the matrix
dimnames(S6.IR30.M) <- list(v, v)             # label the vertices
# We let A_ij = 1 if there is an interaction. 
#A <- 1*(M > 0) #right now we keep the weights as they are
S6.IR30.M #The adjacency Matrix not as an igraph object, easier to manipulate, and to subset for smaller Sensor specific networks
quantile(S6.IR30.M,c(0,0.25,0.5,0.75,1))
#
# 0%   25%   50%   75%  100% 
#0.000 0.000 0.000 0.000 0.991  

#Creating the HeatMap
# make a heatmap image with colors: white below 0,
# grey from 0.4 to 1 with darker closer to 1
breaks <- c(0, seq(0.4, 0.991, length=11)) # breaks in the color bins
cols <- grey(1-c(0,seq(0.5,1,length=10))) # colors

#Red Blue Colors
library(RColorBrewer)
breaks <- c(seq(-1, 1, length=11))
breaks
cols <- gray(1-c(0,seq(0,1,length=10))) # colors
cols2 <- hcl.colors(10 , palette = "Blues3" , alpha = 1 , rev = TRUE)
cols3 <- colorRamp("red" , space = "rgb" , interpolate = "linear")
HMCol<- brewer.pal(10,name = "RdBu")
HMCol<- rev(HMCol)
HMCol


# color the vertex for CCNB1 blue, its neighbors red, and the others white
SensorIDBind <- which(v == "SIRT6") # index for ccnb1
#SensorIDBind<-as.factor(SensorIDBind)
vcols <- rep("white",n)
vcols[SensorIDBind] <- "blue"
#Subsetting S6.IR30.M by doing a loop asking for look into every row and column the index in SensorIDBind worked
for (i in SensorIDBind) {
  vcols[which(S6.IR30.M[,i] > 0 | S6.IR30.M[i,])] <- "red"
}
##Color annotation for functional annotation
MatSet<-subset(Mat,Mat$NAME%in%v)# select only those proteins that are cointained in the node list of our matrixes
match(v,MatSet$NAME) # Check that the indexes of the Annotation Dataframe and the node list are in the same position 
match(rownames(S6.IR30.M),MatSet$NAME) #Check that the row names of the Adj Matrix and the Annotation dataframe are in the same order
match(colnames(S6.IR30.M),MatSet$NAME) #Check that the column names of the Adj Matrix and the Annotation dataframe are in the same order
AnCol<-MatSet$An.color#Vector of colors
AnCol<-as.character(AnCol)#transform it into a character vector, so the colors can be read

#Subsetting S6.IR30.M using the index vector in SensorIDBind do not work
#vcols[which(S6.IR30.M[,SensorIDBind] > 0 | S6.IR30.M[SensorIDBind,])] <- "red"

#Generate a dendogram for keeping the proteins in the sample place
#S6.IR30.dist<-dist(S6.IR30.M,method = "euclidian")
#S6.IR30.Dend<-hclust(S6.IR30.dist, method = "ward.D" )
#plot(S6.IR30.Dend)

#Cut the dendogram
#cut.height <- 10 # try numbers from 4 to 7
#FNN.hclusts<- cutree(FNN.Dend, cut.height)
#plot(FNN.hclusts)

# now actually make the heat map
par(mar=rep(0,4)) # make plot margins 0
heatmap(S6.IR30.M, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=cols, breaks=breaks,frame=T, cexRow = 0.5 , cexCol = 0.5, main = "Sirt6 IR30 AdjMat", ylab="Total Protein Interactors",xlab="Total Protein Interactors")
#
plot.new()
legend("topleft", c("Cell Cycle","Chromatin Org","DNA Repair","Nuclear Pore","p53","Proteasome","Ribosome Biogenesis","RNA Metabolism","Transcription Reg","UPR"), 
       fill=c("coral1","darkviolet","firebrick1","deepskyblue3","tomato","grey","lightgreen","gold","slateblue1","maroon"),  
       bty="n", inset=0, xpd=T,  border=F, cex =2,y.intersp = 0.7,x.intersp = 0.2,pt.cex = 1)


#Now we substract the S6nonIR30 AdjMatrix to the S6IR6 AdjMatrix
S6IR30.AdjMat.Change<-S6.IR30.M-S6.nonIR5.M
range(S6IR30.AdjMat.Change)#Let's check the magnitude of the substracted weights
#[1]-0.813  0.668
any(is.na(S6IR30.AdjMat.Change))
#Generate the breaks and color functions
breaks <- c(seq(-0.5, 0.5, length=10))
breaks
HMCol<- brewer.pal(9,name = "RdBu")
HMCol<- rev(HMCol)
HMCol
[1] "#2166AC" "#4393C3" "#92C5DE" "#D1E5F0" "#F7F7F7" "#FDDBC7" "#F4A582" "#D6604D" "#B2182B"

#Plot the HeatMap
par(mar=rep(0,4)) # make plot margins 0
heatmap(S6IR30.AdjMat.Change, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=HMCol, breaks=breaks,frame=T, cexRow = 0.5 , cexCol = 0.5, main = "Sirt6 IR30 AdjMat (minus nonIR30 edges)", ylab="Total Protein Interactors",xlab="Total Protein Interactors")


#Now we substract the S6IR5 AdjMatrix to the S6IR30 AdjMatrix
S6IR30.AdjMat.Change2<-S6.IR30.M-S6.IR5.M
range(S6IR30.AdjMat.Change)#Let's check the magnitude of the substracted weights
#[1]-0.813  0.668
any(is.na(S6IR30.AdjMat.Change))
#Generate the breaks and color functions
breaks <- c(seq(-0.5, 0.5, length=10))
breaks
HMCol<- brewer.pal(9,name = "RdBu")
HMCol<- rev(HMCol)
HMCol
[1] "#2166AC" "#4393C3" "#92C5DE" "#D1E5F0" "#F7F7F7" "#FDDBC7" "#F4A582" "#D6604D" "#B2182B"

#Plot the HeatMap
par(mar=rep(0,4)) # make plot margins 0
heatmap(S6IR30.AdjMat.Change2, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=HMCol, breaks=breaks,frame=T, cexRow = 0.5 , cexCol = 0.5, main = "Sirt6 IR30 AdjMat (minus IR5 edges)", ylab="Total Protein Interactors",xlab="Total Protein Interactors")

####################### Sirt6 IR2h Networks
###################################################

#Sirt6 IR2 Networks


#Let's keep the Mat as the Node attribute table and load the Neighbor list as the matrixes
Mat# Node attribute matrix for all the networks

#Sirt6ID IR2 nodes Network Edge Lists
S6_IR2<-read.delim("./Sirt6_Files/S6_IR2h.tsv",sep = "\t",col.names = )
str(S6_IR2)
head(S6_IR2)

#Generate the igraph network object
S6.IR2<-graph_from_data_frame(d=S6_IR2,vertices = Mat, directed = F )#Generate the igraph object with the edge and node annotation attributes
str(vertex_attr(S6.IR2))#Check node attributes
str(edge_attr(S6.IR2))#Check for edge attributes
S6.IR2<-igraph::simplify(S6.IR2, remove.multiple=TRUE, remove.loops=TRUE) #Remove repetitive edges and self-connecting edges
S6.IR2
#Small plot
plot(S6.IR2)

#Generate an adjacency Matrix that will be used after
S6.IR2_AdjMat<-as_adjacency_matrix(S6.IR2,type="both",attr = "weight",names=TRUE,sparse = TRUE)
S6.IR2_AdjMat<-as.matrix(S6.IR2_AdjMat)
write.csv(S6.IR2_AdjMat,"./S6.IR2_AdjMat.csv")

#Generate a vertex size based on degree and count
vsize<-igraph::degree(S6.IR2)/vcount(S6.IR2)
#Identify and remove disconnected nodes
isolated<-which(igraph::degree(S6.IR2)==0)
isolated
length(isolated)# [1] 535
#remove isolated nodes
S6.IR2<-igraph::delete_vertices(S6.IR2,isolated)


#Generate layouts for the graph
#Do not run again, to not to change the layouts that I already have
Lyt.Nice<-layout_nicely(S6.IR2,dim = 3)
Lyt.GOpt<-layout_with_graphopt(S6.IR2, niter=5000, mass = 30, charge = E(S6.IR2)$weight)
Lyt.DrL<-layout_with_drl(S6.IR2, weights = E(S6.IR2)$weight)
Lyt.FR<-layout.fruchterman.reingold(S6.IR2, niter=10000, weights=E(S6.IR2)$weight)#This one looks like the best option
Lyt.KK<-layout_with_kk(S6.IR2, maxiter=10000, weights=E(S6.IR2)$weight) #very large separation of peripheral nodes
Lyt.GEM<-layout_with_gem(S6.IR2, maxiter=10000) #Bad separation
#Lyt.DH<-layout_with_dh(S6.IR2, maxiter=150, fineiter=max(10, log2(833)))
Lyt.LGL<-layout.lgl(S6.IR2, maxit=150,maxdelta=833,area=(length(E(S6.IR2)))^2,coolexp=1.5,repulserad=((length(E(S6.IR2)))^2)*(length(E(S6.IR2))))#Do not work


#Star and Tree layouts

V(S6.IR2)[name=="SIRT6"]
which(V(S6.IR2)[name=="SIRT6"])# 228 #Check what is the rowindex for your desired node
which(V(S6.IR2)$name=="SIRT6")# 191 #Check what is the rowindex for your desired node
which(V(S6.IR2)$name=="MRE11A")# 327
which(V(S6.IR2)$name=="XRCC5")# 722
Lyt.Star<-layout_as_star(S6.IR2,center=V(S6.IR2)[name=="SIRT6"]) #it has a special order argument where we can assign an order to the vertex, hence maybe the OPLS-DA:VIP value order may be used as ordering vector

Lyt.Tree<-layout_as_tree(S6.IR2,root=c(251))
Lyt.RT<-layout.reingold.tilford(S6.IR2,root=c(251))# This one looks good, but the indexing do not allow to highlight the correct nodes
Lyt.sugi<-layout_with_sugiyama(S6.IR2, maxiter=100, weights=E(S6.IR22)$weight)#This one looks good for undirected tree layouts

#Check what layout looks fine

plot.igraph(S6.IR2,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(S6.IR2)$An.color, vertex.frame.color=V(S6.IR2)$An.color, vertex.shape="circle",
            vertex.size=V(S6.IR2)$S6.FC.IR2, vertex.label=V(S6.IR2)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(S6.IR2)$weight, edge.lty=1, main="Sirt6ID IR2 Network", sub="FR Tree Layout",
            rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#With ggraph

ggraph(S6.IR2,layout=Lyt.FR[-isolated,])+
  geom_edge_fan(color="gray65", width=1, alpha=0.4)+
  geom_node_point(aes(color=V(S6.IR2)$Annotation, size=V(S6.IR2)$S6.FC.IR2),alpha=0.75)+
  geom_node_text(aes(label=name),color="gray25",size=4)+
  scale_color_manual(limits=V(S6.IR2)$Annotation, values=V(S6.IR2)$An.color)+
  scale_size_continuous(range = c(5, 15)) +
  theme_void()+
  title("Sirt6 IR2 Network")+
  guides(color=guide_legend(override.aes = list(size=7)))+
  theme(legend.text=element_text(size=12,face = "bold"),
        legend.title = element_text(size=13,face = "bold"))



#Community Detection
#clp <- igraph::cluster_optimal(S6.IR2)
class(clp)
V(S6.IR2)$community <- clp$membership
colrs <- adjustcolor( c("gray50", "tomato", "gold", "yellowgreen"), alpha=.6)


plot.igraph(clp,S6.IR2,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(S6.IR2)$An.color, vertex.frame.color=V(S6.IR2)$An.color, vertex.shape="circle",
            vertex.size==V(S6.IR2)$S6.FC.IR2, vertex.label=V(S6.IR2)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(S6.IR2)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#Distances edges highlight
dist.from.Sirt6 <- distances(S6.IR2, v=V(S6.IR2)[name=="SIRT6"],
                             to=V(S6.IR2), weights=E(S6.IR2)$weight)
dist.from.Sirt6<-dist.from.Sirt6[,!is.infinite(colSums(dist.from.Sirt6))]

dist.from.Sirt6
range(dist.from.Sirt6)

#Set colors to plot the distances:
oranges <- colorRampPalette(c("dark red", "gold"))
col <- oranges(max(dist.from.Sirt6))
col <- col[dist.from.Sirt6+1]
plot(S6.IR2, layout=Lyt.FR[-isolated,],vertex.color=col, vertex.label=dist.from.Sirt6,
     vertex.label.color="white",vertex.shape="circle",
     vertex.size=V(S6.IR2)$S6.FC.IR2, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(S6.IR2)$weight, edge.lty=1, main="Sirt6 IR2 Full Nuclear Network", sub="FR Layout",
     rescale=T)


#Shortest path from Sirt6 to Nucleolin
S6Ncl.path <- shortest_paths(S6.IR2,
                             from = V(S6.IR2)[name=="SIRT6"],
                             to = V(S6.IR2)[name=="NCL"],
                             output = "both") # both path nodes and edges

# Generate edge color variable to plot the path:
ecol <- rep("gray80", ecount(S6.IR2))
ecol[unlist(S6Ncl.path$epath)] <- "orange"
# Generate edge width variable to plot the path:
ew <- rep(2, ecount(S6.IR2))
ew[unlist(S6Ncl.path$epath)] <- 4
# Generate node color variable to plot the path:
vcol <- rep("gray40", vcount(S6.IR2))
vcol[unlist(S6Ncl.path$vpath)] <- "gold"
plot(S6.IR2, layout=Lyt.FR[-isolated,], vertex.color=vcol, edge.color=ecol,
     edge.width=ew,vertex.label.color="white",vertex.shape="circle",
     vertex.size=5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(S6.IR2)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)



#Plot based on Sirt6 as incident
inc.edges <- incident(S6.IR2, V(S6.IR2)[name=="SIRT6"], mode="all")
# Set colors to plot the selected edges.
ecol <- rep("gray80", ecount(S6.IR2))
ecol[inc.edges] <- "orange"
vcol <- rep("gray40", vcount(S6.IR2))
vcol[V(S6.IR2)$name=="SIRT6"] <- "gold"

plot(S6.IR2, layout=Lyt.FR[-isolated,], vertex.color=vcol, edge.color=ecol, 
     edge.width=ew,vertex.label.color="white",vertex.shape="circle",
     vertex.size=V(S6.IR2)$S6.FC.IR2, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(S6.IR2)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)


#Plot Degree Distribution
S6.IR2.Deg<-degree_distribution(S6.IR2, cumulative = T, mode="all")
plot( x=0:max(igraph::degree(S6.IR2)), y=1-S6.IR2.Deg, pch=19, cex=1.2, col="orange",
      xlab="Degree", ylab="Cumulative Frequency")

#Network Descriptives
#Number of nodes and edges
vcount(S6.IR2) #[1] 298
ecount(S6.IR2)#[1]  2422
#Density (Number of edges in relation to the total number of edges in a full connected graph)
edge_density(S6.IR2,loops = F) #[1] 0.05473075
#Transitivity (Total number of conected triangles)
transitivity(S6.IR2,type="global") #[1] 0.6488348
#modularity
mod<-V(S6.IR2)$S6.FuzzyCMeans #lets check for fuzzy Cmeans modularity
mod<-mod+1#The non clustered interactors that appears as 0, needs to be transformed to 1, so we add 1 to every cluster
modularity(S6.IR2,membership = mod)
#S6.FuzzyCMeans modularity [1] -0.001870669
#Ku.FuzzyCMeans modularity [1] 0.02570463
#Mre.FuzzyCMeans modularity [1] 0.009508217

#Assorativity
assortativity_nominal(S6.IR2,factor(V(S6.IR2)$Annotation),directed=F)
#Annotation [1] [1] 0.2205071
#S6.FuzzyCMeans [1] -0.002730246
#Ku.FuzzyCMeans [1] 0.0406572
#Mre.FuzzyCMeans [1] 0.01170151

#Diameter (longest geodesic distance), shortest path between two most peripheral nodes
diameter(S6.IR2,directed = F, weights = E(S6.IR2)$weights)#[1] 4.799
diam<-get_diameter(S6.IR2,directed = F,weights = E(S6.IR2)$weights)
diam #+ 9/298 vertices, named, from 79e08fc: [1] NOL9      WDR18     RPL5      EFTUD2    MRE11A    HIST2H2BE DDB1      CCT5      CCT4 
#Plot Diameter
vcol <- rep("gray50", vcount(S6.IR2))
vcol[diam] <- "gold"
ecol <- rep("gray90", ecount(S6.IR2))
ecol[E(S6.IR2, path=diam)] <- "orange"
# E(net, path=diam) finds edges along a path, here 'diam'
plot(S6.IR2, vertex.color=vcol, edge.color=ecol, vertex.label.color="black",edge.arrow.mode=0,layout=Lyt.FR[-isolated,],vertex.size=V(S6.IR2)$S6.FC.IR2)

#Node topological descriptives
degree<-igraph::centr_degree(S6.IR2,mode = "all",normalized = T)
degree<-degree$res
eigen_cent<-igraph::eigen_centrality(S6.IR2)
eigen_cent<-eigen_cent$vector
clos_cent <- igraph::closeness(S6.IR2, normalized=T)
betw_cent <- igraph::betweenness(S6.IR2, normalized=T)
hs <- hub_score(S6.IR2, weights=E(S6.IR2)$weight)$vector
as <- authority_score(S6.IR2, weights=E(S6.IR2)$weight)$vector
Trans<-transitivity(S6.IR2,type="local")
Trans
#write results
S6.IR2.Centrality<-data.frame(NAME=V(S6.IR2)$name, UniprotID=V(S6.IR2)$Uniprot, Degree=degree, Eigen= eigen_cent, Closeness= clos_cent, Betweeness= betw_cent, Hubs=hs, Authority=as,Transitivity=Trans)
S6.IR2.Centrality
write.csv(S6.IR2.Centrality,"./S6.IR2.Centrality.csv")

#Cocitation, Number of shared interactions among nodes
Cocitation<-cocitation(S6.IR2)
Cocitation#This matrix is similar to an adjacency matrix
write.csv(Cocitation,"./S6.IR2.cocitation.csv")
#Plot with degree as node size
ggraph(S6.IR2,layout=Lyt.FR[-isolated,])+
  geom_edge_fan(color="gray85", width=1, alpha=0.4)+
  geom_node_point(aes(color=V(S6.IR2)$Annotation, size=degree),alpha=0.75)+
  geom_node_text(aes(label=name),color="gray25",size=4)+
  scale_color_manual(limits=V(S6.IR2)$Annotation, values=V(S6.IR2)$An.color)+
  scale_size_continuous(range = c(5, 15)) +
  theme_void()+
  title("Sirt6 IR2 Network")+
  guides(color=guide_legend(override.aes = list(size=7)))+
  theme(legend.text=element_text(size=12,face = "bold"),
        legend.title = element_text(size=13,face = "bold"))


#Modularity Analysis

#Edge_Betweeness Modularity
S6.IR2.ebmod=edge.betweenness.community(S6.IR2,weights=E(S6.IR2)$weight,directed = F) 
S6.IR2.ebmod# IGRAPH clustering edge betweenness, groups: 29, mod: 0.46

#
plot.igraph(S6.IR2.ebmod,S6.IR2,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(S6.IR2)$An.color, vertex.frame.color=V(S6.IR2)$An.color, vertex.shape="circle",
            vertex.size=5, vertex.label=V(S6.IR2)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(S6.IR2)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)

#
plot(S6.IR2.ebmod, S6.IR2,layout=Lyt.FR[-isolated,], alpha=0.6, 
     vertex.color=V(S6.IR2)$An.color, vertex.frame.color=V(S6.IR2)$An.color, vertex.shape="circle", 
     vertex.size=5, vertex.label=V(S6.IR2)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black", 
     edge.width=E(S6.IR2)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout", 
     rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)
#Dendogram of the Edge_Betweeness
dendPlot(S6.IR2.ebmod, mode="hclust")

#Louvain modularity
S6.IR2.Louvain.mod<-cluster_louvain(S6.IR2,weights=E(S6.IR2)$weight)
S6.IR2.Louvain.mod #IGRAPH clustering multi level, groups: 11, mod: 0.52

#Plot the Louvain communities
plot(S6.IR2.Louvain.mod, S6.IR2,layout=Lyt.FR[-isolated,], alpha=0.6, 
     vertex.color=V(S6.IR2)$An.color, vertex.frame.color=V(S6.IR2)$An.color, vertex.shape="circle", 
     vertex.size=5, vertex.label=V(S6.IR2)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black", 
     edge.width=E(S6.IR2)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout", 
     rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#Generate HeatMaps of the adjacency matrixes
#Let's use the neighbor list for generating the Adjaency Matrix
SensorID_FNN #Let's use the edge list that contain the total amount of possible interactors for the three sensors
v<-levels(as.factor(unlist(SensorID_FNN[,1:2]))) #get vertex names from the network list
n<-length(v)#get the number of vertex in the network
e<-matrix(match(as.character(unlist(SensorID_FNN[,1:2])),v),ncol = 2)#Edgelist with row index instead of Names
w<-S6_IR2$weight #Edge weights, only from the S6_IR2 network to quit those edges that do not exist at S6_IR2

# M is our co-expression network adjacency matrix.

# distringuish between i,j and j,i we want to make M symmetric (undirected) by
# considering an undirected network.
S6.IR2.M <- matrix(0, n, n)                  # set up a co-expression matrix
S6.IR2.M[e] <- w                             # fill it in with edge weights
S6.IR2.M <- S6.IR2.M + t(S6.IR2.M)                         # make this symmetric
any(is.na(S6.IR2.M))#check no NAs in the matrix
dimnames(S6.IR2.M) <- list(v, v)             # label the vertices
# We let A_ij = 1 if there is an interaction. 
#A <- 1*(M > 0) #right now we keep the weights as they are
S6.IR2.M #The adjacency Matrix not as an igraph object, easier to manipulate, and to subset for smaller Sensor specific networks
quantile(S6.IR2.M,c(0,0.25,0.5,0.75,1))
#
#    0%   25%   50%   75%  100% 
#  0.000 0.000 0.000 0.000 0.991  

#Creating the HeatMap
# make a heatmap image with colors: white below 0,
# grey from 0.4 to 1 with darker closer to 1
breaks <- c(0, seq(0.4, 0.991, length=11)) # breaks in the color bins
cols <- grey(1-c(0,seq(0.5,1,length=10))) # colors

#Red Blue Colors
library(RColorBrewer)
breaks <- c(seq(-1, 1, length=11))
breaks
cols <- gray(1-c(0,seq(0,1,length=10))) # colors
cols2 <- hcl.colors(10 , palette = "Blues3" , alpha = 1 , rev = TRUE)
cols3 <- colorRamp("red" , space = "rgb" , interpolate = "linear")
HMCol<- brewer.pal(10,name = "RdBu")
HMCol<- rev(HMCol)
HMCol


# color the vertex for CCNB1 blue, its neighbors red, and the others white
SensorIDBind <- which(v == "SIRT6") # index for ccnb1
#SensorIDBind<-as.factor(SensorIDBind)
vcols <- rep("white",n)
vcols[SensorIDBind] <- "blue"
#Subsetting S6.IR2.M by doing a loop asking for look into every row and column the index in SensorIDBind worked
for (i in SensorIDBind) {
  vcols[which(S6.IR2.M[,i] > 0 | S6.IR2.M[i,])] <- "red"
}
##Color annotation for functional annotation
MatSet<-subset(Mat,Mat$NAME%in%v)# select only those proteins that are cointained in the node list of our matrixes
match(v,MatSet$NAME) # Check that the indexes of the Annotation Dataframe and the node list are in the same position 
match(rownames(S6.IR2.M),MatSet$NAME) #Check that the row names of the Adj Matrix and the Annotation dataframe are in the same order
match(colnames(S6.IR2.M),MatSet$NAME) #Check that the column names of the Adj Matrix and the Annotation dataframe are in the same order
AnCol<-MatSet$An.color#Vector of colors
AnCol<-as.character(AnCol)#transform it into a character vector, so the colors can be read

#Subsetting S6.IR2.M using the index vector in SensorIDBind do not work
#vcols[which(S6.IR2.M[,SensorIDBind] > 0 | S6.IR2.M[SensorIDBind,])] <- "red"

#Generate a dendogram for keeping the proteins in the sample place
#S6.IR2.dist<-dist(S6.IR2.M,method = "euclidian")
#S6.IR2.Dend<-hclust(S6.IR2.dist, method = "ward.D" )
#plot(S6.IR2.Dend)

#Cut the dendogram
#cut.height <- 10 # try numbers from 4 to 7
#FNN.hclusts<- cutree(FNN.Dend, cut.height)
#plot(FNN.hclusts)

# now actually make the heat map
par(mar=rep(0,4)) # make plot margins 0
heatmap(S6.IR2.M, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=cols, breaks=breaks,frame=T, cexRow = 0.5 , cexCol = 0.5, main = "Sirt6 IR2 AdjMat", ylab="Total Protein Interactors",xlab="Total Protein Interactors")
#
plot.new()
legend("topleft", c("Cell Cycle","Chromatin Org","DNA Repair","Nuclear Pore","p53","Proteasome","Ribosome Biogenesis","RNA Metabolism","Transcription Reg","UPR"), 
       fill=c("coral1","darkviolet","firebrick1","deepskyblue3","tomato","grey","lightgreen","gold","slateblue1","maroon"),  
       bty="n", inset=0, xpd=T,  border=F, cex =2,y.intersp = 0.7,x.intersp = 0.2,pt.cex = 1)


#Now we substract the S6nonIR2 AdjMatrix to the S6IR6 AdjMatrix
S6IR2.AdjMat.Change<-S6.IR2.M-S6.nonIR5.M
range(S6IR2.AdjMat.Change)#Let's check the magnitude of the substracted weights
#[1]-0.813  0.668
any(is.na(S6IR2.AdjMat.Change))
#Generate the breaks and color functions
breaks <- c(seq(-0.5, 0.5, length=10))
breaks
HMCol<- brewer.pal(9,name = "RdBu")
HMCol<- rev(HMCol)
HMCol
[1] "#2166AC" "#4393C3" "#92C5DE" "#D1E5F0" "#F7F7F7" "#FDDBC7" "#F4A582" "#D6604D" "#B2182B"

#Plot the HeatMap
par(mar=rep(0,4)) # make plot margins 0
heatmap(S6IR2.AdjMat.Change, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=HMCol, breaks=breaks,frame=T, cexRow = 0.5 , cexCol = 0.5, main = "Sirt6 IR2 AdjMat (minus nonIR2 edges)", ylab="Total Protein Interactors",xlab="Total Protein Interactors")


#Now we substract the S6IR5 AdjMatrix to the S6IR2 AdjMatrix
S6IR2.AdjMat.Change2<-S6.IR2.M-S6.IR30.M
range(S6IR2.AdjMat.Change)#Let's check the magnitude of the substracted weights
#[1]-0.813  0.668
any(is.na(S6IR2.AdjMat.Change))
#Generate the breaks and color functions
breaks <- c(seq(-0.5, 0.5, length=10))
breaks
HMCol<- brewer.pal(9,name = "RdBu")
HMCol<- rev(HMCol)
HMCol
[1] "#2166AC" "#4393C3" "#92C5DE" "#D1E5F0" "#F7F7F7" "#FDDBC7" "#F4A582" "#D6604D" "#B2182B"

#Plot the HeatMap
par(mar=rep(0,4)) # make plot margins 0
heatmap(S6IR2.AdjMat.Change2, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=HMCol, breaks=breaks,frame=T, cexRow = 0.5 , cexCol = 0.5, main = "Sirt6 IR2 AdjMat (minus IR30 edges)", ylab="Total Protein Interactors",xlab="Total Protein Interactors")

#
#check with the FNN
# now actually make the heat map
par(mar=rep(0,4)) # make plot margins 0
heatmap(FNN.M, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=0:1, frame=T, cexRow = 0.5 , cexCol = 0.5)
#
legend("topleft", c("Neighbors SIRT6","SIRT6"), fill=c("red","blue"),  
       bty="n", inset=0, xpd=T,  border=F)


####################### Sirt6 IR8h Networks
###################################################

#Sirt6 IR8 Networks


#Let's keep the Mat as the Node attribute table and load the Neighbor list as the matrixes
Mat# Node attribute matrix for all the networks

#Sirt6ID IR8 nodes Network Edge Lists
S6_IR8<-read.delim("./Sirt6_Files/S6_IR8h.tsv",sep = "\t",col.names = )
str(S6_IR8)
head(S6_IR8)

#Generate the igraph network object
S6.IR8<-graph_from_data_frame(d=S6_IR8,vertices = Mat, directed = F )#Generate the igraph object with the edge and node annotation attributes
str(vertex_attr(S6.IR8))#Check node attributes
str(edge_attr(S6.IR8))#Check for edge attributes
S6.IR8<-igraph::simplify(S6.IR8, remove.multiple=TRUE, remove.loops=TRUE) #Remove repetitive edges and self-connecting edges
S6.IR8
#Small plot
plot(S6.IR8)

#Generate an adjacency Matrix that will be used after
S6.IR8_AdjMat<-as_adjacency_matrix(S6.IR8,type="both",attr = "weight",names=TRUE,sparse = TRUE)
S6.IR8_AdjMat<-as.matrix(S6.IR8_AdjMat)
write.csv(S6.IR8_AdjMat,"./S6.IR8_AdjMat.csv")

#Generate a vertex size based on degree and count
vsize<-igraph::degree(S6.IR8)/vcount(S6.IR8)
#Identify and remove disconnected nodes
isolated<-which(igraph::degree(S6.IR8)==0)
isolated
length(isolated)# [1] 535
#remove isolated nodes
S6.IR8<-igraph::delete_vertices(S6.IR8,isolated)


#Generate layouts for the graph
#Do not run again, to not to change the layouts that I already have
Lyt.Nice<-layout_nicely(S6.IR8,dim = 3)
Lyt.GOpt<-layout_with_graphopt(S6.IR8, niter=5000, mass = 30, charge = E(S6.IR8)$weight)
Lyt.DrL<-layout_with_drl(S6.IR8, weights = E(S6.IR8)$weight)
Lyt.FR<-layout.fruchterman.reingold(S6.IR8, niter=10000, weights=E(S6.IR8)$weight)#This one looks like the best option
Lyt.KK<-layout_with_kk(S6.IR8, maxiter=10000, weights=E(S6.IR8)$weight) #very large separation of peripheral nodes
Lyt.GEM<-layout_with_gem(S6.IR8, maxiter=10000) #Bad separation
#Lyt.DH<-layout_with_dh(S6.IR8, maxiter=150, fineiter=max(10, log2(833)))
Lyt.LGL<-layout.lgl(S6.IR8, maxit=150,maxdelta=833,area=(length(E(S6.IR8)))^2,coolexp=1.5,repulserad=((length(E(S6.IR8)))^2)*(length(E(S6.IR8))))#Do not work


#Star and Tree layouts

V(S6.IR8)[name=="SIRT6"]
which(V(S6.IR8)[name=="SIRT6"])# 228 #Check what is the rowindex for your desired node
which(V(S6.IR8)$name=="SIRT6")# 191 #Check what is the rowindex for your desired node
which(V(S6.IR8)$name=="MRE11A")# 327
which(V(S6.IR8)$name=="XRCC5")# 722
Lyt.Star<-layout_as_star(S6.IR8,center=V(S6.IR8)[name=="SIRT6"]) #it has a special order argument where we can assign an order to the vertex, hence maybe the OPLS-DA:VIP value order may be used as ordering vector

Lyt.Tree<-layout_as_tree(S6.IR8,root=c(288))
Lyt.RT<-layout.reingold.tilford(S6.IR8,root=c(288))# This one looks good, but the indexing do not allow to highlight the correct nodes
Lyt.sugi<-layout_with_sugiyama(S6.IR8, maxiter=100, weights=E(S6.IR82)$weight)#This one looks good for undirected tree layouts

#Check what layout looks fine

plot.igraph(S6.IR8,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(S6.IR8)$An.color, vertex.frame.color=V(S6.IR8)$An.color, vertex.shape="circle",
            vertex.size=V(S6.IR8)$S6S.FC.IR8, vertex.label=V(S6.IR8)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(S6.IR8)$weight, edge.lty=1, main="Sirt6ID IR8 Network", sub="FR Tree Layout",
            rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#With ggraph

ggraph(S6.IR8,layout=Lyt.RT)+#FR[-isolated,])+
  geom_edge_fan(color="gray65", width=1, alpha=0.4)+
  geom_node_point(aes(color=V(S6.IR8)$Annotation, size=V(S6.IR8)$S6S.FC.IR8),alpha=0.75)+
  geom_node_text(aes(label=name),color="gray25",size=4)+
  scale_color_manual(limits=V(S6.IR8)$Annotation, values=V(S6.IR8)$An.color)+
  scale_size_continuous(range = c(5, 15)) +
  theme_void()+
  title("Sirt6 IR8 Network")+
  guides(color=guide_legend(override.aes = list(size=7)))+
  theme(legend.text=element_text(size=12,face = "bold"),
        legend.title = element_text(size=13,face = "bold"))



#Community Detection
#clp <- igraph::cluster_optimal(S6.IR8)
class(clp)
V(S6.IR8)$community <- clp$membership
colrs <- adjustcolor( c("gray50", "tomato", "gold", "yellowgreen"), alpha=.6)


plot.igraph(clp,S6.IR8,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(S6.IR8)$An.color, vertex.frame.color=V(S6.IR8)$An.color, vertex.shape="circle",
            vertex.size==V(S6.IR8)$S6.FC.IR8, vertex.label=V(S6.IR8)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(S6.IR8)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#Distances edges highlight
dist.from.Sirt6 <- distances(S6.IR8, v=V(S6.IR8)[name=="SIRT6"],
                             to=V(S6.IR8), weights=E(S6.IR8)$weight)
dist.from.Sirt6<-dist.from.Sirt6[,!is.infinite(colSums(dist.from.Sirt6))]

dist.from.Sirt6
range(dist.from.Sirt6)

#Set colors to plot the distances:
oranges <- colorRampPalette(c("dark red", "gold"))
col <- oranges(max(dist.from.Sirt6))
col <- col[dist.from.Sirt6+1]
plot(S6.IR8, layout=Lyt.FR[-isolated,],vertex.color=col, vertex.label=dist.from.Sirt6,
     vertex.label.color="white",vertex.shape="circle",
     vertex.size=V(S6.IR8)$S6S.FC.IR8, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(S6.IR8)$weight, edge.lty=1, main="Sirt6 IR8 Full Nuclear Network", sub="FR Layout",
     rescale=T)


#Shortest path from Sirt6 to Nucleolin
S6Ncl.path <- shortest_paths(S6.IR8,
                             from = V(S6.IR8)[name=="SIRT6"],
                             to = V(S6.IR8)[name=="NCL"],
                             output = "both") # both path nodes and edges
S6Ncl.path
# Generate edge color variable to plot the path:
ecol <- rep("gray80", ecount(S6.IR8))
ecol[unlist(S6Ncl.path$epath)] <- "orange"
# Generate edge width variable to plot the path:
ew <- rep(2, ecount(S6.IR8))
ew[unlist(S6Ncl.path$epath)] <- 4
# Generate node color variable to plot the path:
vcol <- rep("gray40", vcount(S6.IR8))
vcol[unlist(S6Ncl.path$vpath)] <- "gold"
plot(S6.IR8, layout=Lyt.FR[-isolated,], vertex.color=vcol, edge.color=ecol,
     edge.width=ew,vertex.label.color="white",vertex.shape="circle",
     vertex.size=5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(S6.IR8)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)



#Plot based on Sirt6 as incident
inc.edges <- incident(S6.IR8, V(S6.IR8)[name=="SIRT6"], mode="all")
# Set colors to plot the selected edges.
ecol <- rep("gray80", ecount(S6.IR8))
ecol[inc.edges] <- "orange"
vcol <- rep("gray40", vcount(S6.IR8))
vcol[V(S6.IR8)$name=="SIRT6"] <- "gold"

plot(S6.IR8, layout=Lyt.FR[-isolated,], vertex.color=vcol, edge.color=ecol, 
     edge.width=ew,vertex.label.color="white",vertex.shape="circle",
     vertex.size=V(S6.IR8)$S6S.FC.IR8, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(S6.IR8)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)


#Plot Degree Distribution
S6.IR8.Deg<-degree_distribution(S6.IR8, cumulative = T, mode="all")
plot( x=0:max(igraph::degree(S6.IR8)), y=1-S6.IR8.Deg, pch=19, cex=1.2, col="orange",
      xlab="Degree", ylab="Cumulative Frequency")

#Network Descriptives
#Number of nodes and edges
vcount(S6.IR8) #[1] 343
ecount(S6.IR8)#[1]  3216
#Density (Number of edges in relation to the total number of edges in a full connected graph)
edge_density(S6.IR8,loops = F) #[1] 0.05483095
#Transitivity (Total number of conected triangles)
transitivity(S6.IR8,type="global") #[1] 0.6160493
#modularity
mod<-V(S6.IR8)$Mre.FuzzyCMeans #lets check for fuzzy Cmeans modularity
mod<-mod+1#The non clustered interactors that appears as 0, needs to be transformed to 1, so we add 1 to every cluster
modularity(S6.IR8,membership = mod)
#S6.FuzzyCMeans modularity [1] [1] 0.0115187
#Ku.FuzzyCMeans modularity [1] 0.02166503
#Mre.FuzzyCMeans modularity [1] 0.008707531

#Assorativity
assortativity_nominal(S6.IR8,factor(V(S6.IR8)$Mre.FuzzyCMeans),directed=F)
#Annotation [1] [1] 0.2076288
#S6.FuzzyCMeans [1] 0.01595325
#Ku.FuzzyCMeans [1] 0.03519505
#Mre.FuzzyCMeans [1] 0.01078713

#Diameter (longest geodesic distance), shortest path between two most peripheral nodes
diameter(S6.IR8,directed = F, weights = E(S6.IR8)$weights)#[1] 4.008
diam<-get_diameter(S6.IR8,directed = F,weights = E(S6.IR8)$weights)
diam #+ + 9/343 vertices, named, from 8856246:  [1] EMD     BANF1   PARP1   SIRT6   XRCC5   HNRNPA1 SFPQ    CPSF7   NUDT21 
#Plot Diameter
vcol <- rep("gray50", vcount(S6.IR8))
vcol[diam] <- "gold"
ecol <- rep("gray90", ecount(S6.IR8))
ecol[E(S6.IR8, path=diam)] <- "orange"
# E(net, path=diam) finds edges along a path, here 'diam'
plot(S6.IR8, vertex.color=vcol, edge.color=ecol, vertex.label.color="black",edge.arrow.mode=0,layout=Lyt.FR[-isolated,],vertex.size=V(S6.IR8)$S6S.FC.IR8)

#Node topological descriptives
degree<-igraph::centr_degree(S6.IR8,mode = "all",normalized = T)
degree<-degree$res
eigen_cent<-igraph::eigen_centrality(S6.IR8)
eigen_cent<-eigen_cent$vector
clos_cent <- igraph::closeness(S6.IR8, normalized=T)
betw_cent <- igraph::betweenness(S6.IR8, normalized=T)
hs <- hub_score(S6.IR8, weights=E(S6.IR8)$weight)$vector
as <- authority_score(S6.IR8, weights=E(S6.IR8)$weight)$vector
Trans<-transitivity(S6.IR8,type="local")
Trans
#write results
S6.IR8.Centrality<-data.frame(NAME=V(S6.IR8)$name, UniprotID=V(S6.IR8)$Uniprot, Degree=degree, Eigen= eigen_cent, Closeness= clos_cent, Betweeness= betw_cent, Hubs=hs, Authority=as,Transitivity=Trans)
S6.IR8.Centrality
write.csv(S6.IR8.Centrality,"./S6.IR8.Centrality.csv")

#Cocitation, Number of shared interactions among nodes
Cocitation<-cocitation(S6.IR8)
Cocitation#This matrix is similar to an adjacency matrix
write.csv(Cocitation,"./S6.IR8.cocitation.csv")
#Plot with degree as node size
ggraph(S6.IR8,layout=Lyt.FR[-isolated,])+
  geom_edge_fan(color="gray85", width=1, alpha=0.4)+
  geom_node_point(aes(color=V(S6.IR8)$Annotation, size=degree),alpha=0.75)+
  geom_node_text(aes(label=name),color="gray25",size=4)+
  scale_color_manual(limits=V(S6.IR8)$Annotation, values=V(S6.IR8)$An.color)+
  scale_size_continuous(range = c(5, 15)) +
  theme_void()+
  title("Sirt6 IR8 Network")+
  guides(color=guide_legend(override.aes = list(size=7)))+
  theme(legend.text=element_text(size=12,face = "bold"),
        legend.title = element_text(size=13,face = "bold"))


#Modularity Analysis

#Edge_Betweeness Modularity
S6.IR8.ebmod=edge.betweenness.community(S6.IR8,weights=E(S6.IR8)$weight,directed = F) 
S6.IR8.ebmod# IGRAPH clustering edge betweenness, groups: 32, mod: 0.45

#
plot.igraph(S6.IR8.ebmod,S6.IR8,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(S6.IR8)$An.color, vertex.frame.color=V(S6.IR8)$An.color, vertex.shape="circle",
            vertex.size=5, vertex.label=V(S6.IR8)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(S6.IR8)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)

#
plot(S6.IR8.ebmod, S6.IR8,layout=Lyt.FR[-isolated,], alpha=0.6, 
     vertex.color=V(S6.IR8)$An.color, vertex.frame.color=V(S6.IR8)$An.color, vertex.shape="circle", 
     vertex.size=5, vertex.label=V(S6.IR8)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black", 
     edge.width=E(S6.IR8)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout", 
     rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)
#Dendogram of the Edge_Betweeness
dendPlot(S6.IR8.ebmod, mode="hclust")

#Louvain modularity
S6.IR8.Louvain.mod<-cluster_louvain(S6.IR8,weights=E(S6.IR8)$weight)
S6.IR8.Louvain.mod #IGRAPH clustering multi level, groups: 11, mod: 0.5

#Plot the Louvain communities
plot(S6.IR8.Louvain.mod, S6.IR8,layout=Lyt.FR[-isolated,], alpha=0.6, 
     vertex.color=V(S6.IR8)$An.color, vertex.frame.color=V(S6.IR8)$An.color, vertex.shape="circle", 
     vertex.size=5, vertex.label=V(S6.IR8)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black", 
     edge.width=E(S6.IR8)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout", 
     rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#Generate HeatMaps of the adjacency matrixes
#Let's use the neighbor list for generating the Adjaency Matrix
SensorID_FNN #Let's use the edge list that contain the total amount of possible interactors for the three sensors
v<-levels(as.factor(unlist(SensorID_FNN[,1:2]))) #get vertex names from the network list
n<-length(v)#get the number of vertex in the network
e<-matrix(match(as.character(unlist(SensorID_FNN[,1:2])),v),ncol = 2)#Edgelist with row index instead of Names
w<-S6_IR8$weight #Edge weights, only from the S6_IR8 network to quit those edges that do not exist at S6_IR8

# M is our co-expression network adjacency matrix.

# distringuish between i,j and j,i we want to make M symmetric (undirected) by
# considering an undirected network.
S6.IR8.M <- matrix(0, n, n)                  # set up a co-expression matrix
S6.IR8.M[e] <- w                             # fill it in with edge weights
S6.IR8.M <- S6.IR8.M + t(S6.IR8.M)                         # make this symmetric
any(is.na(S6.IR8.M))#check no NAs in the matrix
dimnames(S6.IR8.M) <- list(v, v)             # label the vertices
# We let A_ij = 1 if there is an interaction. 
#A <- 1*(M > 0) #right now we keep the weights as they are
S6.IR8.M #The adjacency Matrix not as an igraph object, easier to manipulate, and to subset for smaller Sensor specific networks
quantile(S6.IR8.M,c(0,0.25,0.5,0.75,1))
#
#       0%   25%   50%   75%  100% 
#     0.000 0.000 0.000 0.000 0.991 

#Creating the HeatMap
# make a heatmap image with colors: white below 0,
# grey from 0.4 to 1 with darker closer to 1
breaks <- c(0, seq(0.4, 0.991, length=11)) # breaks in the color bins
cols <- grey(1-c(0,seq(0.5,1,length=10))) # colors

#Red Blue Colors
library(RColorBrewer)
breaks <- c(seq(-1, 1, length=11))
breaks
cols <- gray(1-c(0,seq(0,1,length=10))) # colors
cols2 <- hcl.colors(10 , palette = "Blues3" , alpha = 1 , rev = TRUE)
cols3 <- colorRamp("red" , space = "rgb" , interpolate = "linear")
HMCol<- brewer.pal(10,name = "RdBu")
HMCol<- rev(HMCol)
HMCol


# color the vertex for CCNB1 blue, its neighbors red, and the others white
SensorIDBind <- which(v == "SIRT6") # index for ccnb1
#SensorIDBind<-as.factor(SensorIDBind)
vcols <- rep("white",n)
vcols[SensorIDBind] <- "blue"
#Subsetting S6.IR8.M by doing a loop asking for look into every row and column the index in SensorIDBind worked
for (i in SensorIDBind) {
  vcols[which(S6.IR8.M[,i] > 0 | S6.IR8.M[i,])] <- "red"
}
##Color annotation for functional annotation
MatSet<-subset(Mat,Mat$NAME%in%v)# select only those proteins that are cointained in the node list of our matrixes
match(v,MatSet$NAME) # Check that the indexes of the Annotation Dataframe and the node list are in the same position 
match(rownames(S6.IR8.M),MatSet$NAME) #Check that the row names of the Adj Matrix and the Annotation dataframe are in the same order
match(colnames(S6.IR8.M),MatSet$NAME) #Check that the column names of the Adj Matrix and the Annotation dataframe are in the same order
AnCol<-MatSet$An.color#Vector of colors
AnCol<-as.character(AnCol)#transform it into a character vector, so the colors can be read

#Subsetting S6.IR8.M using the index vector in SensorIDBind do not work
#vcols[which(S6.IR8.M[,SensorIDBind] > 0 | S6.IR8.M[SensorIDBind,])] <- "red"

#Generate a dendogram for keeping the proteins in the sample place
#S6.IR8.dist<-dist(S6.IR8.M,method = "euclidian")
#S6.IR8.Dend<-hclust(S6.IR8.dist, method = "ward.D" )
#plot(S6.IR8.Dend)

#Cut the dendogram
#cut.height <- 10 # try numbers from 4 to 7
#FNN.hclusts<- cutree(FNN.Dend, cut.height)
#plot(FNN.hclusts)

# now actually make the heat map
par(mar=rep(0,4)) # make plot margins 0
heatmap(S6.IR8.M, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=cols, breaks=breaks,frame=T, cexRow = 0.5 , cexCol = 0.5, main = "Sirt6 IR8 AdjMat", ylab="Total Protein Interactors",xlab="Total Protein Interactors")
#
plot.new()
legend("topleft", c("Cell Cycle","Chromatin Org","DNA Repair","Nuclear Pore","p53","Proteasome","Ribosome Biogenesis","RNA Metabolism","Transcription Reg","UPR"), 
       fill=c("coral1","darkviolet","firebrick1","deepskyblue3","tomato","grey","lightgreen","gold","slateblue1","maroon"),  
       bty="n", inset=0, xpd=T,  border=F, cex =2,y.intersp = 0.7,x.intersp = 0.2,pt.cex = 1)


#Now we substract the S6nonIR8 AdjMatrix to the S6IR6 AdjMatrix
S6IR8.AdjMat.Change<-S6.IR8.M-S6.nonIR5.M
range(S6IR8.AdjMat.Change)#Let's check the magnitude of the substracted weights
#[1]-0.813  0.668
any(is.na(S6IR8.AdjMat.Change))
#Generate the breaks and color functions
breaks <- c(seq(-0.5, 0.5, length=10))
breaks
HMCol<- brewer.pal(9,name = "RdBu")
HMCol<- rev(HMCol)
HMCol
[1] "#2166AC" "#4393C3" "#92C5DE" "#D1E5F0" "#F7F7F7" "#FDDBC7" "#F4A582" "#D6604D" "#B2182B"

#Plot the HeatMap
par(mar=rep(0,4)) # make plot margins 0
heatmap(S6IR8.AdjMat.Change, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=HMCol, breaks=breaks,frame=T, cexRow = 0.5 , cexCol = 0.5, main = "Sirt6 IR8 AdjMat (minus nonIR8 edges)", ylab="Total Protein Interactors",xlab="Total Protein Interactors")


#Now we substract the S6IR5 AdjMatrix to the S6IR8 AdjMatrix
S6IR8.AdjMat.Change2<-S6.IR8.M-S6.IR2.M
range(S6IR8.AdjMat.Change)#Let's check the magnitude of the substracted weights
#[1]-0.813  0.668
any(is.na(S6IR8.AdjMat.Change))
#Generate the breaks and color functions
breaks <- c(seq(-0.5, 0.5, length=10))
breaks
HMCol<- brewer.pal(9,name = "RdBu")
HMCol<- rev(HMCol)
HMCol
[1] "#2166AC" "#4393C3" "#92C5DE" "#D1E5F0" "#F7F7F7" "#FDDBC7" "#F4A582" "#D6604D" "#B2182B"

#Plot the HeatMap
par(mar=rep(0,4)) # make plot margins 0
heatmap(S6IR8.AdjMat.Change2, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=HMCol, breaks=breaks,frame=T, cexRow = 0.5 , cexCol = 0.5, main = "Sirt6 IR8 AdjMat (minus IR2h edges)", ylab="Total Protein Interactors",xlab="Total Protein Interactors")
#
#check with the FNN
# now actually make the heat map
par(mar=rep(0,4)) # make plot margins 0
heatmap(FNN.M, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=0:1, frame=T, cexRow = 0.5 , cexCol = 0.5)
#
legend("topleft", c("Neighbors SIRT6","SIRT6"), fill=c("red","blue"),  
       bty="n", inset=0, xpd=T,  border=F)


####################### Sirt6 IR24h Networks
###################################################

#Sirt6 IR24 Networks


#Let's keep the Mat as the Node attribute table and load the Neighbor list as the matrixes
Mat# Node attribute matrix for all the networks

#Sirt6ID IR24 nodes Network Edge Lists
S6_IR24<-read.delim("./Sirt6_Files/S6_IR24h.tsv",sep = "\t",col.names = )
str(S6_IR24)
head(S6_IR24)

#Generate the igraph network object
S6.IR24<-graph_from_data_frame(d=S6_IR24,vertices = Mat, directed = F )#Generate the igraph object with the edge and node annotation attributes
str(vertex_attr(S6.IR24))#Check node attributes
str(edge_attr(S6.IR24))#Check for edge attributes
S6.IR24<-igraph::simplify(S6.IR24, remove.multiple=TRUE, remove.loops=TRUE) #Remove repetitive edges and self-connecting edges
S6.IR24
#Small plot
plot(S6.IR24)

#Generate an adjacency Matrix that will be used after
S6.IR24_AdjMat<-as_adjacency_matrix(S6.IR24,type="both",attr = "weight",names=TRUE,sparse = TRUE)
S6.IR24_AdjMat<-as.matrix(S6.IR24_AdjMat)
write.csv(S6.IR24_AdjMat,"./S6.IR24_AdjMat.csv")

#Generate a vertex size based on degree and count
vsize<-igraph::degree(S6.IR24)/vcount(S6.IR24)
#Identify and remove disconnected nodes
isolated<-which(igraph::degree(S6.IR24)==0)
isolated
length(isolated)# [1] 621
#remove isolated nodes
S6.IR24<-igraph::delete_vertices(S6.IR24,isolated)


#Generate layouts for the graph
#Do not run again, to not to change the layouts that I already have
Lyt.Nice<-layout_nicely(S6.IR24,dim = 3)
Lyt.GOpt<-layout_with_graphopt(S6.IR24, niter=5000, mass = 30, charge = E(S6.IR24)$weight)
Lyt.DrL<-layout_with_drl(S6.IR24, weights = E(S6.IR24)$weight)
Lyt.FR<-layout.fruchterman.reingold(S6.IR24, niter=10000, weights=E(S6.IR24)$weight)#This one looks like the best option
Lyt.KK<-layout_with_kk(S6.IR24, maxiter=10000, weights=E(S6.IR24)$weight) #very large separation of peripheral nodes
Lyt.GEM<-layout_with_gem(S6.IR24, maxiter=10000) #Bad separation
#Lyt.DH<-layout_with_dh(S6.IR24, maxiter=150, fineiter=max(10, log2(833)))
Lyt.LGL<-layout.lgl(S6.IR24, maxit=150,maxdelta=833,area=(length(E(S6.IR24)))^2,coolexp=1.5,repulserad=((length(E(S6.IR24)))^2)*(length(E(S6.IR24))))#Do not work


#Star and Tree layouts

V(S6.IR24)[name=="SIRT6"]
which(V(S6.IR24)[name=="SIRT6"])# 228 #Check what is the rowindex for your desired node
which(V(S6.IR24)$name=="SIRT6")# 176 #Check what is the rowindex for your desired node
which(V(S6.IR24)$name=="MRE11A")# 327
which(V(S6.IR24)$name=="XRCC5")# 722
Lyt.Star<-layout_as_star(S6.IR24,center=V(S6.IR24)[name=="SIRT6"]) #it has a special order argument where we can assign an order to the vertex, hence maybe the OPLS-DA:VIP value order may be used as ordering vector

Lyt.Tree<-layout_as_tree(S6.IR24,root=c(176))
Lyt.RT<-layout.reingold.tilford(S6.IR24,root=c(176))# This one looks good, but the indexing do not allow to highlight the correct nodes
Lyt.sugi<-layout_with_sugiyama(S6.IR24, maxiter=100, weights=E(S6.IR242)$weight)#This one looks good for undirected tree layouts

#Check what layout looks fine

plot.igraph(S6.IR24,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(S6.IR24)$An.color, vertex.frame.color=V(S6.IR24)$An.color, vertex.shape="circle",
            vertex.size=V(S6.IR24)$S6.FC.IR24, vertex.label=V(S6.IR24)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(S6.IR24)$weight, edge.lty=1, main="Sirt6ID IR24 Network", sub="FR Tree Layout",
            rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#With ggraph

ggraph(S6.IR24,layout=Lyt.FR[-isolated,])+
  geom_edge_fan(color="gray65", width=1, alpha=0.4)+
  geom_node_point(aes(color=V(S6.IR24)$Annotation, size=V(S6.IR24)$S6.FC.IR24),alpha=0.75)+
  geom_node_text(aes(label=name),color="gray25",size=4)+
  scale_color_manual(limits=V(S6.IR24)$Annotation, values=V(S6.IR24)$An.color)+
  scale_size_continuous(range = c(5, 15)) +
  theme_void()+
  title("Sirt6 IR24 Network")+
  guides(color=guide_legend(override.aes = list(size=7)))+
  theme(legend.text=element_text(size=12,face = "bold"),
        legend.title = element_text(size=13,face = "bold"))



#Community Detection
#clp <- igraph::cluster_optimal(S6.IR24)
class(clp)
V(S6.IR24)$community <- clp$membership
colrs <- adjustcolor( c("gray50", "tomato", "gold", "yellowgreen"), alpha=.6)


plot.igraph(clp,S6.IR24,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(S6.IR24)$An.color, vertex.frame.color=V(S6.IR24)$An.color, vertex.shape="circle",
            vertex.size==V(S6.IR24)$S6.FC.IR24, vertex.label=V(S6.IR24)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(S6.IR24)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#Distances edges highlight
dist.from.Sirt6 <- distances(S6.IR24, v=V(S6.IR24)[name=="SIRT6"],
                             to=V(S6.IR24), weights=E(S6.IR24)$weight)
dist.from.Sirt6<-dist.from.Sirt6[,!is.infinite(colSums(dist.from.Sirt6))]

dist.from.Sirt6
range(dist.from.Sirt6)

#Set colors to plot the distances:
oranges <- colorRampPalette(c("dark red", "gold"))
col <- oranges(max(dist.from.Sirt6))
col <- col[dist.from.Sirt6+1]
plot(S6.IR24, layout=Lyt.FR[-isolated,],vertex.color=col, vertex.label=dist.from.Sirt6,
     vertex.label.color="white",vertex.shape="circle",
     vertex.size=V(S6.IR24)$S6.FC.IR24, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(S6.IR24)$weight, edge.lty=1, main="Sirt6 IR24 Full Nuclear Network", sub="FR Layout",
     rescale=T)


#Shortest path from Sirt6 to Nucleolin
S6Ncl.path <- shortest_paths(S6.IR24,
                             from = V(S6.IR24)[name=="SIRT6"],
                             to = V(S6.IR24)[name=="NCL"],
                             output = "both") # both path nodes and edges

# Generate edge color variable to plot the path:
ecol <- rep("gray80", ecount(S6.IR24))
ecol[unlist(S6Ncl.path$epath)] <- "orange"
# Generate edge width variable to plot the path:
ew <- rep(2, ecount(S6.IR24))
ew[unlist(S6Ncl.path$epath)] <- 4
# Generate node color variable to plot the path:
vcol <- rep("gray40", vcount(S6.IR24))
vcol[unlist(S6Ncl.path$vpath)] <- "gold"
plot(S6.IR24, layout=Lyt.FR[-isolated,], vertex.color=vcol, edge.color=ecol,
     edge.width=ew,vertex.label.color="white",vertex.shape="circle",
     vertex.size=5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(S6.IR24)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)



#Plot based on Sirt6 as incident
inc.edges <- incident(S6.IR24, V(S6.IR24)[name=="SIRT6"], mode="all")
# Set colors to plot the selected edges.
ecol <- rep("gray80", ecount(S6.IR24))
ecol[inc.edges] <- "orange"
vcol <- rep("gray40", vcount(S6.IR24))
vcol[V(S6.IR24)$name=="SIRT6"] <- "gold"

plot(S6.IR24, layout=Lyt.FR[-isolated,], vertex.color=vcol, edge.color=ecol, 
     edge.width=ew,vertex.label.color="white",vertex.shape="circle",
     vertex.size=V(S6.IR24)$S6.FC.IR24, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(S6.IR24)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)


#Plot Degree Distribution
S6.IR24.Deg<-degree_distribution(S6.IR24, cumulative = T, mode="all")
plot( x=0:max(igraph::degree(S6.IR24)), y=1-S6.IR24.Deg, pch=19, cex=1.2, col="orange",
      xlab="Degree", ylab="Cumulative Frequency")

#Network Descriptives
#Number of nodes and edges
vcount(S6.IR24) #[1] 212
ecount(S6.IR24)#[1]  1177
#Density (Number of edges in relation to the total number of edges in a full connected graph)
edge_density(S6.IR24,loops = F) #[1] 0.05262452
#Transitivity (Total number of conected triangles)
transitivity(S6.IR24,type="global") #[1] 0.5319826
#modularity
mod<-V(S6.IR24)$Mre.FuzzyCMeans #lets check for fuzzy Cmeans modularity
mod<-mod+1#The non clustered interactors that appears as 0, needs to be transformed to 1, so we add 1 to every cluster
modularity(S6.IR24,membership = mod)
#S6.FuzzyCMeans modularity [1] 0.001586627
#Ku.FuzzyCMeans modularity [1] 0.03026646
#Mre.FuzzyCMeans modularity [1] 0.0112587

#Assorativity
assortativity_nominal(S6.IR24,factor(V(S6.IR24)$Mre.FuzzyCMeans),directed=F)
#Annotation [1] [1] 0.2314902
#S6.FuzzyCMeans [1]0.002410038
#Ku.FuzzyCMeans [1] 0.05135872
#Mre.FuzzyCMeans [1] 0.01452613

#Diameter (longest geodesic distance), shortest path between two most peripheral nodes
diameter(S6.IR24,directed = F, weights = E(S6.IR24)$weights)#[1] 4.799
diam<-get_diameter(S6.IR24,directed = F,weights = E(S6.IR24)$weights)
diam #+ 7/212 vertices, named, from 89c1301:  [1] RPN2   RPN1   TMPO   ILF3   SF3A2  PUF60  SRSF11
#Plot Diameter
vcol <- rep("gray50", vcount(S6.IR24))
vcol[diam] <- "gold"
ecol <- rep("gray90", ecount(S6.IR24))
ecol[E(S6.IR24, path=diam)] <- "orange"
# E(net, path=diam) finds edges along a path, here 'diam'
plot(S6.IR24, vertex.color=vcol, edge.color=ecol, vertex.label.color="black",edge.arrow.mode=0,layout=Lyt.FR[-isolated,],vertex.size=V(S6.IR24)$S6.FC.IR24)

#Node topological descriptives
degree<-igraph::centr_degree(S6.IR24,mode = "all",normalized = T)
degree<-degree$res
eigen_cent<-igraph::eigen_centrality(S6.IR24)
eigen_cent<-eigen_cent$vector
clos_cent <- igraph::closeness(S6.IR24, normalized=T)
betw_cent <- igraph::betweenness(S6.IR24, normalized=T)
hs <- hub_score(S6.IR24, weights=E(S6.IR24)$weight)$vector
as <- authority_score(S6.IR24, weights=E(S6.IR24)$weight)$vector
Trans<-transitivity(S6.IR24,type="local")
Trans
#write results
S6.IR24.Centrality<-data.frame(NAME=V(S6.IR24)$name, UniprotID=V(S6.IR24)$Uniprot, Degree=degree, Eigen= eigen_cent, Closeness= clos_cent, Betweeness= betw_cent, Hubs=hs, Authority=as,Transitivity=Trans)
S6.IR24.Centrality
write.csv(S6.IR24.Centrality,"./S6.IR24.Centrality.csv")

#Cocitation, Number of shared interactions among nodes
Cocitation<-cocitation(S6.IR24)
Cocitation#This matrix is similar to an adjacency matrix
write.csv(Cocitation,"./S6.IR24.cocitation.csv")
#Plot with degree as node size
ggraph(S6.IR24,layout=Lyt.FR[-isolated,])+
  geom_edge_fan(color="gray85", width=1, alpha=0.4)+
  geom_node_point(aes(color=V(S6.IR24)$Annotation, size=degree),alpha=0.75)+
  geom_node_text(aes(label=name),color="gray25",size=4)+
  scale_color_manual(limits=V(S6.IR24)$Annotation, values=V(S6.IR24)$An.color)+
  scale_size_continuous(range = c(5, 15)) +
  theme_void()+
  title("Sirt6 IR24 Network")+
  guides(color=guide_legend(override.aes = list(size=7)))+
  theme(legend.text=element_text(size=12,face = "bold"),
        legend.title = element_text(size=13,face = "bold"))


#Modularity Analysis

#Edge_Betweeness Modularity
S6.IR24.ebmod=edge.betweenness.community(S6.IR24,weights=E(S6.IR24)$weight,directed = F) 
S6.IR24.ebmod# IGRAPH clustering edge betweenness, groups: 22, mod: 0.6

#
plot.igraph(S6.IR24.ebmod,S6.IR24,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(S6.IR24)$An.color, vertex.frame.color=V(S6.IR24)$An.color, vertex.shape="circle",
            vertex.size=5, vertex.label=V(S6.IR24)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(S6.IR24)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)

#
plot(S6.IR24.ebmod, S6.IR24,layout=Lyt.FR[-isolated,], alpha=0.6, 
     vertex.color=V(S6.IR24)$An.color, vertex.frame.color=V(S6.IR24)$An.color, vertex.shape="circle", 
     vertex.size=5, vertex.label=V(S6.IR24)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black", 
     edge.width=E(S6.IR24)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout", 
     rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)
#Dendogram of the Edge_Betweeness
dendPlot(S6.IR24.ebmod, mode="hclust")

#Louvain modularity
S6.IR24.Louvain.mod<-cluster_louvain(S6.IR24,weights=E(S6.IR24)$weight)
S6.IR24.Louvain.mod #IGRAPH clustering multi level, groups: 8, mod: 0.62

#Plot the Louvain communities
plot(S6.IR24.Louvain.mod, S6.IR24,layout=Lyt.FR[-isolated,], alpha=0.6, 
     vertex.color=V(S6.IR24)$An.color, vertex.frame.color=V(S6.IR24)$An.color, vertex.shape="circle", 
     vertex.size=5, vertex.label=V(S6.IR24)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black", 
     edge.width=E(S6.IR24)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout", 
     rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#Generate HeatMaps of the adjacency matrixes
#Let's use the neighbor list for generating the Adjaency Matrix
SensorID_FNN #Let's use the edge list that contain the total amount of possible interactors for the three sensors
v<-levels(as.factor(unlist(SensorID_FNN[,1:2]))) #get vertex names from the network list
n<-length(v)#get the number of vertex in the network
e<-matrix(match(as.character(unlist(SensorID_FNN[,1:2])),v),ncol = 2)#Edgelist with row index instead of Names
w<-S6_IR24$weight #Edge weights, only from the S6_IR24 network to quit those edges that do not exist at S6_IR24

# M is our co-expression network adjacency matrix.

# distringuish between i,j and j,i we want to make M symmetric (undirected) by
# considering an undirected network.
S6.IR24.M <- matrix(0, n, n)                  # set up a co-expression matrix
S6.IR24.M[e] <- w                             # fill it in with edge weights
S6.IR24.M <- S6.IR24.M + t(S6.IR24.M)                         # make this symmetric
any(is.na(S6.IR24.M))#check no NAs in the matrix
dimnames(S6.IR24.M) <- list(v, v)             # label the vertices
# We let A_ij = 1 if there is an interaction. 
#A <- 1*(M > 0) #right now we keep the weights as they are
S6.IR24.M #The adjacency Matrix not as an igraph object, easier to manipulate, and to subset for smaller Sensor specific networks
quantile(S6.IR24.M,c(0,0.25,0.5,0.75,1))
#
#      0%  25%  50%  75% 100% 
#    0.00 0.00 0.00 0.00 0.99   

#Creating the HeatMap
# make a heatmap image with colors: white below 0,
# grey from 0.4 to 1 with darker closer to 1
breaks <- c(0, seq(0.4, 0.991, length=11)) # breaks in the color bins
cols <- grey(1-c(0,seq(0.5,1,length=10))) # colors

#Red Blue Colors
library(RColorBrewer)
breaks <- c(seq(-1, 1, length=11))
breaks
cols <- gray(1-c(0,seq(0,1,length=10))) # colors
cols2 <- hcl.colors(10 , palette = "Blues3" , alpha = 1 , rev = TRUE)
cols3 <- colorRamp("red" , space = "rgb" , interpolate = "linear")
HMCol<- brewer.pal(10,name = "RdBu")
HMCol<- rev(HMCol)
HMCol


# color the vertex for CCNB1 blue, its neighbors red, and the others white
SensorIDBind <- which(v == "SIRT6") # index for ccnb1
#SensorIDBind<-as.factor(SensorIDBind)
vcols <- rep("white",n)
vcols[SensorIDBind] <- "blue"
#Subsetting S6.IR24.M by doing a loop asking for look into every row and column the index in SensorIDBind worked
for (i in SensorIDBind) {
  vcols[which(S6.IR24.M[,i] > 0 | S6.IR24.M[i,])] <- "red"
}
##Color annotation for functional annotation
MatSet<-subset(Mat,Mat$NAME%in%v)# select only those proteins that are cointained in the node list of our matrixes
match(v,MatSet$NAME) # Check that the indexes of the Annotation Dataframe and the node list are in the same position 
match(rownames(S6.IR24.M),MatSet$NAME) #Check that the row names of the Adj Matrix and the Annotation dataframe are in the same order
match(colnames(S6.IR24.M),MatSet$NAME) #Check that the column names of the Adj Matrix and the Annotation dataframe are in the same order
AnCol<-MatSet$An.color#Vector of colors
AnCol<-as.character(AnCol)#transform it into a character vector, so the colors can be read

#Subsetting S6.IR24.M using the index vector in SensorIDBind do not work
#vcols[which(S6.IR24.M[,SensorIDBind] > 0 | S6.IR24.M[SensorIDBind,])] <- "red"

#Generate a dendogram for keeping the proteins in the sample place
#S6.IR24.dist<-dist(S6.IR24.M,method = "euclidian")
#S6.IR24.Dend<-hclust(S6.IR24.dist, method = "ward.D" )
#plot(S6.IR24.Dend)

#Cut the dendogram
#cut.height <- 10 # try numbers from 4 to 7
#FNN.hclusts<- cutree(FNN.Dend, cut.height)
#plot(FNN.hclusts)

# now actually make the heat map
par(mar=rep(0,4)) # make plot margins 0
heatmap(S6.IR24.M, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=cols, breaks=breaks,frame=T, cexRow = 0.5 , cexCol = 0.5, main = "Sirt6 IR24 AdjMat", ylab="Total Protein Interactors",xlab="Total Protein Interactors")
#
plot.new()
legend("topleft", c("Cell Cycle","Chromatin Org","DNA Repair","Nuclear Pore","p53","Proteasome","Ribosome Biogenesis","RNA Metabolism","Transcription Reg","UPR"), 
       fill=c("coral1","darkviolet","firebrick1","deepskyblue3","tomato","grey","lightgreen","gold","slateblue1","maroon"),  
       bty="n", inset=0, xpd=T,  border=F, cex =2,y.intersp = 0.7,x.intersp = 0.2,pt.cex = 1)


#Now we substract the S6nonIR24 AdjMatrix to the S6IR6 AdjMatrix
S6IR24.AdjMat.Change<-S6.IR24.M-S6.nonIR5.M
range(S6IR24.AdjMat.Change)#Let's check the magnitude of the substracted weights
#[1] -0.806  0.805[1] -0.806  0.805
any(is.na(S6IR24.AdjMat.Change))
#Generate the breaks and color functions
breaks <- c(seq(-0.5, 0.5, length=10))
breaks
HMCol<- brewer.pal(9,name = "RdBu")
HMCol<- rev(HMCol)
HMCol
[1] "#2166AC" "#4393C3" "#92C5DE" "#D1E5F0" "#F7F7F7" "#FDDBC7" "#F4A582" "#D6604D" "#B2182B"

#Plot the HeatMap
par(mar=rep(0,4)) # make plot margins 0
heatmap(S6IR24.AdjMat.Change, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=HMCol, breaks=breaks,frame=T, cexRow = 0.5 , cexCol = 0.5, main = "Sirt6 IR24 AdjMat (minus nonIR24 edges)", ylab="Total Protein Interactors",xlab="Total Protein Interactors")


#Now we substract the S6IR5 AdjMatrix to the S6IR24 AdjMatrix
S6IR24.AdjMat.Change2<-S6.IR24.M-S6.IR8.M
range(S6IR24.AdjMat.Change)#Let's check the magnitude of the substracted weights
#[1]-0.813  0.668
any(is.na(S6IR24.AdjMat.Change))
#Generate the breaks and color functions
breaks <- c(seq(-0.5, 0.5, length=10))
breaks
HMCol<- brewer.pal(9,name = "RdBu")
HMCol<- rev(HMCol)
HMCol
[1] "#2166AC" "#4393C3" "#92C5DE" "#D1E5F0" "#F7F7F7" "#FDDBC7" "#F4A582" "#D6604D" "#B2182B"

#Plot the HeatMap
par(mar=rep(0,4)) # make plot margins 0
heatmap(S6IR24.AdjMat.Change2, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=HMCol, breaks=breaks,frame=T, cexRow = 0.5 , cexCol = 0.5, main = "Sirt6 IR24 AdjMat (minus IR8h edges)", ylab="Total Protein Interactors",xlab="Total Protein Interactors")

#
#check with the FNN
# now actually make the heat map
par(mar=rep(0,4)) # make plot margins 0
heatmap(FNN.M, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=0:1, frame=T, cexRow = 0.5 , cexCol = 0.5)
#
legend("topleft", c("Neighbors SIRT6","SIRT6"), fill=c("red","blue"),  
       bty="n", inset=0, xpd=T,  border=F)


####################### Sirt6 nonIR24h Networks
###################################################

#Sirt6 nonIR24 Networks


#Let's keep the Mat as the Node attribute table and load the Neighbor list as the matrixes
Mat# Node attribute matrix for all the networks

#Sirt6ID nonIR24 nodes Network Edge Lists
S6_nonIR24<-read.delim("./Sirt6_Files/S6_nonIR24.tsv",sep = "\t",col.names = )
str(S6_nonIR24)
head(S6_nonIR24)

#Generate the igraph network object
S6.nonIR24<-graph_from_data_frame(d=S6_nonIR24,vertices = Mat, directed = F )#Generate the igraph object with the edge and node annotation attributes
str(vertex_attr(S6.nonIR24))#Check node attributes
str(edge_attr(S6.nonIR24))#Check for edge attributes
S6.nonIR24<-igraph::simplify(S6.nonIR24, remove.multiple=TRUE, remove.loops=TRUE) #Remove repetitive edges and self-connecting edges
S6.nonIR24
#Small plot
plot(S6.nonIR24)

#Generate an adjacency Matrix that will be used after
S6.nonIR24_AdjMat<-as_adjacency_matrix(S6.nonIR24,type="both",attr = "weight",names=TRUE,sparse = TRUE)
S6.nonIR24_AdjMat<-as.matrix(S6.nonIR24_AdjMat)
write.csv(S6.nonIR24_AdjMat,"./S6.nonIR24_AdjMat.csv")

#Generate a vertex size based on degree and count
vsize<-igraph::degree(S6.nonIR24)/vcount(S6.nonIR24)
#Identify and remove disconnected nodes
isolated<-which(igraph::degree(S6.nonIR24)==0)
isolated
length(isolated)# [1] 618
#remove isolated nodes
S6.nonIR24<-igraph::delete_vertices(S6.nonIR24,isolated)


#Generate layouts for the graph
#Do not run again, to not to change the layouts that I already have
Lyt.Nice<-layout_nicely(S6.nonIR24,dim = 3)
Lyt.GOpt<-layout_with_graphopt(S6.nonIR24, niter=5000, mass = 30, charge = E(S6.nonIR24)$weight)
Lyt.DrL<-layout_with_drl(S6.nonIR24, weights = E(S6.nonIR24)$weight)
Lyt.FR<-layout.fruchterman.reingold(S6.nonIR24, niter=10000, weights=E(S6.nonIR24)$weight)#This one looks like the best option
Lyt.KK<-layout_with_kk(S6.nonIR24, maxiter=10000, weights=E(S6.nonIR24)$weight) #very large separation of peripheral nodes
Lyt.GEM<-layout_with_gem(S6.nonIR24, maxiter=10000) #Bad separation
#Lyt.DH<-layout_with_dh(S6.nonIR24, maxiter=150, fineiter=max(10, log2(833)))
Lyt.LGL<-layout.lgl(S6.nonIR24, maxit=150,maxdelta=833,area=(length(E(S6.nonIR24)))^2,coolexp=1.5,repulserad=((length(E(S6.nonIR24)))^2)*(length(E(S6.nonIR24))))#Do not work


#Star and Tree layouts

V(S6.nonIR24)[name=="SIRT6"]
which(V(S6.nonIR24)[name=="SIRT6"])# 228 #Check what is the rowindex for your desired node
which(V(S6.nonIR24)$name=="SIRT6")# 176 #Check what is the rowindex for your desired node
which(V(S6.nonIR24)$name=="MRE11A")# 327
which(V(S6.nonIR24)$name=="XRCC5")# 722
Lyt.Star<-layout_as_star(S6.nonIR24,center=V(S6.nonIR24)[name=="SIRT6"]) #it has a special order argument where we can assign an order to the vertex, hence maybe the OPLS-DA:VIP value order may be used as ordering vector

Lyt.Tree<-layout_as_tree(S6.nonIR24,root=c(180))
Lyt.RT<-layout.reingold.tilford(S6.nonIR24,root=c(180))# This one looks good, but the indexing do not allow to highlight the correct nodes
Lyt.sugi<-layout_with_sugiyama(S6.nonIR24, maxiter=100, weights=E(S6.nonIR242)$weight)#This one looks good for undirected tree layouts

#Check what layout looks fine

plot.igraph(S6.nonIR24,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(S6.nonIR24)$An.color, vertex.frame.color=V(S6.nonIR24)$An.color, vertex.shape="circle",
            vertex.size=V(S6.nonIR24)$S6.FC.nonIR24, vertex.label=V(S6.nonIR24)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(S6.nonIR24)$weight, edge.lty=1, main="Sirt6ID nonIR24 Network", sub="FR Tree Layout",
            rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#With ggraph

ggraph(S6.nonIR24,layout=Lyt.FR[-isolated,])+
  geom_edge_fan(color="gray65", width=1, alpha=0.4)+
  geom_node_point(aes(color=V(S6.nonIR24)$Annotation, size=V(S6.nonIR24)$S6.FC.nonIR24),alpha=0.75)+
  geom_node_text(aes(label=name),color="gray25",size=4)+
  scale_color_manual(limits=V(S6.nonIR24)$Annotation, values=V(S6.nonIR24)$An.color)+
  scale_size_continuous(range = c(5, 15)) +
  theme_void()+
  title("Sirt6 nonIR24 Network")+
  guides(color=guide_legend(override.aes = list(size=7)))+
  theme(legend.text=element_text(size=12,face = "bold"),
        legend.title = element_text(size=13,face = "bold"))



#Community Detection
#clp <- igraph::cluster_optimal(S6.nonIR24)
class(clp)
V(S6.nonIR24)$community <- clp$membership
colrs <- adjustcolor( c("gray50", "tomato", "gold", "yellowgreen"), alpha=.6)


plot.igraph(clp,S6.nonIR24,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(S6.nonIR24)$An.color, vertex.frame.color=V(S6.nonIR24)$An.color, vertex.shape="circle",
            vertex.size==V(S6.nonIR24)$S6.FC.nonIR24, vertex.label=V(S6.nonIR24)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(S6.nonIR24)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#Distances edges highlight
dist.from.Sirt6 <- distances(S6.nonIR24, v=V(S6.nonIR24)[name=="SIRT6"],
                             to=V(S6.nonIR24), weights=E(S6.nonIR24)$weight)
dist.from.Sirt6<-dist.from.Sirt6[,!is.infinite(colSums(dist.from.Sirt6))]

dist.from.Sirt6
range(dist.from.Sirt6)

#Set colors to plot the distances:
oranges <- colorRampPalette(c("dark red", "gold"))
col <- oranges(max(dist.from.Sirt6))
col <- col[dist.from.Sirt6+1]
plot(S6.nonIR24, layout=Lyt.FR[-isolated,],vertex.color=col, vertex.label=dist.from.Sirt6,
     vertex.label.color="white",vertex.shape="circle",
     vertex.size=V(S6.nonIR24)$S6.FC.nonIR24, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(S6.nonIR24)$weight, edge.lty=1, main="Sirt6 nonIR24 Full Nuclear Network", sub="FR Layout",
     rescale=T)


#Shortest path from Sirt6 to Nucleolin
S6Ncl.path <- shortest_paths(S6.nonIR24,
                             from = V(S6.nonIR24)[name=="SIRT6"],
                             to = V(S6.nonIR24)[name=="NCL"],
                             output = "both") # both path nodes and edges

# Generate edge color variable to plot the path:
ecol <- rep("gray80", ecount(S6.nonIR24))
ecol[unlist(S6Ncl.path$epath)] <- "orange"
# Generate edge width variable to plot the path:
ew <- rep(2, ecount(S6.nonIR24))
ew[unlist(S6Ncl.path$epath)] <- 4
# Generate node color variable to plot the path:
vcol <- rep("gray40", vcount(S6.nonIR24))
vcol[unlist(S6Ncl.path$vpath)] <- "gold"
plot(S6.nonIR24, layout=Lyt.FR[-isolated,], vertex.color=vcol, edge.color=ecol,
     edge.width=ew,vertex.label.color="white",vertex.shape="circle",
     vertex.size=5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(S6.nonIR24)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)



#Plot based on Sirt6 as incident
inc.edges <- incident(S6.nonIR24, V(S6.nonIR24)[name=="SIRT6"], mode="all")
# Set colors to plot the selected edges.
ecol <- rep("gray80", ecount(S6.nonIR24))
ecol[inc.edges] <- "orange"
vcol <- rep("gray40", vcount(S6.nonIR24))
vcol[V(S6.nonIR24)$name=="SIRT6"] <- "gold"

plot(S6.nonIR24, layout=Lyt.FR[-isolated,], vertex.color=vcol, edge.color=ecol, 
     edge.width=ew,vertex.label.color="white",vertex.shape="circle",
     vertex.size=V(S6.nonIR24)$S6.FC.nonIR24, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(S6.nonIR24)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)


#Plot Degree Distribution
S6.nonIR24.Deg<-degree_distribution(S6.nonIR24, cumulative = T, mode="all")
plot( x=0:max(igraph::degree(S6.nonIR24)), y=1-S6.nonIR24.Deg, pch=19, cex=1.2, col="orange",
      xlab="Degree", ylab="Cumulative Frequency")

#Network Descriptives
#Number of nodes and edges
vcount(S6.nonIR24) #[1] 215
ecount(S6.nonIR24)#[1] 1398
#Density (Number of edges in relation to the total number of edges in a full connected graph)
edge_density(S6.nonIR24,loops = F) #[1] 0.0607694
#Transitivity (Total number of conected triangles)
transitivity(S6.nonIR24,type="global") #[1] 0.5490721
#modularity
mod<-V(S6.nonIR24)$S6.FuzzyCMeans #lets check for fuzzy Cmeans modularity
mod<-mod+1#The non clustered interactors that appears as 0, needs to be transformed to 1, so we add 1 to every cluster
modularity(S6.nonIR24,membership = mod)
#S6.FuzzyCMeans modularity [1]  0.01404955
#Ku.FuzzyCMeans modularity [1] 0.01062933
#Mre.FuzzyCMeans modularity [1] 0.007612039

#Assorativity
assortativity_nominal(S6.nonIR24,factor(V(S6.nonIR24)$Annotation),directed=F)
#Annotation [1] [1] 0.1898548
#S6.FuzzyCMeans [1]0.0212191
#Ku.FuzzyCMeans [1] 0.01860126
#Mre.FuzzyCMeans [1] 0.009590151

#Diameter (longest geodesic distance), shortest path between two most peripheral nodes
diameter(S6.nonIR24,directed = F, weights = E(S6.nonIR24)$weights)#[1] 4.799
diam<-get_diameter(S6.nonIR24,directed = F,weights = E(S6.nonIR24)$weights)
diam #+ 6/215 vertices, named, from e0161ea:[1] CCT4   ACTG1  RUVBL2 XRCC5  ORC2   MCM3  
#Plot Diameter
vcol <- rep("gray50", vcount(S6.nonIR24))
vcol[diam] <- "gold"
ecol <- rep("gray90", ecount(S6.nonIR24))
ecol[E(S6.nonIR24, path=diam)] <- "orange"
# E(net, path=diam) finds edges along a path, here 'diam'
plot(S6.nonIR24, vertex.color=vcol, edge.color=ecol, vertex.label.color="black",edge.arrow.mode=0,layout=Lyt.FR[-isolated,],vertex.size=V(S6.nonIR24)$S6.FC.nonIR24)

#Node topological descriptives
degree<-igraph::centr_degree(S6.nonIR24,mode = "all",normalized = T)
degree<-degree$res
eigen_cent<-igraph::eigen_centrality(S6.nonIR24)
eigen_cent<-eigen_cent$vector
clos_cent <- igraph::closeness(S6.nonIR24, normalized=T)
betw_cent <- igraph::betweenness(S6.nonIR24, normalized=T)
hs <- hub_score(S6.nonIR24, weights=E(S6.nonIR24)$weight)$vector
as <- authority_score(S6.nonIR24, weights=E(S6.nonIR24)$weight)$vector
Trans<-transitivity(S6.nonIR24,type="local")
Trans
#write results
S6.nonIR24.Centrality<-data.frame(NAME=V(S6.nonIR24)$name, UniprotID=V(S6.nonIR24)$Uniprot, Degree=degree, Eigen= eigen_cent, Closeness= clos_cent, Betweeness= betw_cent, Hubs=hs, Authority=as,Transitivity=Trans)
S6.nonIR24.Centrality
write.csv(S6.nonIR24.Centrality,"./S6.nonIR24.Centrality.csv")

#Cocitation, Number of shared interactions among nodes
Cocitation<-cocitation(S6.nonIR24)
Cocitation#This matrix is similar to an adjacency matrix
write.csv(Cocitation,"./S6.nonIR24.cocitation.csv")
#Plot with degree as node size
ggraph(S6.nonIR24,layout=Lyt.FR[-isolated,])+
  geom_edge_fan(color="gray85", width=1, alpha=0.4)+
  geom_node_point(aes(color=V(S6.nonIR24)$Annotation, size=degree),alpha=0.75)+
  geom_node_text(aes(label=name),color="gray25",size=4)+
  scale_color_manual(limits=V(S6.nonIR24)$Annotation, values=V(S6.nonIR24)$An.color)+
  scale_size_continuous(range = c(5, 15)) +
  theme_void()+
  title("Sirt6 nonIR24 Network")+
  guides(color=guide_legend(override.aes = list(size=7)))+
  theme(legend.text=element_text(size=12,face = "bold"),
        legend.title = element_text(size=13,face = "bold"))


#Modularity Analysis

#Edge_Betweeness Modularity
S6.nonIR24.ebmod=edge.betweenness.community(S6.nonIR24,weights=E(S6.nonIR24)$weight,directed = F) 
S6.nonIR24.ebmod# IGRAPH clustering edge betweenness, groups: 21, mod: 0.47

#
plot.igraph(S6.nonIR24.ebmod,S6.nonIR24,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(S6.nonIR24)$An.color, vertex.frame.color=V(S6.nonIR24)$An.color, vertex.shape="circle",
            vertex.size=5, vertex.label=V(S6.nonIR24)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(S6.nonIR24)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)

#
plot(S6.nonIR24.ebmod, S6.nonIR24,layout=Lyt.FR[-isolated,], alpha=0.6, 
     vertex.color=V(S6.nonIR24)$An.color, vertex.frame.color=V(S6.nonIR24)$An.color, vertex.shape="circle", 
     vertex.size=5, vertex.label=V(S6.nonIR24)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black", 
     edge.width=E(S6.nonIR24)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout", 
     rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)
#Dendogram of the Edge_Betweeness
dendPlot(S6.nonIR24.ebmod, mode="hclust")

#Louvain modularity
S6.nonIR24.Louvain.mod<-cluster_louvain(S6.nonIR24,weights=E(S6.nonIR24)$weight)
S6.nonIR24.Louvain.mod #IGRAPH clustering multi level, groups: 9, mod: 0.51

#Plot the Louvain communities
plot(S6.nonIR24.Louvain.mod, S6.nonIR24,layout=Lyt.FR[-isolated,], alpha=0.6, 
     vertex.color=V(S6.nonIR24)$An.color, vertex.frame.color=V(S6.nonIR24)$An.color, vertex.shape="circle", 
     vertex.size=5, vertex.label=V(S6.nonIR24)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black", 
     edge.width=E(S6.nonIR24)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout", 
     rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#Generate HeatMaps of the adjacency matrixes
#Let's use the neighbor list for generating the Adjaency Matrix
SensorID_FNN #Let's use the edge list that contain the total amount of possible interactors for the three sensors
v<-levels(as.factor(unlist(SensorID_FNN[,1:2]))) #get vertex names from the network list
n<-length(v)#get the number of vertex in the network
e<-matrix(match(as.character(unlist(SensorID_FNN[,1:2])),v),ncol = 2)#Edgelist with row index instead of Names
w<-S6_nonIR24$weight #Edge weights, only from the S6_nonIR24 network to quit those edges that do not exist at S6_nonIR24

# M is our co-expression network adjacency matrix.

# distringuish between i,j and j,i we want to make M symmetric (undirected) by
# considering an undirected network.
S6.nonIR24.M <- matrix(0, n, n)                  # set up a co-expression matrix
S6.nonIR24.M[e] <- w                             # fill it in with edge weights
S6.nonIR24.M <- S6.nonIR24.M + t(S6.nonIR24.M)                         # make this symmetric
any(is.na(S6.nonIR24.M))#check no NAs in the matrix
dimnames(S6.nonIR24.M) <- list(v, v)             # label the vertices
# We let A_ij = 1 if there is an interaction. 
#A <- 1*(M > 0) #right now we keep the weights as they are
S6.nonIR24.M #The adjacency Matrix not as an igraph object, easier to manipulate, and to subset for smaller Sensor specific networks
quantile(S6.nonIR24.M,c(0,0.25,0.5,0.75,1))
#
#      0%  25%  50%  75% 100% 
#    0.00 0.00 0.00 0.00 0.99   

#Creating the HeatMap
# make a heatmap image with colors: white below 0,
# grey from 0.4 to 1 with darker closer to 1
breaks <- c(0, seq(0.4, 0.991, length=11)) # breaks in the color bins
cols <- grey(1-c(0,seq(0.5,1,length=10))) # colors

#Red Blue Colors
library(RColorBrewer)
breaks <- c(seq(-1, 1, length=11))
breaks
cols <- gray(1-c(0,seq(0,1,length=10))) # colors
cols2 <- hcl.colors(10 , palette = "Blues3" , alpha = 1 , rev = TRUE)
cols3 <- colorRamp("red" , space = "rgb" , interpolate = "linear")
HMCol<- brewer.pal(10,name = "RdBu")
HMCol<- rev(HMCol)
HMCol


# color the vertex for CCNB1 blue, its neighbors red, and the others white
SensorIDBind <- which(v == "SIRT6") # index for ccnb1
#SensorIDBind<-as.factor(SensorIDBind)
vcols <- rep("white",n)
vcols[SensorIDBind] <- "blue"
#Subsetting S6.nonIR24.M by doing a loop asking for look into every row and column the index in SensorIDBind worked
for (i in SensorIDBind) {
  vcols[which(S6.nonIR24.M[,i] > 0 | S6.nonIR24.M[i,])] <- "red"
}
##Color annotation for functional annotation
MatSet<-subset(Mat,Mat$NAME%in%v)# select only those proteins that are cointained in the node list of our matrixes
match(v,MatSet$NAME) # Check that the indexes of the Annotation Dataframe and the node list are in the same position 
match(rownames(S6.nonIR24.M),MatSet$NAME) #Check that the row names of the Adj Matrix and the Annotation dataframe are in the same order
match(colnames(S6.nonIR24.M),MatSet$NAME) #Check that the column names of the Adj Matrix and the Annotation dataframe are in the same order
AnCol<-MatSet$An.color#Vector of colors
AnCol<-as.character(AnCol)#transform it into a character vector, so the colors can be read

#Subsetting S6.nonIR24.M using the index vector in SensorIDBind do not work
#vcols[which(S6.nonIR24.M[,SensorIDBind] > 0 | S6.nonIR24.M[SensorIDBind,])] <- "red"

#Generate a dendogram for keeping the proteins in the sample place
#S6.nonIR24.dist<-dist(S6.nonIR24.M,method = "euclidian")
#S6.nonIR24.Dend<-hclust(S6.nonIR24.dist, method = "ward.D" )
#plot(S6.nonIR24.Dend)

#Cut the dendogram
#cut.height <- 10 # try numbers from 4 to 7
#FNN.hclusts<- cutree(FNN.Dend, cut.height)
#plot(FNN.hclusts)

# now actually make the heat map
par(mar=rep(0,4)) # make plot margins 0
heatmap(S6.nonIR24.M, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=cols, breaks=breaks,frame=T, cexRow = 0.5 , cexCol = 0.5, main = "Sirt6 nonIR24 AdjMat", ylab="Total Protein Interactors",xlab="Total Protein Interactors")
#
plot.new()
legend("topleft", c("Cell Cycle","Chromatin Org","DNA Repair","Nuclear Pore","p53","Proteasome","Ribosome Biogenesis","RNA Metabolism","Transcription Reg","UPR"), 
       fill=c("coral1","darkviolet","firebrick1","deepskyblue3","tomato","grey","lightgreen","gold","slateblue1","maroon"),  
       bty="n", inset=0, xpd=T,  border=F, cex =2,y.intersp = 0.7,x.intersp = 0.2,pt.cex = 1)


#Now we substract the S6nonnonIR24 AdjMatrix to the S6IR6 AdjMatrix
S6nonIR24.AdjMat.Change<-S6.nonIR24.M-S6.nonIR5.M
range(S6nonIR24.AdjMat.Change)#Let's check the magnitude of the substracted weights
#[1] -0.806  0.805[1] -0.806  0.805
any(is.na(S6nonIR24.AdjMat.Change))
#Generate the breaks and color functions
breaks <- c(seq(-0.5, 0.5, length=10))
breaks
HMCol<- brewer.pal(9,name = "RdBu")
HMCol<- rev(HMCol)
HMCol
[1] "#2166AC" "#4393C3" "#92C5DE" "#D1E5F0" "#F7F7F7" "#FDDBC7" "#F4A582" "#D6604D" "#B2182B"

#Plot the HeatMap
par(mar=rep(0,4)) # make plot margins 0
heatmap(S6nonIR24.AdjMat.Change, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=HMCol, breaks=breaks,frame=T, cexRow = 0.5 , cexCol = 0.5, main = "Sirt6 nonIR24 AdjMat (minus nonnonIR24 edges)", ylab="Total Protein Interactors",xlab="Total Protein Interactors")


#Now we substract the S6IR5 AdjMatrix to the S6nonIR24 AdjMatrix
S6nonIR24.AdjMat.Change2<-S6.nonIR24.M-S6.IR24.M
range(S6nonIR24.AdjMat.Change)#Let's check the magnitude of the substracted weights
#[1]-0.813  0.668
any(is.na(S6nonIR24.AdjMat.Change))
#Generate the breaks and color functions
breaks <- c(seq(-0.5, 0.5, length=10))
breaks
HMCol<- brewer.pal(9,name = "RdBu")
HMCol<- rev(HMCol)
HMCol
[1] "#2166AC" "#4393C3" "#92C5DE" "#D1E5F0" "#F7F7F7" "#FDDBC7" "#F4A582" "#D6604D" "#B2182B"

#Plot the HeatMap
par(mar=rep(0,4)) # make plot margins 0
heatmap(S6nonIR24.AdjMat.Change2, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=HMCol, breaks=breaks,frame=T, cexRow = 0.5 , cexCol = 0.5, main = "Sirt6 nonIR24 AdjMat (minus IR24 edges)", ylab="Total Protein Interactors",xlab="Total Protein Interactors")

#
#check with the FNN
# now actually make the heat map
par(mar=rep(0,4)) # make plot margins 0
heatmap(FNN.M, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=0:1, frame=T, cexRow = 0.5 , cexCol = 0.5)
#
legend("topleft", c("Neighbors SIRT6","SIRT6"), fill=c("red","blue"),  
       bty="n", inset=0, xpd=T,  border=F)

#
dev.off()
plot.new()


###############################################
#Sirt6ID FNN nodes Network Edge Lists: This is the Sirt6ID Full Nuclear network, contrasting with the previous FNN, this one only contains the Sirt6ID interactors
S6_FNN<-read.delim("./Sirt6_Files/S6_FNN.tsv",sep = "\t",col.names = )
str(S6_FNN)
head(S6_FNN)

#Generate the igraph network object
S6.FNN<-graph_from_data_frame(d=S6_FNN,vertices = Mat, directed = F )#Generate the igraph object with the edge and node annotation attributes
str(vertex_attr(S6.FNN))#Check node attributes
str(edge_attr(S6.FNN))#Check for edge attributes
S6.FNN<-igraph::simplify(S6.FNN, remove.multiple=TRUE, remove.loops=TRUE) #Remove repetitive edges and self-connecting edges
S6.FNN
#Small plot
plot(S6.FNN)

#Generate an adjacency Matrix that will be used after
S6.FNN_AdjMat<-as_adjacency_matrix(S6.FNN,type="both",attr = "weight",names=TRUE,sparse = TRUE)
S6.FNN_AdjMat<-as.matrix(S6.FNN_AdjMat)
write.csv(S6.FNN_AdjMat,"./S6.FNN_AdjMat.csv")

#Generate a vertex size based on degree and count
vsize<-igraph::degree(S6.FNN)/vcount(S6.FNN)
#Identify and remove disconnected nodes
isolated<-which(igraph::degree(S6.FNN)==0)
isolated
length(isolated)# [1] 599
#remove isolated nodes
S6.FNN<-igraph::delete_vertices(S6.FNN,isolated)


#Generate layouts for the graph
#Do not run again, to not to change the layouts that I already have
Lyt.Nice<-layout_nicely(S6.FNN,dim = 3)
Lyt.GOpt<-layout_with_graphopt(S6.FNN, niter=5000, mass = 30, charge = E(S6.FNN)$weight)
Lyt.DrL<-layout_with_drl(S6.FNN, weights = E(S6.FNN)$weight)
Lyt.FR<-layout.fruchterman.reingold(S6.FNN, niter=10000, weights=E(S6.FNN)$weight)#This one looks like the best option
Lyt.KK<-layout_with_kk(S6.FNN, maxiter=10000, weights=E(S6.FNN)$weight) #very large separation of peripheral nodes
Lyt.GEM<-layout_with_gem(S6.FNN, maxiter=10000) #Bad separation
#Lyt.DH<-layout_with_dh(S6.FNN, maxiter=150, fineiter=max(10, log2(833)))
Lyt.LGL<-layout.lgl(S6.FNN, maxit=150,maxdelta=833,area=(length(E(S6.FNN)))^2,coolexp=1.5,repulserad=((length(E(S6.FNN)))^2)*(length(E(S6.FNN))))#Do not work


#Star and Tree layouts

V(S6.FNN)[name=="SIRT6"]
which(V(S6.FNN)[name=="SIRT6"])# 610 #Check what is the rowindex for your desired node
which(V(S6.FNN)$name=="SIRT6")# 610 #Check what is the rowindex for your desired node
which(V(S6.FNN)$name=="MRE11A")# 327
which(V(S6.FNN)$name=="XRCC5")# 722
Lyt.Star<-layout_as_star(S6.FNN,center=V(S6.FNN)[name=="SIRT6"]) #it has a special order argument where we can assign an order to the vertex, hence maybe the OPLS-DA:VIP value order may be used as ordering vector

Lyt.Tree<-layout_as_tree(S6.FNN,root=c(190))
Lyt.RT<-layout.reingold.tilford(S6.FNN,root=c(190))# This one looks good, but the indexing do not allow to highlight the correct nodes
Lyt.sugi<-layout_with_sugiyama(S6.FNN, maxiter=100, weights=E(S6.FNN2)$weight)#This one looks good for undirected tree layouts

#Check what layout looks fine

plot.igraph(S6.FNN,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(S6.FNN)$An.color, vertex.frame.color=V(S6.FNN)$An.color, vertex.shape="circle",
            vertex.size=5, vertex.label=V(S6.FNN)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(S6.FNN)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)

#With ggraph

ggraph(S6.FNN,layout=Lyt.FR[-isolated,])+
  geom_edge_fan(color="gray65", width=1, alpha=0.4)+
  geom_node_point(aes(color=V(S6.FNN)$Annotation, size=V(S6.FNN)$S6.FC.nonIR5),alpha=0.75)+
  geom_node_text(aes(label=name),color="gray25",size=4)+
  scale_color_manual(limits=V(S6.FNN)$Annotation, values=V(S6.FNN)$An.color)+
  scale_size_continuous(range = c(5, 15)) +
  theme_void()+
  title("Sirt6 nonIR5 Network")+
  guides(color=guide_legend(override.aes = list(size=7)))+
  theme(legend.text=element_text(size=12,face = "bold"),
        legend.title = element_text(size=13,face = "bold"))



#Community Detection
#clp <- igraph::cluster_optimal(S6.FNN)
class(clp)
V(S6.FNN)$community <- clp$membership
colrs <- adjustcolor( c("gray50", "tomato", "gold", "yellowgreen"), alpha=.6)


plot.igraph(clp,S6.FNN,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(S6.FNN)$An.color, vertex.frame.color=V(S6.FNN)$An.color, vertex.shape="circle",
            vertex.size==V(S6.FNN)$S6.FC.nonIR5, vertex.label=V(S6.FNN)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(S6.FNN)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#Distances edges highlight
dist.from.Sirt6 <- distances(S6.FNN, v=V(S6.FNN)[name=="SIRT6"],
                             to=V(S6.FNN), weights=E(S6.FNN)$weight)
dist.from.Sirt6<-dist.from.Sirt6[,!is.infinite(colSums(dist.from.Sirt6))]

dist.from.Sirt6
range(dist.from.Sirt6)

#Set colors to plot the distances:
oranges <- colorRampPalette(c("dark red", "gold"))
col <- oranges(max(dist.from.Sirt6))
col <- col[dist.from.Sirt6+1]
plot(S6.FNN, layout=Lyt.FR[-isolated,],vertex.color=col, vertex.label=dist.from.Sirt6,
     vertex.label.color="white",vertex.shape="circle",
     vertex.size=V(S6.FNN)$S6.FC.nonIR5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(S6.FNN)$weight, edge.lty=1, main="Sirt6 nonIR5 Full Nuclear Network", sub="FR Layout",
     rescale=T)


#Shortest path from Sirt6 to Nucleolin
S6Ncl.path <- shortest_paths(S6.FNN,
                             from = V(S6.FNN)[name=="SIRT6"],
                             to = V(S6.FNN)[name=="NCL"],
                             output = "both") # both path nodes and edges
# Generate edge color variable to plot the path:
ecol <- rep("gray80", ecount(S6.FNN))
ecol[unlist(S6Ncl.path$epath)] <- "orange"
# Generate edge width variable to plot the path:
ew <- rep(2, ecount(S6.FNN))
ew[unlist(S6Ncl.path$epath)] <- 4
# Generate node color variable to plot the path:
vcol <- rep("gray40", vcount(S6.FNN))
vcol[unlist(S6Ncl.path$vpath)] <- "gold"
plot(S6.FNN, layout=Lyt.FR[-isolated,], vertex.color=vcol, edge.color=ecol,
     edge.width=ew,vertex.label.color="white",vertex.shape="circle",
     vertex.size=5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(S6.FNN)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)



#Plot based on Sirt6 as incident
inc.edges <- incident(S6.FNN, V(S6.FNN)[name=="SIRT6"], mode="all")
# Set colors to plot the selected edges.
ecol <- rep("gray80", ecount(S6.FNN))
ecol[inc.edges] <- "orange"
vcol <- rep("gray40", vcount(S6.FNN))
vcol[V(S6.FNN)$name=="SIRT6"] <- "gold"

plot(S6.FNN, layout=Lyt.FR[-isolated,], vertex.color=vcol, edge.color=ecol, 
     edge.width=ew,vertex.label.color="white",vertex.shape="circle",
     vertex.size=V(S6.FNN)$S6.FC.nonIR5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(S6.FNN)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)


#Plot Degree Distribution
S6.FNN.Deg<-degree_distribution(S6.FNN, cumulative = T, mode="all")
plot( x=0:max(igraph::degree(S6.FNN)), y=1-S6.FNN.Deg, pch=19, cex=1.2, col="orange",
      xlab="Degree", ylab="Cumulative Frequency")

#Network Descriptives
#Number of nodes and edges
vcount(S6.FNN) #[1] 463
ecount(S6.FNN)#[1] 5002
#Density (Number of edges in relation to the total number of edges in a full connected graph)
edge_density(S6.FNN,loops = F) #[1] 0.04676821
#Transitivity (Total number of conected triangles)
transitivity(S6.FNN,type="global") #[1] 0.5871439
#modularity
mod<-V(S6.FNN)$S6.FuzzyCMeans #lets check for fuzzy Cmeans modularity
mod<-mod+1#The non clustered interactors that appears as 0, needs to be transformed to 1, so we add 1 to every cluster
modularity(S6.FNN,membership = mod)
#S6.FuzzyCMeans modularity [1] [1] 0.007007653
#Ku.FuzzyCMeans modularity [1] 0.0208517
#Mre.FuzzyCMeans modularity [1]0.01368445

#Assorativity
assortativity_nominal(S6.FNN,factor(V(S6.FNN)$Annotation),directed=F)
#Annotation [1] [1]0.2215214
#S6.FuzzyCMeans [1] 0.009339672
#Ku.FuzzyCMeans [1] 0.0354847
#Mre.FuzzyCMeans [1]0.01682863

#Diameter (longest geodesic distance), shortest path between two most peripheral nodes
diameter(S6.FNN,directed = F, weights = E(S6.FNN)$weights)#3.361
diam<-get_diameter(S6.FNN,directed = F,weights = E(S6.FNN)$weights)
diam #+ 6/234 vertices, named, from 74b8ef8:[1] DDX46    SNRPD3   HIST1H4F TMPO     RPN1     RPN2  
#Plot Diameter
vcol <- rep("gray50", vcount(S6.FNN))
vcol[diam] <- "gold"
ecol <- rep("gray90", ecount(S6.FNN))
ecol[E(S6.FNN, path=diam)] <- "orange"
# E(net, path=diam) finds edges along a path, here 'diam'
plot(S6.FNN, vertex.color=vcol, edge.color=ecol, vertex.label.color="black",edge.arrow.mode=0,layout=Lyt.FR[-isolated,],vertex.size=V(S6.FNN)$S6.FC.nonIR5)

#Node topological descriptives
degree<-igraph::centr_degree(S6.FNN,mode = "all",normalized = T)
degree<-degree$res
eigen_cent<-igraph::eigen_centrality(S6.FNN)
eigen_cent<-eigen_cent$vector
clos_cent <- igraph::closeness(S6.FNN, normalized=T)
betw_cent <- igraph::betweenness(S6.FNN, normalized=T)
hs <- hub_score(S6.FNN, weights=E(S6.FNN)$weight)$vector
as <- authority_score(S6.FNN, weights=E(S6.FNN)$weight)$vector
Trans<-transitivity(S6.FNN,type="local")
Trans
#write results
S6.FNN.Centrality<-data.frame(NAME=V(S6.FNN)$name, UniprotID=V(S6.FNN)$Uniprot, Degree=degree, Eigen= eigen_cent, Closeness= clos_cent, Betweeness= betw_cent, Hubs=hs, Authority=as,Transitivity=Trans)
S6.FNN.Centrality
write.csv(S6.FNN.Centrality,"./S6.FNN.Centrality.csv")

#Cocitation, Number of shared interactions among nodes
Cocitation<-cocitation(S6.FNN)
Cocitation#This matrix is similar to an adjacency matrix
write.csv(Cocitation,"./S6.FNN.cocitation.csv")
#Plot with degree as node size
ggraph(S6.FNN,layout=Lyt.FR[-isolated,])+
  geom_edge_fan(color="gray85", width=1, alpha=0.4)+
  geom_node_point(aes(color=V(S6.FNN)$Annotation, size=degree),alpha=0.75)+
  geom_node_text(aes(label=name),color="gray25",size=4)+
  scale_color_manual(limits=V(S6.FNN)$Annotation, values=V(S6.FNN)$An.color)+
  scale_size_continuous(range = c(5, 15)) +
  theme_void()+
  title("Sirt6 nonIR5 Network")+
  guides(color=guide_legend(override.aes = list(size=7)))+
  theme(legend.text=element_text(size=12,face = "bold"),
        legend.title = element_text(size=13,face = "bold"))


#Modularity Analysis

#Edge_Betweeness Modularity
S6.FNN.ebmod=edge.betweenness.community(S6.FNN,weights=E(S6.FNN)$weight,directed = F) 
S6.FNN.ebmod# groups: 15, mod: 0.48
S6.FNN.ebmod=edge.betweenness.community(S6.FNN,weights=E(S6.FNN)$weight,directed = F) 
S6.FNN.ebmod #37 groups, good enough mod 0.47
#
plot.igraph(S6.FNN.ebmod,S6.FNN,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(S6.FNN)$An.color, vertex.frame.color=V(S6.FNN)$An.color, vertex.shape="circle",
            vertex.size=5, vertex.label=V(S6.FNN)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(S6.FNN)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)

#
plot(S6.FNN.ebmod, S6.FNN,layout=Lyt.FR[-isolated,], alpha=0.6, 
     vertex.color=V(S6.FNN)$An.color, vertex.frame.color=V(S6.FNN)$An.color, vertex.shape="circle", 
     vertex.size=5, vertex.label=V(S6.FNN)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black", 
     edge.width=E(S6.FNN)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout", 
     rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)
#Dendogram of the Edge_Betweeness
dendPlot(S6.FNN.ebmod, mode="hclust")

#Louvain modularity
S6.FNN.Louvain.mod<-cluster_louvain(S6.FNN,weights=E(S6.FNN)$weight)
S6.FNN.Louvain.mod #IGRAPH clustering multi level, groups: 13, mod: 0.53

#Plot the Louvain communities
plot(S6.FNN.Louvain.mod, S6.FNN,layout=Lyt.FR[-isolated,], alpha=0.6, 
     vertex.color=V(S6.FNN)$An.color, vertex.frame.color=V(S6.FNN)$An.color, vertex.shape="circle", 
     vertex.size=5, vertex.label=V(S6.FNN)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black", 
     edge.width=E(S6.FNN)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout", 
     rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#Generate HeatMaps of the adjacency matrixes
#Let's use the neighbor list for generating the Adjaency Matrix
SensorID_FNN #Let's use the edge list that contain the total amount of possible interactors for the three sensors
v<-levels(as.factor(unlist(SensorID_FNN[,1:2]))) #get vertex names from the network list
n<-length(v)#get the number of vertex in the network
e<-matrix(match(as.character(unlist(SensorID_FNN[,1:2])),v),ncol = 2)#Edgelist with row index instead of Names
w<-S6_FNN$weight #Edge weights, only from the S6_FNN network to quit those edges that do not exist at S6_FNN

# M is our co-expression network adjacency matrix.

# distringuish between i,j and j,i we want to make M symmetric (undirected) by
# considering an undirected network.
S6.FNN.M <- matrix(0, n, n)                  # set up a co-expression matrix
S6.FNN.M[e] <- w                             # fill it in with edge weights
S6.FNN.M <- S6.FNN.M + t(S6.FNN.M)                         # make this symmetric
any(is.na(S6.FNN.M))#check no NAs in the matrix
dimnames(S6.FNN.M) <- list(v, v)             # label the vertices
# We let A_ij = 1 if there is an interaction. 
#A <- 1*(M > 0) #right now we keep the weights as they are
S6.FNN.M #The adjacency Matrix not as an igraph object, easier to manipulate, and to subset for smaller Sensor specific networks
quantile(S6.FNN.M,c(0,0.25,0.5,0.75,1))
#
#  0%   25%   50%   75%  100% 
#0.000 0.000 0.000 0.000 0.991 

#Creating the HeatMap
# make a heatmap image with colors: white below 0,
# grey from 0.4 to 1 with darker closer to 1
breaks <- c(0, seq(0.4, 0.991, length=11)) # breaks in the color bins
cols <- grey(1-c(0,seq(0.5,1,length=10))) # colors

#Red Blue Colors
library(RColorBrewer)
breaks <- c(seq(-1, 1, length=10))
breaks
cols <- gray(1-c(0,seq(0,1,length=10))) # colors
cols2 <- hcl.colors(10 , palette = "Blues3" , alpha = 1 , rev = TRUE)
cols3 <- colorRamp("red" , space = "rgb" , interpolate = "linear")
HMCol<- brewer.pal(10,name = "RdBu")
HMCol<- rev(HMCol)
HMCol


# color the vertex for CCNB1 blue, its neighbors red, and the others white
SensorIDBind <- which(v == "SIRT6") # index for ccnb1
#SensorIDBind<-as.factor(SensorIDBind)
vcols <- rep("white",n)
vcols[SensorIDBind] <- "blue"
#Subsetting S6.FNN.M by doing a loop asking for look into every row and column the index in SensorIDBind worked
for (i in SensorIDBind) {
  vcols[which(S6.FNN.M[,i] > 0 | S6.FNN.M[i,])] <- "red"
}
##Color annotation for functional annotation
MatSet<-subset(Mat,Mat$NAME%in%v)# select only those proteins that are cointained in the node list of our matrixes
match(v,MatSet$NAME) # Check that the indexes of the Annotation Dataframe and the node list are in the same position 
match(rownames(S6.FNN.M),MatSet$NAME) #Check that the row names of the Adj Matrix and the Annotation dataframe are in the same order
match(colnames(S6.FNN.M),MatSet$NAME) #Check that the column names of the Adj Matrix and the Annotation dataframe are in the same order
AnCol<-MatSet$An.color#Vector of colors
AnCol<-as.character(AnCol)#transform it into a character vector, so the colors can be read

#Subsetting S6.FNN.M using the index vector in SensorIDBind do not work
#vcols[which(S6.FNN.M[,SensorIDBind] > 0 | S6.FNN.M[SensorIDBind,])] <- "red"

#Generate a dendogram for keeping the proteins in the sample place
#S6.FNN.dist<-dist(S6.FNN.M,method = "euclidian")
#S6.FNN.Dend<-hclust(S6.FNN.dist, method = "ward.D" )
#plot(S6.FNN.Dend)

#Cut the dendogram
#cut.height <- 10 # try numbers from 4 to 7
#FNN.hclusts<- cutree(FNN.Dend, cut.height)
#plot(FNN.hclusts)

# now actually make the heat map
par(mar=rep(0,4)) # make plot margins 0
heatmap(S6.FNN.M, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=cols, breaks=breaks,frame=T, cexRow = 0.5 , cexCol = 0.5, main = "Sirt6 nonIR5 AdjMat", ylab="Total Protein Interactors",xlab="Total Protein Interactors")

#
plot.new()
legend("topleft", c("Cell Cycle","Chromatin Org","DNA Repair","Nuclear Pore","p53","Proteasome","Ribosome Biogenesis","RNA Metabolism","Transcription Reg","UPR"), 
       fill=c("coral1","darkviolet","firebrick1","deepskyblue3","tomato","grey","lightgreen","gold","slateblue1","maroon"),  
       bty="n", inset=0, xpd=T,  border=F, cex =2,y.intersp = 0.7,x.intersp = 0.2,pt.cex = 1)

#check with the FNN
# now actually make the heat map
par(mar=rep(0,4)) # make plot margins 0
heatmap(FNN.M, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=0:1, frame=T, cexRow = 0.5 , cexCol = 0.5)
#
legend("topleft", c("Neighbors SIRT6","SIRT6"), fill=c("red","blue"),  
       bty="n", inset=0, xpd=T,  border=F)



################################################Ku80 Networks
###########################################################
#Ku80 Networks


#Let's keep the Mat as the Node attribute table and load the Neighbor list as the matrixes
Mat# Node attribute matrix for all the networks

#Ku80ID nonIR5 nodes Network Edge Lists
Ku_nonIR5<-read.delim("./Ku80_Files/Ku_nonIR5.tsv",sep = "\t",col.names = )
str(Ku_nonIR5)
head(Ku_nonIR5)

#Generate the igraph network object
Ku.nonIR5<-graph_from_data_frame(d=Ku_nonIR5,vertices = Mat, directed = F )#Generate the igraph object with the edge and node annotation attributes
str(vertex_attr(Ku.nonIR5))#Check node attributes
str(edge_attr(Ku.nonIR5))#Check for edge attributes
Ku.nonIR5<-igraph::simplify(Ku.nonIR5, remove.multiple=TRUE, remove.loops=TRUE) #Remove repetitive edges and self-connecting edges
Ku.nonIR5
#Small plot
plot(Ku.nonIR5)

#Generate an adjacency Matrix that will be used after
Ku.nonIR5_AdjMat<-as_adjacency_matrix(Ku.nonIR5,type="both",attr = "weight",names=TRUE,sparse = TRUE)
Ku.nonIR5_AdjMat<-as.matrix(Ku.nonIR5_AdjMat)
write.csv(Ku.nonIR5_AdjMat,"./Ku.nonIR5_AdjMat.csv")

#Generate a vertex size based on degree and count
vsize<-igraph::degree(Ku.nonIR5)/vcount(Ku.nonIR5)
#Identify and remove disconnected nodes
isolated<-which(igraph::degree(Ku.nonIR5)==0)
isolated
length(isolated)# [1] 690
#remove isolated nodes
Ku.nonIR5<-igraph::delete_vertices(Ku.nonIR5,isolated)


#Generate layouts for the graph
#Do not run again, to not to change the layouts that I already have
Lyt.Nice<-layout_nicely(Ku.nonIR5,dim = 3)
Lyt.GOpt<-layout_with_graphopt(Ku.nonIR5, niter=5000, mass = 30, charge = E(Ku.nonIR5)$weight)
Lyt.DrL<-layout_with_drl(Ku.nonIR5, weights = E(Ku.nonIR5)$weight)
Lyt.FR<-layout.fruchterman.reingold(Ku.nonIR5, niter=10000, weights=E(Ku.nonIR5)$weight)#This one looks like the best option
Lyt.KK<-layout_with_kk(Ku.nonIR5, maxiter=10000, weights=E(Ku.nonIR5)$weight) #very large separation of peripheral nodes
Lyt.GEM<-layout_with_gem(Ku.nonIR5, maxiter=10000) #Bad separation
#Lyt.DH<-layout_with_dh(Ku.nonIR5, maxiter=150, fineiter=max(10, log2(833)))
Lyt.LGL<-layout.lgl(Ku.nonIR5, maxit=150,maxdelta=833,area=(length(E(Ku.nonIR5)))^2,coolexp=1.5,repulserad=((length(E(Ku.nonIR5)))^2)*(length(E(Ku.nonIR5))))#Do not work


#Star and Tree layouts

V(Ku.nonIR5)[name=="XRCC5"]
which(V(Ku.nonIR5)[name=="Ku80"])# 610 #Check what is the rowindex for your desired node
which(V(Ku.nonIR5)$name=="XRCC5")# 610 #Check what is the rowindex for your desired node
which(V(Ku.nonIR5)$name=="MRE11A")# 327
which(V(Ku.nonIR5)$name=="XRCC5")# 141
Lyt.Star<-layout_as_star(Ku.nonIR5,center=V(Ku.nonIR5)[name=="XRCC5"]) #it has a special order argument where we can assign an order to the vertex, hence maybe the OPLS-DA:VIP value order may be used as ordering vector

Lyt.Tree<-layout_as_tree(Ku.nonIR5,root=c(141))
Lyt.RT<-layout.reingold.tilford(Ku.nonIR5,root=c(141))# This one looks good, but the indexing do not allow to highlight the correct nodes
Lyt.sugi<-layout_with_sugiyama(Ku.nonIR5, maxiter=100, weights=E(Ku.nonIR52)$weight)#This one looks good for undirected tree layouts

#Check what layout looks fine

plot.igraph(Ku.nonIR5,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(Ku.nonIR5)$An.color, vertex.frame.color=V(Ku.nonIR5)$An.color, vertex.shape="circle",
            vertex.size=5, vertex.label=V(Ku.nonIR5)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(Ku.nonIR5)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)

#With ggraph

ggraph(Ku.nonIR5,layout=Lyt.FR[-isolated,])+
  geom_edge_fan(color="gray65", width=1, alpha=0.4)+
  geom_node_point(aes(color=V(Ku.nonIR5)$Annotation, size=V(Ku.nonIR5)$Ku.FC.nonIR5),alpha=0.75)+
  geom_node_text(aes(label=name),color="gray25",size=4)+
  scale_color_manual(limits=V(Ku.nonIR5)$Annotation, values=V(Ku.nonIR5)$An.color)+
  scale_size_continuous(range = c(5, 15)) +
  theme_void()+
  title("Ku80 nonIR5 Network")+
  guides(color=guide_legend(override.aes = list(size=7)))+
  theme(legend.text=element_text(size=12,face = "bold"),
        legend.title = element_text(size=13,face = "bold"))



#Community Detection
#clp <- igraph::cluster_optimal(Ku.nonIR5)
class(clp)
V(Ku.nonIR5)$community <- clp$membership
colrs <- adjustcolor( c("gray50", "tomato", "gold", "yellowgreen"), alpha=.6)


plot.igraph(clp,Ku.nonIR5,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(Ku.nonIR5)$An.color, vertex.frame.color=V(Ku.nonIR5)$An.color, vertex.shape="circle",
            vertex.size==V(Ku.nonIR5)$Ku.FC.nonIR5, vertex.label=V(Ku.nonIR5)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(Ku.nonIR5)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#Distances edges highlight
dist.from.Ku80 <- distances(Ku.nonIR5, v=V(Ku.nonIR5)[name=="XRCC5"],
                            to=V(Ku.nonIR5), weights=E(Ku.nonIR5)$weight)
dist.from.Ku80<-dist.from.Ku80[,!is.infinite(colSums(dist.from.Ku80))]

dist.from.Ku80
range(dist.from.Ku80)

#Set colors to plot the distances:
oranges <- colorRampPalette(c("dark red", "gold"))
col <- oranges(max(dist.from.Ku80))
col <- col[dist.from.Ku80+1]
plot(Ku.nonIR5, layout=Lyt.FR[-isolated,],vertex.color=col, vertex.label=dist.from.Ku80,
     vertex.label.color="white",vertex.shape="circle",
     vertex.size=V(Ku.nonIR5)$Ku.FC.nonIR5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(Ku.nonIR5)$weight, edge.lty=1, main="Ku80 nonIR5 Full Nuclear Network", sub="FR Layout",
     rescale=T)


#Shortest path from Ku80 to Nucleolin
KuNcl.path <- shortest_paths(Ku.nonIR5,
                             from = V(Ku.nonIR5)[name=="XRCC5"],
                             to = V(Ku.nonIR5)[name=="NCL"],
                             output = "both") # both path nodes and edges
# Generate edge color variable to plot the path:
ecol <- rep("gray80", ecount(Ku.nonIR5))
ecol[unlist(KuNcl.path$epath)] <- "orange"
# Generate edge width variable to plot the path:
ew <- rep(2, ecount(Ku.nonIR5))
ew[unlist(KuNcl.path$epath)] <- 4
# Generate node color variable to plot the path:
vcol <- rep("gray40", vcount(Ku.nonIR5))
vcol[unlist(KuNcl.path$vpath)] <- "gold"
plot(Ku.nonIR5, layout=Lyt.FR[-isolated,], vertex.color=vcol, edge.color=ecol,
     edge.width=ew,vertex.label.color="white",vertex.shape="circle",
     vertex.size=5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(Ku.nonIR5)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)



#Plot based on Ku80 as incident
inc.edges <- incident(Ku.nonIR5, V(Ku.nonIR5)[name=="XRCC5"], mode="all")
# Set colors to plot the selected edges.
ecol <- rep("gray80", ecount(Ku.nonIR5))
ecol[inc.edges] <- "orange"
vcol <- rep("gray40", vcount(Ku.nonIR5))
vcol[V(Ku.nonIR5)$name=="Ku80"] <- "gold"

plot(Ku.nonIR5, layout=Lyt.FR[-isolated,], vertex.color=vcol, edge.color=ecol, 
     edge.width=ew,vertex.label.color="white",vertex.shape="circle",
     vertex.size=V(Ku.nonIR5)$Ku.FC.nonIR5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(Ku.nonIR5)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)


#Plot Degree Distribution
Ku.nonIR5.Deg<-degree_distribution(Ku.nonIR5, cumulative = T, mode="all")
plot( x=0:max(igraph::degree(Ku.nonIR5)), y=1-Ku.nonIR5.Deg, pch=19, cex=1.2, col="orange",
      xlab="Degree", ylab="Cumulative Frequency")

#Network Descriptives
#Number of nodes and edges
vcount(Ku.nonIR5) #[1] 143
ecount(Ku.nonIR5)#[1] 895
#Density (Number of edges in relation to the total number of edges in a full connected graph)
edge_density(Ku.nonIR5,loops = F) #[1] 0.08815129
#Transitivity (Total number of conected triangles)
transitivity(Ku.nonIR5,type="global") #[1] 0.6377146
#modularity
mod<-V(Ku.nonIR5)$S6.FuzzyCMeans #lets check for fuzzy Cmeans modularity
mod<-mod+1#The non clustered interactors that appears as 0, needs to be transformed to 1, so we add 1 to every cluster
modularity(Ku.nonIR5,membership = mod)
#Ku.FuzzyCMeans modularity [1] -0.001878843
#Ku.FuzzyCMeans modularity [1] 0.02609219
#Mre.FuzzyCMeans modularity [1]-0.001919545

#Assorativity
assortativity_nominal(Ku.nonIR5,factor(V(Ku.nonIR5)$Annotation),directed=F)
#Annotation [1] [1] 0.1605138
#Ku.FuzzyCMeans [1] -0.003066766
#Ku.FuzzyCMeans [1]  0.04549021
#Mre.FuzzyCMeans [1] -0.001919545

#Diameter (longest geodesic distance), shortest path between two most peripheral nodes
diameter(Ku.nonIR5,directed = F, weights = E(Ku.nonIR5)$weights)#4.506
diam<-get_diameter(Ku.nonIR5,directed = F,weights = E(Ku.nonIR5)$weights)
diam #+ 9/143 vertices, named, from c2e359b: [1] ACTL6A    HIST1H2AC WDR5      DDB1      XRCC5     HNRNPA1   RALY      RANBP2    RAN  
#Plot Diameter
vcol <- rep("gray50", vcount(Ku.nonIR5))
vcol[diam] <- "gold"
ecol <- rep("gray90", ecount(Ku.nonIR5))
ecol[E(Ku.nonIR5, path=diam)] <- "orange"
# E(net, path=diam) finds edges along a path, here 'diam'
plot(Ku.nonIR5, vertex.color=vcol, edge.color=ecol, vertex.label.color="black",edge.arrow.mode=0,layout=Lyt.FR[-isolated,],vertex.size=V(Ku.nonIR5)$Ku.FC.nonIR5)

#Node topological descriptives
degree<-igraph::centr_degree(Ku.nonIR5,mode = "all",normalized = T)
degree<-degree$res
eigen_cent<-igraph::eigen_centrality(Ku.nonIR5)
eigen_cent<-eigen_cent$vector
clos_cent <- igraph::closeness(Ku.nonIR5, normalized=T)
betw_cent <- igraph::betweenness(Ku.nonIR5, normalized=T)
hs <- hub_score(Ku.nonIR5, weights=E(Ku.nonIR5)$weight)$vector
as <- authority_score(Ku.nonIR5, weights=E(Ku.nonIR5)$weight)$vector
Trans<-transitivity(Ku.nonIR5,type="local")
Trans
#write results
Ku.nonIR5.Centrality<-data.frame(NAME=V(Ku.nonIR5)$name, UniprotID=V(Ku.nonIR5)$Uniprot, Degree=degree, Eigen= eigen_cent, Closeness= clos_cent, Betweeness= betw_cent, Hubs=hs, Authority=as,Transitivity=Trans)
Ku.nonIR5.Centrality
write.csv(Ku.nonIR5.Centrality,"./Ku.nonIR5.Centrality.csv")

#Cocitation, Number of shared interactions among nodes
Cocitation<-cocitation(Ku.nonIR5)
Cocitation#This matrix is similar to an adjacency matrix
write.csv(Cocitation,"./Ku.nonIR5.cocitation.csv")
#Plot with degree as node size
ggraph(Ku.nonIR5,layout=Lyt.FR[-isolated,])+
  geom_edge_fan(color="gray85", width=1, alpha=0.4)+
  geom_node_point(aes(color=V(Ku.nonIR5)$Annotation, size=degree),alpha=0.75)+
  geom_node_text(aes(label=name),color="gray25",size=4)+
  scale_color_manual(limits=V(Ku.nonIR5)$Annotation, values=V(Ku.nonIR5)$An.color)+
  scale_size_continuous(range = c(5, 15)) +
  theme_void()+
  title("Ku80 nonIR5 Network")+
  guides(color=guide_legend(override.aes = list(size=7)))+
  theme(legend.text=element_text(size=12,face = "bold"),
        legend.title = element_text(size=13,face = "bold"))


#Modularity Analysis

#Edge_Betweeness Modularity
Ku.nonIR5.ebmod=edge.betweenness.community(Ku.nonIR5,weights=E(Ku.nonIR5)$weight,directed = F) 
Ku.nonIR5.ebmod# IGRAPH clustering edge betweenness, groups: 13, mod: 0.39
#
plot.igraph(Ku.nonIR5.ebmod,Ku.nonIR5,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(Ku.nonIR5)$An.color, vertex.frame.color=V(Ku.nonIR5)$An.color, vertex.shape="circle",
            vertex.size=5, vertex.label=V(Ku.nonIR5)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(Ku.nonIR5)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)

#
plot(Ku.nonIR5.ebmod, Ku.nonIR5,layout=Lyt.FR[-isolated,], alpha=0.6, 
     vertex.color=V(Ku.nonIR5)$An.color, vertex.frame.color=V(Ku.nonIR5)$An.color, vertex.shape="circle", 
     vertex.size=5, vertex.label=V(Ku.nonIR5)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black", 
     edge.width=E(Ku.nonIR5)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout", 
     rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)
#Dendogram of the Edge_Betweeness
dendPlot(Ku.nonIR5.ebmod, mode="hclust")

#Louvain modularity
Ku.nonIR5.Louvain.mod<-cluster_louvain(Ku.nonIR5,weights=E(Ku.nonIR5)$weight)
Ku.nonIR5.Louvain.mod #IGRAPH clustering multi level, groups: 7, mod: 0.42

#Plot the Louvain communities
plot(Ku.nonIR5.Louvain.mod, Ku.nonIR5,layout=Lyt.FR[-isolated,], alpha=0.6, 
     vertex.color=V(Ku.nonIR5)$An.color, vertex.frame.color=V(Ku.nonIR5)$An.color, vertex.shape="circle", 
     vertex.size=5, vertex.label=V(Ku.nonIR5)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black", 
     edge.width=E(Ku.nonIR5)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout", 
     rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#Generate HeatMaps of the adjacency matrixes
#Let's use the neighbor list for generating the Adjaency Matrix
SensorID_FNN #Let's use the edge list that contain the total amount of possible interactors for the three sensors
v<-levels(as.factor(unlist(SensorID_FNN[,1:2]))) #get vertex names from the network list
n<-length(v)#get the number of vertex in the network
e<-matrix(match(as.character(unlist(SensorID_FNN[,1:2])),v),ncol = 2)#Edgelist with row index instead of Names
w<-Ku_nonIR5$weight #Edge weights, only from the Ku_nonIR5 network to quit those edges that do not exist at Ku_nonIR5

# M is our co-expression network adjacency matrix.

# distringuish between i,j and j,i we want to make M symmetric (undirected) by
# considering an undirected network.
Ku.nonIR5.M <- matrix(0, n, n)                  # set up a co-expression matrix
Ku.nonIR5.M[e] <- w                             # fill it in with edge weights
Ku.nonIR5.M <- Ku.nonIR5.M + t(Ku.nonIR5.M)                         # make this symmetric
any(is.na(Ku.nonIR5.M))#check no NAs in the matrix
dimnames(Ku.nonIR5.M) <- list(v, v)             # label the vertices
# We let A_ij = 1 if there is an interaction. 
#A <- 1*(M > 0) #right now we keep the weights as they are
Ku.nonIR5.M #The adjacency Matrix not as an igraph object, easier to manipulate, and to subset for smaller Sensor specific networks
quantile(Ku.nonIR5.M,c(0,0.25,0.5,0.75,1))
#
#  0%   25%   50%   75%  100% 
#0.000 0.000 0.000 0.000 0.991 

#Creating the HeatMap
# make a heatmap image with colors: white below 0,
# grey from 0.4 to 1 with darker closer to 1
breaks <- c(0, seq(0.4, 0.991, length=11)) # breaks in the color bins
cols <- grey(1-c(0,seq(0.5,1,length=10))) # colors

#Red Blue Colors
library(RColorBrewer)
breaks <- c(seq(-1, 1, length=10))
breaks
cols <- gray(1-c(0,seq(0,1,length=10))) # colors
cols2 <- hcl.colors(10 , palette = "Blues3" , alpha = 1 , rev = TRUE)
cols3 <- colorRamp("red" , space = "rgb" , interpolate = "linear")
HMCol<- brewer.pal(10,name = "RdBu")
HMCol<- rev(HMCol)
HMCol


# color the vertex for CCNB1 blue, its neighbors red, and the others white
SensorIDBind <- which(v == "Ku80") # index for ccnb1
#SensorIDBind<-as.factor(SensorIDBind)
vcols <- rep("white",n)
vcols[SensorIDBind] <- "blue"
#Subsetting Ku.nonIR5.M by doing a loop asking for look into every row and column the index in SensorIDBind worked
for (i in SensorIDBind) {
  vcols[which(Ku.nonIR5.M[,i] > 0 | Ku.nonIR5.M[i,])] <- "red"
}
##Color annotation for functional annotation
MatSet<-subset(Mat,Mat$NAME%in%v)# select only those proteins that are cointained in the node list of our matrixes
match(v,MatSet$NAME) # Check that the indexes of the Annotation Dataframe and the node list are in the same position 
match(rownames(Ku.nonIR5.M),MatSet$NAME) #Check that the row names of the Adj Matrix and the Annotation dataframe are in the same order
match(colnames(Ku.nonIR5.M),MatSet$NAME) #Check that the column names of the Adj Matrix and the Annotation dataframe are in the same order
AnCol<-MatSet$An.color#Vector of colors
AnCol<-as.character(AnCol)#transform it into a character vector, so the colors can be read

#Subsetting Ku.nonIR5.M using the index vector in SensorIDBind do not work
#vcols[which(Ku.nonIR5.M[,SensorIDBind] > 0 | Ku.nonIR5.M[SensorIDBind,])] <- "red"

#Generate a dendogram for keeping the proteins in the sample place
#Ku.nonIR5.dist<-dist(Ku.nonIR5.M,method = "euclidian")
#Ku.nonIR5.Dend<-hclust(Ku.nonIR5.dist, method = "ward.D" )
#plot(Ku.nonIR5.Dend)

#Cut the dendogram
#cut.height <- 10 # try numbers from 4 to 7
#FNN.hclusts<- cutree(FNN.Dend, cut.height)
#plot(FNN.hclusts)

# now actually make the heat map
par(mar=rep(0,4)) # make plot margins 0
heatmap(Ku.nonIR5.M, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=cols, breaks=breaks,frame=T, cexRow = 0.5 , cexCol = 0.5, main = "Ku80 nonIR5 AdjMat", ylab="Total Protein Interactors",xlab="Total Protein Interactors")

#
plot.new()
legend("topleft", c("Cell Cycle","Chromatin Org","DNA Repair","Nuclear Pore","p53","Proteasome","Ribosome Biogenesis","RNA Metabolism","Transcription Reg","UPR"), 
       fill=c("coral1","darkviolet","firebrick1","deepskyblue3","tomato","grey","lightgreen","gold","slateblue1","maroon"),  
       bty="n", inset=0, xpd=T,  border=F, cex =2,y.intersp = 0.7,x.intersp = 0.2,pt.cex = 1)

#check with the FNN
# now actually make the heat map
par(mar=rep(0,4)) # make plot margins 0
heatmap(FNN.M, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=0:1, frame=T, cexRow = 0.5 , cexCol = 0.5)
#
legend("topleft", c("Neighbors Ku80","Ku80"), fill=c("red","blue"),  
       bty="n", inset=0, xpd=T,  border=F)


############################################Ku80 IR5 Networks
##############################################
#Ku80 IR5 Networks


#Let's keep the Mat as the Node attribute table and load the Neighbor list as the matrixes
Mat# Node attribute matrix for all the networks

#Ku80ID IR5 nodes Network Edge Lists
Ku_IR5<-read.delim("./Ku80_Files/Ku_IR5.tsv",sep = "\t",col.names = )
str(Ku_IR5)
head(Ku_IR5)

#Generate the igraph network object
Ku.IR5<-graph_from_data_frame(d=Ku_IR5,vertices = Mat, directed = F )#Generate the igraph object with the edge and node annotation attributes
str(vertex_attr(Ku.IR5))#Check node attributes
str(edge_attr(Ku.IR5))#Check for edge attributes
Ku.IR5<-igraph::simplify(Ku.IR5, remove.multiple=TRUE, remove.loops=TRUE) #Remove repetitive edges and self-connecting edges
Ku.IR5
#Small plot
plot(Ku.IR5)

#Generate an adjacency Matrix that will be used after
Ku.IR5_AdjMat<-as_adjacency_matrix(Ku.IR5,type="both",attr = "weight",names=TRUE,sparse = TRUE)
Ku.IR5_AdjMat<-as.matrix(Ku.IR5_AdjMat)
write.csv(Ku.IR5_AdjMat,"./Ku.IR5_AdjMat.csv")

#Generate a vertex size based on degree and count
vsize<-igraph::degree(Ku.IR5)/vcount(Ku.IR5)
#Identify and remove disconnected nodes
isolated<-which(igraph::degree(Ku.IR5)==0)
isolated
length(isolated)# [1] 687
#remove isolated nodes
Ku.IR5<-igraph::delete_vertices(Ku.IR5,isolated)


#Generate layouts for the graph
#Do not run again, to not to change the layouts that I already have
Lyt.Nice<-layout_nicely(Ku.IR5,dim = 3)
Lyt.GOpt<-layout_with_graphopt(Ku.IR5, niter=5000, mass = 30, charge = E(Ku.IR5)$weight)
Lyt.DrL<-layout_with_drl(Ku.IR5, weights = E(Ku.IR5)$weight)
Lyt.FR<-layout.fruchterman.reingold(Ku.IR5, niter=10000, weights=E(Ku.IR5)$weight)#This one looks like the best option
Lyt.KK<-layout_with_kk(Ku.IR5, maxiter=10000, weights=E(Ku.IR5)$weight) #very large separation of peripheral nodes
Lyt.GEM<-layout_with_gem(Ku.IR5, maxiter=10000) #Bad separation
#Lyt.DH<-layout_with_dh(Ku.IR5, maxiter=150, fineiter=max(10, log2(833)))
Lyt.LGL<-layout.lgl(Ku.IR5, maxit=150,maxdelta=833,area=(length(E(Ku.IR5)))^2,coolexp=1.5,repulserad=((length(E(Ku.IR5)))^2)*(length(E(Ku.IR5))))#Do not work


#Star and Tree layouts

V(Ku.IR5)[name=="XRCC5"]
which(V(Ku.IR5)[name=="Ku80"])# 228 #Check what is the rowindex for your desired node
which(V(Ku.IR5)$name=="Ku80")# 191 #Check what is the rowindex for your desired node
which(V(Ku.IR5)$name=="MRE11A")# 327
which(V(Ku.IR5)$name=="XRCC5")# 143
Lyt.Star<-layout_as_star(Ku.IR5,center=V(Ku.IR5)[name=="XRCC5"]) #it has a special order argument where we can assign an order to the vertex, hence maybe the OPLS-DA:VIP value order may be used as ordering vector

Lyt.Tree<-layout_as_tree(Ku.IR5,root=c(143))
Lyt.RT<-layout.reingold.tilford(Ku.IR5,root=c(143))# This one looks good, but the indexing do not allow to highlight the correct nodes
Lyt.sugi<-layout_with_sugiyama(Ku.IR5, maxiter=100, weights=E(Ku.IR52)$weight)#This one looks good for undirected tree layouts

#Check what layout looks fine

plot.igraph(Ku.IR5,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(Ku.IR5)$An.color, vertex.frame.color=V(Ku.IR5)$An.color, vertex.shape="circle",
            vertex.size=V(Ku.IR5)$Ku.FC.IR5, vertex.label=V(Ku.IR5)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(Ku.IR5)$weight, edge.lty=1, main="Ku80ID IR5 Network", sub="FR Tree Layout",
            rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#With ggraph

ggraph(Ku.IR5,layout=Lyt.FR[-isolated,])+
  geom_edge_fan(color="gray65", width=1, alpha=0.4)+
  geom_node_point(aes(color=V(Ku.IR5)$Annotation, size=V(Ku.IR5)$Ku.FC.IR5),alpha=0.75)+
  geom_node_text(aes(label=name),color="gray25",size=4)+
  scale_color_manual(limits=V(Ku.IR5)$Annotation, values=V(Ku.IR5)$An.color)+
  scale_size_continuous(range = c(5, 15)) +
  theme_void()+
  title("Ku80 IR5 Network")+
  guides(color=guide_legend(override.aes = list(size=7)))+
  theme(legend.text=element_text(size=12,face = "bold"),
        legend.title = element_text(size=13,face = "bold"))



#Community Detection
#clp <- igraph::cluster_optimal(Ku.IR5)
class(clp)
V(Ku.IR5)$community <- clp$membership
colrs <- adjustcolor( c("gray50", "tomato", "gold", "yellowgreen"), alpha=.6)


plot.igraph(clp,Ku.IR5,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(Ku.IR5)$An.color, vertex.frame.color=V(Ku.IR5)$An.color, vertex.shape="circle",
            vertex.size==V(Ku.IR5)$Ku.FC.IR5, vertex.label=V(Ku.IR5)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(Ku.IR5)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#Distances edges highlight
dist.from.Ku80 <- distances(Ku.IR5, v=V(Ku.IR5)[name=="XRCC5"],
                            to=V(Ku.IR5), weights=E(Ku.IR5)$weight)
dist.from.Ku80<-dist.from.Ku80[,!is.infinite(colSums(dist.from.Ku80))]

dist.from.Ku80
range(dist.from.Ku80)

#Set colors to plot the distances:
oranges <- colorRampPalette(c("dark red", "gold"))
col <- oranges(max(dist.from.Ku80))
col <- col[dist.from.Ku80+1]
plot(Ku.IR5, layout=Lyt.FR[-isolated,],vertex.color=col, vertex.label=dist.from.Ku80,
     vertex.label.color="white",vertex.shape="circle",
     vertex.size=V(Ku.IR5)$Ku.FC.IR5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(Ku.IR5)$weight, edge.lty=1, main="Ku80 IR5 Full Nuclear Network", sub="FR Layout",
     rescale=T)


#Shortest path from Ku80 to Nucleolin
KuNcl.path <- shortest_paths(Ku.IR5,
                             from = V(Ku.IR5)[name=="XRCC5"],
                             to = V(Ku.IR5)[name=="NCL"],
                             output = "both") # both path nodes and edges
# Generate edge color variable to plot the path:
ecol <- rep("gray80", ecount(Ku.IR5))
ecol[unlist(KuNcl.path$epath)] <- "orange"
# Generate edge width variable to plot the path:
ew <- rep(2, ecount(Ku.IR5))
ew[unlist(KuNcl.path$epath)] <- 4
# Generate node color variable to plot the path:
vcol <- rep("gray40", vcount(Ku.IR5))
vcol[unlist(KuNcl.path$vpath)] <- "gold"
plot(Ku.IR5, layout=Lyt.FR[-isolated,], vertex.color=vcol, edge.color=ecol,
     edge.width=ew,vertex.label.color="white",vertex.shape="circle",
     vertex.size=5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(Ku.IR5)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)



#Plot based on Ku80 as incident
inc.edges <- incident(Ku.IR5, V(Ku.IR5)[name=="XRCC5"], mode="all")
# Set colors to plot the selected edges.
ecol <- rep("gray80", ecount(Ku.IR5))
ecol[inc.edges] <- "orange"
vcol <- rep("gray40", vcount(Ku.IR5))
vcol[V(Ku.IR5)$name=="Ku80"] <- "gold"

plot(Ku.IR5, layout=Lyt.FR[-isolated,], vertex.color=vcol, edge.color=ecol, 
     edge.width=ew,vertex.label.color="white",vertex.shape="circle",
     vertex.size=V(Ku.IR5)$Ku.FC.IR5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(Ku.IR5)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)


#Plot Degree Distribution
Ku.IR5.Deg<-degree_distribution(Ku.IR5, cumulative = T, mode="all")
plot( x=0:max(igraph::degree(Ku.IR5)), y=1-Ku.IR5.Deg, pch=19, cex=1.2, col="orange",
      xlab="Degree", ylab="Cumulative Frequency")

#Network Descriptives
#Number of nodes and edges
vcount(Ku.IR5) #[1] 146
ecount(Ku.IR5)#[1]  712
#Density (Number of edges in relation to the total number of edges in a full connected graph)
edge_density(Ku.IR5,loops = F) #[1] 0.067265
#Transitivity (Total number of conected triangles)
transitivity(Ku.IR5,type="global") #[1] 0.547224
#modularity
mod<-V(Ku.IR5)$Mre.FuzzyCMeans #lets check for fuzzy Cmeans modularity
mod<-mod+1#The non clustered interactors that appears as 0, needs to be transformed to 1, so we add 1 to every cluster
modularity(Ku.IR5,membership = mod)
#Ku.FuzzyCMeans modularity [1] 0.01465053
#Ku.FuzzyCMeans modularity [1] 0.01458051
#Mre.FuzzyCMeans modularity [1] -0.001212165

#Assorativity
assortativity_nominal(Ku.IR5,factor(V(Ku.IR5)$Mre.FuzzyCMeans),directed=F)
#Annotation [1] [1] 0.1927308
#Ku.FuzzyCMeans [1] 0.02541517
#Ku.FuzzyCMeans [1] 0.02573575
#Mre.FuzzyCMeans [1] -0.001538059

#Diameter (longest geodesic distance), shortest path between two most peripheral nodes
diameter(Ku.IR5,directed = F, weights = E(Ku.IR5)$weights)#[1] 3.732
diam<-get_diameter(Ku.IR5,directed = F,weights = E(Ku.IR5)$weights)
diam #+ 7/146 vertices, named, from 3e48a7b: [1] CFL1   ACTG1  RUVBL2 XRCC5  TMPO   RPN1   RPN2  
#Plot Diameter
vcol <- rep("gray50", vcount(Ku.IR5))
vcol[diam] <- "gold"
ecol <- rep("gray90", ecount(Ku.IR5))
ecol[E(Ku.IR5, path=diam)] <- "orange"
# E(net, path=diam) finds edges along a path, here 'diam'
plot(Ku.IR5, vertex.color=vcol, edge.color=ecol, vertex.label.color="black",edge.arrow.mode=0,layout=Lyt.FR[-isolated,],vertex.size=V(Ku.IR5)$Ku.FC.IR5)

#Node topological descriptives
degree<-igraph::centr_degree(Ku.IR5,mode = "all",normalized = T)
degree<-degree$res
eigen_cent<-igraph::eigen_centrality(Ku.IR5)
eigen_cent<-eigen_cent$vector
clos_cent <- igraph::closeness(Ku.IR5, normalized=T)
betw_cent <- igraph::betweenness(Ku.IR5, normalized=T)
hs <- hub_score(Ku.IR5, weights=E(Ku.IR5)$weight)$vector
as <- authority_score(Ku.IR5, weights=E(Ku.IR5)$weight)$vector
Trans<-transitivity(Ku.IR5,type="local")
Trans
#write results
Ku.IR5.Centrality<-data.frame(NAME=V(Ku.IR5)$name, UniprotID=V(Ku.IR5)$Uniprot, Degree=degree, Eigen= eigen_cent, Closeness= clos_cent, Betweeness= betw_cent, Hubs=hs, Authority=as,Transitivity=Trans)
Ku.IR5.Centrality
write.csv(Ku.IR5.Centrality,"./Ku.IR5.Centrality.csv")

#Cocitation, Number of shared interactions among nodes
Cocitation<-cocitation(Ku.IR5)
Cocitation#This matrix is similar to an adjacency matrix
write.csv(Cocitation,"./Ku.IR5.cocitation.csv")
#Plot with degree as node size
ggraph(Ku.IR5,layout=Lyt.FR[-isolated,])+
  geom_edge_fan(color="gray85", width=1, alpha=0.4)+
  geom_node_point(aes(color=V(Ku.IR5)$Annotation, size=degree),alpha=0.75)+
  geom_node_text(aes(label=name),color="gray25",size=4)+
  scale_color_manual(limits=V(Ku.IR5)$Annotation, values=V(Ku.IR5)$An.color)+
  scale_size_continuous(range = c(5, 15)) +
  theme_void()+
  title("Ku80 IR5 Network")+
  guides(color=guide_legend(override.aes = list(size=7)))+
  theme(legend.text=element_text(size=12,face = "bold"),
        legend.title = element_text(size=13,face = "bold"))


#Modularity Analysis

#Edge_Betweeness Modularity
Ku.IR5.ebmod=edge.betweenness.community(Ku.IR5,weights=E(Ku.IR5)$weight,directed = F) 
Ku.IR5.ebmod# IGRAPH clustering edge betweenness, groups: 13, mod: 0.47

#
plot.igraph(Ku.IR5.ebmod,Ku.IR5,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(Ku.IR5)$An.color, vertex.frame.color=V(Ku.IR5)$An.color, vertex.shape="circle",
            vertex.size=5, vertex.label=V(Ku.IR5)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(Ku.IR5)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)

#
plot(Ku.IR5.ebmod, Ku.IR5,layout=Lyt.FR[-isolated,], alpha=0.6, 
     vertex.color=V(Ku.IR5)$An.color, vertex.frame.color=V(Ku.IR5)$An.color, vertex.shape="circle", 
     vertex.size=5, vertex.label=V(Ku.IR5)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black", 
     edge.width=E(Ku.IR5)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout", 
     rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)
#Dendogram of the Edge_Betweeness
dendPlot(Ku.IR5.ebmod, mode="hclust")

#Louvain modularity
Ku.IR5.Louvain.mod<-cluster_louvain(Ku.IR5,weights=E(Ku.IR5)$weight)
Ku.IR5.Louvain.mod #IGRAPH clustering multi level, groups: 7, mod: 0.51

#Plot the Louvain communities
plot(Ku.IR5.Louvain.mod, Ku.IR5,layout=Lyt.FR[-isolated,], alpha=0.6, 
     vertex.color=V(Ku.IR5)$An.color, vertex.frame.color=V(Ku.IR5)$An.color, vertex.shape="circle", 
     vertex.size=5, vertex.label=V(Ku.IR5)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black", 
     edge.width=E(Ku.IR5)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout", 
     rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#Generate HeatMaps of the adjacency matrixes
#Let's use the neighbor list for generating the Adjaency Matrix
SensorID_FNN #Let's use the edge list that contain the total amount of possible interactors for the three sensors
v<-levels(as.factor(unlist(SensorID_FNN[,1:2]))) #get vertex names from the network list
n<-length(v)#get the number of vertex in the network
e<-matrix(match(as.character(unlist(SensorID_FNN[,1:2])),v),ncol = 2)#Edgelist with row index instead of Names
w<-Ku_IR5$weight #Edge weights, only from the Ku_IR5 network to quit those edges that do not exist at Ku_IR5

# M is our co-expression network adjacency matrix.

# distringuish between i,j and j,i we want to make M symmetric (undirected) by
# considering an undirected network.
Ku.IR5.M <- matrix(0, n, n)                  # set up a co-expression matrix
Ku.IR5.M[e] <- w                             # fill it in with edge weights
Ku.IR5.M <- Ku.IR5.M + t(Ku.IR5.M)                         # make this symmetric
any(is.na(Ku.IR5.M))#check no NAs in the matrix
dimnames(Ku.IR5.M) <- list(v, v)             # label the vertices
# We let A_ij = 1 if there is an interaction. 
#A <- 1*(M > 0) #right now we keep the weights as they are
Ku.IR5.M #The adjacency Matrix not as an igraph object, easier to manipulate, and to subset for smaller Sensor specific networks
quantile(Ku.IR5.M,c(0,0.25,0.5,0.75,1))
#
# 0%   25%   50%   75%  100% 
#0.000 0.000 0.000 0.000 0.991  

#Creating the HeatMap
# make a heatmap image with colors: white below 0,
# grey from 0.4 to 1 with darker closer to 1
breaks <- c(0, seq(0.4, 0.991, length=11)) # breaks in the color bins
cols <- grey(1-c(0,seq(0.5,1,length=10))) # colors

#Red Blue Colors
library(RColorBrewer)
breaks <- c(seq(-1, 1, length=11))
breaks
cols <- gray(1-c(0,seq(0,1,length=10))) # colors
cols2 <- hcl.colors(10 , palette = "Blues3" , alpha = 1 , rev = TRUE)
cols3 <- colorRamp("red" , space = "rgb" , interpolate = "linear")
HMCol<- brewer.pal(10,name = "RdBu")
HMCol<- rev(HMCol)
HMCol


# color the vertex for CCNB1 blue, its neighbors red, and the others white
SensorIDBind <- which(v == "Ku80") # index for ccnb1
#SensorIDBind<-as.factor(SensorIDBind)
vcols <- rep("white",n)
vcols[SensorIDBind] <- "blue"
#Subsetting Ku.IR5.M by doing a loop asking for look into every row and column the index in SensorIDBind worked
for (i in SensorIDBind) {
  vcols[which(Ku.IR5.M[,i] > 0 | Ku.IR5.M[i,])] <- "red"
}
##Color annotation for functional annotation
MatSet<-subset(Mat,Mat$NAME%in%v)# select only those proteins that are cointained in the node list of our matrixes
match(v,MatSet$NAME) # Check that the indexes of the Annotation Dataframe and the node list are in the same position 
match(rownames(Ku.IR5.M),MatSet$NAME) #Check that the row names of the Adj Matrix and the Annotation dataframe are in the same order
match(colnames(Ku.IR5.M),MatSet$NAME) #Check that the column names of the Adj Matrix and the Annotation dataframe are in the same order
AnCol<-MatSet$An.color#Vector of colors
AnCol<-as.character(AnCol)#transform it into a character vector, so the colors can be read

#Subsetting Ku.IR5.M using the index vector in SensorIDBind do not work
#vcols[which(Ku.IR5.M[,SensorIDBind] > 0 | Ku.IR5.M[SensorIDBind,])] <- "red"

#Generate a dendogram for keeping the proteins in the sample place
#Ku.IR5.dist<-dist(Ku.IR5.M,method = "euclidian")
#Ku.IR5.Dend<-hclust(Ku.IR5.dist, method = "ward.D" )
#plot(Ku.IR5.Dend)

#Cut the dendogram
#cut.height <- 10 # try numbers from 4 to 7
#FNN.hclusts<- cutree(FNN.Dend, cut.height)
#plot(FNN.hclusts)

# now actually make the heat map
par(mar=rep(0,4)) # make plot margins 0
heatmap(Ku.IR5.M, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=cols, breaks=breaks,frame=T, cexRow = 0.5 , cexCol = 0.5, main = "Ku80 IR5 AdjMat", ylab="Total Protein Interactors",xlab="Total Protein Interactors")
#
plot.new()
legend("topleft", c("Cell Cycle","Chromatin Org","DNA Repair","Nuclear Pore","p53","Proteasome","Ribosome Biogenesis","RNA Metabolism","Transcription Reg","UPR"), 
       fill=c("coral1","darkviolet","firebrick1","deepskyblue3","tomato","grey","lightgreen","gold","slateblue1","maroon"),  
       bty="n", inset=0, xpd=T,  border=F, cex =2,y.intersp = 0.7,x.intersp = 0.2,pt.cex = 1)


#Now we substract the KunonIR5 AdjMatrix to the KuIR6 AdjMatrix
KuIR5.AdjMat.Change<-Ku.IR5.M-Ku.nonIR5.M
range(KuIR5.AdjMat.Change)#Let's check the magnitude of the substracted weights
#[1] -0.807  0.805
any(is.na(KuIR5.AdjMat.Change))
#Generate the breaks and color functions
breaks <- c(seq(-0.5, 0.5, length=10))
breaks
HMCol<- brewer.pal(9,name = "RdBu")
HMCol<- rev(HMCol)
HMCol
[1] "#2166AC" "#4393C3" "#92C5DE" "#D1E5F0" "#F7F7F7" "#FDDBC7" "#F4A582" "#D6604D" "#B2182B"

#Plot the HeatMap
par(mar=rep(0,4)) # make plot margins 0
heatmap(KuIR5.AdjMat.Change, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=HMCol, breaks=breaks,frame=T, cexRow = 0.5 , cexCol = 0.5, main = "Ku80 IR5 AdjMat (minus nonIR5 edges)", ylab="Total Protein Interactors",xlab="Total Protein Interactors")
#
#check with the FNN
# now actually make the heat map
par(mar=rep(0,4)) # make plot margins 0
heatmap(FNN.M, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=0:1, frame=T, cexRow = 0.5 , cexCol = 0.5)
#
legend("topleft", c("Neighbors Ku80","Ku80"), fill=c("red","blue"),  
       bty="n", inset=0, xpd=T,  border=F)

##########################################Ku80 IR30
##############################################
#Ku80 IR30 Networks


#Let's keep the Mat as the Node attribute table and load the Neighbor list as the matrixes
Mat# Node attribute matrix for all the networks

#Ku80ID IR30 nodes Network Edge Lists
Ku_IR30<-read.delim("./Ku80_Files/Ku_IR30.tsv",sep = "\t",col.names = )
str(Ku_IR30)
head(Ku_IR30)

#Generate the igraph network object
Ku.IR30<-graph_from_data_frame(d=Ku_IR30,vertices = Mat, directed = F )#Generate the igraph object with the edge and node annotation attributes
str(vertex_attr(Ku.IR30))#Check node attributes
str(edge_attr(Ku.IR30))#Check for edge attributes
Ku.IR30<-igraph::simplify(Ku.IR30, remove.multiple=TRUE, remove.loops=TRUE) #Remove repetitive edges and self-connecting edges
Ku.IR30
#Small plot
plot(Ku.IR30)

#Generate an adjacency Matrix that will be used after
Ku.IR30_AdjMat<-as_adjacency_matrix(Ku.IR30,type="both",attr = "weight",names=TRUE,sparse = TRUE)
Ku.IR30_AdjMat<-as.matrix(Ku.IR30_AdjMat)
write.csv(Ku.IR30_AdjMat,"./Ku.IR30_AdjMat.csv")

#Generate a vertex size based on degree and count
vsize<-igraph::degree(Ku.IR30)/vcount(Ku.IR30)
#Identify and remove disconnected nodes
isolated<-which(igraph::degree(Ku.IR30)==0)
isolated
length(isolated)# [1] 618
#remove isolated nodes
Ku.IR30<-igraph::delete_vertices(Ku.IR30,isolated)


#Generate layouts for the graph
#Do not run again, to not to change the layouts that I already have
Lyt.Nice<-layout_nicely(Ku.IR30,dim = 3)
Lyt.GOpt<-layout_with_graphopt(Ku.IR30, niter=5000, mass = 30, charge = E(Ku.IR30)$weight)
Lyt.DrL<-layout_with_drl(Ku.IR30, weights = E(Ku.IR30)$weight)
Lyt.FR<-layout.fruchterman.reingold(Ku.IR30, niter=10000, weights=E(Ku.IR30)$weight)#This one looks like the best option
Lyt.KK<-layout_with_kk(Ku.IR30, maxiter=10000, weights=E(Ku.IR30)$weight) #very large separation of peripheral nodes
Lyt.GEM<-layout_with_gem(Ku.IR30, maxiter=10000) #Bad separation
#Lyt.DH<-layout_with_dh(Ku.IR30, maxiter=150, fineiter=max(10, log2(833)))
Lyt.LGL<-layout.lgl(Ku.IR30, maxit=150,maxdelta=833,area=(length(E(Ku.IR30)))^2,coolexp=1.5,repulserad=((length(E(Ku.IR30)))^2)*(length(E(Ku.IR30))))#Do not work


#Star and Tree layouts

V(Ku.IR30)[name=="XRCC5"]
which(V(Ku.IR30)[name=="Ku80"])# 228 #Check what is the rowindex for your desired node
which(V(Ku.IR30)$name=="XRCC5")# 191 #Check what is the rowindex for your desired node
which(V(Ku.IR30)$name=="MRE11A")# 327
which(V(Ku.IR30)$name=="XRCC5")# 213
Lyt.Star<-layout_as_star(Ku.IR30,center=V(Ku.IR30)[name=="XRCC5"]) #it has a special order argument where we can assign an order to the vertex, hence maybe the OPLS-DA:VIP value order may be used as ordering vector

Lyt.Tree<-layout_as_tree(Ku.IR30,root=c(213))
Lyt.RT<-layout.reingold.tilford(Ku.IR30,root=c(213))# This one looks good, but the indexing do not allow to highlight the correct nodes
Lyt.sugi<-layout_with_sugiyama(Ku.IR30, maxiter=100, weights=E(Ku.IR302)$weight)#This one looks good for undirected tree layouts

#Check what layout looks fine

plot.igraph(Ku.IR30,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(Ku.IR30)$An.color, vertex.frame.color=V(Ku.IR30)$An.color, vertex.shape="circle",
            vertex.size=V(Ku.IR30)$Ku.FC.IR30, vertex.label=V(Ku.IR30)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(Ku.IR30)$weight, edge.lty=1, main="Ku80ID IR30 Network", sub="FR Tree Layout",
            rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#With ggraph

ggraph(Ku.IR30,layout=Lyt.FR[-isolated,])+
  geom_edge_fan(color="gray65", width=1, alpha=0.4)+
  geom_node_point(aes(color=V(Ku.IR30)$Annotation, size=V(Ku.IR30)$Ku.FC.IR30),alpha=0.75)+
  geom_node_text(aes(label=name),color="gray25",size=4)+
  scale_color_manual(limits=V(Ku.IR30)$Annotation, values=V(Ku.IR30)$An.color)+
  scale_size_continuous(range = c(5, 15)) +
  theme_void()+
  title("Ku80 IR30 Network")+
  guides(color=guide_legend(override.aes = list(size=7)))+
  theme(legend.text=element_text(size=12,face = "bold"),
        legend.title = element_text(size=13,face = "bold"))



#Community Detection
#clp <- igraph::cluster_optimal(Ku.IR30)
class(clp)
V(Ku.IR30)$community <- clp$membership
colrs <- adjustcolor( c("gray50", "tomato", "gold", "yellowgreen"), alpha=.6)


plot.igraph(clp,Ku.IR30,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(Ku.IR30)$An.color, vertex.frame.color=V(Ku.IR30)$An.color, vertex.shape="circle",
            vertex.size==V(Ku.IR30)$Ku.FC.IR30, vertex.label=V(Ku.IR30)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(Ku.IR30)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#Distances edges highlight
dist.from.Ku80 <- distances(Ku.IR30, v=V(Ku.IR30)[name=="XRCC5"],
                            to=V(Ku.IR30), weights=E(Ku.IR30)$weight)
dist.from.Ku80<-dist.from.Ku80[,!is.infinite(colSums(dist.from.Ku80))]

dist.from.Ku80
range(dist.from.Ku80)

#Set colors to plot the distances:
oranges <- colorRampPalette(c("dark red", "gold"))
col <- oranges(max(dist.from.Ku80))
col <- col[dist.from.Ku80+1]
plot(Ku.IR30, layout=Lyt.FR[-isolated,],vertex.color=col, vertex.label=dist.from.Ku80,
     vertex.label.color="white",vertex.shape="circle",
     vertex.size=V(Ku.IR30)$Ku.FC.IR30, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(Ku.IR30)$weight, edge.lty=1, main="Ku80 IR30 Full Nuclear Network", sub="FR Layout",
     rescale=T)


#Shortest path from Ku80 to Nucleolin
KuNcl.path <- shortest_paths(Ku.IR30,
                             from = V(Ku.IR30)[name=="XRCC5"],
                             to = V(Ku.IR30)[name=="NCL"],
                             output = "both") # both path nodes and edges
# Generate edge color variable to plot the path:
ecol <- rep("gray80", ecount(Ku.IR30))
ecol[unlist(KuNcl.path$epath)] <- "orange"
# Generate edge width variable to plot the path:
ew <- rep(2, ecount(Ku.IR30))
ew[unlist(KuNcl.path$epath)] <- 4
# Generate node color variable to plot the path:
vcol <- rep("gray40", vcount(Ku.IR30))
vcol[unlist(KuNcl.path$vpath)] <- "gold"
plot(Ku.IR30, layout=Lyt.FR[-isolated,], vertex.color=vcol, edge.color=ecol,
     edge.width=ew,vertex.label.color="white",vertex.shape="circle",
     vertex.size=5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(Ku.IR30)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)



#Plot based on Ku80 as incident
inc.edges <- incident(Ku.IR30, V(Ku.IR30)[name=="XRCC5"], mode="all")
# Set colors to plot the selected edges.
ecol <- rep("gray80", ecount(Ku.IR30))
ecol[inc.edges] <- "orange"
vcol <- rep("gray40", vcount(Ku.IR30))
vcol[V(Ku.IR30)$name=="Ku80"] <- "gold"

plot(Ku.IR30, layout=Lyt.FR[-isolated,], vertex.color=vcol, edge.color=ecol, 
     edge.width=ew,vertex.label.color="white",vertex.shape="circle",
     vertex.size=V(Ku.IR30)$Ku.FC.IR30, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(Ku.IR30)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)


#Plot Degree Distribution
Ku.IR30.Deg<-degree_distribution(Ku.IR30, cumulative = T, mode="all")
plot( x=0:max(igraph::degree(Ku.IR30)), y=1-Ku.IR30.Deg, pch=19, cex=1.2, col="orange",
      xlab="Degree", ylab="Cumulative Frequency")

#Network Descriptives
#Number of nodes and edges
vcount(Ku.IR30) #[1] 215
ecount(Ku.IR30)#[1]  1488
#Density (Number of edges in relation to the total number of edges in a full connected graph)
edge_density(Ku.IR30,loops = F) #[1] 0.06468159
#Transitivity (Total number of conected triangles)
transitivity(Ku.IR30,type="global") #[1] 0.6140637
#modularity
mod<-V(Ku.IR30)$Ku.FuzzyCMeans #lets check for fuzzy Cmeans modularity
mod<-mod+1#The non clustered interactors that appears as 0, needs to be transformed to 1, so we add 1 to every cluster
modularity(Ku.IR30,membership = mod)
#Ku.FuzzyCMeans modularity [1] 0.0176059
#Ku.FuzzyCMeans modularity [1] 0.0244555
#Mre.FuzzyCMeans modularity [1] 0.02041376

#Assorativity
assortativity_nominal(Ku.IR30,factor(V(Ku.IR30)$Annotation),directed=F)
#Annotation [1] [1] 0.1897503
#S6.FuzzyCMeans [1] 0.02560364
#Ku.FuzzyCMeans [1] 0.0408695
#Mre.FuzzyCMeans [1] 0.02494721

#Diameter (longest geodesic distance), shortest path between two most peripheral nodes
diameter(Ku.IR30,directed = F, weights = E(Ku.IR30)$weights)#[1] 4.799
diam<-get_diameter(Ku.IR30,directed = F,weights = E(Ku.IR30)$weights)
diam #9/215 vertices, named, from e733417: [1] HIST1H2AJ HIST2H2AC HIST2H2BE XRCC5     RBM39     SF3B3     LMNA      LMNB1     VIM 
#Plot Diameter
vcol <- rep("gray50", vcount(Ku.IR30))
vcol[diam] <- "gold"
ecol <- rep("gray90", ecount(Ku.IR30))
ecol[E(Ku.IR30, path=diam)] <- "orange"
# E(net, path=diam) finds edges along a path, here 'diam'
plot(Ku.IR30, vertex.color=vcol, edge.color=ecol, vertex.label.color="black",edge.arrow.mode=0,layout=Lyt.FR[-isolated,],vertex.size=V(Ku.IR30)$Ku.FC.IR30)

#Node topological descriptives
degree<-igraph::centr_degree(Ku.IR30,mode = "all",normalized = T)
degree<-degree$res
eigen_cent<-igraph::eigen_centrality(Ku.IR30)
eigen_cent<-eigen_cent$vector
clos_cent <- igraph::closeness(Ku.IR30, normalized=T)
betw_cent <- igraph::betweenness(Ku.IR30, normalized=T)
hs <- hub_score(Ku.IR30, weights=E(Ku.IR30)$weight)$vector
as <- authority_score(Ku.IR30, weights=E(Ku.IR30)$weight)$vector
Trans<-transitivity(Ku.IR30,type="local")
Trans
#write results
Ku.IR30.Centrality<-data.frame(NAME=V(Ku.IR30)$name, UniprotID=V(Ku.IR30)$Uniprot, Degree=degree, Eigen= eigen_cent, Closeness= clos_cent, Betweeness= betw_cent, Hubs=hs, Authority=as,Transitivity=Trans)
Ku.IR30.Centrality
write.csv(Ku.IR30.Centrality,"./Ku.IR30.Centrality.csv")

#Cocitation, Number of shared interactions among nodes
Cocitation<-cocitation(Ku.IR30)
Cocitation#This matrix is similar to an adjacency matrix
write.csv(Cocitation,"./Ku.IR30.cocitation.csv")
#Plot with degree as node size
ggraph(Ku.IR30,layout=Lyt.FR[-isolated,])+
  geom_edge_fan(color="gray85", width=1, alpha=0.4)+
  geom_node_point(aes(color=V(Ku.IR30)$Annotation, size=degree),alpha=0.75)+
  geom_node_text(aes(label=name),color="gray25",size=4)+
  scale_color_manual(limits=V(Ku.IR30)$Annotation, values=V(Ku.IR30)$An.color)+
  scale_size_continuous(range = c(5, 15)) +
  theme_void()+
  title("Ku80 IR30 Network")+
  guides(color=guide_legend(override.aes = list(size=7)))+
  theme(legend.text=element_text(size=12,face = "bold"),
        legend.title = element_text(size=13,face = "bold"))


#Modularity Analysis

#Edge_Betweeness Modularity
Ku.IR30.ebmod=edge.betweenness.community(Ku.IR30,weights=E(Ku.IR30)$weight,directed = F) 
Ku.IR30.ebmod# IGRAPH clustering edge betweenness, groups: 25, mod: 0.41

#
plot.igraph(Ku.IR30.ebmod,Ku.IR30,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(Ku.IR30)$An.color, vertex.frame.color=V(Ku.IR30)$An.color, vertex.shape="circle",
            vertex.size=5, vertex.label=V(Ku.IR30)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(Ku.IR30)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)

#
plot(Ku.IR30.ebmod, Ku.IR30,layout=Lyt.FR[-isolated,], alpha=0.6, 
     vertex.color=V(Ku.IR30)$An.color, vertex.frame.color=V(Ku.IR30)$An.color, vertex.shape="circle", 
     vertex.size=5, vertex.label=V(Ku.IR30)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black", 
     edge.width=E(Ku.IR30)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout", 
     rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)
#Dendogram of the Edge_Betweeness
dendPlot(Ku.IR30.ebmod, mode="hclust")

#Louvain modularity
Ku.IR30.Louvain.mod<-cluster_louvain(Ku.IR30,weights=E(Ku.IR30)$weight)
Ku.IR30.Louvain.mod #IGRAPH clustering multi level, groups: 13, mod: 0.46

#Plot the Louvain communities
plot(Ku.IR30.Louvain.mod, Ku.IR30,layout=Lyt.FR[-isolated,], alpha=0.6, 
     vertex.color=V(Ku.IR30)$An.color, vertex.frame.color=V(Ku.IR30)$An.color, vertex.shape="circle", 
     vertex.size=5, vertex.label=V(Ku.IR30)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black", 
     edge.width=E(Ku.IR30)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout", 
     rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#Generate HeatMaps of the adjacency matrixes
#Let's use the neighbor list for generating the Adjaency Matrix
SensorID_FNN #Let's use the edge list that contain the total amount of possible interactors for the three sensors
v<-levels(as.factor(unlist(SensorID_FNN[,1:2]))) #get vertex names from the network list
n<-length(v)#get the number of vertex in the network
e<-matrix(match(as.character(unlist(SensorID_FNN[,1:2])),v),ncol = 2)#Edgelist with row index instead of Names
w<-Ku_IR30$weight #Edge weights, only from the Ku_IR30 network to quit those edges that do not exist at Ku_IR30

# M is our co-expression network adjacency matrix.

# distringuish between i,j and j,i we want to make M symmetric (undirected) by
# considering an undirected network.
Ku.IR30.M <- matrix(0, n, n)                  # set up a co-expression matrix
Ku.IR30.M[e] <- w                             # fill it in with edge weights
Ku.IR30.M <- Ku.IR30.M + t(Ku.IR30.M)                         # make this symmetric
any(is.na(Ku.IR30.M))#check no NAs in the matrix
dimnames(Ku.IR30.M) <- list(v, v)             # label the vertices
# We let A_ij = 1 if there is an interaction. 
#A <- 1*(M > 0) #right now we keep the weights as they are
Ku.IR30.M #The adjacency Matrix not as an igraph object, easier to manipulate, and to subset for smaller Sensor specific networks
quantile(Ku.IR30.M,c(0,0.25,0.5,0.75,1))
#
# 0%   25%   50%   75%  100% 
#0.000 0.000 0.000 0.000 0.991  

#Creating the HeatMap
# make a heatmap image with colors: white below 0,
# grey from 0.4 to 1 with darker closer to 1
breaks <- c(0, seq(0.4, 0.991, length=11)) # breaks in the color bins
cols <- grey(1-c(0,seq(0.5,1,length=10))) # colors

#Red Blue Colors
library(RColorBrewer)
breaks <- c(seq(-1, 1, length=11))
breaks
cols <- gray(1-c(0,seq(0,1,length=10))) # colors
cols2 <- hcl.colors(10 , palette = "Blues3" , alpha = 1 , rev = TRUE)
cols3 <- colorRamp("red" , space = "rgb" , interpolate = "linear")
HMCol<- brewer.pal(10,name = "RdBu")
HMCol<- rev(HMCol)
HMCol


# color the vertex for CCNB1 blue, its neighbors red, and the others white
SensorIDBind <- which(v == "Ku80") # index for ccnb1
#SensorIDBind<-as.factor(SensorIDBind)
vcols <- rep("white",n)
vcols[SensorIDBind] <- "blue"
#Subsetting Ku.IR30.M by doing a loop asking for look into every row and column the index in SensorIDBind worked
for (i in SensorIDBind) {
  vcols[which(Ku.IR30.M[,i] > 0 | Ku.IR30.M[i,])] <- "red"
}
##Color annotation for functional annotation
MatSet<-subset(Mat,Mat$NAME%in%v)# select only those proteins that are cointained in the node list of our matrixes
match(v,MatSet$NAME) # Check that the indexes of the Annotation Dataframe and the node list are in the same position 
match(rownames(Ku.IR30.M),MatSet$NAME) #Check that the row names of the Adj Matrix and the Annotation dataframe are in the same order
match(colnames(Ku.IR30.M),MatSet$NAME) #Check that the column names of the Adj Matrix and the Annotation dataframe are in the same order
AnCol<-MatSet$An.color#Vector of colors
AnCol<-as.character(AnCol)#transform it into a character vector, so the colors can be read

#Subsetting Ku.IR30.M using the index vector in SensorIDBind do not work
#vcols[which(Ku.IR30.M[,SensorIDBind] > 0 | Ku.IR30.M[SensorIDBind,])] <- "red"

#Generate a dendogram for keeping the proteins in the sample place
#Ku.IR30.dist<-dist(Ku.IR30.M,method = "euclidian")
#Ku.IR30.Dend<-hclust(Ku.IR30.dist, method = "ward.D" )
#plot(Ku.IR30.Dend)

#Cut the dendogram
#cut.height <- 10 # try numbers from 4 to 7
#FNN.hclusts<- cutree(FNN.Dend, cut.height)
#plot(FNN.hclusts)

# now actually make the heat map
par(mar=rep(0,4)) # make plot margins 0
heatmap(Ku.IR30.M, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=cols, breaks=breaks,frame=T, cexRow = 0.5 , cexCol = 0.5, main = "Ku80 IR30 AdjMat", ylab="Total Protein Interactors",xlab="Total Protein Interactors")
#
plot.new()
legend("topleft", c("Cell Cycle","Chromatin Org","DNA Repair","Nuclear Pore","p53","Proteasome","Ribosome Biogenesis","RNA Metabolism","Transcription Reg","UPR"), 
       fill=c("coral1","darkviolet","firebrick1","deepskyblue3","tomato","grey","lightgreen","gold","slateblue1","maroon"),  
       bty="n", inset=0, xpd=T,  border=F, cex =2,y.intersp = 0.7,x.intersp = 0.2,pt.cex = 1)


#Now we substract the KunonIR30 AdjMatrix to the KuIR6 AdjMatrix
KuIR30.AdjMat.Change<-Ku.IR30.M-Ku.nonIR5.M
range(KuIR30.AdjMat.Change)#Let's check the magnitude of the substracted weights
#[1]-0.813  0.668
any(is.na(KuIR30.AdjMat.Change))
#Generate the breaks and color functions
breaks <- c(seq(-0.5, 0.5, length=10))
breaks
HMCol<- brewer.pal(9,name = "RdBu")
HMCol<- rev(HMCol)
HMCol
[1] "#2166AC" "#4393C3" "#92C5DE" "#D1E5F0" "#F7F7F7" "#FDDBC7" "#F4A582" "#D6604D" "#B2182B"

#Plot the HeatMap
par(mar=rep(0,4)) # make plot margins 0
heatmap(KuIR30.AdjMat.Change, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=HMCol, breaks=breaks,frame=T, cexRow = 0.5 , cexCol = 0.5, main = "Ku80 IR30 AdjMat (minus nonIR30 edges)", ylab="Total Protein Interactors",xlab="Total Protein Interactors")


#Now we substract the KuIR5 AdjMatrix to the KuIR30 AdjMatrix
KuIR30.AdjMat.Change2<-Ku.IR30.M-Ku.IR5.M
range(KuIR30.AdjMat.Change)#Let's check the magnitude of the substracted weights
#[1]-0.813  0.668
any(is.na(KuIR30.AdjMat.Change))
#Generate the breaks and color functions
breaks <- c(seq(-0.5, 0.5, length=10))
breaks
HMCol<- brewer.pal(9,name = "RdBu")
HMCol<- rev(HMCol)
HMCol
[1] "#2166AC" "#4393C3" "#92C5DE" "#D1E5F0" "#F7F7F7" "#FDDBC7" "#F4A582" "#D6604D" "#B2182B"

#Plot the HeatMap
par(mar=rep(0,4)) # make plot margins 0
heatmap(KuIR30.AdjMat.Change2, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=HMCol, breaks=breaks,frame=T, cexRow = 0.5 , cexCol = 0.5, main = "Ku80 IR30 AdjMat (minus IR5 edges)", ylab="Total Protein Interactors",xlab="Total Protein Interactors")

####################### Ku80 IR2h Networks
###################################################

#Ku80 IR2 Networks


#Let's keep the Mat as the Node attribute table and load the Neighbor list as the matrixes
Mat# Node attribute matrix for all the networks

#Ku80ID IR2 nodes Network Edge Lists
Ku_IR2<-read.delim("./Ku80_Files/Ku_IR2h.tsv",sep = "\t",col.names = )
str(Ku_IR2)
head(Ku_IR2)

#Generate the igraph network object
Ku.IR2<-graph_from_data_frame(d=Ku_IR2,vertices = Mat, directed = F )#Generate the igraph object with the edge and node annotation attributes
str(vertex_attr(Ku.IR2))#Check node attributes
str(edge_attr(Ku.IR2))#Check for edge attributes
Ku.IR2<-igraph::simplify(Ku.IR2, remove.multiple=TRUE, remove.loops=TRUE) #Remove repetitive edges and self-connecting edges
Ku.IR2
#Small plot
plot(Ku.IR2)

#Generate an adjacency Matrix that will be used after
Ku.IR2_AdjMat<-as_adjacency_matrix(Ku.IR2,type="both",attr = "weight",names=TRUE,sparse = TRUE)
Ku.IR2_AdjMat<-as.matrix(Ku.IR2_AdjMat)
write.csv(Ku.IR2_AdjMat,"./Ku.IR2_AdjMat.csv")

#Generate a vertex size based on degree and count
vsize<-igraph::degree(Ku.IR2)/vcount(Ku.IR2)
#Identify and remove disconnected nodes
isolated<-which(igraph::degree(Ku.IR2)==0)
isolated
length(isolated)# [1] 538
#remove isolated nodes
Ku.IR2<-igraph::delete_vertices(Ku.IR2,isolated)


#Generate layouts for the graph
#Do not run again, to not to change the layouts that I already have
Lyt.Nice<-layout_nicely(Ku.IR2,dim = 3)
Lyt.GOpt<-layout_with_graphopt(Ku.IR2, niter=5000, mass = 30, charge = E(Ku.IR2)$weight)
Lyt.DrL<-layout_with_drl(Ku.IR2, weights = E(Ku.IR2)$weight)
Lyt.FR<-layout.fruchterman.reingold(Ku.IR2, niter=10000, weights=E(Ku.IR2)$weight)#This one looks like the best option
Lyt.KK<-layout_with_kk(Ku.IR2, maxiter=10000, weights=E(Ku.IR2)$weight) #very large separation of peripheral nodes
Lyt.GEM<-layout_with_gem(Ku.IR2, maxiter=10000) #Bad separation
#Lyt.DH<-layout_with_dh(Ku.IR2, maxiter=150, fineiter=max(10, log2(833)))
Lyt.LGL<-layout.lgl(Ku.IR2, maxit=150,maxdelta=833,area=(length(E(Ku.IR2)))^2,coolexp=1.5,repulserad=((length(E(Ku.IR2)))^2)*(length(E(Ku.IR2))))#Do not work


#Star and Tree layouts

V(Ku.IR2)[name=="XRCC5"]
which(V(Ku.IR2)[name=="Ku80"])# 228 #Check what is the rowindex for your desired node
which(V(Ku.IR2)$name=="Ku80")# 191 #Check what is the rowindex for your desired node
which(V(Ku.IR2)$name=="MRE11A")# 327
which(V(Ku.IR2)$name=="XRCC5")# 290
Lyt.Star<-layout_as_star(Ku.IR2,center=V(Ku.IR2)[name=="XRCC5"]) #it has a special order argument where we can assign an order to the vertex, hence maybe the OPLS-DA:VIP value order may be used as ordering vector

Lyt.Tree<-layout_as_tree(Ku.IR2,root=c(290))
Lyt.RT<-layout.reingold.tilford(Ku.IR2,root=c(290))# This one looks good, but the indexing do not allow to highlight the correct nodes
Lyt.sugi<-layout_with_sugiyama(Ku.IR2, maxiter=100, weights=E(Ku.IR22)$weight)#This one looks good for undirected tree layouts

#Check what layout looks fine

plot.igraph(Ku.IR2,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(Ku.IR2)$An.color, vertex.frame.color=V(Ku.IR2)$An.color, vertex.shape="circle",
            vertex.size=V(Ku.IR2)$Ku.FC.IR2, vertex.label=V(Ku.IR2)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(Ku.IR2)$weight, edge.lty=1, main="Ku80ID IR2 Network", sub="FR Tree Layout",
            rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#With ggraph

ggraph(Ku.IR2,layout=Lyt.RT)+#FR[-isolated,])+
  geom_edge_fan(color="gray65", width=1, alpha=0.4)+
  geom_node_point(aes(color=V(Ku.IR2)$Annotation, size=V(Ku.IR2)$Ku.FC.IR2),alpha=0.75)+
  geom_node_text(aes(label=name),color="gray25",size=4)+
  scale_color_manual(limits=V(Ku.IR2)$Annotation, values=V(Ku.IR2)$An.color)+
  scale_size_continuous(range = c(5, 15)) +
  theme_void()+
  title("Ku80 IR2 Network")+
  guides(color=guide_legend(override.aes = list(size=7)))+
  theme(legend.text=element_text(size=12,face = "bold"),
        legend.title = element_text(size=13,face = "bold"))



#Community Detection
#clp <- igraph::cluster_optimal(Ku.IR2)
class(clp)
V(Ku.IR2)$community <- clp$membership
colrs <- adjustcolor( c("gray50", "tomato", "gold", "yellowgreen"), alpha=.6)


plot.igraph(clp,Ku.IR2,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(Ku.IR2)$An.color, vertex.frame.color=V(Ku.IR2)$An.color, vertex.shape="circle",
            vertex.size==V(Ku.IR2)$Ku.FC.IR2, vertex.label=V(Ku.IR2)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(Ku.IR2)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#Distances edges highlight
dist.from.Ku80 <- distances(Ku.IR2, v=V(Ku.IR2)[name=="XRCC5"],
                            to=V(Ku.IR2), weights=E(Ku.IR2)$weight)
dist.from.Ku80<-dist.from.Ku80[,!is.infinite(colSums(dist.from.Ku80))]

dist.from.Ku80
range(dist.from.Ku80)

#Set colors to plot the distances:
oranges <- colorRampPalette(c("dark red", "gold"))
col <- oranges(max(dist.from.Ku80))
col <- col[dist.from.Ku80+1]
plot(Ku.IR2, layout=Lyt.FR[-isolated,],vertex.color=col, vertex.label=dist.from.Ku80,
     vertex.label.color="white",vertex.shape="circle",
     vertex.size=V(Ku.IR2)$Ku.FC.IR2, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(Ku.IR2)$weight, edge.lty=1, main="Ku80 IR2 Full Nuclear Network", sub="FR Layout",
     rescale=T)


#Shortest path from Ku80 to Nucleolin
KuNcl.path <- shortest_paths(Ku.IR2,
                             from = V(Ku.IR2)[name=="XRCC5"],
                             to = V(Ku.IR2)[name=="NCL"],
                             output = "both") # both path nodes and edges
KuNcl.path
# Generate edge color variable to plot the path:
ecol <- rep("gray80", ecount(Ku.IR2))
ecol[unlist(KuNcl.path$epath)] <- "orange"
# Generate edge width variable to plot the path:
ew <- rep(2, ecount(Ku.IR2))
ew[unlist(KuNcl.path$epath)] <- 4
# Generate node color variable to plot the path:
vcol <- rep("gray40", vcount(Ku.IR2))
vcol[unlist(KuNcl.path$vpath)] <- "gold"
plot(Ku.IR2, layout=Lyt.FR[-isolated,], vertex.color=vcol, edge.color=ecol,
     edge.width=ew,vertex.label.color="white",vertex.shape="circle",
     vertex.size=5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(Ku.IR2)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)



#Plot based on Ku80 as incident
inc.edges <- incident(Ku.IR2, V(Ku.IR2)[name=="XRCC5"], mode="all")
inc.edges
# Set colors to plot the selected edges.
ecol <- rep("gray80", ecount(Ku.IR2))
ecol[inc.edges] <- "orange"
vcol <- rep("gray40", vcount(Ku.IR2))
vcol[V(Ku.IR2)$name=="XRCC5"] <- "gold"

plot(Ku.IR2, layout=Lyt.FR[-isolated,], vertex.color=vcol, edge.color=ecol, 
     edge.width=ew,vertex.label.color="white",vertex.shape="circle",
     vertex.size=V(Ku.IR2)$Ku.FC.IR2, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(Ku.IR2)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)


#Plot Degree Distribution
Ku.IR2.Deg<-degree_distribution(Ku.IR2, cumulative = T, mode="all")
plot( x=0:max(igraph::degree(Ku.IR2)), y=1-Ku.IR2.Deg, pch=19, cex=1.2, col="orange",
      xlab="Degree", ylab="Cumulative Frequency")

#Network Descriptives
#Number of nodes and edges
vcount(Ku.IR2) #[1] 295
ecount(Ku.IR2)#[1]  2442
#Density (Number of edges in relation to the total number of edges in a full connected graph)
edge_density(Ku.IR2,loops = F) #[1]0.05631269
#Transitivity (Total number of conected triangles)
transitivity(Ku.IR2,type="global") #[1] 0.6437669
#modularity
mod<-V(Ku.IR2)$S6.FuzzyCMeans #lets check for fuzzy Cmeans modularity
mod<-mod+1#The non clustered interactors that appears as 0, needs to be transformed to 1, so we add 1 to every cluster
modularity(Ku.IR2,membership = mod)
#S6.FuzzyCMeans modularity [1] 0.01398724
#Ku.FuzzyCMeans modularity [1] 0.04469868
#Mre.FuzzyCMeans modularity [1] 0.02629766

#Assorativity
assortativity_nominal(Ku.IR2,factor(V(Ku.IR2)$Mre.FuzzyCMeans),directed=F)
#Annotation [1] [1] 0.2412916
#S6.FuzzyCMeans [1] 0.02085077
#Ku.FuzzyCMeans [1] 0.06838574
#Mre.FuzzyCMeans [1]  0.03297979

#Diameter (longest geodesic distance), shortest path between two most peripheral nodes
diameter(Ku.IR2,directed = F, weights = E(Ku.IR2)$weights)#[1] 4.362
diam<-get_diameter(Ku.IR2,directed = F,weights = E(Ku.IR2)$weights)
diam #+ + 10/295 vertices, named, from 22aebf2:[1] CALR    HSP90B1 HSPA5   HSPA1B  APEX1   XRCC5   HNRNPA1 SFPQ    CPSF7   NUDT21 
#Plot Diameter
vcol <- rep("gray50", vcount(Ku.IR2))
vcol[diam] <- "gold"
ecol <- rep("gray90", ecount(Ku.IR2))
ecol[E(Ku.IR2, path=diam)] <- "orange"
# E(net, path=diam) finds edges along a path, here 'diam'
plot(Ku.IR2, vertex.color=vcol, edge.color=ecol, vertex.label.color="black",edge.arrow.mode=0,layout=Lyt.FR[-isolated,],vertex.size=V(Ku.IR2)$Ku.FC.IR2)

#Node topological descriptives
degree<-igraph::centr_degree(Ku.IR2,mode = "all",normalized = T)
degree<-degree$res
eigen_cent<-igraph::eigen_centrality(Ku.IR2)
eigen_cent<-eigen_cent$vector
clos_cent <- igraph::closeness(Ku.IR2, normalized=T)
betw_cent <- igraph::betweenness(Ku.IR2, normalized=T)
hs <- hub_score(Ku.IR2, weights=E(Ku.IR2)$weight)$vector
as <- authority_score(Ku.IR2, weights=E(Ku.IR2)$weight)$vector
Trans<-transitivity(Ku.IR2,type="local")
Trans
#write results
Ku.IR2.Centrality<-data.frame(NAME=V(Ku.IR2)$name, UniprotID=V(Ku.IR2)$Uniprot, Degree=degree, Eigen= eigen_cent, Closeness= clos_cent, Betweeness= betw_cent, Hubs=hs, Authority=as,Transitivity=Trans)
Ku.IR2.Centrality
write.csv(Ku.IR2.Centrality,"./Ku.IR2.Centrality.csv")

#Cocitation, Number of shared interactions among nodes
Cocitation<-cocitation(Ku.IR2)
Cocitation#This matrix is similar to an adjacency matrix
write.csv(Cocitation,"./Ku.IR2.cocitation.csv")
#Plot with degree as node size
ggraph(Ku.IR2,layout=Lyt.FR[-isolated,])+
  geom_edge_fan(color="gray85", width=1, alpha=0.4)+
  geom_node_point(aes(color=V(Ku.IR2)$Annotation, size=degree),alpha=0.75)+
  geom_node_text(aes(label=name),color="gray25",size=4)+
  scale_color_manual(limits=V(Ku.IR2)$Annotation, values=V(Ku.IR2)$An.color)+
  scale_size_continuous(range = c(5, 15)) +
  theme_void()+
  title("Ku80 IR2 Network")+
  guides(color=guide_legend(override.aes = list(size=7)))+
  theme(legend.text=element_text(size=12,face = "bold"),
        legend.title = element_text(size=13,face = "bold"))


#Modularity Analysis

#Edge_Betweeness Modularity
Ku.IR2.ebmod=edge.betweenness.community(Ku.IR2,weights=E(Ku.IR2)$weight,directed = F) 
Ku.IR2.ebmod# IGRAPH clustering edge betweenness, groups: 30, mod: 0.45

#
plot.igraph(Ku.IR2.ebmod,Ku.IR2,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(Ku.IR2)$An.color, vertex.frame.color=V(Ku.IR2)$An.color, vertex.shape="circle",
            vertex.size=5, vertex.label=V(Ku.IR2)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(Ku.IR2)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)

#
plot(Ku.IR2.ebmod, Ku.IR2,layout=Lyt.FR[-isolated,], alpha=0.6, 
     vertex.color=V(Ku.IR2)$An.color, vertex.frame.color=V(Ku.IR2)$An.color, vertex.shape="circle", 
     vertex.size=5, vertex.label=V(Ku.IR2)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black", 
     edge.width=E(Ku.IR2)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout", 
     rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)
#Dendogram of the Edge_Betweeness
dendPlot(Ku.IR2.ebmod, mode="hclust")

#Louvain modularity
Ku.IR2.Louvain.mod<-cluster_louvain(Ku.IR2,weights=E(Ku.IR2)$weight)
Ku.IR2.Louvain.mod #IGRAPH clustering multi level, groups: 12, mod: 0.48

#Plot the Louvain communities
plot(Ku.IR2.Louvain.mod, Ku.IR2,layout=Lyt.FR[-isolated,], alpha=0.6, 
     vertex.color=V(Ku.IR2)$An.color, vertex.frame.color=V(Ku.IR2)$An.color, vertex.shape="circle", 
     vertex.size=5, vertex.label=V(Ku.IR2)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black", 
     edge.width=E(Ku.IR2)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout", 
     rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#Generate HeatMaps of the adjacency matrixes
#Let's use the neighbor list for generating the Adjaency Matrix
SensorID_FNN #Let's use the edge list that contain the total amount of possible interactors for the three sensors
v<-levels(as.factor(unlist(SensorID_FNN[,1:2]))) #get vertex names from the network list
n<-length(v)#get the number of vertex in the network
e<-matrix(match(as.character(unlist(SensorID_FNN[,1:2])),v),ncol = 2)#Edgelist with row index instead of Names
w<-Ku_IR2$weight #Edge weights, only from the Ku_IR2 network to quit those edges that do not exist at Ku_IR2

# M is our co-expression network adjacency matrix.

# distringuish between i,j and j,i we want to make M symmetric (undirected) by
# considering an undirected network.
Ku.IR2.M <- matrix(0, n, n)                  # set up a co-expression matrix
Ku.IR2.M[e] <- w                             # fill it in with edge weights
Ku.IR2.M <- Ku.IR2.M + t(Ku.IR2.M)                         # make this symmetric
any(is.na(Ku.IR2.M))#check no NAs in the matrix
dimnames(Ku.IR2.M) <- list(v, v)             # label the vertices
# We let A_ij = 1 if there is an interaction. 
#A <- 1*(M > 0) #right now we keep the weights as they are
Ku.IR2.M #The adjacency Matrix not as an igraph object, easier to manipulate, and to subset for smaller Sensor specific networks
quantile(Ku.IR2.M,c(0,0.25,0.5,0.75,1))
#
#    0%   25%   50%   75%  100% 
#  0.000 0.000 0.000 0.000 0.991  

#Creating the HeatMap
# make a heatmap image with colors: white below 0,
# grey from 0.4 to 1 with darker closer to 1
breaks <- c(0, seq(0.4, 0.991, length=11)) # breaks in the color bins
cols <- grey(1-c(0,seq(0.5,1,length=10))) # colors

#Red Blue Colors
library(RColorBrewer)
breaks <- c(seq(-1, 1, length=11))
breaks
cols <- gray(1-c(0,seq(0,1,length=10))) # colors
cols2 <- hcl.colors(10 , palette = "Blues3" , alpha = 1 , rev = TRUE)
cols3 <- colorRamp("red" , space = "rgb" , interpolate = "linear")
HMCol<- brewer.pal(10,name = "RdBu")
HMCol<- rev(HMCol)
HMCol


# color the vertex for CCNB1 blue, its neighbors red, and the others white
SensorIDBind <- which(v == "Ku80") # index for ccnb1
#SensorIDBind<-as.factor(SensorIDBind)
vcols <- rep("white",n)
vcols[SensorIDBind] <- "blue"
#Subsetting Ku.IR2.M by doing a loop asking for look into every row and column the index in SensorIDBind worked
for (i in SensorIDBind) {
  vcols[which(Ku.IR2.M[,i] > 0 | Ku.IR2.M[i,])] <- "red"
}
##Color annotation for functional annotation
MatSet<-subset(Mat,Mat$NAME%in%v)# select only those proteins that are cointained in the node list of our matrixes
match(v,MatSet$NAME) # Check that the indexes of the Annotation Dataframe and the node list are in the same position 
match(rownames(Ku.IR2.M),MatSet$NAME) #Check that the row names of the Adj Matrix and the Annotation dataframe are in the same order
match(colnames(Ku.IR2.M),MatSet$NAME) #Check that the column names of the Adj Matrix and the Annotation dataframe are in the same order
AnCol<-MatSet$An.color#Vector of colors
AnCol<-as.character(AnCol)#transform it into a character vector, so the colors can be read
#Subsetting Ku.IR2.M using the index vector in SensorIDBind do not work
#vcols[which(Ku.IR2.M[,SensorIDBind] > 0 | Ku.IR2.M[SensorIDBind,])] <- "red"

#Generate a dendogram for keeping the proteins in the sample place
#Ku.IR2.dist<-dist(Ku.IR2.M,method = "euclidian")
#Ku.IR2.Dend<-hclust(Ku.IR2.dist, method = "ward.D" )
#plot(Ku.IR2.Dend)

#Cut the dendogram
#cut.height <- 10 # try numbers from 4 to 7
#FNN.hclusts<- cutree(FNN.Dend, cut.height)
#plot(FNN.hclusts)

# now actually make the heat map
par(mar=rep(0,4)) # make plot margins 0
heatmap(Ku.IR2.M, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=cols, breaks=breaks,frame=T, cexRow = 0.5 , cexCol = 0.5, main = "Ku80 IR2 AdjMat", ylab="Total Protein Interactors",xlab="Total Protein Interactors")
#
plot.new()
legend("topleft", c("Cell Cycle","Chromatin Org","DNA Repair","Nuclear Pore","p53","Proteasome","Ribosome Biogenesis","RNA Metabolism","Transcription Reg","UPR"), 
       fill=c("coral1","darkviolet","firebrick1","deepskyblue3","tomato","grey","lightgreen","gold","slateblue1","maroon"),  
       bty="n", inset=0, xpd=T,  border=F, cex =2,y.intersp = 0.7,x.intersp = 0.2,pt.cex = 1)


#Now we substract the KunonIR2 AdjMatrix to the KuIR6 AdjMatrix
KuIR2.AdjMat.Change<-Ku.IR2.M-Ku.nonIR5.M
range(KuIR2.AdjMat.Change)#Let's check the magnitude of the substracted weights
#[1]-0.813  0.668
any(is.na(KuIR2.AdjMat.Change))
#Generate the breaks and color functions
breaks <- c(seq(-0.5, 0.5, length=10))
breaks
HMCol<- brewer.pal(9,name = "RdBu")
HMCol<- rev(HMCol)
HMCol
[1] "#2166AC" "#4393C3" "#92C5DE" "#D1E5F0" "#F7F7F7" "#FDDBC7" "#F4A582" "#D6604D" "#B2182B"

#Plot the HeatMap
par(mar=rep(0,4)) # make plot margins 0
heatmap(KuIR2.AdjMat.Change, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=HMCol, breaks=breaks,frame=T, cexRow = 0.5 , cexCol = 0.5, main = "Ku80 IR2 AdjMat (minus nonIR2 edges)", ylab="Total Protein Interactors",xlab="Total Protein Interactors")


#Now we substract the KuIR5 AdjMatrix to the KuIR2 AdjMatrix
KuIR2.AdjMat.Change2<-Ku.IR2.M-Ku.IR30.M
range(KuIR2.AdjMat.Change)#Let's check the magnitude of the substracted weights
#[1]-0.813  0.668
any(is.na(KuIR2.AdjMat.Change))
#Generate the breaks and color functions
breaks <- c(seq(-0.5, 0.5, length=10))
breaks
HMCol<- brewer.pal(9,name = "RdBu")
HMCol<- rev(HMCol)
HMCol
[1] "#2166AC" "#4393C3" "#92C5DE" "#D1E5F0" "#F7F7F7" "#FDDBC7" "#F4A582" "#D6604D" "#B2182B"

#Plot the HeatMap
par(mar=rep(0,4)) # make plot margins 0
heatmap(KuIR2.AdjMat.Change2, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=HMCol, breaks=breaks,frame=T, cexRow = 0.5 , cexCol = 0.5, main = "Ku80 IR2 AdjMat (minus IR30 edges)", ylab="Total Protein Interactors",xlab="Total Protein Interactors")

#
#check with the FNN
# now actually make the heat map
par(mar=rep(0,4)) # make plot margins 0
heatmap(FNN.M, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=0:1, frame=T, cexRow = 0.5 , cexCol = 0.5)
#
legend("topleft", c("Neighbors Ku80","Ku80"), fill=c("red","blue"),  
       bty="n", inset=0, xpd=T,  border=F)


####################### Ku80 IR8h Networks
###################################################

#Ku80 IR8 Networks


#Let's keep the Mat as the Node attribute table and load the Neighbor list as the matrixes
Mat# Node attribute matrix for all the networks

#Ku80ID IR8 nodes Network Edge Lists
Ku_IR8<-read.delim("./Ku80_Files/Ku_IR8h.tsv",sep = "\t",col.names = )
str(Ku_IR8)
head(Ku_IR8)

#Generate the igraph network object
Ku.IR8<-graph_from_data_frame(d=Ku_IR8,vertices = Mat, directed = F )#Generate the igraph object with the edge and node annotation attributes
str(vertex_attr(Ku.IR8))#Check node attributes
str(edge_attr(Ku.IR8))#Check for edge attributes
Ku.IR8<-igraph::simplify(Ku.IR8, remove.multiple=TRUE, remove.loops=TRUE) #Remove repetitive edges and self-connecting edges
Ku.IR8
#Small plot
plot(Ku.IR8)

#Generate an adjacency Matrix that will be used after
Ku.IR8_AdjMat<-as_adjacency_matrix(Ku.IR8,type="both",attr = "weight",names=TRUE,sparse = TRUE)
Ku.IR8_AdjMat<-as.matrix(Ku.IR8_AdjMat)
write.csv(Ku.IR8_AdjMat,"./Ku.IR8_AdjMat.csv")

#Generate a vertex size based on degree and count
vsize<-igraph::degree(Ku.IR8)/vcount(Ku.IR8)
#Identify and remove disconnected nodes
isolated<-which(igraph::degree(Ku.IR8)==0)
isolated
length(isolated)# [1] 535
#remove isolated nodes
Ku.IR8<-igraph::delete_vertices(Ku.IR8,isolated)


#Generate layouts for the graph
#Do not run again, to not to change the layouts that I already have
Lyt.Nice<-layout_nicely(Ku.IR8,dim = 3)
Lyt.GOpt<-layout_with_graphopt(Ku.IR8, niter=5000, mass = 30, charge = E(Ku.IR8)$weight)
Lyt.DrL<-layout_with_drl(Ku.IR8, weights = E(Ku.IR8)$weight)
Lyt.FR<-layout.fruchterman.reingold(Ku.IR8, niter=10000, weights=E(Ku.IR8)$weight)#This one looks like the best option
Lyt.KK<-layout_with_kk(Ku.IR8, maxiter=10000, weights=E(Ku.IR8)$weight) #very large separation of peripheral nodes
Lyt.GEM<-layout_with_gem(Ku.IR8, maxiter=10000) #Bad separation
#Lyt.DH<-layout_with_dh(Ku.IR8, maxiter=150, fineiter=max(10, log2(833)))
Lyt.LGL<-layout.lgl(Ku.IR8, maxit=150,maxdelta=833,area=(length(E(Ku.IR8)))^2,coolexp=1.5,repulserad=((length(E(Ku.IR8)))^2)*(length(E(Ku.IR8))))#Do not work


#Star and Tree layouts

V(Ku.IR8)[name=="XRCC5"]
which(V(Ku.IR8)[name=="Ku80"])# 228 #Check what is the rowindex for your desired node
which(V(Ku.IR8)$name=="XRCC5")# 191 #Check what is the rowindex for your desired node
which(V(Ku.IR8)$name=="MRE11A")# 327
which(V(Ku.IR8)$name=="XRCC5")# 350
Lyt.Star<-layout_as_star(Ku.IR8,center=V(Ku.IR8)[name=="XRCC5"]) #it has a special order argument where we can assign an order to the vertex, hence maybe the OPLS-DA:VIP value order may be used as ordering vector

Lyt.Tree<-layout_as_tree(Ku.IR8,root=c(350))
Lyt.RT<-layout.reingold.tilford(Ku.IR8,root=c(350))# This one looks good, but the indexing do not allow to highlight the correct nodes
Lyt.sugi<-layout_with_sugiyama(Ku.IR8, maxiter=100, weights=E(Ku.IR82)$weight)#This one looks good for undirected tree layouts

#Check what layout looks fine

plot.igraph(Ku.IR8,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(Ku.IR8)$An.color, vertex.frame.color=V(Ku.IR8)$An.color, vertex.shape="circle",
            vertex.size=V(Ku.IR8)$Ku.FC.IR8, vertex.label=V(Ku.IR8)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(Ku.IR8)$weight, edge.lty=1, main="Ku80ID IR8 Network", sub="FR Tree Layout",
            rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#With ggraph

ggraph(Ku.IR8,layout=Lyt.FR[-isolated,])+
  geom_edge_fan(color="gray65", width=1, alpha=0.4)+
  geom_node_point(aes(color=V(Ku.IR8)$Annotation, size=V(Ku.IR8)$Ku.FC.IR8),alpha=0.75)+
  geom_node_text(aes(label=name),color="gray25",size=4)+
  scale_color_manual(limits=V(Ku.IR8)$Annotation, values=V(Ku.IR8)$An.color)+
  scale_size_continuous(range = c(5, 15)) +
  theme_void()+
  title("Ku80 IR8 Network")+
  guides(color=guide_legend(override.aes = list(size=7)))+
  theme(legend.text=element_text(size=12,face = "bold"),
        legend.title = element_text(size=13,face = "bold"))



#Community Detection
#clp <- igraph::cluster_optimal(Ku.IR8)
class(clp)
V(Ku.IR8)$community <- clp$membership
colrs <- adjustcolor( c("gray50", "tomato", "gold", "yellowgreen"), alpha=.6)


plot.igraph(clp,Ku.IR8,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(Ku.IR8)$An.color, vertex.frame.color=V(Ku.IR8)$An.color, vertex.shape="circle",
            vertex.size==V(Ku.IR8)$Ku.FC.IR8, vertex.label=V(Ku.IR8)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(Ku.IR8)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#Distances edges highlight
dist.from.Ku80 <- distances(Ku.IR8, v=V(Ku.IR8)[name=="XRCC5"],
                            to=V(Ku.IR8), weights=E(Ku.IR8)$weight)
dist.from.Ku80<-dist.from.Ku80[,!is.infinite(colSums(dist.from.Ku80))]

dist.from.Ku80
range(dist.from.Ku80)

#Set colors to plot the distances:
oranges <- colorRampPalette(c("dark red", "gold"))
col <- oranges(max(dist.from.Ku80))
col <- col[dist.from.Ku80+1]
plot(Ku.IR8, layout=Lyt.FR[-isolated,],vertex.color=col, vertex.label=dist.from.Ku80,
     vertex.label.color="white",vertex.shape="circle",
     vertex.size=V(Ku.IR8)$Ku.FC.IR8, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(Ku.IR8)$weight, edge.lty=1, main="Ku80 IR8 Full Nuclear Network", sub="FR Layout",
     rescale=T)


#Shortest path from Ku80 to Nucleolin
KuNcl.path <- shortest_paths(Ku.IR8,
                             from = V(Ku.IR8)[name=="XRCC5"],
                             to = V(Ku.IR8)[name=="NCL"],
                             output = "both") # both path nodes and edges
KuNcl.path
# Generate edge color variable to plot the path:
ecol <- rep("gray80", ecount(Ku.IR8))
ecol[unlist(KuNcl.path$epath)] <- "orange"
# Generate edge width variable to plot the path:
ew <- rep(2, ecount(Ku.IR8))
ew[unlist(KuNcl.path$epath)] <- 4
# Generate node color variable to plot the path:
vcol <- rep("gray40", vcount(Ku.IR8))
vcol[unlist(KuNcl.path$vpath)] <- "gold"
plot(Ku.IR8, layout=Lyt.FR[-isolated,], vertex.color=vcol, edge.color=ecol,
     edge.width=ew,vertex.label.color="white",vertex.shape="circle",
     vertex.size=5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(Ku.IR8)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)



#Plot based on Ku80 as incident
inc.edges <- incident(Ku.IR8, V(Ku.IR8)[name=="XRCC5"], mode="all")
# Set colors to plot the selected edges.
ecol <- rep("gray80", ecount(Ku.IR8))
ecol[inc.edges] <- "orange"
vcol <- rep("gray40", vcount(Ku.IR8))
vcol[V(Ku.IR8)$name=="Ku80"] <- "gold"

plot(Ku.IR8, layout=Lyt.FR[-isolated,], vertex.color=vcol, edge.color=ecol, 
     edge.width=ew,vertex.label.color="white",vertex.shape="circle",
     vertex.size=V(Ku.IR8)$Ku.FC.IR8, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(Ku.IR8)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)


#Plot Degree Distribution
Ku.IR8.Deg<-degree_distribution(Ku.IR8, cumulative = T, mode="all")
plot( x=0:max(igraph::degree(Ku.IR8)), y=1-Ku.IR8.Deg, pch=19, cex=1.2, col="orange",
      xlab="Degree", ylab="Cumulative Frequency")

#Network Descriptives
#Number of nodes and edges
vcount(Ku.IR8) #[1] 353
ecount(Ku.IR8)#[1]  3738
#Density (Number of edges in relation to the total number of edges in a full connected graph)
edge_density(Ku.IR8,loops = F) #[1] 0.06016611
#Transitivity (Total number of conected triangles)
transitivity(Ku.IR8,type="global") #[1] 0.6511962
#modularity
mod<-V(Ku.IR8)$S6.FuzzyCMeans #lets check for fuzzy Cmeans modularity
mod<-mod+1#The non clustered interactors that appears as 0, needs to be transformed to 1, so we add 1 to every cluster
modularity(Ku.IR8,membership = mod)
#Ku.FuzzyCMeans modularity [1] 0.01537486
#Ku.FuzzyCMeans modularity [1] 0.01047776
#Mre.FuzzyCMeans modularity [1] 0.01609502

#Assorativity
assortativity_nominal(Ku.IR8,factor(V(Ku.IR8)$Annotation),directed=F)
#Annotation [1] [1] 0.1798564
#Ku.FuzzyCMeans [1] 0.02216852
#Ku.FuzzyCMeans [1]  0.01573453
#Mre.FuzzyCMeans [1] 0.01989416

#Diameter (longest geodesic distance), shortest path between two most peripheral nodes
diameter(Ku.IR8,directed = F, weights = E(Ku.IR8)$weights)#[1] 3.821
diam<-get_diameter(Ku.IR8,directed = F,weights = E(Ku.IR8)$weights)
diam #+ 9/353 vertices, named, from 4b04e55:[1] DDX1     RTCB     HIST1H4F PRPF8    XRCC5    EGFR     ACTN4    ACTN1    TTN 
#Plot Diameter
vcol <- rep("gray50", vcount(Ku.IR8))
vcol[diam] <- "gold"
ecol <- rep("gray90", ecount(Ku.IR8))
ecol[E(Ku.IR8, path=diam)] <- "orange"
# E(net, path=diam) finds edges along a path, here 'diam'
plot(Ku.IR8, vertex.color=vcol, edge.color=ecol, vertex.label.color="black",edge.arrow.mode=0,layout=Lyt.FR[-isolated,],vertex.size=V(Ku.IR8)$Ku.FC.IR8)

#Node topological descriptives
degree<-igraph::centr_degree(Ku.IR8,mode = "all",normalized = T)
degree<-degree$res
eigen_cent<-igraph::eigen_centrality(Ku.IR8)
eigen_cent<-eigen_cent$vector
clos_cent <- igraph::closeness(Ku.IR8, normalized=T)
betw_cent <- igraph::betweenness(Ku.IR8, normalized=T)
hs <- hub_score(Ku.IR8, weights=E(Ku.IR8)$weight)$vector
as <- authority_score(Ku.IR8, weights=E(Ku.IR8)$weight)$vector
Trans<-transitivity(Ku.IR8,type="local")
Trans
#write results
Ku.IR8.Centrality<-data.frame(NAME=V(Ku.IR8)$name, UniprotID=V(Ku.IR8)$Uniprot, Degree=degree, Eigen= eigen_cent, Closeness= clos_cent, Betweeness= betw_cent, Hubs=hs, Authority=as,Transitivity=Trans)
Ku.IR8.Centrality
write.csv(Ku.IR8.Centrality,"./Ku.IR8.Centrality.csv")

#Cocitation, Number of shared interactions among nodes
Cocitation<-cocitation(Ku.IR8)
Cocitation#This matrix is similar to an adjacency matrix
write.csv(Cocitation,"./Ku.IR8.cocitation.csv")
#Plot with degree as node size
ggraph(Ku.IR8,layout=Lyt.FR[-isolated,])+
  geom_edge_fan(color="gray85", width=1, alpha=0.4)+
  geom_node_point(aes(color=V(Ku.IR8)$Annotation, size=degree),alpha=0.75)+
  geom_node_text(aes(label=name),color="gray25",size=4)+
  scale_color_manual(limits=V(Ku.IR8)$Annotation, values=V(Ku.IR8)$An.color)+
  scale_size_continuous(range = c(5, 15)) +
  theme_void()+
  title("Ku80 IR8 Network")+
  guides(color=guide_legend(override.aes = list(size=7)))+
  theme(legend.text=element_text(size=12,face = "bold"),
        legend.title = element_text(size=13,face = "bold"))


#Modularity Analysis

#Edge_Betweeness Modularity
Ku.IR8.ebmod=edge.betweenness.community(Ku.IR8,weights=E(Ku.IR8)$weight,directed = F) 
Ku.IR8.ebmod# GRAPH clustering edge betweenness, groups: 56, mod: 0.33

#
plot.igraph(Ku.IR8.ebmod,Ku.IR8,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(Ku.IR8)$An.color, vertex.frame.color=V(Ku.IR8)$An.color, vertex.shape="circle",
            vertex.size=5, vertex.label=V(Ku.IR8)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(Ku.IR8)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)

#
plot(Ku.IR8.ebmod, Ku.IR8,layout=Lyt.FR[-isolated,], alpha=0.6, 
     vertex.color=V(Ku.IR8)$An.color, vertex.frame.color=V(Ku.IR8)$An.color, vertex.shape="circle", 
     vertex.size=5, vertex.label=V(Ku.IR8)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black", 
     edge.width=E(Ku.IR8)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout", 
     rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)
#Dendogram of the Edge_Betweeness
dendPlot(Ku.IR8.ebmod, mode="hclust")

#Louvain modularity
Ku.IR8.Louvain.mod<-cluster_louvain(Ku.IR8,weights=E(Ku.IR8)$weight)
Ku.IR8.Louvain.mod #IGRAPH clustering multi level, groups: 12, mod: 0.38

#Plot the Louvain communities
plot(Ku.IR8.Louvain.mod, Ku.IR8,layout=Lyt.FR[-isolated,], alpha=0.6, 
     vertex.color=V(Ku.IR8)$An.color, vertex.frame.color=V(Ku.IR8)$An.color, vertex.shape="circle", 
     vertex.size=5, vertex.label=V(Ku.IR8)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black", 
     edge.width=E(Ku.IR8)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout", 
     rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#Generate HeatMaps of the adjacency matrixes
#Let's use the neighbor list for generating the Adjaency Matrix
SensorID_FNN #Let's use the edge list that contain the total amount of possible interactors for the three sensors
v<-levels(as.factor(unlist(SensorID_FNN[,1:2]))) #get vertex names from the network list
n<-length(v)#get the number of vertex in the network
e<-matrix(match(as.character(unlist(SensorID_FNN[,1:2])),v),ncol = 2)#Edgelist with row index instead of Names
w<-Ku_IR8$weight #Edge weights, only from the Ku_IR8 network to quit those edges that do not exist at Ku_IR8

# M is our co-expression network adjacency matrix.

# distringuish between i,j and j,i we want to make M symmetric (undirected) by
# considering an undirected network.
Ku.IR8.M <- matrix(0, n, n)                  # set up a co-expression matrix
Ku.IR8.M[e] <- w                             # fill it in with edge weights
Ku.IR8.M <- Ku.IR8.M + t(Ku.IR8.M)                         # make this symmetric
any(is.na(Ku.IR8.M))#check no NAs in the matrix
dimnames(Ku.IR8.M) <- list(v, v)             # label the vertices
# We let A_ij = 1 if there is an interaction. 
#A <- 1*(M > 0) #right now we keep the weights as they are
Ku.IR8.M #The adjacency Matrix not as an igraph object, easier to manipulate, and to subset for smaller Sensor specific networks
quantile(Ku.IR8.M,c(0,0.25,0.5,0.75,1))
#
#       0%   25%   50%   75%  100% 
#     0.000 0.000 0.000 0.000 0.991 

#Creating the HeatMap
# make a heatmap image with colors: white below 0,
# grey from 0.4 to 1 with darker closer to 1
breaks <- c(0, seq(0.4, 0.991, length=11)) # breaks in the color bins
cols <- grey(1-c(0,seq(0.5,1,length=10))) # colors

#Red Blue Colors
library(RColorBrewer)
breaks <- c(seq(-1, 1, length=11))
breaks
cols <- gray(1-c(0,seq(0,1,length=10))) # colors
cols2 <- hcl.colors(10 , palette = "Blues3" , alpha = 1 , rev = TRUE)
cols3 <- colorRamp("red" , space = "rgb" , interpolate = "linear")
HMCol<- brewer.pal(10,name = "RdBu")
HMCol<- rev(HMCol)
HMCol


# color the vertex for CCNB1 blue, its neighbors red, and the others white
SensorIDBind <- which(v == "Ku80") # index for ccnb1
#SensorIDBind<-as.factor(SensorIDBind)
vcols <- rep("white",n)
vcols[SensorIDBind] <- "blue"
#Subsetting Ku.IR8.M by doing a loop asking for look into every row and column the index in SensorIDBind worked
for (i in SensorIDBind) {
  vcols[which(Ku.IR8.M[,i] > 0 | Ku.IR8.M[i,])] <- "red"
}
##Color annotation for functional annotation
MatSet<-subset(Mat,Mat$NAME%in%v)# select only those proteins that are cointained in the node list of our matrixes
match(v,MatSet$NAME) # Check that the indexes of the Annotation Dataframe and the node list are in the same position 
match(rownames(Ku.IR8.M),MatSet$NAME) #Check that the row names of the Adj Matrix and the Annotation dataframe are in the same order
match(colnames(Ku.IR8.M),MatSet$NAME) #Check that the column names of the Adj Matrix and the Annotation dataframe are in the same order
AnCol<-MatSet$An.color#Vector of colors
AnCol<-as.character(AnCol)#transform it into a character vector, so the colors can be read

#Subsetting Ku.IR8.M using the index vector in SensorIDBind do not work
#vcols[which(Ku.IR8.M[,SensorIDBind] > 0 | Ku.IR8.M[SensorIDBind,])] <- "red"

#Generate a dendogram for keeping the proteins in the sample place
#Ku.IR8.dist<-dist(Ku.IR8.M,method = "euclidian")
#Ku.IR8.Dend<-hclust(Ku.IR8.dist, method = "ward.D" )
#plot(Ku.IR8.Dend)

#Cut the dendogram
#cut.height <- 10 # try numbers from 4 to 7
#FNN.hclusts<- cutree(FNN.Dend, cut.height)
#plot(FNN.hclusts)

# now actually make the heat map
par(mar=rep(0,4)) # make plot margins 0
heatmap(Ku.IR8.M, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=cols, breaks=breaks,frame=T, cexRow = 0.5 , cexCol = 0.5, main = "Ku80 IR8 AdjMat", ylab="Total Protein Interactors",xlab="Total Protein Interactors")
#
plot.new()
legend("topleft", c("Cell Cycle","Chromatin Org","DNA Repair","Nuclear Pore","p53","Proteasome","Ribosome Biogenesis","RNA Metabolism","Transcription Reg","UPR"), 
       fill=c("coral1","darkviolet","firebrick1","deepskyblue3","tomato","grey","lightgreen","gold","slateblue1","maroon"),  
       bty="n", inset=0, xpd=T,  border=F, cex =2,y.intersp = 0.7,x.intersp = 0.2,pt.cex = 1)


#Now we substract the KunonIR8 AdjMatrix to the KuIR6 AdjMatrix
KuIR8.AdjMat.Change<-Ku.IR8.M-Ku.nonIR5.M
range(KuIR8.AdjMat.Change)#Let's check the magnitude of the substracted weights
#[1]-0.813  0.668
any(is.na(KuIR8.AdjMat.Change))
#Generate the breaks and color functions
breaks <- c(seq(-0.5, 0.5, length=10))
breaks
HMCol<- brewer.pal(9,name = "RdBu")
HMCol<- rev(HMCol)
HMCol
[1] "#2166AC" "#4393C3" "#92C5DE" "#D1E5F0" "#F7F7F7" "#FDDBC7" "#F4A582" "#D6604D" "#B2182B"

#Plot the HeatMap
par(mar=rep(0,4)) # make plot margins 0
heatmap(KuIR8.AdjMat.Change, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=HMCol, breaks=breaks,frame=T, cexRow = 0.5 , cexCol = 0.5, main = "Ku80 IR8 AdjMat (minus nonIR5 edges)", ylab="Total Protein Interactors",xlab="Total Protein Interactors")


#Now we substract the KuIR5 AdjMatrix to the KuIR8 AdjMatrix
KuIR8.AdjMat.Change2<-Ku.IR8.M-Ku.IR2.M
range(KuIR8.AdjMat.Change)#Let's check the magnitude of the substracted weights
#[1]-0.813  0.668
any(is.na(KuIR8.AdjMat.Change))
#Generate the breaks and color functions
breaks <- c(seq(-0.5, 0.5, length=10))
breaks
HMCol<- brewer.pal(9,name = "RdBu")
HMCol<- rev(HMCol)
HMCol
[1] "#2166AC" "#4393C3" "#92C5DE" "#D1E5F0" "#F7F7F7" "#FDDBC7" "#F4A582" "#D6604D" "#B2182B"

#Plot the HeatMap
par(mar=rep(0,4)) # make plot margins 0
heatmap(KuIR8.AdjMat.Change2, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=HMCol, breaks=breaks,frame=T, cexRow = 0.5 , cexCol = 0.5, main = "Ku80 IR8 AdjMat (minus IR2 edges)", ylab="Total Protein Interactors",xlab="Total Protein Interactors")
#
#check with the FNN
# now actually make the heat map
par(mar=rep(0,4)) # make plot margins 0
heatmap(FNN.M, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=0:1, frame=T, cexRow = 0.5 , cexCol = 0.5)
#
legend("topleft", c("Neighbors Ku80","Ku80"), fill=c("red","blue"),  
       bty="n", inset=0, xpd=T,  border=F)


####################### Ku80 IR24h Networks
###################################################

#Ku80 IR24 Networks


#Let's keep the Mat as the Node attribute table and load the Neighbor list as the matrixes
Mat# Node attribute matrix for all the networks

#Ku80ID IR24 nodes Network Edge Lists
Ku_IR24<-read.delim("./Ku80_Files/Ku_IR24h.tsv",sep = "\t",col.names = )
str(Ku_IR24)
head(Ku_IR24)

#Generate the igraph network object
Ku.IR24<-graph_from_data_frame(d=Ku_IR24,vertices = Mat, directed = F )#Generate the igraph object with the edge and node annotation attributes
str(vertex_attr(Ku.IR24))#Check node attributes
str(edge_attr(Ku.IR24))#Check for edge attributes
Ku.IR24<-igraph::simplify(Ku.IR24, remove.multiple=TRUE, remove.loops=TRUE) #Remove repetitive edges and self-connecting edges
Ku.IR24
#Small plot
plot(Ku.IR24)

#Generate an adjacency Matrix that will be used after
Ku.IR24_AdjMat<-as_adjacency_matrix(Ku.IR24,type="both",attr = "weight",names=TRUE,sparse = TRUE)
Ku.IR24_AdjMat<-as.matrix(Ku.IR24_AdjMat)
write.csv(Ku.IR24_AdjMat,"./Ku.IR24_AdjMat.csv")

#Generate a vertex size based on degree and count
vsize<-igraph::degree(Ku.IR24)/vcount(Ku.IR24)
#Identify and remove disconnected nodes
isolated<-which(igraph::degree(Ku.IR24)==0)
isolated
length(isolated)# [1] 685
#remove isolated nodes
Ku.IR24<-igraph::delete_vertices(Ku.IR24,isolated)


#Generate layouts for the graph
#Do not run again, to not to change the layouts that I already have
Lyt.Nice<-layout_nicely(Ku.IR24,dim = 3)
Lyt.GOpt<-layout_with_graphopt(Ku.IR24, niter=5000, mass = 30, charge = E(Ku.IR24)$weight)
Lyt.DrL<-layout_with_drl(Ku.IR24, weights = E(Ku.IR24)$weight)
Lyt.FR<-layout.fruchterman.reingold(Ku.IR24, niter=10000, weights=E(Ku.IR24)$weight)#This one looks like the best option
Lyt.KK<-layout_with_kk(Ku.IR24, maxiter=10000, weights=E(Ku.IR24)$weight) #very large separation of peripheral nodes
Lyt.GEM<-layout_with_gem(Ku.IR24, maxiter=10000) #Bad separation
#Lyt.DH<-layout_with_dh(Ku.IR24, maxiter=150, fineiter=max(10, log2(833)))
Lyt.LGL<-layout.lgl(Ku.IR24, maxit=150,maxdelta=833,area=(length(E(Ku.IR24)))^2,coolexp=1.5,repulserad=((length(E(Ku.IR24)))^2)*(length(E(Ku.IR24))))#Do not work


#Star and Tree layouts

V(Ku.IR24)[name=="XRCC5"]
which(V(Ku.IR24)[name=="Ku80"])# 228 #Check what is the rowindex for your desired node
which(V(Ku.IR24)$name=="XRCC5")# 147 #Check what is the rowindex for your desired node
which(V(Ku.IR24)$name=="MRE11A")# 327
which(V(Ku.IR24)$name=="XRCC5")# 147
Lyt.Star<-layout_as_star(Ku.IR24,center=V(Ku.IR24)[name=="XRCC5"]) #it has a special order argument where we can assign an order to the vertex, hence maybe the OPLS-DA:VIP value order may be used as ordering vector

Lyt.Tree<-layout_as_tree(Ku.IR24,root=c(147))
Lyt.RT<-layout.reingold.tilford(Ku.IR24,root=c(147))# This one looks good, but the indexing do not allow to highlight the correct nodes
Lyt.sugi<-layout_with_sugiyama(Ku.IR24, maxiter=100, weights=E(Ku.IR242)$weight)#This one looks good for undirected tree layouts

#Check what layout looks fine

plot.igraph(Ku.IR24,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(Ku.IR24)$An.color, vertex.frame.color=V(Ku.IR24)$An.color, vertex.shape="circle",
            vertex.size=V(Ku.IR24)$Ku.FC.IR24, vertex.label=V(Ku.IR24)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(Ku.IR24)$weight, edge.lty=1, main="Ku80ID IR24 Network", sub="FR Tree Layout",
            rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#With ggraph

ggraph(Ku.IR24,layout=Lyt.FR[-isolated,])+
  geom_edge_fan(color="gray65", width=1, alpha=0.4)+
  geom_node_point(aes(color=V(Ku.IR24)$Annotation, size=V(Ku.IR24)$Ku.FC.IR24),alpha=0.75)+
  geom_node_text(aes(label=name),color="gray25",size=4)+
  scale_color_manual(limits=V(Ku.IR24)$Annotation, values=V(Ku.IR24)$An.color)+
  scale_size_continuous(range = c(5, 15)) +
  theme_void()+
  title("Ku80 IR24 Network")+
  guides(color=guide_legend(override.aes = list(size=7)))+
  theme(legend.text=element_text(size=12,face = "bold"),
        legend.title = element_text(size=13,face = "bold"))



#Community Detection
#clp <- igraph::cluster_optimal(Ku.IR24)
class(clp)
V(Ku.IR24)$community <- clp$membership
colrs <- adjustcolor( c("gray50", "tomato", "gold", "yellowgreen"), alpha=.6)


plot.igraph(clp,Ku.IR24,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(Ku.IR24)$An.color, vertex.frame.color=V(Ku.IR24)$An.color, vertex.shape="circle",
            vertex.size==V(Ku.IR24)$Ku.FC.IR24, vertex.label=V(Ku.IR24)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(Ku.IR24)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#Distances edges highlight
dist.from.Ku80 <- distances(Ku.IR24, v=V(Ku.IR24)[name=="XRCC5"],
                            to=V(Ku.IR24), weights=E(Ku.IR24)$weight)
dist.from.Ku80<-dist.from.Ku80[,!is.infinite(colSums(dist.from.Ku80))]

dist.from.Ku80
range(dist.from.Ku80)

#Set colors to plot the distances:
oranges <- colorRampPalette(c("dark red", "gold"))
col <- oranges(max(dist.from.Ku80))
col <- col[dist.from.Ku80+1]
plot(Ku.IR24, layout=Lyt.FR[-isolated,],vertex.color=col, vertex.label=dist.from.Ku80,
     vertex.label.color="white",vertex.shape="circle",
     vertex.size=V(Ku.IR24)$Ku.FC.IR24, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(Ku.IR24)$weight, edge.lty=1, main="Ku80 IR24 Full Nuclear Network", sub="FR Layout",
     rescale=T)


#Shortest path from Ku80 to Nucleolin
KuNcl.path <- shortest_paths(Ku.IR24,
                             from = V(Ku.IR24)[name=="XRCC5"],
                             to = V(Ku.IR24)[name=="NCL"],
                             output = "both") # both path nodes and edges
KuNcl.path
# Generate edge color variable to plot the path:
ecol <- rep("gray80", ecount(Ku.IR24))
ecol[unlist(KuNcl.path$epath)] <- "orange"
# Generate edge width variable to plot the path:
ew <- rep(2, ecount(Ku.IR24))
ew[unlist(KuNcl.path$epath)] <- 4
# Generate node color variable to plot the path:
vcol <- rep("gray40", vcount(Ku.IR24))
vcol[unlist(KuNcl.path$vpath)] <- "gold"
plot(Ku.IR24, layout=Lyt.FR[-isolated,], vertex.color=vcol, edge.color=ecol,
     edge.width=ew,vertex.label.color="white",vertex.shape="circle",
     vertex.size=5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(Ku.IR24)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)



#Plot based on Ku80 as incident
inc.edges <- incident(Ku.IR24, V(Ku.IR24)[name=="XRCC5"], mode="all")
# Set colors to plot the selected edges.
ecol <- rep("gray80", ecount(Ku.IR24))
ecol[inc.edges] <- "orange"
vcol <- rep("gray40", vcount(Ku.IR24))
vcol[V(Ku.IR24)$name=="Ku80"] <- "gold"

plot(Ku.IR24, layout=Lyt.FR[-isolated,], vertex.color=vcol, edge.color=ecol, 
     edge.width=ew,vertex.label.color="white",vertex.shape="circle",
     vertex.size=V(Ku.IR24)$Ku.FC.IR24, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(Ku.IR24)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)


#Plot Degree Distribution
Ku.IR24.Deg<-degree_distribution(Ku.IR24, cumulative = T, mode="all")
plot( x=0:max(igraph::degree(Ku.IR24)), y=1-Ku.IR24.Deg, pch=19, cex=1.2, col="orange",
      xlab="Degree", ylab="Cumulative Frequency")

#Network Descriptives
#Number of nodes and edges
vcount(Ku.IR24) #[1] 148
ecount(Ku.IR24)#[1]  758
#Density (Number of edges in relation to the total number of edges in a full connected graph)
edge_density(Ku.IR24,loops = F) #[1] 0.06968193
#Transitivity (Total number of conected triangles)
transitivity(Ku.IR24,type="global") #[1] 0.5731848
#modularity
mod<-V(Ku.IR24)$S6.FuzzyCMeans #lets check for fuzzy Cmeans modularity
mod<-mod+1#The non clustered interactors that appears as 0, needs to be transformed to 1, so we add 1 to every cluster
modularity(Ku.IR24,membership = mod)
#Ku.FuzzyCMeans modularity [1]-0.01577805
#Ku.FuzzyCMeans modularity [1] -0.003096261
#Mre.FuzzyCMeans modularity [1] 0.009537667

#Assorativity
assortativity_nominal(Ku.IR24,factor(V(Ku.IR24)$Annotation),directed=F)
#Annotation [1] [1] 0.2038583
#Ku.FuzzyCMeans [1] -0.02507076
#Ku.FuzzyCMeans [1] -0.00570132
#Mre.FuzzyCMeans [1] 0.0123745

#Diameter (longest geodesic distance), shortest path between two most peripheral nodes
diameter(Ku.IR24,directed = F, weights = E(Ku.IR24)$weights)#[1] 4.799
diam<-get_diameter(Ku.IR24,directed = F,weights = E(Ku.IR24)$weights)
diam #+ + 8/148 vertices, named, from 1258251: [1] CCT4  TUBB  RPL35 RBM39 SF3B3 LMNA  LMNB1 VIM
#Plot Diameter
vcol <- rep("gray50", vcount(Ku.IR24))
vcol[diam] <- "gold"
ecol <- rep("gray90", ecount(Ku.IR24))
ecol[E(Ku.IR24, path=diam)] <- "orange"
# E(net, path=diam) finds edges along a path, here 'diam'
plot(Ku.IR24, vertex.color=vcol, edge.color=ecol, vertex.label.color="black",edge.arrow.mode=0,layout=Lyt.FR[-isolated,],vertex.size=V(Ku.IR24)$Ku.FC.IR24)

#Node topological descriptives
degree<-igraph::centr_degree(Ku.IR24,mode = "all",normalized = T)
degree<-degree$res
eigen_cent<-igraph::eigen_centrality(Ku.IR24)
eigen_cent<-eigen_cent$vector
clos_cent <- igraph::closeness(Ku.IR24, normalized=T)
betw_cent <- igraph::betweenness(Ku.IR24, normalized=T)
hs <- hub_score(Ku.IR24, weights=E(Ku.IR24)$weight)$vector
as <- authority_score(Ku.IR24, weights=E(Ku.IR24)$weight)$vector
Trans<-transitivity(Ku.IR24,type="local")
Trans
#write results
Ku.IR24.Centrality<-data.frame(NAME=V(Ku.IR24)$name, UniprotID=V(Ku.IR24)$Uniprot, Degree=degree, Eigen= eigen_cent, Closeness= clos_cent, Betweeness= betw_cent, Hubs=hs, Authority=as,Transitivity=Trans)
Ku.IR24.Centrality
write.csv(Ku.IR24.Centrality,"./Ku.IR24.Centrality.csv")

#Cocitation, Number of shared interactions among nodes
Cocitation<-cocitation(Ku.IR24)
Cocitation#This matrix is similar to an adjacency matrix
write.csv(Cocitation,"./Ku.IR24.cocitation.csv")
#Plot with degree as node size
ggraph(Ku.IR24,layout=Lyt.FR[-isolated,])+
  geom_edge_fan(color="gray85", width=1, alpha=0.4)+
  geom_node_point(aes(color=V(Ku.IR24)$Annotation, size=degree),alpha=0.75)+
  geom_node_text(aes(label=name),color="gray25",size=4)+
  scale_color_manual(limits=V(Ku.IR24)$Annotation, values=V(Ku.IR24)$An.color)+
  scale_size_continuous(range = c(5, 15)) +
  theme_void()+
  title("Ku80 IR24 Network")+
  guides(color=guide_legend(override.aes = list(size=7)))+
  theme(legend.text=element_text(size=12,face = "bold"),
        legend.title = element_text(size=13,face = "bold"))


#Modularity Analysis

#Edge_Betweeness Modularity
Ku.IR24.ebmod=edge.betweenness.community(Ku.IR24,weights=E(Ku.IR24)$weight,directed = F) 
Ku.IR24.ebmod# IGRAPH clustering edge betweenness, groups: 16, mod: 0.48

#
plot.igraph(Ku.IR24.ebmod,Ku.IR24,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(Ku.IR24)$An.color, vertex.frame.color=V(Ku.IR24)$An.color, vertex.shape="circle",
            vertex.size=5, vertex.label=V(Ku.IR24)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(Ku.IR24)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)

#
plot(Ku.IR24.ebmod, Ku.IR24,layout=Lyt.FR[-isolated,], alpha=0.6, 
     vertex.color=V(Ku.IR24)$An.color, vertex.frame.color=V(Ku.IR24)$An.color, vertex.shape="circle", 
     vertex.size=5, vertex.label=V(Ku.IR24)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black", 
     edge.width=E(Ku.IR24)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout", 
     rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)
#Dendogram of the Edge_Betweeness
dendPlot(Ku.IR24.ebmod, mode="hclust")

#Louvain modularity
Ku.IR24.Louvain.mod<-cluster_louvain(Ku.IR24,weights=E(Ku.IR24)$weight)
Ku.IR24.Louvain.mod #IGRAPH clustering multi level, groups: 12, mod: 0.52

#Plot the Louvain communities
plot(Ku.IR24.Louvain.mod, Ku.IR24,layout=Lyt.FR[-isolated,], alpha=0.6, 
     vertex.color=V(Ku.IR24)$An.color, vertex.frame.color=V(Ku.IR24)$An.color, vertex.shape="circle", 
     vertex.size=5, vertex.label=V(Ku.IR24)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black", 
     edge.width=E(Ku.IR24)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout", 
     rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#Generate HeatMaps of the adjacency matrixes
#Let's use the neighbor list for generating the Adjaency Matrix
SensorID_FNN #Let's use the edge list that contain the total amount of possible interactors for the three sensors
v<-levels(as.factor(unlist(SensorID_FNN[,1:2]))) #get vertex names from the network list
n<-length(v)#get the number of vertex in the network
e<-matrix(match(as.character(unlist(SensorID_FNN[,1:2])),v),ncol = 2)#Edgelist with row index instead of Names
w<-Ku_IR24$weight #Edge weights, only from the Ku_IR24 network to quit those edges that do not exist at Ku_IR24

# M is our co-expression network adjacency matrix.

# distringuish between i,j and j,i we want to make M symmetric (undirected) by
# considering an undirected network.
Ku.IR24.M <- matrix(0, n, n)                  # set up a co-expression matrix
Ku.IR24.M[e] <- w                             # fill it in with edge weights
Ku.IR24.M <- Ku.IR24.M + t(Ku.IR24.M)                         # make this symmetric
any(is.na(Ku.IR24.M))#check no NAs in the matrix
dimnames(Ku.IR24.M) <- list(v, v)             # label the vertices
# We let A_ij = 1 if there is an interaction. 
#A <- 1*(M > 0) #right now we keep the weights as they are
Ku.IR24.M #The adjacency Matrix not as an igraph object, easier to manipulate, and to subset for smaller Sensor specific networks
quantile(Ku.IR24.M,c(0,0.25,0.5,0.75,1))
#
#      0%  25%  50%  75% 100% 
#    0.00 0.00 0.00 0.00 0.99   

#Creating the HeatMap
# make a heatmap image with colors: white below 0,
# grey from 0.4 to 1 with darker closer to 1
breaks <- c(0, seq(0.4, 0.991, length=11)) # breaks in the color bins
cols <- grey(1-c(0,seq(0.5,1,length=10))) # colors

#Red Blue Colors
library(RColorBrewer)
breaks <- c(seq(-1, 1, length=11))
breaks
cols <- gray(1-c(0,seq(0,1,length=10))) # colors
cols2 <- hcl.colors(10 , palette = "Blues3" , alpha = 1 , rev = TRUE)
cols3 <- colorRamp("red" , space = "rgb" , interpolate = "linear")
HMCol<- brewer.pal(10,name = "RdBu")
HMCol<- rev(HMCol)
HMCol


# color the vertex for CCNB1 blue, its neighbors red, and the others white
SensorIDBind <- which(v == "Ku80") # index for ccnb1
#SensorIDBind<-as.factor(SensorIDBind)
vcols <- rep("white",n)
vcols[SensorIDBind] <- "blue"
#Subsetting Ku.IR24.M by doing a loop asking for look into every row and column the index in SensorIDBind worked
for (i in SensorIDBind) {
  vcols[which(Ku.IR24.M[,i] > 0 | Ku.IR24.M[i,])] <- "red"
}
##Color annotation for functional annotation
MatSet<-subset(Mat,Mat$NAME%in%v)# select only those proteins that are cointained in the node list of our matrixes
match(v,MatSet$NAME) # Check that the indexes of the Annotation Dataframe and the node list are in the same position 
match(rownames(Ku.IR24.M),MatSet$NAME) #Check that the row names of the Adj Matrix and the Annotation dataframe are in the same order
match(colnames(Ku.IR24.M),MatSet$NAME) #Check that the column names of the Adj Matrix and the Annotation dataframe are in the same order
AnCol<-MatSet$An.color#Vector of colors
AnCol<-as.character(AnCol)#transform it into a character vector, so the colors can be read

#Subsetting Ku.IR24.M using the index vector in SensorIDBind do not work
#vcols[which(Ku.IR24.M[,SensorIDBind] > 0 | Ku.IR24.M[SensorIDBind,])] <- "red"

#Generate a dendogram for keeping the proteins in the sample place
#Ku.IR24.dist<-dist(Ku.IR24.M,method = "euclidian")
#Ku.IR24.Dend<-hclust(Ku.IR24.dist, method = "ward.D" )
#plot(Ku.IR24.Dend)

#Cut the dendogram
#cut.height <- 10 # try numbers from 4 to 7
#FNN.hclusts<- cutree(FNN.Dend, cut.height)
#plot(FNN.hclusts)

# now actually make the heat map
par(mar=rep(0,4)) # make plot margins 0
heatmap(Ku.IR24.M, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=cols, breaks=breaks,frame=T, cexRow = 0.5 , cexCol = 0.5, main = "Ku80 IR24 AdjMat", ylab="Total Protein Interactors",xlab="Total Protein Interactors")
#
plot.new()
legend("topleft", c("Cell Cycle","Chromatin Org","DNA Repair","Nuclear Pore","p53","Proteasome","Ribosome Biogenesis","RNA Metabolism","Transcription Reg","UPR"), 
       fill=c("coral1","darkviolet","firebrick1","deepskyblue3","tomato","grey","lightgreen","gold","slateblue1","maroon"),  
       bty="n", inset=0, xpd=T,  border=F, cex =2,y.intersp = 0.7,x.intersp = 0.2,pt.cex = 1)


#Now we substract the KunonIR24 AdjMatrix to the KuIR6 AdjMatrix
KuIR24.AdjMat.Change<-Ku.IR24.M-Ku.nonIR5.M
range(KuIR24.AdjMat.Change)#Let's check the magnitude of the substracted weights
#[1] -0.806  0.805[1] -0.806  0.805
any(is.na(KuIR24.AdjMat.Change))
#Generate the breaks and color functions
breaks <- c(seq(-0.5, 0.5, length=10))
breaks
HMCol<- brewer.pal(9,name = "RdBu")
HMCol<- rev(HMCol)
HMCol
[1] "#2166AC" "#4393C3" "#92C5DE" "#D1E5F0" "#F7F7F7" "#FDDBC7" "#F4A582" "#D6604D" "#B2182B"

#Plot the HeatMap
par(mar=rep(0,4)) # make plot margins 0
heatmap(KuIR24.AdjMat.Change, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=HMCol, breaks=breaks,frame=T, cexRow = 0.5 , cexCol = 0.5, main = "Ku80 IR24 AdjMat (minus nonIR24 edges)", ylab="Total Protein Interactors",xlab="Total Protein Interactors")


#Now we substract the KuIR5 AdjMatrix to the KuIR24 AdjMatrix
KuIR24.AdjMat.Change2<-Ku.IR24.M-Ku.IR8.M
range(KuIR24.AdjMat.Change)#Let's check the magnitude of the substracted weights
#[1]-0.813  0.668
any(is.na(KuIR24.AdjMat.Change))
#Generate the breaks and color functions
breaks <- c(seq(-0.5, 0.5, length=10))
breaks
HMCol<- brewer.pal(9,name = "RdBu")
HMCol<- rev(HMCol)
HMCol
[1] "#2166AC" "#4393C3" "#92C5DE" "#D1E5F0" "#F7F7F7" "#FDDBC7" "#F4A582" "#D6604D" "#B2182B"

#Plot the HeatMap
par(mar=rep(0,4)) # make plot margins 0
heatmap(KuIR24.AdjMat.Change2, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=HMCol, breaks=breaks,frame=T, cexRow = 0.5 , cexCol = 0.5, main = "Ku80 IR24 AdjMat (minus IR8h edges)", ylab="Total Protein Interactors",xlab="Total Protein Interactors")

#
#check with the FNN
# now actually make the heat map
par(mar=rep(0,4)) # make plot margins 0
heatmap(FNN.M, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=0:1, frame=T, cexRow = 0.5 , cexCol = 0.5)
#
legend("topleft", c("Neighbors Ku80","Ku80"), fill=c("red","blue"),  
       bty="n", inset=0, xpd=T,  border=F)


####################### Ku80 nonIR24h Networks
###################################################

#Ku80 nonIR24 Networks


#Let's keep the Mat as the Node attribute table and load the Neighbor list as the matrixes
Mat# Node attribute matrix for all the networks

#Ku80ID nonIR24 nodes Network Edge Lists
Ku_nonIR24<-read.delim("./Ku80_Files/Ku_nonIR24h.tsv",sep = "\t",col.names = )
str(Ku_nonIR24)
head(Ku_nonIR24)

#Generate the igraph network object
Ku.nonIR24<-graph_from_data_frame(d=Ku_nonIR24,vertices = Mat, directed = F )#Generate the igraph object with the edge and node annotation attributes
str(vertex_attr(Ku.nonIR24))#Check node attributes
str(edge_attr(Ku.nonIR24))#Check for edge attributes
Ku.nonIR24<-igraph::simplify(Ku.nonIR24, remove.multiple=TRUE, remove.loops=TRUE) #Remove repetitive edges and self-connecting edges
Ku.nonIR24
#Small plot
plot(Ku.nonIR24)

#Generate an adjacency Matrix that will be used after
Ku.nonIR24_AdjMat<-as_adjacency_matrix(Ku.nonIR24,type="both",attr = "weight",names=TRUE,sparse = TRUE)
Ku.nonIR24_AdjMat<-as.matrix(Ku.nonIR24_AdjMat)
write.csv(Ku.nonIR24_AdjMat,"./Ku.nonIR24_AdjMat.csv")

#Generate a vertex size based on degree and count
vsize<-igraph::degree(Ku.nonIR24)/vcount(Ku.nonIR24)
#Identify and remove disconnected nodes
isolated<-which(igraph::degree(Ku.nonIR24)==0)
isolated
length(isolated)# [1] 673
#remove isolated nodes
Ku.nonIR24<-igraph::delete_vertices(Ku.nonIR24,isolated)


#Generate layouts for the graph
#Do not run again, to not to change the layouts that I already have
Lyt.Nice<-layout_nicely(Ku.nonIR24,dim = 3)
Lyt.GOpt<-layout_with_graphopt(Ku.nonIR24, niter=5000, mass = 30, charge = E(Ku.nonIR24)$weight)
Lyt.DrL<-layout_with_drl(Ku.nonIR24, weights = E(Ku.nonIR24)$weight)
Lyt.FR<-layout.fruchterman.reingold(Ku.nonIR24, niter=10000, weights=E(Ku.nonIR24)$weight)#This one looks like the best option
Lyt.KK<-layout_with_kk(Ku.nonIR24, maxiter=10000, weights=E(Ku.nonIR24)$weight) #very large separation of peripheral nodes
Lyt.GEM<-layout_with_gem(Ku.nonIR24, maxiter=10000) #Bad separation
#Lyt.DH<-layout_with_dh(Ku.nonIR24, maxiter=150, fineiter=max(10, log2(833)))
Lyt.LGL<-layout.lgl(Ku.nonIR24, maxit=150,maxdelta=833,area=(length(E(Ku.nonIR24)))^2,coolexp=1.5,repulserad=((length(E(Ku.nonIR24)))^2)*(length(E(Ku.nonIR24))))#Do not work


#Star and Tree layouts

V(Ku.nonIR24)[name=="XRCC5"]
which(V(Ku.nonIR24)[name=="XRCC5"])# 228 #Check what is the rowindex for your desired node
which(V(Ku.nonIR24)$name=="XRCC5")# 176 #Check what is the rowindex for your desired node
which(V(Ku.nonIR24)$name=="MRE11A")# 327
which(V(Ku.nonIR24)$name=="XRCC5")# 722
Lyt.Star<-layout_as_star(Ku.nonIR24,center=V(Ku.nonIR24)[name=="XRCC5"]) #it has a special order argument where we can assign an order to the vertex, hence maybe the OPLS-DA:VIP value order may be used as ordering vector

Lyt.Tree<-layout_as_tree(Ku.nonIR24,root=c(159))
Lyt.RT<-layout.reingold.tilford(Ku.nonIR24,root=c(159))# This one looks good, but the indexing do not allow to highlight the correct nodes
Lyt.sugi<-layout_with_sugiyama(Ku.nonIR24, maxiter=100, weights=E(Ku.nonIR242)$weight)#This one looks good for undirected tree layouts

#Check what layout looks fine

plot.igraph(Ku.nonIR24,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(Ku.nonIR24)$An.color, vertex.frame.color=V(Ku.nonIR24)$An.color, vertex.shape="circle",
            vertex.size=V(Ku.nonIR24)$Ku.FC.nonIR24, vertex.label=V(Ku.nonIR24)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(Ku.nonIR24)$weight, edge.lty=1, main="Ku80ID nonIR24 Network", sub="FR Tree Layout",
            rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#With ggraph

ggraph(Ku.nonIR24,layout=Lyt.FR[-isolated,])+
  geom_edge_fan(color="gray65", width=1, alpha=0.4)+
  geom_node_point(aes(color=V(Ku.nonIR24)$Annotation, size=V(Ku.nonIR24)$Ku.FC.nonIR24),alpha=0.75)+
  geom_node_text(aes(label=name),color="gray25",size=4)+
  scale_color_manual(limits=V(Ku.nonIR24)$Annotation, values=V(Ku.nonIR24)$An.color)+
  scale_size_continuous(range = c(5, 15)) +
  theme_void()+
  title("Ku80 nonIR24 Network")+
  guides(color=guide_legend(override.aes = list(size=7)))+
  theme(legend.text=element_text(size=12,face = "bold"),
        legend.title = element_text(size=13,face = "bold"))



#Community Detection
#clp <- igraph::cluster_optimal(Ku.nonIR24)
class(clp)
V(Ku.nonIR24)$community <- clp$membership
colrs <- adjustcolor( c("gray50", "tomato", "gold", "yellowgreen"), alpha=.6)


plot.igraph(clp,Ku.nonIR24,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(Ku.nonIR24)$An.color, vertex.frame.color=V(Ku.nonIR24)$An.color, vertex.shape="circle",
            vertex.size==V(Ku.nonIR24)$Ku.FC.nonIR24, vertex.label=V(Ku.nonIR24)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(Ku.nonIR24)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#Distances edges highlight
dist.from.Ku80 <- distances(Ku.nonIR24, v=V(Ku.nonIR24)[name=="XRCC5"],
                            to=V(Ku.nonIR24), weights=E(Ku.nonIR24)$weight)
dist.from.Ku80<-dist.from.Ku80[,!is.infinite(colSums(dist.from.Ku80))]

dist.from.Ku80
range(dist.from.Ku80)

#Set colors to plot the distances:
oranges <- colorRampPalette(c("dark red", "gold"))
col <- oranges(max(dist.from.Ku80))
col <- col[dist.from.Ku80+1]
plot(Ku.nonIR24, layout=Lyt.FR[-isolated,],vertex.color=col, vertex.label=dist.from.Ku80,
     vertex.label.color="white",vertex.shape="circle",
     vertex.size=V(Ku.nonIR24)$Ku.FC.nonIR24, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(Ku.nonIR24)$weight, edge.lty=1, main="Ku80 nonIR24 Full Nuclear Network", sub="FR Layout",
     rescale=T)


#Shortest path from Ku80 to Nucleolin
KuNcl.path <- shortest_paths(Ku.nonIR24,
                             from = V(Ku.nonIR24)[name=="XRCC5"],
                             to = V(Ku.nonIR24)[name=="NCL"],
                             output = "both") # both path nodes and edges
KuNcl.path
# Generate edge color variable to plot the path:
ecol <- rep("gray80", ecount(Ku.nonIR24))
ecol[unlist(KuNcl.path$epath)] <- "orange"
# Generate edge width variable to plot the path:
ew <- rep(2, ecount(Ku.nonIR24))
ew[unlist(KuNcl.path$epath)] <- 4
# Generate node color variable to plot the path:
vcol <- rep("gray40", vcount(Ku.nonIR24))
vcol[unlist(KuNcl.path$vpath)] <- "gold"
plot(Ku.nonIR24, layout=Lyt.FR[-isolated,], vertex.color=vcol, edge.color=ecol,
     edge.width=ew,vertex.label.color="white",vertex.shape="circle",
     vertex.size=5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(Ku.nonIR24)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)



#Plot based on Ku80 as incident
inc.edges <- incident(Ku.nonIR24, V(Ku.nonIR24)[name=="XRCC5"], mode="all")
# Set colors to plot the selected edges.
ecol <- rep("gray80", ecount(Ku.nonIR24))
ecol[inc.edges] <- "orange"
vcol <- rep("gray40", vcount(Ku.nonIR24))
vcol[V(Ku.nonIR24)$name=="Ku80"] <- "gold"

plot(Ku.nonIR24, layout=Lyt.FR[-isolated,], vertex.color=vcol, edge.color=ecol, 
     edge.width=ew,vertex.label.color="white",vertex.shape="circle",
     vertex.size=V(Ku.nonIR24)$Ku.FC.nonIR24, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(Ku.nonIR24)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)


#Plot Degree Distribution
Ku.nonIR24.Deg<-degree_distribution(Ku.nonIR24, cumulative = T, mode="all")
plot( x=0:max(igraph::degree(Ku.nonIR24)), y=1-Ku.nonIR24.Deg, pch=19, cex=1.2, col="orange",
      xlab="Degree", ylab="Cumulative Frequency")

#Network Descriptives
#Number of nodes and edges
vcount(Ku.nonIR24) #[1] 160
ecount(Ku.nonIR24)#[1] 913
#Density (Number of edges in relation to the total number of edges in a full connected graph)
edge_density(Ku.nonIR24,loops = F) #[1] 0.07177673
#Transitivity (Total number of conected triangles)
transitivity(Ku.nonIR24,type="global") #[1] 0.6277133
#modularity
mod<-V(Ku.nonIR24)$Mre.FuzzyCMeans #lets check for fuzzy Cmeans modularity
mod<-mod+1#The non clustered interactors that appears as 0, needs to be transformed to 1, so we add 1 to every cluster
modularity(Ku.nonIR24,membership = mod)
#Ku.FuzzyCMeans modularity [1]  0.009224191
#Ku.FuzzyCMeans modularity [1] 0.00604569
#Mre.FuzzyCMeans modularity [1] 0.01608025

#Assorativity
assortativity_nominal(Ku.nonIR24,factor(V(Ku.nonIR24)$Mre.FuzzyCMeans),directed=F)
#Annotation [1] [1] 0.1663444
#S6.FuzzyCMeans [1]0.01476397
#Ku.FuzzyCMeans [1] 0.01127578
#Mre.FuzzyCMeans [1] 0.02028693

#Diameter (longest geodesic distance), shortest path between two most peripheral nodes
diameter(Ku.nonIR24,directed = F, weights = E(Ku.nonIR24)$weights)#[1] 4.241
diam<-get_diameter(Ku.nonIR24,directed = F,weights = E(Ku.nonIR24)$weights)
diam #+ 9/160 vertices, named, from 6cd078b: [1] PSMA1    PSMA2    HSP90AB1 PHB      RBM39    SF3B3    LMNA     LMNB1    VIM   
#Plot Diameter
vcol <- rep("gray50", vcount(Ku.nonIR24))
vcol[diam] <- "gold"
ecol <- rep("gray90", ecount(Ku.nonIR24))
ecol[E(Ku.nonIR24, path=diam)] <- "orange"
# E(net, path=diam) finds edges along a path, here 'diam'
plot(Ku.nonIR24, vertex.color=vcol, edge.color=ecol, vertex.label.color="black",edge.arrow.mode=0,layout=Lyt.FR[-isolated,],vertex.size=V(Ku.nonIR24)$Ku.FC.nonIR24)

#Node topological descriptives
degree<-igraph::centr_degree(Ku.nonIR24,mode = "all",normalized = T)
degree<-degree$res
eigen_cent<-igraph::eigen_centrality(Ku.nonIR24)
eigen_cent<-eigen_cent$vector
clos_cent <- igraph::closeness(Ku.nonIR24, normalized=T)
betw_cent <- igraph::betweenness(Ku.nonIR24, normalized=T)
hs <- hub_score(Ku.nonIR24, weights=E(Ku.nonIR24)$weight)$vector
as <- authority_score(Ku.nonIR24, weights=E(Ku.nonIR24)$weight)$vector
Trans<-transitivity(Ku.nonIR24,type="local")
Trans
#write results
Ku.nonIR24.Centrality<-data.frame(NAME=V(Ku.nonIR24)$name, UniprotID=V(Ku.nonIR24)$Uniprot, Degree=degree, Eigen= eigen_cent, Closeness= clos_cent, Betweeness= betw_cent, Hubs=hs, Authority=as,Transitivity=Trans)
Ku.nonIR24.Centrality
write.csv(Ku.nonIR24.Centrality,"./Ku.nonIR24.Centrality.csv")

#Cocitation, Number of shared interactions among nodes
Cocitation<-cocitation(Ku.nonIR24)
Cocitation#This matrix is similar to an adjacency matrix
write.csv(Cocitation,"./Ku.nonIR24.cocitation.csv")
#Plot with degree as node size
ggraph(Ku.nonIR24,layout=Lyt.FR[-isolated,])+
  geom_edge_fan(color="gray85", width=1, alpha=0.4)+
  geom_node_point(aes(color=V(Ku.nonIR24)$Annotation, size=degree),alpha=0.75)+
  geom_node_text(aes(label=name),color="gray25",size=4)+
  scale_color_manual(limits=V(Ku.nonIR24)$Annotation, values=V(Ku.nonIR24)$An.color)+
  scale_size_continuous(range = c(5, 15)) +
  theme_void()+
  title("Ku80 nonIR24 Network")+
  guides(color=guide_legend(override.aes = list(size=7)))+
  theme(legend.text=element_text(size=12,face = "bold"),
        legend.title = element_text(size=13,face = "bold"))


#Modularity Analysis

#Edge_Betweeness Modularity
Ku.nonIR24.ebmod=edge.betweenness.community(Ku.nonIR24,weights=E(Ku.nonIR24)$weight,directed = F) 
Ku.nonIR24.ebmod# IGRAPH clustering edge betweenness, groups: 9, mod: 0.41

#
plot.igraph(Ku.nonIR24.ebmod,Ku.nonIR24,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(Ku.nonIR24)$An.color, vertex.frame.color=V(Ku.nonIR24)$An.color, vertex.shape="circle",
            vertex.size=5, vertex.label=V(Ku.nonIR24)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(Ku.nonIR24)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)

#
plot(Ku.nonIR24.ebmod, Ku.nonIR24,layout=Lyt.FR[-isolated,], alpha=0.6, 
     vertex.color=V(Ku.nonIR24)$An.color, vertex.frame.color=V(Ku.nonIR24)$An.color, vertex.shape="circle", 
     vertex.size=5, vertex.label=V(Ku.nonIR24)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black", 
     edge.width=E(Ku.nonIR24)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout", 
     rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)
#Dendogram of the Edge_Betweeness
dendPlot(Ku.nonIR24.ebmod, mode="hclust")

#Louvain modularity
Ku.nonIR24.Louvain.mod<-cluster_louvain(Ku.nonIR24,weights=E(Ku.nonIR24)$weight)
Ku.nonIR24.Louvain.mod #IGRAPH clustering multi level, groups: 9, mod: 0.45

#Plot the Louvain communities
plot(Ku.nonIR24.Louvain.mod, Ku.nonIR24,layout=Lyt.FR[-isolated,], alpha=0.6, 
     vertex.color=V(Ku.nonIR24)$An.color, vertex.frame.color=V(Ku.nonIR24)$An.color, vertex.shape="circle", 
     vertex.size=5, vertex.label=V(Ku.nonIR24)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black", 
     edge.width=E(Ku.nonIR24)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout", 
     rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#Generate HeatMaps of the adjacency matrixes
#Let's use the neighbor list for generating the Adjaency Matrix
SensorID_FNN #Let's use the edge list that contain the total amount of possible interactors for the three sensors
v<-levels(as.factor(unlist(SensorID_FNN[,1:2]))) #get vertex names from the network list
n<-length(v)#get the number of vertex in the network
e<-matrix(match(as.character(unlist(SensorID_FNN[,1:2])),v),ncol = 2)#Edgelist with row index instead of Names
w<-Ku_nonIR24$weight #Edge weights, only from the Ku_nonIR24 network to quit those edges that do not exist at Ku_nonIR24

# M is our co-expression network adjacency matrix.

# distringuish between i,j and j,i we want to make M symmetric (undirected) by
# considering an undirected network.
Ku.nonIR24.M <- matrix(0, n, n)                  # set up a co-expression matrix
Ku.nonIR24.M[e] <- w                             # fill it in with edge weights
Ku.nonIR24.M <- Ku.nonIR24.M + t(Ku.nonIR24.M)                         # make this symmetric
any(is.na(Ku.nonIR24.M))#check no NAs in the matrix
dimnames(Ku.nonIR24.M) <- list(v, v)             # label the vertices
# We let A_ij = 1 if there is an interaction. 
#A <- 1*(M > 0) #right now we keep the weights as they are
Ku.nonIR24.M #The adjacency Matrix not as an igraph object, easier to manipulate, and to subset for smaller Sensor specific networks
quantile(Ku.nonIR24.M,c(0,0.25,0.5,0.75,1))
#
#      0%  25%  50%  75% 100% 
#    0.00 0.00 0.00 0.00 0.99   

#Creating the HeatMap
# make a heatmap image with colors: white below 0,
# grey from 0.4 to 1 with darker closer to 1
breaks <- c(0, seq(0.4, 0.991, length=11)) # breaks in the color bins
cols <- grey(1-c(0,seq(0.5,1,length=10))) # colors

#Red Blue Colors
library(RColorBrewer)
breaks <- c(seq(-1, 1, length=11))
breaks
cols <- gray(1-c(0,seq(0,1,length=10))) # colors
cols2 <- hcl.colors(10 , palette = "Blues3" , alpha = 1 , rev = TRUE)
cols3 <- colorRamp("red" , space = "rgb" , interpolate = "linear")
HMCol<- brewer.pal(10,name = "RdBu")
HMCol<- rev(HMCol)
HMCol


# color the vertex for CCNB1 blue, its neighbors red, and the others white
SensorIDBind <- which(v == "Ku80") # index for ccnb1
#SensorIDBind<-as.factor(SensorIDBind)
vcols <- rep("white",n)
vcols[SensorIDBind] <- "blue"
#Subsetting Ku.nonIR24.M by doing a loop asking for look into every row and column the index in SensorIDBind worked
for (i in SensorIDBind) {
  vcols[which(Ku.nonIR24.M[,i] > 0 | Ku.nonIR24.M[i,])] <- "red"
}
##Color annotation for functional annotation
MatSet<-subset(Mat,Mat$NAME%in%v)# select only those proteins that are cointained in the node list of our matrixes
match(v,MatSet$NAME) # Check that the indexes of the Annotation Dataframe and the node list are in the same position 
match(rownames(Ku.nonIR24.M),MatSet$NAME) #Check that the row names of the Adj Matrix and the Annotation dataframe are in the same order
match(colnames(Ku.nonIR24.M),MatSet$NAME) #Check that the column names of the Adj Matrix and the Annotation dataframe are in the same order
AnCol<-MatSet$An.color#Vector of colors
AnCol<-as.character(AnCol)#transform it into a character vector, so the colors can be read

#Subsetting Ku.nonIR24.M using the index vector in SensorIDBind do not work
#vcols[which(Ku.nonIR24.M[,SensorIDBind] > 0 | Ku.nonIR24.M[SensorIDBind,])] <- "red"

#Generate a dendogram for keeping the proteins in the sample place
#Ku.nonIR24.dist<-dist(Ku.nonIR24.M,method = "euclidian")
#Ku.nonIR24.Dend<-hclust(Ku.nonIR24.dist, method = "ward.D" )
#plot(Ku.nonIR24.Dend)

#Cut the dendogram
#cut.height <- 10 # try numbers from 4 to 7
#FNN.hclusts<- cutree(FNN.Dend, cut.height)
#plot(FNN.hclusts)

# now actually make the heat map
par(mar=rep(0,4)) # make plot margins 0
heatmap(Ku.nonIR24.M, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=cols, breaks=breaks,frame=T, cexRow = 0.5 , cexCol = 0.5, main = "Ku80 nonIR24 AdjMat", ylab="Total Protein Interactors",xlab="Total Protein Interactors")
#
plot.new()
legend("topleft", c("Cell Cycle","Chromatin Org","DNA Repair","Nuclear Pore","p53","Proteasome","Ribosome Biogenesis","RNA Metabolism","Transcription Reg","UPR"), 
       fill=c("coral1","darkviolet","firebrick1","deepskyblue3","tomato","grey","lightgreen","gold","slateblue1","maroon"),  
       bty="n", inset=0, xpd=T,  border=F, cex =2,y.intersp = 0.7,x.intersp = 0.2,pt.cex = 1)


#Now we substract the KunonnonIR24 AdjMatrix to the KuIR6 AdjMatrix
KunonIR24.AdjMat.Change<-Ku.nonIR24.M-Ku.nonIR5.M
range(KunonIR24.AdjMat.Change)#Let's check the magnitude of the substracted weights
#[1] -0.806  0.805[1] -0.806  0.805
any(is.na(KunonIR24.AdjMat.Change))
#Generate the breaks and color functions
breaks <- c(seq(-0.5, 0.5, length=6))
breaks
HMCol<- brewer.pal(5,name = "RdBu")
HMCol<- rev(HMCol)
HMCol
[1] "#2166AC" "#4393C3" "#92C5DE" "#D1E5F0" "#F7F7F7" "#FDDBC7" "#F4A582" "#D6604D" "#B2182B"
[1] "#0571B0" "#92C5DE" "#F7F7F7" "#F4A582" "#CA0020"
#Plot the HeatMap
par(mar=rep(0,4)) # make plot margins 0
heatmap(KunonIR24.AdjMat.Change, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=HMCol, breaks=breaks,frame=T, cexRow = 0.5 , cexCol = 0.5, main = "Ku80 nonIR24 AdjMat (minus nonIR5 edges)", ylab="Total Protein Interactors",xlab="Total Protein Interactors")


#Now we substract the KuIR5 AdjMatrix to the KunonIR24 AdjMatrix
KunonIR24.AdjMat.Change2<-Ku.nonIR24.M-Ku.IR24.M
range(KunonIR24.AdjMat.Change)#Let's check the magnitude of the substracted weights
#[1]-0.813  0.668
any(is.na(KunonIR24.AdjMat.Change))
#Generate the breaks and color functions
breaks <- c(seq(-0.5, 0.5, length=10))
breaks
HMCol<- brewer.pal(9,name = "RdBu")
HMCol<- rev(HMCol)
HMCol
[1] "#2166AC" "#4393C3" "#92C5DE" "#D1E5F0" "#F7F7F7" "#FDDBC7" "#F4A582" "#D6604D" "#B2182B"

#Plot the HeatMap
par(mar=rep(0,4)) # make plot margins 0
heatmap(KunonIR24.AdjMat.Change2, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=HMCol, breaks=breaks,frame=T, cexRow = 0.5 , cexCol = 0.5, main = "Ku80 IR24 AdjMat (minus IR24 edges)", ylab="Total Protein Interactors",xlab="Total Protein Interactors")

#
#check with the FNN
# now actually make the heat map
par(mar=rep(0,4)) # make plot margins 0
heatmap(FNN.M, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=0:1, frame=T, cexRow = 0.5 , cexCol = 0.5)
#
legend("topleft", c("Neighbors Ku80","Ku80"), fill=c("red","blue"),  
       bty="n", inset=0, xpd=T,  border=F)

#

#
dev.off()
plot.new()


##################################
#Ku80ID FNN nodes Network Edge Lists: This is the Ku80ID FNN, this only contains the Ku80ID interactors in contrast with the first FNN
Ku_FNN<-read.delim("./Ku80_Files/Ku_FNN.tsv",sep = "\t",col.names = )
str(Ku_FNN)
head(Ku_FNN)

#Generate the igraph network object
Ku.FNN<-graph_from_data_frame(d=Ku_FNN,vertices = Mat, directed = F )#Generate the igraph object with the edge and node annotation attributes
str(vertex_attr(Ku.FNN))#Check node attributes
str(edge_attr(Ku.FNN))#Check for edge attributes
Ku.FNN<-igraph::simplify(Ku.FNN, remove.multiple=TRUE, remove.loops=TRUE) #Remove repetitive edges and self-connecting edges
Ku.FNN
#Small plot
plot(Ku.FNN)

#Generate an adjacency Matrix that will be used after
Ku.FNN_AdjMat<-as_adjacency_matrix(Ku.FNN,type="both",attr = "weight",names=TRUE,sparse = TRUE)
Ku.FNN_AdjMat<-as.matrix(Ku.FNN_AdjMat)
write.csv(Ku.FNN_AdjMat,"./Ku.FNN_AdjMat.csv")

#Generate a vertex size based on degree and count
vsize<-igraph::degree(Ku.FNN)/vcount(Ku.FNN)
#Identify and remove disconnected nodes
isolated<-which(igraph::degree(Ku.FNN)==0)
isolated
length(isolated)# [1] 599
#remove isolated nodes
Ku.FNN<-igraph::delete_vertices(Ku.FNN,isolated)


#Generate layouts for the graph
#Do not run again, to not to change the layouts that I already have
Lyt.Nice<-layout_nicely(Ku.FNN,dim = 3)
Lyt.GOpt<-layout_with_graphopt(Ku.FNN, niter=5000, mass = 30, charge = E(Ku.FNN)$weight)
Lyt.DrL<-layout_with_drl(Ku.FNN, weights = E(Ku.FNN)$weight)
Lyt.FR<-layout.fruchterman.reingold(Ku.FNN, niter=10000, weights=E(Ku.FNN)$weight)#This one looks like the best option
Lyt.KK<-layout_with_kk(Ku.FNN, maxiter=10000, weights=E(Ku.FNN)$weight) #very large separation of peripheral nodes
Lyt.GEM<-layout_with_gem(Ku.FNN, maxiter=10000) #Bad separation
#Lyt.DH<-layout_with_dh(Ku.FNN, maxiter=150, fineiter=max(10, log2(833)))
Lyt.LGL<-layout.lgl(Ku.FNN, maxit=150,maxdelta=833,area=(length(E(Ku.FNN)))^2,coolexp=1.5,repulserad=((length(E(Ku.FNN)))^2)*(length(E(Ku.FNN))))#Do not work


#Star and Tree layouts

V(Ku.FNN)[name=="SIRT6"]
which(V(Ku.FNN)[name=="SIRT6"])# 610 #Check what is the rowindex for your desired node
which(V(Ku.FNN)$name=="SIRT6")# 610 #Check what is the rowindex for your desired node
which(V(Ku.FNN)$name=="MRE11A")# 327
which(V(Ku.FNN)$name=="XRCC5")# 722
Lyt.Star<-layout_as_star(Ku.FNN,center=V(Ku.FNN)[name=="SIRT6"]) #it has a special order argument where we can assign an order to the vertex, hence maybe the OPLS-DA:VIP value order may be used as ordering vector

Lyt.Tree<-layout_as_tree(Ku.FNN,root=c(190))
Lyt.RT<-layout.reingold.tilford(Ku.FNN,root=c(190))# This one looks good, but the indexing do not allow to highlight the correct nodes
Lyt.sugi<-layout_with_sugiyama(Ku.FNN, maxiter=100, weights=E(Ku.FNN2)$weight)#This one looks good for undirected tree layouts

#Check what layout looks fine

plot.igraph(Ku.FNN,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(Ku.FNN)$An.color, vertex.frame.color=V(Ku.FNN)$An.color, vertex.shape="circle",
            vertex.size=5, vertex.label=V(Ku.FNN)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(Ku.FNN)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)

#With ggraph

ggraph(Ku.FNN,layout=Lyt.FR[-isolated,])+
  geom_edge_fan(color="gray65", width=1, alpha=0.4)+
  geom_node_point(aes(color=V(Ku.FNN)$Annotation, size=V(Ku.FNN)$Ku.FC.nonIR5),alpha=0.75)+
  geom_node_text(aes(label=name),color="gray25",size=4)+
  scale_color_manual(limits=V(Ku.FNN)$Annotation, values=V(Ku.FNN)$An.color)+
  scale_size_continuous(range = c(5, 15)) +
  theme_void()+
  title("Sirt6 nonIR5 Network")+
  guides(color=guide_legend(override.aes = list(size=7)))+
  theme(legend.text=element_text(size=12,face = "bold"),
        legend.title = element_text(size=13,face = "bold"))



#Community Detection
#clp <- igraph::cluster_optimal(Ku.FNN)
class(clp)
V(Ku.FNN)$community <- clp$membership
colrs <- adjustcolor( c("gray50", "tomato", "gold", "yellowgreen"), alpha=.6)


plot.igraph(clp,Ku.FNN,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(Ku.FNN)$An.color, vertex.frame.color=V(Ku.FNN)$An.color, vertex.shape="circle",
            vertex.size==V(Ku.FNN)$Ku.FC.nonIR5, vertex.label=V(Ku.FNN)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(Ku.FNN)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#Distances edges highlight
dist.from.Sirt6 <- distances(Ku.FNN, v=V(Ku.FNN)[name=="SIRT6"],
                             to=V(Ku.FNN), weights=E(Ku.FNN)$weight)
dist.from.Sirt6<-dist.from.Sirt6[,!is.infinite(colSums(dist.from.Sirt6))]

dist.from.Sirt6
range(dist.from.Sirt6)

#Set colors to plot the distances:
oranges <- colorRampPalette(c("dark red", "gold"))
col <- oranges(max(dist.from.Sirt6))
col <- col[dist.from.Sirt6+1]
plot(Ku.FNN, layout=Lyt.FR[-isolated,],vertex.color=col, vertex.label=dist.from.Sirt6,
     vertex.label.color="white",vertex.shape="circle",
     vertex.size=V(Ku.FNN)$Ku.FC.nonIR5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(Ku.FNN)$weight, edge.lty=1, main="Sirt6 nonIR5 Full Nuclear Network", sub="FR Layout",
     rescale=T)


#Shortest path from Sirt6 to Nucleolin
KuNcl.path <- shortest_paths(Ku.FNN,
                             from = V(Ku.FNN)[name=="SIRT6"],
                             to = V(Ku.FNN)[name=="NCL"],
                             output = "both") # both path nodes and edges
# Generate edge color variable to plot the path:
ecol <- rep("gray80", ecount(Ku.FNN))
ecol[unlist(KuNcl.path$epath)] <- "orange"
# Generate edge width variable to plot the path:
ew <- rep(2, ecount(Ku.FNN))
ew[unlist(KuNcl.path$epath)] <- 4
# Generate node color variable to plot the path:
vcol <- rep("gray40", vcount(Ku.FNN))
vcol[unlist(KuNcl.path$vpath)] <- "gold"
plot(Ku.FNN, layout=Lyt.FR[-isolated,], vertex.color=vcol, edge.color=ecol,
     edge.width=ew,vertex.label.color="white",vertex.shape="circle",
     vertex.size=5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(Ku.FNN)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)



#Plot based on Sirt6 as incident
inc.edges <- incident(Ku.FNN, V(Ku.FNN)[name=="SIRT6"], mode="all")
# Set colors to plot the selected edges.
ecol <- rep("gray80", ecount(Ku.FNN))
ecol[inc.edges] <- "orange"
vcol <- rep("gray40", vcount(Ku.FNN))
vcol[V(Ku.FNN)$name=="SIRT6"] <- "gold"

plot(Ku.FNN, layout=Lyt.FR[-isolated,], vertex.color=vcol, edge.color=ecol, 
     edge.width=ew,vertex.label.color="white",vertex.shape="circle",
     vertex.size=V(Ku.FNN)$Ku.FC.nonIR5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(Ku.FNN)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)


#Plot Degree Distribution
Ku.FNN.Deg<-degree_distribution(Ku.FNN, cumulative = T, mode="all")
plot( x=0:max(igraph::degree(Ku.FNN)), y=1-Ku.FNN.Deg, pch=19, cex=1.2, col="orange",
      xlab="Degree", ylab="Cumulative Frequency")

#Network Descriptives
#Number of nodes and edges
vcount(Ku.FNN) #[1] 405
ecount(Ku.FNN)#[1] 4402
#Density (Number of edges in relation to the total number of edges in a full connected graph)
edge_density(Ku.FNN,loops = F) #[1]0.0538076
#Transitivity (Total number of conected triangles)
transitivity(Ku.FNN,type="global") #[1] 0.6395053
#modularity
mod<-V(Ku.FNN)$Mre.FuzzyCMeans #lets check for fuzzy Cmeans modularity
mod<-mod+1#The non clustered interactors that appears as 0, needs to be transformed to 1, so we add 1 to every cluster
modularity(Ku.FNN,membership = mod)
#S6.FuzzyCMeans modularity [1] [1] 0.01375601
#Ku.FuzzyCMeans modularity [1] 0.0208517
#Mre.FuzzyCMeans modularity [1]0.02361523

#Assorativity
assortativity_nominal(Ku.FNN,factor(V(Ku.FNN)$Mre.FuzzyCMeans),directed=F)
#Annotation [1] [1]0.2244201
#S6.FuzzyCMeans [1] 0.02020786
#Ku.FuzzyCMeans [1] 0.03150289
#Mre.FuzzyCMeans [1]0.02929151

#Diameter (longest geodesic distance), shortest path between two most peripheral nodes
diameter(Ku.FNN,directed = F, weights = E(Ku.FNN)$weights)#4.293
diam<-get_diameter(Ku.FNN,directed = F,weights = E(Ku.FNN)$weights)
diam #+ 6/234 vertices, named, from 74b8ef8:[1] DDX46    SNRPD3   HIST1H4F TMPO     RPN1     RPN2  
#Plot Diameter
vcol <- rep("gray50", vcount(Ku.FNN))
vcol[diam] <- "gold"
ecol <- rep("gray90", ecount(Ku.FNN))
ecol[E(Ku.FNN, path=diam)] <- "orange"
# E(net, path=diam) finds edges along a path, here 'diam'
plot(Ku.FNN, vertex.color=vcol, edge.color=ecol, vertex.label.color="black",edge.arrow.mode=0,layout=Lyt.FR[-isolated,],vertex.size=V(Ku.FNN)$Ku.FC.nonIR5)

#Node topological descriptives
degree<-igraph::centr_degree(Ku.FNN,mode = "all",normalized = T)
degree<-degree$res
eigen_cent<-igraph::eigen_centrality(Ku.FNN)
eigen_cent<-eigen_cent$vector
clos_cent <- igraph::closeness(Ku.FNN, normalized=T)
betw_cent <- igraph::betweenness(Ku.FNN, normalized=T)
hs <- hub_score(Ku.FNN, weights=E(Ku.FNN)$weight)$vector
as <- authority_score(Ku.FNN, weights=E(Ku.FNN)$weight)$vector
Trans<-transitivity(Ku.FNN,type="local")
Trans
#write results
Ku.FNN.Centrality<-data.frame(NAME=V(Ku.FNN)$name, UniprotID=V(Ku.FNN)$Uniprot, Degree=degree, Eigen= eigen_cent, Closeness= clos_cent, Betweeness= betw_cent, Hubs=hs, Authority=as,Transitivity=Trans)
Ku.FNN.Centrality
write.csv(Ku.FNN.Centrality,"./Ku.FNN.Centrality.csv")

#Cocitation, Number of shared interactions among nodes
Cocitation<-cocitation(Ku.FNN)
Cocitation#This matrix is similar to an adjacency matrix
write.csv(Cocitation,"./Ku.FNN.cocitation.csv")
#Plot with degree as node size
ggraph(Ku.FNN,layout=Lyt.FR[-isolated,])+
  geom_edge_fan(color="gray85", width=1, alpha=0.4)+
  geom_node_point(aes(color=V(Ku.FNN)$Annotation, size=degree),alpha=0.75)+
  geom_node_text(aes(label=name),color="gray25",size=4)+
  scale_color_manual(limits=V(Ku.FNN)$Annotation, values=V(Ku.FNN)$An.color)+
  scale_size_continuous(range = c(5, 15)) +
  theme_void()+
  title("Sirt6 nonIR5 Network")+
  guides(color=guide_legend(override.aes = list(size=7)))+
  theme(legend.text=element_text(size=12,face = "bold"),
        legend.title = element_text(size=13,face = "bold"))


#Modularity Analysis

#Edge_Betweeness Modularity
Ku.FNN.ebmod=edge.betweenness.community(Ku.FNN,weights=E(Ku.FNN)$weight,directed = F) 
Ku.FNN.ebmod# groups: 15, mod: 0.48
Ku.FNN.ebmod=edge.betweenness.community(Ku.FNN,weights=E(Ku.FNN)$weight,directed = F) 
Ku.FNN.ebmod #37 groups, good enough mod 0.47
#
plot.igraph(Ku.FNN.ebmod,Ku.FNN,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(Ku.FNN)$An.color, vertex.frame.color=V(Ku.FNN)$An.color, vertex.shape="circle",
            vertex.size=5, vertex.label=V(Ku.FNN)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(Ku.FNN)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)

#
plot(Ku.FNN.ebmod, Ku.FNN,layout=Lyt.FR[-isolated,], alpha=0.6, 
     vertex.color=V(Ku.FNN)$An.color, vertex.frame.color=V(Ku.FNN)$An.color, vertex.shape="circle", 
     vertex.size=5, vertex.label=V(Ku.FNN)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black", 
     edge.width=E(Ku.FNN)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout", 
     rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)
#Dendogram of the Edge_Betweeness
dendPlot(Ku.FNN.ebmod, mode="hclust")

#Louvain modularity
Ku.FNN.Louvain.mod<-cluster_louvain(Ku.FNN,weights=E(Ku.FNN)$weight)
Ku.FNN.Louvain.mod #IGRAPH clustering multi level, groups: 11, mod: 0.43

#Plot the Louvain communities
plot(Ku.FNN.Louvain.mod, Ku.FNN,layout=Lyt.FR[-isolated,], alpha=0.6, 
     vertex.color=V(Ku.FNN)$An.color, vertex.frame.color=V(Ku.FNN)$An.color, vertex.shape="circle", 
     vertex.size=5, vertex.label=V(Ku.FNN)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black", 
     edge.width=E(Ku.FNN)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout", 
     rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#Generate HeatMaps of the adjacency matrixes
#Let's use the neighbor list for generating the Adjaency Matrix
SensorID_FNN #Let's use the edge list that contain the total amount of possible interactors for the three sensors
v<-levels(as.factor(unlist(SensorID_FNN[,1:2]))) #get vertex names from the network list
n<-length(v)#get the number of vertex in the network
e<-matrix(match(as.character(unlist(SensorID_FNN[,1:2])),v),ncol = 2)#Edgelist with row index instead of Names
w<-Ku_FNN$weight #Edge weights, only from the Ku_FNN network to quit those edges that do not exist at Ku_FNN

# M is our co-expression network adjacency matrix.

# distringuish between i,j and j,i we want to make M symmetric (undirected) by
# considering an undirected network.
Ku.FNN.M <- matrix(0, n, n)                  # set up a co-expression matrix
Ku.FNN.M[e] <- w                             # fill it in with edge weights
Ku.FNN.M <- Ku.FNN.M + t(Ku.FNN.M)                         # make this symmetric
any(is.na(Ku.FNN.M))#check no NAs in the matrix
dimnames(Ku.FNN.M) <- list(v, v)             # label the vertices
# We let A_ij = 1 if there is an interaction. 
#A <- 1*(M > 0) #right now we keep the weights as they are
Ku.FNN.M #The adjacency Matrix not as an igraph object, easier to manipulate, and to subset for smaller Sensor specific networks
quantile(Ku.FNN.M,c(0,0.25,0.5,0.75,1))
#
#  0%   25%   50%   75%  100% 
#0.000 0.000 0.000 0.000 0.991 

#Creating the HeatMap
# make a heatmap image with colors: white below 0,
# grey from 0.4 to 1 with darker closer to 1
breaks <- c(0, seq(0.4, 0.991, length=11)) # breaks in the color bins
cols <- grey(1-c(0,seq(0.5,1,length=10))) # colors

#Red Blue Colors
library(RColorBrewer)
breaks <- c(seq(-1, 1, length=10))
breaks
cols <- gray(1-c(0,seq(0,1,length=10))) # colors
cols2 <- hcl.colors(10 , palette = "Blues3" , alpha = 1 , rev = TRUE)
cols3 <- colorRamp("red" , space = "rgb" , interpolate = "linear")
HMCol<- brewer.pal(10,name = "RdBu")
HMCol<- rev(HMCol)
HMCol


# color the vertex for CCNB1 blue, its neighbors red, and the others white
SensorIDBind <- which(v == "SIRT6") # index for ccnb1
#SensorIDBind<-as.factor(SensorIDBind)
vcols <- rep("white",n)
vcols[SensorIDBind] <- "blue"
#Subsetting Ku.FNN.M by doing a loop asking for look into every row and column the index in SensorIDBind worked
for (i in SensorIDBind) {
  vcols[which(Ku.FNN.M[,i] > 0 | Ku.FNN.M[i,])] <- "red"
}
##Color annotation for functional annotation
MatSet<-subset(Mat,Mat$NAME%in%v)# select only those proteins that are cointained in the node list of our matrixes
match(v,MatSet$NAME) # Check that the indexes of the Annotation Dataframe and the node list are in the same position 
match(rownames(Ku.FNN.M),MatSet$NAME) #Check that the row names of the Adj Matrix and the Annotation dataframe are in the same order
match(colnames(Ku.FNN.M),MatSet$NAME) #Check that the column names of the Adj Matrix and the Annotation dataframe are in the same order
AnCol<-MatSet$An.color#Vector of colors
AnCol<-as.character(AnCol)#transform it into a character vector, so the colors can be read

#Subsetting Ku.FNN.M using the index vector in SensorIDBind do not work
#vcols[which(Ku.FNN.M[,SensorIDBind] > 0 | Ku.FNN.M[SensorIDBind,])] <- "red"

#Generate a dendogram for keeping the proteins in the sample place
#Ku.FNN.dist<-dist(Ku.FNN.M,method = "euclidian")
#Ku.FNN.Dend<-hclust(Ku.FNN.dist, method = "ward.D" )
#plot(Ku.FNN.Dend)

#Cut the dendogram
#cut.height <- 10 # try numbers from 4 to 7
#FNN.hclusts<- cutree(FNN.Dend, cut.height)
#plot(FNN.hclusts)

# now actually make the heat map
par(mar=rep(0,4)) # make plot margins 0
heatmap(Ku.FNN.M, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=cols, breaks=breaks,frame=T, cexRow = 0.5 , cexCol = 0.5, main = "Sirt6 nonIR5 AdjMat", ylab="Total Protein Interactors",xlab="Total Protein Interactors")

#
plot.new()
legend("topleft", c("Cell Cycle","Chromatin Org","DNA Repair","Nuclear Pore","p53","Proteasome","Ribosome Biogenesis","RNA Metabolism","Transcription Reg","UPR"), 
       fill=c("coral1","darkviolet","firebrick1","deepskyblue3","tomato","grey","lightgreen","gold","slateblue1","maroon"),  
       bty="n", inset=0, xpd=T,  border=F, cex =2,y.intersp = 0.7,x.intersp = 0.2,pt.cex = 1)

#check with the FNN
# now actually make the heat map
par(mar=rep(0,4)) # make plot margins 0
heatmap(FNN.M, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=0:1, frame=T, cexRow = 0.5 , cexCol = 0.5)
#
legend("topleft", c("Neighbors SIRT6","SIRT6"), fill=c("red","blue"),  
       bty="n", inset=0, xpd=T,  border=F)





################################################Mre11 Networks
###########################################################
#Mre11 Networks


#Let's keep the Mat as the Node attribute table and load the Neighbor list as the matrixes
Mat# Node attribute matrix for all the networks

#Mre11ID nonIR5 nodes Network Edge Lists
Mre_nonIR5<-read.delim("./Mre11_Files/Mre_nonIR5.tsv",sep = "\t",col.names = )
str(Mre_nonIR5)
head(Mre_nonIR5)

#Generate the igraph network object
Mre.nonIR5<-graph_from_data_frame(d=Mre_nonIR5,vertices = Mat, directed = F )#Generate the igraph object with the edge and node annotation attributes
str(vertex_attr(Mre.nonIR5))#Check node attributes
str(edge_attr(Mre.nonIR5))#Check for edge attributes
Mre.nonIR5<-igraph::simplify(Mre.nonIR5, remove.multiple=TRUE, remove.loops=TRUE) #Remove repetitive edges and self-connecting edges
Mre.nonIR5
#Small plot
plot(Mre.nonIR5)

#Generate an adjacency Matrix that will be used after
Mre.nonIR5_AdjMat<-as_adjacency_matrix(Mre.nonIR5,type="both",attr = "weight",names=TRUE,sparse = TRUE)
Mre.nonIR5_AdjMat<-as.matrix(Mre.nonIR5_AdjMat)
write.csv(Mre.nonIR5_AdjMat,"./Mre.nonIR5_AdjMat.csv")

#Generate a vertex size based on degree and count
vsize<-igraph::degree(Mre.nonIR5)/vcount(Mre.nonIR5)
#Identify and remove disconnected nodes
isolated<-which(igraph::degree(Mre.nonIR5)==0)
isolated
length(isolated)# [1] 748
#remove isolated nodes
Mre.nonIR5<-igraph::delete_vertices(Mre.nonIR5,isolated)


#Generate layouts for the graph
#Do not run again, to not to change the layouts that I already have
Lyt.Nice<-layout_nicely(Mre.nonIR5,dim = 3)
Lyt.GOpt<-layout_with_graphopt(Mre.nonIR5, niter=5000, mass = 30, charge = E(Mre.nonIR5)$weight)
Lyt.DrL<-layout_with_drl(Mre.nonIR5, weights = E(Mre.nonIR5)$weight)
Lyt.FR<-layout.fruchterman.reingold(Mre.nonIR5, niter=10000, weights=E(Mre.nonIR5)$weight)#This one looks like the best option
Lyt.KK<-layout_with_kk(Mre.nonIR5, maxiter=10000, weights=E(Mre.nonIR5)$weight) #very large separation of peripheral nodes
Lyt.GEM<-layout_with_gem(Mre.nonIR5, maxiter=10000) #Bad separation
#Lyt.DH<-layout_with_dh(Mre.nonIR5, maxiter=150, fineiter=max(10, log2(833)))
Lyt.LGL<-layout.lgl(Mre.nonIR5, maxit=150,maxdelta=833,area=(length(E(Mre.nonIR5)))^2,coolexp=1.5,repulserad=((length(E(Mre.nonIR5)))^2)*(length(E(Mre.nonIR5))))#Do not work


#Star and Tree layouts

V(Mre.nonIR5)[name=="MRE11A"]
which(V(Mre.nonIR5)[name=="Mre11"])# 610 #Check what is the rowindex for your desired node
which(V(Mre.nonIR5)$name=="MRE11A")# 610 #Check what is the rowindex for your desired node
which(V(Mre.nonIR5)$name=="MRE11A")# 35
which(V(Mre.nonIR5)$name=="MRE11A")# 35
Lyt.Star<-layout_as_star(Mre.nonIR5,center=V(Mre.nonIR5)[name=="MRE11A"]) #it has a special order argument where we can assign an order to the vertex, hence maybe the OPLS-DA:VIP value order may be used as ordering vector

Lyt.Tree<-layout_as_tree(Mre.nonIR5,root=c(35))
Lyt.RT<-layout.reingold.tilford(Mre.nonIR5,root=c(35))# This one looks good, but the indexing do not allow to highlight the correct nodes
Lyt.sugi<-layout_with_sugiyama(Mre.nonIR5, maxiter=100, weights=E(Mre.nonIR52)$weight)#This one looks good for undirected tree layouts

#Check what layout looks fine

plot.igraph(Mre.nonIR5,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(Mre.nonIR5)$An.color, vertex.frame.color=V(Mre.nonIR5)$An.color, vertex.shape="circle",
            vertex.size=5, vertex.label=V(Mre.nonIR5)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(Mre.nonIR5)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)

#With ggraph

ggraph(Mre.nonIR5,layout=Lyt.FR[-isolated,])+
  geom_edge_fan(color="gray65", width=1, alpha=0.4)+
  geom_node_point(aes(color=V(Mre.nonIR5)$Annotation, size=V(Mre.nonIR5)$Mre.FC.nonIR5),alpha=0.75)+
  geom_node_text(aes(label=name),color="gray25",size=4)+
  scale_color_manual(limits=V(Mre.nonIR5)$Annotation, values=V(Mre.nonIR5)$An.color)+
  scale_size_continuous(range = c(5, 15)) +
  theme_void()+
  title("Mre11 nonIR5 Network")+
  guides(color=guide_legend(override.aes = list(size=7)))+
  theme(legend.text=element_text(size=12,face = "bold"),
        legend.title = element_text(size=13,face = "bold"))



#Community Detection
#clp <- igraph::cluster_optimal(Mre.nonIR5)
class(clp)
V(Mre.nonIR5)$community <- clp$membership
colrs <- adjustcolor( c("gray50", "tomato", "gold", "yellowgreen"), alpha=.6)


plot.igraph(clp,Mre.nonIR5,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(Mre.nonIR5)$An.color, vertex.frame.color=V(Mre.nonIR5)$An.color, vertex.shape="circle",
            vertex.size==V(Mre.nonIR5)$Mre.FC.nonIR5, vertex.label=V(Mre.nonIR5)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(Mre.nonIR5)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#Distances edges highlight
dist.from.Mre11 <- distances(Mre.nonIR5, v=V(Mre.nonIR5)[name=="MRE11A"],
                             to=V(Mre.nonIR5), weights=E(Mre.nonIR5)$weight)
dist.from.Mre11<-dist.from.Mre11[,!is.infinite(colSums(dist.from.Mre11))]

dist.from.Mre11
range(dist.from.Mre11)

#Set colors to plot the distances:
oranges <- colorRampPalette(c("dark red", "gold"))
col <- oranges(max(dist.from.Mre11))
col <- col[dist.from.Mre11+1]
plot(Mre.nonIR5, layout=Lyt.FR[-isolated,],vertex.color=col, vertex.label=dist.from.Mre11,
     vertex.label.color="white",vertex.shape="circle",
     vertex.size=5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(Mre.nonIR5)$weight, edge.lty=1, main="Mre11 nonIR5 Full Nuclear Network", sub="FR Layout",
     rescale=T)


#Shortest path from Mre11 to Nucleolin
MreNcl.path <- shortest_paths(Mre.nonIR5,
                              from = V(Mre.nonIR5)[name=="MRE11A"],
                              to = V(Mre.nonIR5)[name=="NCL"],
                              output = "both") # both path nodes and edges
# Generate edge color variable to plot the path:
ecol <- rep("gray80", ecount(Mre.nonIR5))
ecol[unlist(MreNcl.path$epath)] <- "orange"
# Generate edge width variable to plot the path:
ew <- rep(2, ecount(Mre.nonIR5))
ew[unlist(MreNcl.path$epath)] <- 4
# Generate node color variable to plot the path:
vcol <- rep("gray40", vcount(Mre.nonIR5))
vcol[unlist(MreNcl.path$vpath)] <- "gold"
plot(Mre.nonIR5, layout=Lyt.FR[-isolated,], vertex.color=vcol, edge.color=ecol,
     edge.width=ew,vertex.label.color="white",vertex.shape="circle",
     vertex.size=5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(Mre.nonIR5)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)



#Plot based on Mre11 as incident
inc.edges <- incident(Mre.nonIR5, V(Mre.nonIR5)[name=="MRE11A"], mode="all")
# Set colors to plot the selected edges.
ecol <- rep("gray80", ecount(Mre.nonIR5))
ecol[inc.edges] <- "orange"
vcol <- rep("gray40", vcount(Mre.nonIR5))
vcol[V(Mre.nonIR5)$name=="Mre11"] <- "gold"

plot(Mre.nonIR5, layout=Lyt.FR[-isolated,], vertex.color=vcol, edge.color=ecol, 
     edge.width=ew,vertex.label.color="white",vertex.shape="circle",
     vertex.size=5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(Mre.nonIR5)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)


#Plot Degree Distribution
Mre.nonIR5.Deg<-degree_distribution(Mre.nonIR5, cumulative = T, mode="all")
plot( x=0:max(igraph::degree(Mre.nonIR5)), y=1-Mre.nonIR5.Deg, pch=19, cex=1.2, col="orange",
      xlab="Degree", ylab="Cumulative Frequency")

#Network Descriptives
#Number of nodes and edges
vcount(Mre.nonIR5) #[1] 85
ecount(Mre.nonIR5)#[1] 242
#Density (Number of edges in relation to the total number of edges in a full connected graph)
edge_density(Mre.nonIR5,loops = F) #[1] 0.06778711
#Transitivity (Total number of conected triangles)
transitivity(Mre.nonIR5,type="global") #[1] 0.6456522
#modularity
mod<-V(Mre.nonIR5)$Mre.FuzzyCMeans #lets check for fuzzy Cmeans modularity
mod<-mod+1#The non clustered interactors that appears as 0, needs to be transformed to 1, so we add 1 to every cluster
modularity(Mre.nonIR5,membership = mod)
#S6.FuzzyCMeans modularity [1] -0.001348952
#Ku.FuzzyCMeans modularity [1] -0.002680828
#Mre.FuzzyCMeans modularity [1]-0.01143194

#Assorativity
assortativity_nominal(Mre.nonIR5,factor(V(Mre.nonIR5)$Mre.FuzzyCMeans ),directed=F)
#Annotation [1] [1]  0.2481616
#S6.FuzzyCMeans [1] -0.002320595
#Ku.FuzzyCMeans [1]  -0.00509426
#Mre.FuzzyCMeans [1] -0.01782173

#Diameter (longest geodesic distance), shortest path between two most peripheral nodes
diameter(Mre.nonIR5,directed = F, weights = E(Mre.nonIR5)$weights)#4.29
diam<-get_diameter(Mre.nonIR5,directed = F,weights = E(Mre.nonIR5)$weights)
diam #+ 9/85 vertices, named, from 30db4c6: [1] DNAJA2   HSPA1B   HSP90AA1 VCP      ATAD3A   PHB      RPL7     EIF4A3   DDX39B    
#Plot Diameter
vcol <- rep("gray50", vcount(Mre.nonIR5))
vcol[diam] <- "gold"
ecol <- rep("gray90", ecount(Mre.nonIR5))
ecol[E(Mre.nonIR5, path=diam)] <- "orange"
# E(net, path=diam) finds edges along a path, here 'diam'
plot(Mre.nonIR5, vertex.color=vcol, edge.color=ecol, vertex.label.color="black",edge.arrow.mode=0,layout=Lyt.FR[-isolated,],vertex.size=5)

#Node topological descriptives
degree<-igraph::centr_degree(Mre.nonIR5,mode = "all",normalized = T)
degree<-degree$res
eigen_cent<-igraph::eigen_centrality(Mre.nonIR5)
eigen_cent<-eigen_cent$vector
clos_cent <- igraph::closeness(Mre.nonIR5, normalized=T)
betw_cent <- igraph::betweenness(Mre.nonIR5, normalized=T)
hs <- hub_score(Mre.nonIR5, weights=E(Mre.nonIR5)$weight)$vector
as <- authority_score(Mre.nonIR5, weights=E(Mre.nonIR5)$weight)$vector
Trans<-transitivity(Mre.nonIR5,type="local")
Trans
#write results
Mre.nonIR5.Centrality<-data.frame(NAME=V(Mre.nonIR5)$name, UniprotID=V(Mre.nonIR5)$Uniprot, Degree=degree, Eigen= eigen_cent, Closeness= clos_cent, Betweeness= betw_cent, Hubs=hs, Authority=as,Transitivity=Trans)
Mre.nonIR5.Centrality
write.csv(Mre.nonIR5.Centrality,"./Mre.nonIR5.Centrality.csv")

#Cocitation, Number of shared interactions among nodes
Cocitation<-cocitation(Mre.nonIR5)
Cocitation#This matrix is similar to an adjacency matrix
write.csv(Cocitation,"./Mre.nonIR5.cocitation.csv")
#Plot with degree as node size
ggraph(Mre.nonIR5,layout=Lyt.FR[-isolated,])+
  geom_edge_fan(color="gray85", width=1, alpha=0.4)+
  geom_node_point(aes(color=V(Mre.nonIR5)$Annotation, size=degree),alpha=0.75)+
  geom_node_text(aes(label=name),color="gray25",size=4)+
  scale_color_manual(limits=V(Mre.nonIR5)$Annotation, values=V(Mre.nonIR5)$An.color)+
  scale_size_continuous(range = c(5, 15)) +
  theme_void()+
  title("Mre11 nonIR5 Network")+
  guides(color=guide_legend(override.aes = list(size=7)))+
  theme(legend.text=element_text(size=12,face = "bold"),
        legend.title = element_text(size=13,face = "bold"))


#Modularity Analysis

#Edge_Betweeness Modularity
Mre.nonIR5.ebmod=edge.betweenness.community(Mre.nonIR5,weights=E(Mre.nonIR5)$weight,directed = F) 
Mre.nonIR5.ebmod# IGRAPH clustering edge betweenness, groups: 11, mod: 0.46
#
plot.igraph(Mre.nonIR5.ebmod,Mre.nonIR5,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(Mre.nonIR5)$An.color, vertex.frame.color=V(Mre.nonIR5)$An.color, vertex.shape="circle",
            vertex.size=5, vertex.label=V(Mre.nonIR5)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(Mre.nonIR5)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)

#
plot(Mre.nonIR5.ebmod, Mre.nonIR5,layout=Lyt.FR[-isolated,], alpha=0.6, 
     vertex.color=V(Mre.nonIR5)$An.color, vertex.frame.color=V(Mre.nonIR5)$An.color, vertex.shape="circle", 
     vertex.size=5, vertex.label=V(Mre.nonIR5)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black", 
     edge.width=E(Mre.nonIR5)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout", 
     rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)
#Dendogram of the Edge_Betweeness
dendPlot(Mre.nonIR5.ebmod, mode="hclust")

#Louvain modularity
Mre.nonIR5.Louvain.mod<-cluster_louvain(Mre.nonIR5,weights=E(Mre.nonIR5)$weight)
Mre.nonIR5.Louvain.mod #IGRAPH clustering multi level, groups: 8, mod: 0.48

#Plot the Louvain communities
plot(Mre.nonIR5.Louvain.mod, Mre.nonIR5,layout=Lyt.FR[-isolated,], alpha=0.6, 
     vertex.color=V(Mre.nonIR5)$An.color, vertex.frame.color=V(Mre.nonIR5)$An.color, vertex.shape="circle", 
     vertex.size=5, vertex.label=V(Mre.nonIR5)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black", 
     edge.width=E(Mre.nonIR5)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout", 
     rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#Generate HeatMaps of the adjacency matrixes
#Let's use the neighbor list for generating the Adjaency Matrix
SensorID_FNN #Let's use the edge list that contain the total amount of possible interactors for the three sensors
v<-levels(as.factor(unlist(SensorID_FNN[,1:2]))) #get vertex names from the network list
n<-length(v)#get the number of vertex in the network
e<-matrix(match(as.character(unlist(SensorID_FNN[,1:2])),v),ncol = 2)#Edgelist with row index instead of Names
w<-Mre_nonIR5$weight #Edge weights, only from the Mre_nonIR5 network to quit those edges that do not exist at Mre_nonIR5

# M is our co-expression network adjacency matrix.

# distringuish between i,j and j,i we want to make M symmetric (undirected) by
# considering an undirected network.
Mre.nonIR5.M <- matrix(0, n, n)                  # set up a co-expression matrix
Mre.nonIR5.M[e] <- w                             # fill it in with edge weights
Mre.nonIR5.M <- Mre.nonIR5.M + t(Mre.nonIR5.M)                         # make this symmetric
any(is.na(Mre.nonIR5.M))#check no NAs in the matrix
dimnames(Mre.nonIR5.M) <- list(v, v)             # label the vertices
# We let A_ij = 1 if there is an interaction. 
#A <- 1*(M > 0) #right now we keep the weights as they are
Mre.nonIR5.M #The adjacency Matrix not as an igraph object, easier to manipulate, and to subset for smaller Sensor specific networks
quantile(Mre.nonIR5.M,c(0,0.25,0.5,0.75,1))
#
#  0%   25%   50%   75%  100% 
#0.000 0.000 0.000 0.000 0.991 

#Creating the HeatMap
# make a heatmap image with colors: white below 0,
# grey from 0.4 to 1 with darker closer to 1
breaks <- c(0, seq(0.4, 0.991, length=11)) # breaks in the color bins
cols <- grey(1-c(0,seq(0.5,1,length=10))) # colors

#Red Blue Colors
library(RColorBrewer)
breaks <- c(seq(-1, 1, length=10))
breaks
cols <- gray(1-c(0,seq(0,1,length=10))) # colors
cols2 <- hcl.colors(10 , palette = "Blues3" , alpha = 1 , rev = TRUE)
cols3 <- colorRamp("red" , space = "rgb" , interpolate = "linear")
HMCol<- brewer.pal(10,name = "RdBu")
HMCol<- rev(HMCol)
HMCol


# color the vertex for CCNB1 blue, its neighbors red, and the others white
SensorIDBind <- which(v == "Mre11") # index for ccnb1
#SensorIDBind<-as.factor(SensorIDBind)
vcols <- rep("white",n)
vcols[SensorIDBind] <- "blue"
#Subsetting Mre.nonIR5.M by doing a loop asking for look into every row and column the index in SensorIDBind worked
for (i in SensorIDBind) {
  vcols[which(Mre.nonIR5.M[,i] > 0 | Mre.nonIR5.M[i,])] <- "red"
}
##Color annotation for functional annotation
MatSet<-subset(Mat,Mat$NAME%in%v)# select only those proteins that are cointained in the node list of our matrixes
match(v,MatSet$NAME) # Check that the indexes of the Annotation Dataframe and the node list are in the same position 
match(rownames(Mre.nonIR5.M),MatSet$NAME) #Check that the row names of the Adj Matrix and the Annotation dataframe are in the same order
match(colnames(Mre.nonIR5.M),MatSet$NAME) #Check that the column names of the Adj Matrix and the Annotation dataframe are in the same order
AnCol<-MatSet$An.color#Vector of colors
AnCol<-as.character(AnCol)#transform it into a character vector, so the colors can be read

#Subsetting Mre.nonIR5.M using the index vector in SensorIDBind do not work
#vcols[which(Mre.nonIR5.M[,SensorIDBind] > 0 | Mre.nonIR5.M[SensorIDBind,])] <- "red"

#Generate a dendogram for keeping the proteins in the sample place
#Mre.nonIR5.dist<-dist(Mre.nonIR5.M,method = "euclidian")
#Mre.nonIR5.Dend<-hclust(Mre.nonIR5.dist, method = "ward.D" )
#plot(Mre.nonIR5.Dend)

#Cut the dendogram
#cut.height <- 10 # try numbers from 4 to 7
#FNN.hclusts<- cutree(FNN.Dend, cut.height)
#plot(FNN.hclusts)

# now actually make the heat map
par(mar=rep(0,4)) # make plot margins 0
heatmap(Mre.nonIR5.M, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=cols, breaks=breaks,frame=T, cexRow = 0.5 , cexCol = 0.5, main = "Mre11 nonIR5 AdjMat", ylab="Total Protein Interactors",xlab="Total Protein Interactors")

#
plot.new()
legend("topleft", c("Cell Cycle","Chromatin Org","DNA Repair","Nuclear Pore","p53","Proteasome","Ribosome Biogenesis","RNA Metabolism","Transcription Reg","UPR"), 
       fill=c("coral1","darkviolet","firebrick1","deepskyblue3","tomato","grey","lightgreen","gold","slateblue1","maroon"),  
       bty="n", inset=0, xpd=T,  border=F, cex =2,y.intersp = 0.7,x.intersp = 0.2,pt.cex = 1)

#check with the FNN
# now actually make the heat map
par(mar=rep(0,4)) # make plot margins 0
heatmap(FNN.M, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=0:1, frame=T, cexRow = 0.5 , cexCol = 0.5)
#
legend("topleft", c("Neighbors Mre11","Mre11"), fill=c("red","blue"),  
       bty="n", inset=0, xpd=T,  border=F)


############################################Mre11 IR5 Networks
##############################################
#Mre11 IR5 Networks


#Let's keep the Mat as the Node attribute table and load the Neighbor list as the matrixes
Mat# Node attribute matrix for all the networks

#Mre11ID IR5 nodes Network Edge Lists
Mre_IR5<-read.delim("./Mre11_Files/Mre_IR5.tsv",sep = "\t",col.names = )
str(Mre_IR5)
head(Mre_IR5)

#Generate the igraph network object
Mre.IR5<-graph_from_data_frame(d=Mre_IR5,vertices = Mat, directed = F )#Generate the igraph object with the edge and node annotation attributes
str(vertex_attr(Mre.IR5))#Check node attributes
str(edge_attr(Mre.IR5))#Check for edge attributes
Mre.IR5<-igraph::simplify(Mre.IR5, remove.multiple=TRUE, remove.loops=TRUE) #Remove repetitive edges and self-connecting edges
Mre.IR5
#Small plot
plot(Mre.IR5)

#Generate an adjacency Matrix that will be used after
Mre.IR5_AdjMat<-as_adjacency_matrix(Mre.IR5,type="both",attr = "weight",names=TRUE,sparse = TRUE)
Mre.IR5_AdjMat<-as.matrix(Mre.IR5_AdjMat)
write.csv(Mre.IR5_AdjMat,"./Mre.IR5_AdjMat.csv")

#Generate a vertex size based on degree and count
vsize<-igraph::degree(Mre.IR5)/vcount(Mre.IR5)
#Identify and remove disconnected nodes
isolated<-which(igraph::degree(Mre.IR5)==0)
isolated
length(isolated)# [1] 771
#remove isolated nodes
Mre.IR5<-igraph::delete_vertices(Mre.IR5,isolated)


#Generate layouts for the graph
#Do not run again, to not to change the layouts that I already have
Lyt.Nice<-layout_nicely(Mre.IR5,dim = 3)
Lyt.GOpt<-layout_with_graphopt(Mre.IR5, niter=5000, mass = 30, charge = E(Mre.IR5)$weight)
Lyt.DrL<-layout_with_drl(Mre.IR5, weights = E(Mre.IR5)$weight)
Lyt.FR<-layout.fruchterman.reingold(Mre.IR5, niter=10000, weights=E(Mre.IR5)$weight)#This one looks like the best option
Lyt.KK<-layout_with_kk(Mre.IR5, maxiter=10000, weights=E(Mre.IR5)$weight) #very large separation of peripheral nodes
Lyt.GEM<-layout_with_gem(Mre.IR5, maxiter=10000) #Bad separation
#Lyt.DH<-layout_with_dh(Mre.IR5, maxiter=150, fineiter=max(10, log2(833)))
Lyt.LGL<-layout.lgl(Mre.IR5, maxit=150,maxdelta=833,area=(length(E(Mre.IR5)))^2,coolexp=1.5,repulserad=((length(E(Mre.IR5)))^2)*(length(E(Mre.IR5))))#Do not work


#Star and Tree layouts

V(Mre.IR5)[name=="MRE11A"]
which(V(Mre.IR5)[name=="Mre11"])# 228 #Check what is the rowindex for your desired node
which(V(Mre.IR5)$name=="Mre11")# 191 #Check what is the rowindex for your desired node
which(V(Mre.IR5)$name=="MRE11A")# 327
which(V(Mre.IR5)$name=="MRE11A")# 143
Lyt.Star<-layout_as_star(Mre.IR5,center=V(Mre.IR5)[name=="MRE11A"]) #it has a special order argument where we can assign an order to the vertex, hence maybe the OPLS-DA:VIP value order may be used as ordering vector

Lyt.Tree<-layout_as_tree(Mre.IR5,root=c(30))
Lyt.RT<-layout.reingold.tilford(Mre.IR5,root=c(30))# This one looks good, but the indexing do not allow to highlight the correct nodes
Lyt.sugi<-layout_with_sugiyama(Mre.IR5, maxiter=100, weights=E(Mre.IR52)$weight)#This one looks good for undirected tree layouts

#Check what layout looks fine

plot.igraph(Mre.IR5,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(Mre.IR5)$An.color, vertex.frame.color=V(Mre.IR5)$An.color, vertex.shape="circle",
            vertex.size=5, vertex.label=V(Mre.IR5)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(Mre.IR5)$weight, edge.lty=1, main="Mre11ID IR5 Network", sub="FR Tree Layout",
            rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#With ggraph

ggraph(Mre.IR5,layout=LytFR[-isolated,])+
  geom_edge_fan(color="gray65", width=1, alpha=0.4)+
  geom_node_point(aes(color=V(Mre.IR5)$Annotation, size=V(Mre.IR5)$Mre.FC.IR5),alpha=0.75)+
  geom_node_text(aes(label=name),color="gray25",size=4)+
  scale_color_manual(limits=V(Mre.IR5)$Annotation, values=V(Mre.IR5)$An.color)+
  scale_size_continuous(range = c(5, 15)) +
  theme_void()+
  title("Mre11 IR5 Network")+
  guides(color=guide_legend(override.aes = list(size=7)))+
  theme(legend.text=element_text(size=12,face = "bold"),
        legend.title = element_text(size=13,face = "bold"))



#Community Detection
#clp <- igraph::cluster_optimal(Mre.IR5)
class(clp)
V(Mre.IR5)$community <- clp$membership
colrs <- adjustcolor( c("gray50", "tomato", "gold", "yellowgreen"), alpha=.6)


plot.igraph(clp,Mre.IR5,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(Mre.IR5)$An.color, vertex.frame.color=V(Mre.IR5)$An.color, vertex.shape="circle",
            vertex.size==V(Mre.IR5)$Mre.FC.IR5, vertex.label=V(Mre.IR5)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(Mre.IR5)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#Distances edges highlight
dist.from.Mre11 <- distances(Mre.IR5, v=V(Mre.IR5)[name=="MRE11A"],
                             to=V(Mre.IR5), weights=E(Mre.IR5)$weight)
dist.from.Mre11<-dist.from.Mre11[,!is.infinite(colSums(dist.from.Mre11))]

dist.from.Mre11
range(dist.from.Mre11)

#Set colors to plot the distances:
oranges <- colorRampPalette(c("dark red", "gold"))
col <- oranges(max(dist.from.Mre11))
col <- col[dist.from.Mre11+1]
plot(Mre.IR5, layout=Lyt.FR[-isolated,],vertex.color=col, vertex.label=dist.from.Mre11,
     vertex.label.color="white",vertex.shape="circle",
     vertex.size=5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(Mre.IR5)$weight, edge.lty=1, main="Mre11 IR5 Full Nuclear Network", sub="FR Layout",
     rescale=T)


#Shortest path from Mre11 to Nucleolin
MreNcl.path <- shortest_paths(Mre.IR5,
                              from = V(Mre.IR5)[name=="MRE11A"],
                              to = V(Mre.IR5)[name=="NCL"],
                              output = "both") # both path nodes and edges
MreNcl.path
# Generate edge color variable to plot the path:
ecol <- rep("gray80", ecount(Mre.IR5))
ecol[unlist(MreNcl.path$epath)] <- "orange"
# Generate edge width variable to plot the path:
ew <- rep(2, ecount(Mre.IR5))
ew[unlist(MreNcl.path$epath)] <- 4
# Generate node color variable to plot the path:
vcol <- rep("gray40", vcount(Mre.IR5))
vcol[unlist(MreNcl.path$vpath)] <- "gold"
plot(Mre.IR5, layout=Lyt.FR[-isolated,], vertex.color=vcol, edge.color=ecol,
     edge.width=ew,vertex.label.color="white",vertex.shape="circle",
     vertex.size=5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(Mre.IR5)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)



#Plot based on Mre11 as incident
inc.edges <- incident(Mre.IR5, V(Mre.IR5)[name=="MRE11A"], mode="all")
# Set colors to plot the selected edges.
ecol <- rep("gray80", ecount(Mre.IR5))
ecol[inc.edges] <- "orange"
vcol <- rep("gray40", vcount(Mre.IR5))
vcol[V(Mre.IR5)$name=="Mre11"] <- "gold"

plot(Mre.IR5, layout=Lyt.FR[-isolated,], vertex.color=vcol, edge.color=ecol, 
     edge.width=ew,vertex.label.color="white",vertex.shape="circle",
     vertex.size=10, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(Mre.IR5)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)


#Plot Degree Distribution
Mre.IR5.Deg<-degree_distribution(Mre.IR5, cumulative = T, mode="all")
plot( x=0:max(igraph::degree(Mre.IR5)), y=1-Mre.IR5.Deg, pch=19, cex=1.2, col="orange",
      xlab="Degree", ylab="Cumulative Frequency")

#Network Descriptives
#Number of nodes and edges
vcount(Mre.IR5) #[1] 62
ecount(Mre.IR5)#[1]  160
#Density (Number of edges in relation to the total number of edges in a full connected graph)
edge_density(Mre.IR5,loops = F) #[1] 0.08461132
#Transitivity (Total number of conected triangles)
transitivity(Mre.IR5,type="global") #[1] 0.5107973
#modularity
mod<-V(Mre.IR5)$S6.FuzzyCMeans #lets check for fuzzy Cmeans modularity
mod<-mod+1#The non clustered interactors that appears as 0, needs to be transformed to 1, so we add 1 to every cluster
modularity(Mre.IR5,membership = mod)
#S6.FuzzyCMeans modularity [1] -0.0009375
#Ku.FuzzyCMeans modularity [1] 0.03058594
#Mre.FuzzyCMeans modularity [1] -0.02068359

#Assorativity
assortativity_nominal(Mre.IR5,factor(V(Mre.IR5)$Annotation),directed=F)
#Annotation [1] [1] 0.1529952
#S6.FuzzyCMeans [1] -0.001902346
#Ku.FuzzyCMeans [1] 0.05631878
#Mre.FuzzyCMeans [1] -0.031611

#Diameter (longest geodesic distance), shortest path between two most peripheral nodes
diameter(Mre.IR5,directed = F, weights = E(Mre.IR5)$weights)#[1] 3.54
diam<-get_diameter(Mre.IR5,directed = F,weights = E(Mre.IR5)$weights)
diam #+ + 7/62 vertices, named, from 5e2fa6e: [1] NOLC1  FBL    RPL4   EIF4A3 DDX17  SFPQ   NONO   
#Plot Diameter
vcol <- rep("gray50", vcount(Mre.IR5))
vcol[diam] <- "gold"
ecol <- rep("gray90", ecount(Mre.IR5))
ecol[E(Mre.IR5, path=diam)] <- "orange"
# E(net, path=diam) finds edges along a path, here 'diam'
plot(Mre.IR5, vertex.color=vcol, edge.color=ecol, vertex.label.color="black",edge.arrow.mode=0,layout=Lyt.FR[-isolated,],vertex.size=5)

#Node topological descriptives
degree<-igraph::centr_degree(Mre.IR5,mode = "all",normalized = T)
degree<-degree$res
eigen_cent<-igraph::eigen_centrality(Mre.IR5)
eigen_cent<-eigen_cent$vector
clos_cent <- igraph::closeness(Mre.IR5, normalized=T)
betw_cent <- igraph::betweenness(Mre.IR5, normalized=T)
hs <- hub_score(Mre.IR5, weights=E(Mre.IR5)$weight)$vector
as <- authority_score(Mre.IR5, weights=E(Mre.IR5)$weight)$vector
Trans<-transitivity(Mre.IR5,type="local")
Trans
#write results
Mre.IR5.Centrality<-data.frame(NAME=V(Mre.IR5)$name, UniprotID=V(Mre.IR5)$Uniprot, Degree=degree, Eigen= eigen_cent, Closeness= clos_cent, Betweeness= betw_cent, Hubs=hs, Authority=as,Transitivity=Trans)
Mre.IR5.Centrality
write.csv(Mre.IR5.Centrality,"./Mre.IR5.Centrality.csv")

#Cocitation, Number of shared interactions among nodes
Cocitation<-cocitation(Mre.IR5)
Cocitation#This matrix is similar to an adjacency matrix
write.csv(Cocitation,"./Mre.IR5.cocitation.csv")
#Plot with degree as node size
ggraph(Mre.IR5,layout=Lyt.FR[-isolated,])+
  geom_edge_fan(color="gray85", width=1, alpha=0.4)+
  geom_node_point(aes(color=V(Mre.IR5)$Annotation, size=degree),alpha=0.75)+
  geom_node_text(aes(label=name),color="gray25",size=4)+
  scale_color_manual(limits=V(Mre.IR5)$Annotation, values=V(Mre.IR5)$An.color)+
  scale_size_continuous(range = c(5, 15)) +
  theme_void()+
  title("Mre11 IR5 Network")+
  guides(color=guide_legend(override.aes = list(size=7)))+
  theme(legend.text=element_text(size=12,face = "bold"),
        legend.title = element_text(size=13,face = "bold"))


#Modularity Analysis

#Edge_Betweeness Modularity
Mre.IR5.ebmod=edge.betweenness.community(Mre.IR5,weights=E(Mre.IR5)$weight,directed = F) 
Mre.IR5.ebmod# IGRAPH clustering edge betweenness, groups: 8, mod: 0.43

#
plot.igraph(Mre.IR5.ebmod,Mre.IR5,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(Mre.IR5)$An.color, vertex.frame.color=V(Mre.IR5)$An.color, vertex.shape="circle",
            vertex.size=5, vertex.label=V(Mre.IR5)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(Mre.IR5)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)

#
plot(Mre.IR5.ebmod, Mre.IR5,layout=Lyt.FR[-isolated,], alpha=0.6, 
     vertex.color=V(Mre.IR5)$An.color, vertex.frame.color=V(Mre.IR5)$An.color, vertex.shape="circle", 
     vertex.size=5, vertex.label=V(Mre.IR5)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black", 
     edge.width=E(Mre.IR5)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout", 
     rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)
#Dendogram of the Edge_Betweeness
dendPlot(Mre.IR5.ebmod, mode="hclust")

#Louvain modularity
Mre.IR5.Louvain.mod<-cluster_louvain(Mre.IR5,weights=E(Mre.IR5)$weight)
Mre.IR5.Louvain.mod #IGRAPH clustering multi level, groups: 9, mod: 0.48

#Plot the Louvain communities
plot(Mre.IR5.Louvain.mod, Mre.IR5,layout=Lyt.FR[-isolated,], alpha=0.6, 
     vertex.color=V(Mre.IR5)$An.color, vertex.frame.color=V(Mre.IR5)$An.color, vertex.shape="circle", 
     vertex.size=5, vertex.label=V(Mre.IR5)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black", 
     edge.width=E(Mre.IR5)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout", 
     rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#Generate HeatMaps of the adjacency matrixes
#Let's use the neighbor list for generating the Adjaency Matrix
SensorID_FNN #Let's use the edge list that contain the total amount of possible interactors for the three sensors
v<-levels(as.factor(unlist(SensorID_FNN[,1:2]))) #get vertex names from the network list
n<-length(v)#get the number of vertex in the network
e<-matrix(match(as.character(unlist(SensorID_FNN[,1:2])),v),ncol = 2)#Edgelist with row index instead of Names
w<-Mre_IR5$weight #Edge weights, only from the Mre_IR5 network to quit those edges that do not exist at Mre_IR5

# M is our co-expression network adjacency matrix.

# distringuish between i,j and j,i we want to make M symmetric (undirected) by
# considering an undirected network.
Mre.IR5.M <- matrix(0, n, n)                  # set up a co-expression matrix
Mre.IR5.M[e] <- w                             # fill it in with edge weights
Mre.IR5.M <- Mre.IR5.M + t(Mre.IR5.M)                         # make this symmetric
any(is.na(Mre.IR5.M))#check no NAs in the matrix
dimnames(Mre.IR5.M) <- list(v, v)             # label the vertices
# We let A_ij = 1 if there is an interaction. 
#A <- 1*(M > 0) #right now we keep the weights as they are
Mre.IR5.M #The adjacency Matrix not as an igraph object, easier to manipulate, and to subset for smaller Sensor specific networks
quantile(Mre.IR5.M,c(0,0.25,0.5,0.75,1))
#
# 0%   25%   50%   75%  100% 
#0.000 0.000 0.000 0.000 0.991  

#Creating the HeatMap
# make a heatmap image with colors: white below 0,
# grey from 0.4 to 1 with darker closer to 1
breaks <- c(0, seq(0.4, 0.991, length=11)) # breaks in the color bins
cols <- grey(1-c(0,seq(0.5,1,length=10))) # colors

#Red Blue Colors
library(RColorBrewer)
breaks <- c(seq(-1, 1, length=11))
breaks
cols <- gray(1-c(0,seq(0,1,length=10))) # colors
cols2 <- hcl.colors(10 , palette = "Blues3" , alpha = 1 , rev = TRUE)
cols3 <- colorRamp("red" , space = "rgb" , interpolate = "linear")
HMCol<- brewer.pal(10,name = "RdBu")
HMCol<- rev(HMCol)
HMCol


# color the vertex for CCNB1 blue, its neighbors red, and the others white
SensorIDBind <- which(v == "Mre11") # index for ccnb1
#SensorIDBind<-as.factor(SensorIDBind)
vcols <- rep("white",n)
vcols[SensorIDBind] <- "blue"
#Subsetting Mre.IR5.M by doing a loop asking for look into every row and column the index in SensorIDBind worked
for (i in SensorIDBind) {
  vcols[which(Mre.IR5.M[,i] > 0 | Mre.IR5.M[i,])] <- "red"
}
##Color annotation for functional annotation
MatSet<-subset(Mat,Mat$NAME%in%v)# select only those proteins that are cointained in the node list of our matrixes
match(v,MatSet$NAME) # Check that the indexes of the Annotation Dataframe and the node list are in the same position 
match(rownames(Mre.IR5.M),MatSet$NAME) #Check that the row names of the Adj Matrix and the Annotation dataframe are in the same order
match(colnames(Mre.IR5.M),MatSet$NAME) #Check that the column names of the Adj Matrix and the Annotation dataframe are in the same order
AnCol<-MatSet$An.color#Vector of colors
AnCol<-as.character(AnCol)#transform it into a character vector, so the colors can be read

#Subsetting Mre.IR5.M using the index vector in SensorIDBind do not work
#vcols[which(Mre.IR5.M[,SensorIDBind] > 0 | Mre.IR5.M[SensorIDBind,])] <- "red"

#Generate a dendogram for keeping the proteins in the sample place
#Mre.IR5.dist<-dist(Mre.IR5.M,method = "euclidian")
#Mre.IR5.Dend<-hclust(Mre.IR5.dist, method = "ward.D" )
#plot(Mre.IR5.Dend)

#Cut the dendogram
#cut.height <- 10 # try numbers from 4 to 7
#FNN.hclusts<- cutree(FNN.Dend, cut.height)
#plot(FNN.hclusts)

# now actually make the heat map
par(mar=rep(0,4)) # make plot margins 0
heatmap(Mre.IR5.M, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=cols, breaks=breaks,frame=T, cexRow = 0.5 , cexCol = 0.5, main = "Mre11 IR5 AdjMat", ylab="Total Protein Interactors",xlab="Total Protein Interactors")
#
plot.new()
legend("topleft", c("Cell Cycle","Chromatin Org","DNA Repair","Nuclear Pore","p53","Proteasome","Ribosome Biogenesis","RNA Metabolism","Transcription Reg","UPR"), 
       fill=c("coral1","darkviolet","firebrick1","deepskyblue3","tomato","grey","lightgreen","gold","slateblue1","maroon"),  
       bty="n", inset=0, xpd=T,  border=F, cex =2,y.intersp = 0.7,x.intersp = 0.2,pt.cex = 1)


#Now we substract the MrenonIR5 AdjMatrix to the MreIR6 AdjMatrix
MreIR5.AdjMat.Change<-Mre.IR5.M-Mre.nonIR5.M
range(MreIR5.AdjMat.Change)#Let's check the magnitude of the substracted weights
#[1] -0.807  0.805
any(is.na(MreIR5.AdjMat.Change))
#Generate the breaks and color functions
breaks <- c(seq(-0.5, 0.5, length=6))
breaks
HMCol<- brewer.pal(5,name = "RdBu")
HMCol<- rev(HMCol)
HMCol
[1] "#2166AC" "#4393C3" "#92C5DE" "#D1E5F0" "#F7F7F7" "#FDDBC7" "#F4A582" "#D6604D" "#B2182B"

#Plot the HeatMap
par(mar=rep(0,4)) # make plot margins 0
heatmap(MreIR5.AdjMat.Change, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=HMCol, breaks=breaks,frame=T, cexRow = 0.5 , cexCol = 0.5, main = "Mre11 IR5 AdjMat (minus nonIR5 edges)", ylab="Total Protein Interactors",xlab="Total Protein Interactors")
#
#check with the FNN
# now actually make the heat map
par(mar=rep(0,4)) # make plot margins 0
heatmap(FNN.M, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=0:1, frame=T, cexRow = 0.5 , cexCol = 0.5)
#
legend("topleft", c("Neighbors Mre11","Mre11"), fill=c("red","blue"),  
       bty="n", inset=0, xpd=T,  border=F)

##########################################Mre11 IR30
##############################################
#Mre11 IR30 Networks


#Let's keep the Mat as the Node attribute table and load the Neighbor list as the matrixes
Mat# Node attribute matrix for all the networks

#Mre11ID IR30 nodes Network Edge Lists
Mre_IR30<-read.delim("./Mre11_Files/Mre_IR30.tsv",sep = "\t",col.names = )
str(Mre_IR30)
head(Mre_IR30)

#Generate the igraph network object
Mre.IR30<-graph_from_data_frame(d=Mre_IR30,vertices = Mat, directed = F )#Generate the igraph object with the edge and node annotation attributes
str(vertex_attr(Mre.IR30))#Check node attributes
str(edge_attr(Mre.IR30))#Check for edge attributes
Mre.IR30<-igraph::simplify(Mre.IR30, remove.multiple=TRUE, remove.loops=TRUE) #Remove repetitive edges and self-connecting edges
Mre.IR30
#Small plot
plot(Mre.IR30)

#Generate an adjacency Matrix that will be used after
Mre.IR30_AdjMat<-as_adjacency_matrix(Mre.IR30,type="both",attr = "weight",names=TRUE,sparse = TRUE)
Mre.IR30_AdjMat<-as.matrix(Mre.IR30_AdjMat)
write.csv(Mre.IR30_AdjMat,"./Mre.IR30_AdjMat.csv")

#Generate a vertex size based on degree and count
vsize<-igraph::degree(Mre.IR30)/vcount(Mre.IR30)
#Identify and remove disconnected nodes
isolated<-which(igraph::degree(Mre.IR30)==0)
isolated
length(isolated)# [1] 623
#remove isolated nodes
Mre.IR30<-igraph::delete_vertices(Mre.IR30,isolated)


#Generate layouts for the graph
#Do not run again, to not to change the layouts that I already have
Lyt.Nice<-layout_nicely(Mre.IR30,dim = 3)
Lyt.GOpt<-layout_with_graphopt(Mre.IR30, niter=5000, mass = 30, charge = E(Mre.IR30)$weight)
Lyt.DrL<-layout_with_drl(Mre.IR30, weights = E(Mre.IR30)$weight)
Lyt.FR<-layout.fruchterman.reingold(Mre.IR30, niter=10000, weights=E(Mre.IR30)$weight)#This one looks like the best option
Lyt.KK<-layout_with_kk(Mre.IR30, maxiter=10000, weights=E(Mre.IR30)$weight) #very large separation of peripheral nodes
Lyt.GEM<-layout_with_gem(Mre.IR30, maxiter=10000) #Bad separation
#Lyt.DH<-layout_with_dh(Mre.IR30, maxiter=150, fineiter=max(10, log2(833)))
Lyt.LGL<-layout.lgl(Mre.IR30, maxit=150,maxdelta=833,area=(length(E(Mre.IR30)))^2,coolexp=1.5,repulserad=((length(E(Mre.IR30)))^2)*(length(E(Mre.IR30))))#Do not work


#Star and Tree layouts

V(Mre.IR30)[name=="MRE11A"]
which(V(Mre.IR30)[name=="Mre11"])# 228 #Check what is the rowindex for your desired node
which(V(Mre.IR30)$name=="MRE11A")# 191 #Check what is the rowindex for your desired node
which(V(Mre.IR30)$name=="MRE11A")# 327
which(V(Mre.IR30)$name=="MRE11A")# 213
Lyt.Star<-layout_as_star(Mre.IR30,center=V(Mre.IR30)[name=="MRE11A"]) #it has a special order argument where we can assign an order to the vertex, hence maybe the OPLS-DA:VIP value order may be used as ordering vector

Lyt.Tree<-layout_as_tree(Mre.IR30,root=c(84))
Lyt.RT<-layout.reingold.tilford(Mre.IR30,root=c(84))# This one looks good, but the indexing do not allow to highlight the correct nodes
Lyt.sugi<-layout_with_sugiyama(Mre.IR30, maxiter=100, weights=E(Mre.IR302)$weight)#This one looks good for undirected tree layouts

#Check what layout looks fine

plot.igraph(Mre.IR30,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(Mre.IR30)$An.color, vertex.frame.color=V(Mre.IR30)$An.color, vertex.shape="circle",
            vertex.size=5, vertex.label=V(Mre.IR30)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(Mre.IR30)$weight, edge.lty=1, main="Mre11ID IR30 Network", sub="FR Tree Layout",
            rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#With ggraph

ggraph(Mre.IR30,layout=Lyt.FR[-isolated,])+
  geom_edge_fan(color="gray65", width=1, alpha=0.4)+
  geom_node_point(aes(color=V(Mre.IR30)$Annotation, size=V(Mre.IR30)$Mre.FC.IR30),alpha=0.75)+
  geom_node_text(aes(label=name),color="gray25",size=4)+
  scale_color_manual(limits=V(Mre.IR30)$Annotation, values=V(Mre.IR30)$An.color)+
  scale_size_continuous(range = c(5, 15)) +
  theme_void()+
  title("Mre11 IR30 Network")+
  guides(color=guide_legend(override.aes = list(size=7)))+
  theme(legend.text=element_text(size=12,face = "bold"),
        legend.title = element_text(size=13,face = "bold"))



#Community Detection
#clp <- igraph::cluster_optimal(Mre.IR30)
class(clp)
V(Mre.IR30)$community <- clp$membership
colrs <- adjustcolor( c("gray50", "tomato", "gold", "yellowgreen"), alpha=.6)


plot.igraph(clp,Mre.IR30,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(Mre.IR30)$An.color, vertex.frame.color=V(Mre.IR30)$An.color, vertex.shape="circle",
            vertex.size==V(Mre.IR30)$Mre.FC.IR30, vertex.label=V(Mre.IR30)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(Mre.IR30)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#Distances edges highlight
dist.from.Mre11 <- distances(Mre.IR30, v=V(Mre.IR30)[name=="MRE11A"],
                             to=V(Mre.IR30), weights=E(Mre.IR30)$weight)
dist.from.Mre11<-dist.from.Mre11[,!is.infinite(colSums(dist.from.Mre11))]

dist.from.Mre11
range(dist.from.Mre11)

#Set colors to plot the distances:
oranges <- colorRampPalette(c("dark red", "gold"))
col <- oranges(max(dist.from.Mre11))
col <- col[dist.from.Mre11+1]
plot(Mre.IR30, layout=Lyt.FR[-isolated,],vertex.color=col, vertex.label=dist.from.Mre11,
     vertex.label.color="white",vertex.shape="circle",
     vertex.size=10, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(Mre.IR30)$weight, edge.lty=1, main="Mre11 IR30 Full Nuclear Network", sub="FR Layout",
     rescale=T)


#Shortest path from Mre11 to Nucleolin
MreNcl.path <- shortest_paths(Mre.IR30,
                              from = V(Mre.IR30)[name=="MRE11A"],
                              to = V(Mre.IR30)[name=="NCL"],
                              output = "both") # both path nodes and edges

MreNcl.path
# Generate edge color variable to plot the path:
ecol <- rep("gray80", ecount(Mre.IR30))
ecol[unlist(MreNcl.path$epath)] <- "orange"
# Generate edge width variable to plot the path:
ew <- rep(2, ecount(Mre.IR30))
ew[unlist(MreNcl.path$epath)] <- 4
# Generate node color variable to plot the path:
vcol <- rep("gray40", vcount(Mre.IR30))
vcol[unlist(MreNcl.path$vpath)] <- "gold"
plot(Mre.IR30, layout=Lyt.FR[-isolated,], vertex.color=vcol, edge.color=ecol,
     edge.width=ew,vertex.label.color="white",vertex.shape="circle",
     vertex.size=5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(Mre.IR30)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)



#Plot based on Mre11 as incident
inc.edges <- incident(Mre.IR30, V(Mre.IR30)[name=="MRE11A"], mode="all")
# Set colors to plot the selected edges.
ecol <- rep("gray80", ecount(Mre.IR30))
ecol[inc.edges] <- "orange"
vcol <- rep("gray40", vcount(Mre.IR30))
vcol[V(Mre.IR30)$name=="Mre11"] <- "gold"

plot(Mre.IR30, layout=Lyt.FR[-isolated,], vertex.color=vcol, edge.color=ecol, 
     edge.width=ew,vertex.label.color="white",vertex.shape="circle",
     vertex.size=5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(Mre.IR30)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)


#Plot Degree Distribution
Mre.IR30.Deg<-degree_distribution(Mre.IR30, cumulative = T, mode="all")
plot( x=0:max(igraph::degree(Mre.IR30)), y=1-Mre.IR30.Deg, pch=19, cex=1.2, col="orange",
      xlab="Degree", ylab="Cumulative Frequency")

#Network Descriptives
#Number of nodes and edges
vcount(Mre.IR30) #[1] 210
ecount(Mre.IR30)#[1]  1293
#Density (Number of edges in relation to the total number of edges in a full connected graph)
edge_density(Mre.IR30,loops = F) #[1] 0.05892003
#Transitivity (Total number of conected triangles)
transitivity(Mre.IR30,type="global") #[1] 0.6413934
#modularity
mod<-V(Mre.IR30)$S6.FuzzyCMeans #lets check for fuzzy Cmeans modularity
mod<-mod+1#The non clustered interactors that appears as 0, needs to be transformed to 1, so we add 1 to every cluster
modularity(Mre.IR30,membership = mod)
#S6.FuzzyCMeans modularity [1] -0.007128933
#Ku.FuzzyCMeans modularity [1] 0.001670007
#Mre.FuzzyCMeans modularity [1] -0.00287825

#Assorativity
assortativity_nominal(Mre.IR30,factor(V(Mre.IR30)$Mre.FuzzyCMeans),directed=F)
#Annotation [1] [1] 0.1956823
#S6.FuzzyCMeans [1] -0.01089844
#Ku.FuzzyCMeans [1] 0.002973617
#Mre.FuzzyCMeans [1] -0.003949552

#Diameter (longest geodesic distance), shortest path between two most peripheral nodes
diameter(Mre.IR30,directed = F, weights = E(Mre.IR30)$weights)#[1] 4.422  
diam<-get_diameter(Mre.IR30,directed = F,weights = E(Mre.IR30)$weights)
diam #+ 9/210 vertices, named, from d1227ee: [1] ACTN1  ACTN4  CTNNB1 RUVBL2 XRCC5  TMPO   CANX   PDIA3  PPIB  
#Plot Diameter
vcol <- rep("gray50", vcount(Mre.IR30))
vcol[diam] <- "gold"
ecol <- rep("gray90", ecount(Mre.IR30))
ecol[E(Mre.IR30, path=diam)] <- "orange"
# E(net, path=diam) finds edges along a path, here 'diam'
plot(Mre.IR30, vertex.color=vcol, edge.color=ecol, vertex.label.color="black",edge.arrow.mode=0,layout=Lyt.FR[-isolated,],vertex.size=10)

#Node topological descriptives
degree<-igraph::centr_degree(Mre.IR30,mode = "all",normalized = T)
degree<-degree$res
eigen_cent<-igraph::eigen_centrality(Mre.IR30)
eigen_cent<-eigen_cent$vector
clos_cent <- igraph::closeness(Mre.IR30, normalized=T)
betw_cent <- igraph::betweenness(Mre.IR30, normalized=T)
hs <- hub_score(Mre.IR30, weights=E(Mre.IR30)$weight)$vector
as <- authority_score(Mre.IR30, weights=E(Mre.IR30)$weight)$vector
Trans<-transitivity(Mre.IR30,type="local")
Trans
#write results
Mre.IR30.Centrality<-data.frame(NAME=V(Mre.IR30)$name, UniprotID=V(Mre.IR30)$Uniprot, Degree=degree, Eigen= eigen_cent, Closeness= clos_cent, Betweeness= betw_cent, Hubs=hs, Authority=as,Transitivity=Trans)
Mre.IR30.Centrality
write.csv(Mre.IR30.Centrality,"./Mre.IR30.Centrality.csv")

#Cocitation, Number of shared interactions among nodes
Cocitation<-cocitation(Mre.IR30)
Cocitation#This matrix is similar to an adjacency matrix
write.csv(Cocitation,"./Mre.IR30.cocitation.csv")
#Plot with degree as node size
ggraph(Mre.IR30,layout=Lyt.FR[-isolated,])+
  geom_edge_fan(color="gray85", width=1, alpha=0.4)+
  geom_node_point(aes(color=V(Mre.IR30)$Annotation, size=degree),alpha=0.75)+
  geom_node_text(aes(label=name),color="gray25",size=4)+
  scale_color_manual(limits=V(Mre.IR30)$Annotation, values=V(Mre.IR30)$An.color)+
  scale_size_continuous(range = c(5, 15)) +
  theme_void()+
  title("Mre11 IR30 Network")+
  guides(color=guide_legend(override.aes = list(size=7)))+
  theme(legend.text=element_text(size=12,face = "bold"),
        legend.title = element_text(size=13,face = "bold"))


#Modularity Analysis

#Edge_Betweeness Modularity
Mre.IR30.ebmod=edge.betweenness.community(Mre.IR30,weights=E(Mre.IR30)$weight,directed = F) 
Mre.IR30.ebmod# IGRAPH clustering edge betweenness, groups: 29, mod: 0.37

#
plot.igraph(Mre.IR30.ebmod,Mre.IR30,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(Mre.IR30)$An.color, vertex.frame.color=V(Mre.IR30)$An.color, vertex.shape="circle",
            vertex.size=5, vertex.label=V(Mre.IR30)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(Mre.IR30)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)

#
plot(Mre.IR30.ebmod, Mre.IR30,layout=Lyt.FR[-isolated,], alpha=0.6, 
     vertex.color=V(Mre.IR30)$An.color, vertex.frame.color=V(Mre.IR30)$An.color, vertex.shape="circle", 
     vertex.size=5, vertex.label=V(Mre.IR30)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black", 
     edge.width=E(Mre.IR30)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout", 
     rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)
#Dendogram of the Edge_Betweeness
dendPlot(Mre.IR30.ebmod, mode="hclust")

#Louvain modularity
Mre.IR30.Louvain.mod<-cluster_louvain(Mre.IR30,weights=E(Mre.IR30)$weight)
Mre.IR30.Louvain.mod #IGRAPH clustering multi level, groups: 12, mod: 0.42

#Plot the Louvain communities
plot(Mre.IR30.Louvain.mod, Mre.IR30,layout=Lyt.FR[-isolated,], alpha=0.6, 
     vertex.color=V(Mre.IR30)$An.color, vertex.frame.color=V(Mre.IR30)$An.color, vertex.shape="circle", 
     vertex.size=5, vertex.label=V(Mre.IR30)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black", 
     edge.width=E(Mre.IR30)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout", 
     rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#Generate HeatMaps of the adjacency matrixes
#Let's use the neighbor list for generating the Adjaency Matrix
SensorID_FNN #Let's use the edge list that contain the total amount of possible interactors for the three sensors
v<-levels(as.factor(unlist(SensorID_FNN[,1:2]))) #get vertex names from the network list
n<-length(v)#get the number of vertex in the network
e<-matrix(match(as.character(unlist(SensorID_FNN[,1:2])),v),ncol = 2)#Edgelist with row index instead of Names
w<-Mre_IR30$weight #Edge weights, only from the Mre_IR30 network to quit those edges that do not exist at Mre_IR30

# M is our co-expression network adjacency matrix.

# distringuish between i,j and j,i we want to make M symmetric (undirected) by
# considering an undirected network.
Mre.IR30.M <- matrix(0, n, n)                  # set up a co-expression matrix
Mre.IR30.M[e] <- w                             # fill it in with edge weights
Mre.IR30.M <- Mre.IR30.M + t(Mre.IR30.M)                         # make this symmetric
any(is.na(Mre.IR30.M))#check no NAs in the matrix
dimnames(Mre.IR30.M) <- list(v, v)             # label the vertices
# We let A_ij = 1 if there is an interaction. 
#A <- 1*(M > 0) #right now we keep the weights as they are
Mre.IR30.M #The adjacency Matrix not as an igraph object, easier to manipulate, and to subset for smaller Sensor specific networks
quantile(Mre.IR30.M,c(0,0.25,0.5,0.75,1))
#
# 0%   25%   50%   75%  100% 
#0.000 0.000 0.000 0.000 0.991  

#Creating the HeatMap
# make a heatmap image with colors: white below 0,
# grey from 0.4 to 1 with darker closer to 1
breaks <- c(0, seq(0.4, 0.991, length=11)) # breaks in the color bins
cols <- grey(1-c(0,seq(0.5,1,length=10))) # colors

#Red Blue Colors
library(RColorBrewer)
breaks <- c(seq(-1, 1, length=11))
breaks
cols <- gray(1-c(0,seq(0,1,length=10))) # colors
cols2 <- hcl.colors(10 , palette = "Blues3" , alpha = 1 , rev = TRUE)
cols3 <- colorRamp("red" , space = "rgb" , interpolate = "linear")
HMCol<- brewer.pal(10,name = "RdBu")
HMCol<- rev(HMCol)
HMCol


# color the vertex for CCNB1 blue, its neighbors red, and the others white
SensorIDBind <- which(v == "Mre11") # index for ccnb1
#SensorIDBind<-as.factor(SensorIDBind)
vcols <- rep("white",n)
vcols[SensorIDBind] <- "blue"
#Subsetting Mre.IR30.M by doing a loop asking for look into every row and column the index in SensorIDBind worked
for (i in SensorIDBind) {
  vcols[which(Mre.IR30.M[,i] > 0 | Mre.IR30.M[i,])] <- "red"
}
##Color annotation for functional annotation
MatSet<-subset(Mat,Mat$NAME%in%v)# select only those proteins that are cointained in the node list of our matrixes
match(v,MatSet$NAME) # Check that the indexes of the Annotation Dataframe and the node list are in the same position 
match(rownames(Mre.IR30.M),MatSet$NAME) #Check that the row names of the Adj Matrix and the Annotation dataframe are in the same order
match(colnames(Mre.IR30.M),MatSet$NAME) #Check that the column names of the Adj Matrix and the Annotation dataframe are in the same order
AnCol<-MatSet$An.color#Vector of colors
AnCol<-as.character(AnCol)#transform it into a character vector, so the colors can be read

#Subsetting Mre.IR30.M using the index vector in SensorIDBind do not work
#vcols[which(Mre.IR30.M[,SensorIDBind] > 0 | Mre.IR30.M[SensorIDBind,])] <- "red"

#Generate a dendogram for keeping the proteins in the sample place
#Mre.IR30.dist<-dist(Mre.IR30.M,method = "euclidian")
#Mre.IR30.Dend<-hclust(Mre.IR30.dist, method = "ward.D" )
#plot(Mre.IR30.Dend)

#Cut the dendogram
#cut.height <- 10 # try numbers from 4 to 7
#FNN.hclusts<- cutree(FNN.Dend, cut.height)
#plot(FNN.hclusts)

# now actually make the heat map
par(mar=rep(0,4)) # make plot margins 0
heatmap(Mre.IR30.M, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=cols, breaks=breaks,frame=T, cexRow = 0.5 , cexCol = 0.5, main = "Mre11 IR30 AdjMat", ylab="Total Protein Interactors",xlab="Total Protein Interactors")
#
plot.new()
legend("topleft", c("Cell Cycle","Chromatin Org","DNA Repair","Nuclear Pore","p53","Proteasome","Ribosome Biogenesis","RNA Metabolism","Transcription Reg","UPR"), 
       fill=c("coral1","darkviolet","firebrick1","deepskyblue3","tomato","grey","lightgreen","gold","slateblue1","maroon"),  
       bty="n", inset=0, xpd=T,  border=F, cex =2,y.intersp = 0.7,x.intersp = 0.2,pt.cex = 1)


#Now we substract the MrenonIR30 AdjMatrix to the MreIR6 AdjMatrix
MreIR30.AdjMat.Change<-Mre.IR30.M-Mre.nonIR5.M
range(MreIR30.AdjMat.Change)#Let's check the magnitude of the substracted weights
#[1]-0.813  0.668
any(is.na(MreIR30.AdjMat.Change))
#Generate the breaks and color functions
breaks <- c(seq(-0.5, 0.5, length=6))
breaks
HMCol<- brewer.pal(5,name = "RdBu")
HMCol<- rev(HMCol)
HMCol
[1] "#2166AC" "#4393C3" "#92C5DE" "#D1E5F0" "#F7F7F7" "#FDDBC7" "#F4A582" "#D6604D" "#B2182B"

#Plot the HeatMap
par(mar=rep(0,4)) # make plot margins 0
heatmap(MreIR30.AdjMat.Change, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=HMCol, breaks=breaks,frame=T, cexRow = 0.5 , cexCol = 0.5, main = "Mre11 IR30 AdjMat (minus nonIR30 edges)", ylab="Total Protein Interactors",xlab="Total Protein Interactors")


#Now we substract the MreIR5 AdjMatrix to the MreIR30 AdjMatrix
MreIR30.AdjMat.Change2<-Mre.IR30.M-Mre.IR5.M
range(MreIR30.AdjMat.Change)#Let's check the magnitude of the substracted weights
#[1]-0.813  0.668
any(is.na(MreIR30.AdjMat.Change))
#Generate the breaks and color functions
breaks <- c(seq(-0.5, 0.5, length=6))
breaks
HMCol<- brewer.pal(5,name = "RdBu")
HMCol<- rev(HMCol)
HMCol
[1] "#2166AC" "#4393C3" "#92C5DE" "#D1E5F0" "#F7F7F7" "#FDDBC7" "#F4A582" "#D6604D" "#B2182B"

#Plot the HeatMap
par(mar=rep(0,4)) # make plot margins 0
heatmap(MreIR30.AdjMat.Change2, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=HMCol, breaks=breaks,frame=T, cexRow = 0.5 , cexCol = 0.5, main = "Mre11 IR30 AdjMat (minus nonIR30 edges)", ylab="Total Protein Interactors",xlab="Total Protein Interactors")

####################### Mre11 IR2h Networks
###################################################

#Mre11 IR2 Networks


#Let's keep the Mat as the Node attribute table and load the Neighbor list as the matrixes
Mat# Node attribute matrix for all the networks

#Mre11ID IR2 nodes Network Edge Lists
Mre_IR2<-read.delim("./Mre11_Files/Mre_IR2h.tsv",sep = "\t",col.names = )
str(Mre_IR2)
head(Mre_IR2)

#Generate the igraph network object
Mre.IR2<-graph_from_data_frame(d=Mre_IR2,vertices = Mat, directed = F )#Generate the igraph object with the edge and node annotation attributes
str(vertex_attr(Mre.IR2))#Check node attributes
str(edge_attr(Mre.IR2))#Check for edge attributes
Mre.IR2<-igraph::simplify(Mre.IR2, remove.multiple=TRUE, remove.loops=TRUE) #Remove repetitive edges and self-connecting edges
Mre.IR2
#Small plot
plot(Mre.IR2)

#Generate an adjacency Matrix that will be used after
Mre.IR2_AdjMat<-as_adjacency_matrix(Mre.IR2,type="both",attr = "weight",names=TRUE,sparse = TRUE)
Mre.IR2_AdjMat<-as.matrix(Mre.IR2_AdjMat)
write.csv(Mre.IR2_AdjMat,"./Mre.IR2_AdjMat.csv")

#Generate a vertex size based on degree and count
vsize<-igraph::degree(Mre.IR2)/vcount(Mre.IR2)
#Identify and remove disconnected nodes
isolated<-which(igraph::degree(Mre.IR2)==0)
isolated
length(isolated)# [1] 696
#remove isolated nodes
Mre.IR2<-igraph::delete_vertices(Mre.IR2,isolated)


#Generate layouts for the graph
#Do not run again, to not to change the layouts that I already have
Lyt.Nice<-layout_nicely(Mre.IR2,dim = 3)
Lyt.GOpt<-layout_with_graphopt(Mre.IR2, niter=5000, mass = 30, charge = E(Mre.IR2)$weight)
Lyt.DrL<-layout_with_drl(Mre.IR2, weights = E(Mre.IR2)$weight)
Lyt.FR<-layout.fruchterman.reingold(Mre.IR2, niter=10000, weights=E(Mre.IR2)$weight)#This one looks like the best option
Lyt.KK<-layout_with_kk(Mre.IR2, maxiter=10000, weights=E(Mre.IR2)$weight) #very large separation of peripheral nodes
Lyt.GEM<-layout_with_gem(Mre.IR2, maxiter=10000) #Bad separation
#Lyt.DH<-layout_with_dh(Mre.IR2, maxiter=150, fineiter=max(10, log2(833)))
Lyt.LGL<-layout.lgl(Mre.IR2, maxit=150,maxdelta=833,area=(length(E(Mre.IR2)))^2,coolexp=1.5,repulserad=((length(E(Mre.IR2)))^2)*(length(E(Mre.IR2))))#Do not work


#Star and Tree layouts

V(Mre.IR2)[name=="MRE11A"]
which(V(Mre.IR2)[name=="Mre11"])# 68 #Check what is the rowindex for your desired node
which(V(Mre.IR2)$name=="Mre11")# 191 #Check what is the rowindex for your desired node
which(V(Mre.IR2)$name=="NCL")# 0
which(V(Mre.IR2)$name=="MRE11A")# 68
Lyt.Star<-layout_as_star(Mre.IR2,center=V(Mre.IR2)[name=="MRE11A"]) #it has a special order argument where we can assign an order to the vertex, hence maybe the OPLS-DA:VIP value order may be used as ordering vector

Lyt.Tree<-layout_as_tree(Mre.IR2,root=c(68))
Lyt.RT<-layout.reingold.tilford(Mre.IR2,root=c(68))# This one looks good, but the indexing do not allow to highlight the correct nodes
Lyt.sugi<-layout_with_sugiyama(Mre.IR2, maxiter=100, weights=E(Mre.IR22)$weight)#This one looks good for undirected tree layouts

#Check what layout looks fine

plot.igraph(Mre.IR2,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(Mre.IR2)$An.color, vertex.frame.color=V(Mre.IR2)$An.color, vertex.shape="circle",
            vertex.size=5, vertex.label=V(Mre.IR2)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(Mre.IR2)$weight, edge.lty=1, main="Mre11ID IR2 Network", sub="FR Tree Layout",
            rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#With ggraph

ggraph(Mre.IR2,layout=Lyt.FR[-isolated,])+
  geom_edge_fan(color="gray65", width=1, alpha=0.4)+
  geom_node_point(aes(color=V(Mre.IR2)$Annotation, size=V(Mre.IR2)$Mre.FC.IR2),alpha=0.75)+
  geom_node_text(aes(label=name),color="gray25",size=4)+
  scale_color_manual(limits=V(Mre.IR2)$Annotation, values=V(Mre.IR2)$An.color)+
  scale_size_continuous(range = c(5, 15)) +
  theme_void()+
  title("Mre11 IR2 Network")+
  guides(color=guide_legend(override.aes = list(size=7)))+
  theme(legend.text=element_text(size=12,face = "bold"),
        legend.title = element_text(size=13,face = "bold"))



#Community Detection
#clp <- igraph::cluster_optimal(Mre.IR2)
class(clp)
V(Mre.IR2)$community <- clp$membership
colrs <- adjustcolor( c("gray50", "tomato", "gold", "yellowgreen"), alpha=.6)


plot.igraph(clp,Mre.IR2,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(Mre.IR2)$An.color, vertex.frame.color=V(Mre.IR2)$An.color, vertex.shape="circle",
            vertex.size==V(Mre.IR2)$Mre.FC.IR2, vertex.label=V(Mre.IR2)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(Mre.IR2)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#Distances edges highlight
dist.from.Mre11 <- distances(Mre.IR2, v=V(Mre.IR2)[name=="MRE11A"],
                             to=V(Mre.IR2), weights=E(Mre.IR2)$weight)
dist.from.Mre11<-dist.from.Mre11[,!is.infinite(colSums(dist.from.Mre11))]

dist.from.Mre11
range(dist.from.Mre11)

#Set colors to plot the distances:
oranges <- colorRampPalette(c("dark red", "gold"))
col <- oranges(max(dist.from.Mre11))
col <- col[dist.from.Mre11+1]
plot(Mre.IR2, layout=Lyt.FR[-isolated,],vertex.color=col, vertex.label=dist.from.Mre11,
     vertex.label.color="white",vertex.shape="circle",
     vertex.size=5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(Mre.IR2)$weight, edge.lty=1, main="Mre11 IR2 Full Nuclear Network", sub="FR Layout",
     rescale=T)


#Shortest path from Mre11 to Nucleolin
MreNcl.path <- shortest_paths(Mre.IR2,
                              from = V(Mre.IR2)[name=="MRE11A"],
                              to = V(Mre.IR2)[name=="NCL"],
                              output = "both") # both path nodes and edges
MreNcl.path
# Generate edge color variable to plot the path:
ecol <- rep("gray80", ecount(Mre.IR2))
ecol[unlist(MreNcl.path$epath)] <- "orange"
# Generate edge width variable to plot the path:
ew <- rep(2, ecount(Mre.IR2))
ew[unlist(MreNcl.path$epath)] <- 4
# Generate node color variable to plot the path:
vcol <- rep("gray40", vcount(Mre.IR2))
vcol[unlist(MreNcl.path$vpath)] <- "gold"
plot(Mre.IR2, layout=Lyt.FR[-isolated,], vertex.color=vcol, edge.color=ecol,
     edge.width=ew,vertex.label.color="white",vertex.shape="circle",
     vertex.size=5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(Mre.IR2)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)



#Plot based on Mre11 as incident
inc.edges <- incident(Mre.IR2, V(Mre.IR2)[name=="MRE11A"], mode="all")
inc.edges
# Set colors to plot the selected edges.
ecol <- rep("gray80", ecount(Mre.IR2))
ecol[inc.edges] <- "orange"
vcol <- rep("gray40", vcount(Mre.IR2))
vcol[V(Mre.IR2)$name=="MRE11A"] <- "gold"

plot(Mre.IR2, layout=Lyt.FR[-isolated,], vertex.color=vcol, edge.color=ecol, 
     edge.width=ew,vertex.label.color="white",vertex.shape="circle",
     vertex.size=5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(Mre.IR2)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)


#Plot Degree Distribution
Mre.IR2.Deg<-degree_distribution(Mre.IR2, cumulative = T, mode="all")
plot( x=0:max(igraph::degree(Mre.IR2)), y=1-Mre.IR2.Deg, pch=19, cex=1.2, col="orange",
      xlab="Degree", ylab="Cumulative Frequency")

#Network Descriptives
#Number of nodes and edges
vcount(Mre.IR2) #[1] 295
ecount(Mre.IR2)#[1]  423
#Density (Number of edges in relation to the total number of edges in a full connected graph)
edge_density(Mre.IR2,loops = F) #[1]0.04540575
#Transitivity (Total number of conected triangles)
transitivity(Mre.IR2,type="global") #[1] 0.3811537
#modularity
mod<-V(Mre.IR2)$Mre.FuzzyCMeans #lets check for fuzzy Cmeans modularity
mod<-mod+1#The non clustered interactors that appears as 0, needs to be transformed to 1, so we add 1 to every cluster
modularity(Mre.IR2,membership = mod)
#S6.FuzzyCMeans modularity [1]-0.01500875
#Ku.FuzzyCMeans modularity [1] 0.00413013
#Mre.FuzzyCMeans modularity [1] 0.01087862

#Assorativity
assortativity_nominal(Mre.IR2,factor(V(Mre.IR2)$Annotation),directed=F)
#Annotation [1] [1] 0.3191639
#S6.FuzzyCMeans [1] -0.02417159
#Ku.FuzzyCMeans [1] 0.009355141
#Mre.FuzzyCMeans [1]  0.0160001

#Diameter (longest geodesic distance), shortest path between two most peripheral nodes
diameter(Mre.IR2,directed = F, weights = E(Mre.IR2)$weights)#[1] 4.221
diam<-get_diameter(Mre.IR2,directed = F,weights = E(Mre.IR2)$weights)
diam #++ 8/137 vertices, named, from d3c6453: [1] UTP14A DCAF13 IFI16  XRCC5  YBX1   PRPF19 SNRPF  PPIH  
#Plot Diameter
vcol <- rep("gray50", vcount(Mre.IR2))
vcol[diam] <- "gold"
ecol <- rep("gray90", ecount(Mre.IR2))
ecol[E(Mre.IR2, path=diam)] <- "orange"
# E(net, path=diam) finds edges along a path, here 'diam'
plot(Mre.IR2, vertex.color=vcol, edge.color=ecol, vertex.label.color="black",edge.arrow.mode=0,layout=Lyt.FR[-isolated,],vertex.size=5)

#Node topological descriptives
degree<-igraph::centr_degree(Mre.IR2,mode = "all",normalized = T)
degree<-degree$res
eigen_cent<-igraph::eigen_centrality(Mre.IR2)
eigen_cent<-eigen_cent$vector
clos_cent <- igraph::closeness(Mre.IR2, normalized=T)
betw_cent <- igraph::betweenness(Mre.IR2, normalized=T)
hs <- hub_score(Mre.IR2, weights=E(Mre.IR2)$weight)$vector
as <- authority_score(Mre.IR2, weights=E(Mre.IR2)$weight)$vector
Trans<-transitivity(Mre.IR2,type="local")
Trans
#write results
Mre.IR2.Centrality<-data.frame(NAME=V(Mre.IR2)$name, UniprotID=V(Mre.IR2)$Uniprot, Degree=degree, Eigen= eigen_cent, Closeness= clos_cent, Betweeness= betw_cent, Hubs=hs, Authority=as,Transitivity=Trans)
Mre.IR2.Centrality
write.csv(Mre.IR2.Centrality,"./Mre.IR2.Centrality.csv")

#Cocitation, Number of shared interactions among nodes
Cocitation<-cocitation(Mre.IR2)
Cocitation#This matrix is similar to an adjacency matrix
write.csv(Cocitation,"./Mre.IR2.cocitation.csv")
#Plot with degree as node size
ggraph(Mre.IR2,layout=Lyt.FR[-isolated,])+
  geom_edge_fan(color="gray85", width=1, alpha=0.4)+
  geom_node_point(aes(color=V(Mre.IR2)$Annotation, size=degree),alpha=0.75)+
  geom_node_text(aes(label=name),color="gray25",size=4)+
  scale_color_manual(limits=V(Mre.IR2)$Annotation, values=V(Mre.IR2)$An.color)+
  scale_size_continuous(range = c(5, 15)) +
  theme_void()+
  title("Mre11 IR2 Network")+
  guides(color=guide_legend(override.aes = list(size=7)))+
  theme(legend.text=element_text(size=12,face = "bold"),
        legend.title = element_text(size=13,face = "bold"))


#Modularity Analysis

#Edge_Betweeness Modularity
Mre.IR2.ebmod=edge.betweenness.community(Mre.IR2,weights=E(Mre.IR2)$weight,directed = F) 
Mre.IR2.ebmod# IGRAPH clustering edge betweenness, groups: 13, mod: 0.63

#
plot.igraph(Mre.IR2.ebmod,Mre.IR2,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(Mre.IR2)$An.color, vertex.frame.color=V(Mre.IR2)$An.color, vertex.shape="circle",
            vertex.size=5, vertex.label=V(Mre.IR2)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(Mre.IR2)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)

#
plot(Mre.IR2.ebmod, Mre.IR2,layout=Lyt.FR[-isolated,], alpha=0.6, 
     vertex.color=V(Mre.IR2)$An.color, vertex.frame.color=V(Mre.IR2)$An.color, vertex.shape="circle", 
     vertex.size=5, vertex.label=V(Mre.IR2)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black", 
     edge.width=E(Mre.IR2)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout", 
     rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)
#Dendogram of the Edge_Betweeness
dendPlot(Mre.IR2.ebmod, mode="hclust")

#Louvain modularity
Mre.IR2.Louvain.mod<-cluster_louvain(Mre.IR2,weights=E(Mre.IR2)$weight)
Mre.IR2.Louvain.mod #IGRAPH clustering multi level, groups: 12, mod: 0.66

#Plot the Louvain communities
plot(Mre.IR2.Louvain.mod, Mre.IR2,layout=Lyt.FR[-isolated,], alpha=0.6, 
     vertex.color=V(Mre.IR2)$An.color, vertex.frame.color=V(Mre.IR2)$An.color, vertex.shape="circle", 
     vertex.size=5, vertex.label=V(Mre.IR2)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black", 
     edge.width=E(Mre.IR2)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout", 
     rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#Generate HeatMaps of the adjacency matrixes
#Let's use the neighbor list for generating the Adjaency Matrix
SensorID_FNN #Let's use the edge list that contain the total amount of possible interactors for the three sensors
v<-levels(as.factor(unlist(SensorID_FNN[,1:2]))) #get vertex names from the network list
n<-length(v)#get the number of vertex in the network
e<-matrix(match(as.character(unlist(SensorID_FNN[,1:2])),v),ncol = 2)#Edgelist with row index instead of Names
w<-Mre_IR2$weight #Edge weights, only from the Mre_IR2 network to quit those edges that do not exist at Mre_IR2

# M is our co-expression network adjacency matrix.

# distringuish between i,j and j,i we want to make M symmetric (undirected) by
# considering an undirected network.
Mre.IR2.M <- matrix(0, n, n)                  # set up a co-expression matrix
Mre.IR2.M[e] <- w                             # fill it in with edge weights
Mre.IR2.M <- Mre.IR2.M + t(Mre.IR2.M)                         # make this symmetric
any(is.na(Mre.IR2.M))#check no NAs in the matrix
dimnames(Mre.IR2.M) <- list(v, v)             # label the vertices
# We let A_ij = 1 if there is an interaction. 
#A <- 1*(M > 0) #right now we keep the weights as they are
Mre.IR2.M #The adjacency Matrix not as an igraph object, easier to manipulate, and to subset for smaller Sensor specific networks
quantile(Mre.IR2.M,c(0,0.25,0.5,0.75,1))
#
#    0%   25%   50%   75%  100% 
#  0.000 0.000 0.000 0.000 0.991  

#Creating the HeatMap
# make a heatmap image with colors: white below 0,
# grey from 0.4 to 1 with darker closer to 1
breaks <- c(0, seq(0.4, 0.991, length=11)) # breaks in the color bins
cols <- grey(1-c(0,seq(0.5,1,length=10))) # colors

#Red Blue Colors
library(RColorBrewer)
breaks <- c(seq(-1, 1, length=11))
breaks
cols <- gray(1-c(0,seq(0,1,length=10))) # colors
cols2 <- hcl.colors(10 , palette = "Blues3" , alpha = 1 , rev = TRUE)
cols3 <- colorRamp("red" , space = "rgb" , interpolate = "linear")
HMCol<- brewer.pal(10,name = "RdBu")
HMCol<- rev(HMCol)
HMCol


# color the vertex for CCNB1 blue, its neighbors red, and the others white
SensorIDBind <- which(v == "Mre11") # index for ccnb1
#SensorIDBind<-as.factor(SensorIDBind)
vcols <- rep("white",n)
vcols[SensorIDBind] <- "blue"
#Subsetting Mre.IR2.M by doing a loop asking for look into every row and column the index in SensorIDBind worked
for (i in SensorIDBind) {
  vcols[which(Mre.IR2.M[,i] > 0 | Mre.IR2.M[i,])] <- "red"
}
##Color annotation for functional annotation
MatSet<-subset(Mat,Mat$NAME%in%v)# select only those proteins that are cointained in the node list of our matrixes
match(v,MatSet$NAME) # Check that the indexes of the Annotation Dataframe and the node list are in the same position 
match(rownames(Mre.IR2.M),MatSet$NAME) #Check that the row names of the Adj Matrix and the Annotation dataframe are in the same order
match(colnames(Mre.IR2.M),MatSet$NAME) #Check that the column names of the Adj Matrix and the Annotation dataframe are in the same order
AnCol<-MatSet$An.color#Vector of colors
AnCol<-as.character(AnCol)#transform it into a character vector, so the colors can be read
#Subsetting Mre.IR2.M using the index vector in SensorIDBind do not work
#vcols[which(Mre.IR2.M[,SensorIDBind] > 0 | Mre.IR2.M[SensorIDBind,])] <- "red"

#Generate a dendogram for keeping the proteins in the sample place
#Mre.IR2.dist<-dist(Mre.IR2.M,method = "euclidian")
#Mre.IR2.Dend<-hclust(Mre.IR2.dist, method = "ward.D" )
#plot(Mre.IR2.Dend)

#Cut the dendogram
#cut.height <- 10 # try numbers from 4 to 7
#FNN.hclusts<- cutree(FNN.Dend, cut.height)
#plot(FNN.hclusts)

# now actually make the heat map
par(mar=rep(0,4)) # make plot margins 0
heatmap(Mre.IR2.M, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=cols, breaks=breaks,frame=T, cexRow = 0.5 , cexCol = 0.5, main = "Mre11 IR2 AdjMat", ylab="Total Protein Interactors",xlab="Total Protein Interactors")
#
plot.new()
legend("topleft", c("Cell Cycle","Chromatin Org","DNA Repair","Nuclear Pore","p53","Proteasome","Ribosome Biogenesis","RNA Metabolism","Transcription Reg","UPR"), 
       fill=c("coral1","darkviolet","firebrick1","deepskyblue3","tomato","grey","lightgreen","gold","slateblue1","maroon"),  
       bty="n", inset=0, xpd=T,  border=F, cex =2,y.intersp = 0.7,x.intersp = 0.2,pt.cex = 1)


#Now we substract the MrenonIR2 AdjMatrix to the MreIR6 AdjMatrix
MreIR2.AdjMat.Change<-Mre.IR2.M-Mre.nonIR5.M
range(MreIR2.AdjMat.Change)#Let's check the magnitude of the substracted weights
#[1]-0.813  0.668
any(is.na(MreIR2.AdjMat.Change))
#Generate the breaks and color functions
breaks <- c(seq(-0.5, 0.5, length=6))
breaks
HMCol<- brewer.pal(5,name = "RdBu")
HMCol<- rev(HMCol)
HMCol
[1] "#2166AC" "#4393C3" "#92C5DE" "#D1E5F0" "#F7F7F7" "#FDDBC7" "#F4A582" "#D6604D" "#B2182B"

#Plot the HeatMap
par(mar=rep(0,4)) # make plot margins 0
heatmap(MreIR2.AdjMat.Change, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=HMCol, breaks=breaks,frame=T, cexRow = 0.5 , cexCol = 0.5, main = "Mre11 IR2 AdjMat (minus nonIR5 edges)", ylab="Total Protein Interactors",xlab="Total Protein Interactors")


#Now we substract the MreIR5 AdjMatrix to the MreIR2 AdjMatrix
MreIR2.AdjMat.Change2<-Mre.IR2.M-Mre.IR30.M
range(MreIR2.AdjMat.Change)#Let's check the magnitude of the substracted weights
#[1]-0.813  0.668
any(is.na(MreIR2.AdjMat.Change))
#Generate the breaks and color functions
breaks <- c(seq(-0.5, 0.5, length=6))
breaks
HMCol<- brewer.pal(5,name = "RdBu")
HMCol<- rev(HMCol)
HMCol
[1] "#2166AC" "#4393C3" "#92C5DE" "#D1E5F0" "#F7F7F7" "#FDDBC7" "#F4A582" "#D6604D" "#B2182B"

#Plot the HeatMap
par(mar=rep(0,4)) # make plot margins 0
heatmap(MreIR2.AdjMat.Change2, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=HMCol, breaks=breaks,frame=T, cexRow = 0.5 , cexCol = 0.5, main = "Mre11 IR2 AdjMat (minus IR30 edges)", ylab="Total Protein Interactors",xlab="Total Protein Interactors")

#
#check with the FNN
# now actually make the heat map
par(mar=rep(0,4)) # make plot margins 0
heatmap(FNN.M, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=0:1, frame=T, cexRow = 0.5 , cexCol = 0.5)
#
legend("topleft", c("Neighbors Mre11","Mre11"), fill=c("red","blue"),  
       bty="n", inset=0, xpd=T,  border=F)


####################### Mre11 IR8h Networks
###################################################

#Mre11 IR8 Networks


#Let's keep the Mat as the Node attribute table and load the Neighbor list as the matrixes
Mat# Node attribute matrix for all the networks

#Mre11ID IR8 nodes Network Edge Lists
Mre_IR8<-read.delim("./Mre11_Files/Mre_IR8h.tsv",sep = "\t",col.names = )
str(Mre_IR8)
head(Mre_IR8)

#Generate the igraph network object
Mre.IR8<-graph_from_data_frame(d=Mre_IR8,vertices = Mat, directed = F )#Generate the igraph object with the edge and node annotation attributes
str(vertex_attr(Mre.IR8))#Check node attributes
str(edge_attr(Mre.IR8))#Check for edge attributes
Mre.IR8<-igraph::simplify(Mre.IR8, remove.multiple=TRUE, remove.loops=TRUE) #Remove repetitive edges and self-connecting edges
Mre.IR8
#Small plot
plot(Mre.IR8)

#Generate an adjacency Matrix that will be used after
Mre.IR8_AdjMat<-as_adjacency_matrix(Mre.IR8,type="both",attr = "weight",names=TRUE,sparse = TRUE)
Mre.IR8_AdjMat<-as.matrix(Mre.IR8_AdjMat)
write.csv(Mre.IR8_AdjMat,"./Mre.IR8_AdjMat.csv")

#Generate a vertex size based on degree and count
vsize<-igraph::degree(Mre.IR8)/vcount(Mre.IR8)
#Identify and remove disconnected nodes
isolated<-which(igraph::degree(Mre.IR8)==0)
isolated
length(isolated)# [1] 546
#remove isolated nodes
Mre.IR8<-igraph::delete_vertices(Mre.IR8,isolated)


#Generate layouts for the graph
#Do not run again, to not to change the layouts that I already have
Lyt.Nice<-layout_nicely(Mre.IR8,dim = 3)
Lyt.GOpt<-layout_with_graphopt(Mre.IR8, niter=5000, mass = 30, charge = E(Mre.IR8)$weight)
Lyt.DrL<-layout_with_drl(Mre.IR8, weights = E(Mre.IR8)$weight)
Lyt.FR<-layout.fruchterman.reingold(Mre.IR8, niter=10000, weights=E(Mre.IR8)$weight)#This one looks like the best option
Lyt.KK<-layout_with_kk(Mre.IR8, maxiter=10000, weights=E(Mre.IR8)$weight) #very large separation of peripheral nodes
Lyt.GEM<-layout_with_gem(Mre.IR8, maxiter=10000) #Bad separation
#Lyt.DH<-layout_with_dh(Mre.IR8, maxiter=150, fineiter=max(10, log2(833)))
Lyt.LGL<-layout.lgl(Mre.IR8, maxit=150,maxdelta=833,area=(length(E(Mre.IR8)))^2,coolexp=1.5,repulserad=((length(E(Mre.IR8)))^2)*(length(E(Mre.IR8))))#Do not work


#Star and Tree layouts

V(Mre.IR8)[name=="MRE11A"]
which(V(Mre.IR8)[name=="Mre11"])# 228 #Check what is the rowindex for your desired node
which(V(Mre.IR8)$name=="MRE11A")# 191 #Check what is the rowindex for your desired node
which(V(Mre.IR8)$name=="MRE11A")# 125
which(V(Mre.IR8)$name=="NCL")# 133
Lyt.Star<-layout_as_star(Mre.IR8,center=V(Mre.IR8)[name=="MRE11A"]) #it has a special order argument where we can assign an order to the vertex, hence maybe the OPLS-DA:VIP value order may be used as ordering vector

Lyt.Tree<-layout_as_tree(Mre.IR8,root=c(125))
Lyt.RT<-layout.reingold.tilford(Mre.IR8,root=c(125))# This one looks good, but the indexing do not allow to highlight the correct nodes
Lyt.sugi<-layout_with_sugiyama(Mre.IR8, maxiter=100, weights=E(Mre.IR82)$weight)#This one looks good for undirected tree layouts

#Check what layout looks fine

plot.igraph(Mre.IR8,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(Mre.IR8)$An.color, vertex.frame.color=V(Mre.IR8)$An.color, vertex.shape="circle",
            vertex.size=5, vertex.label=V(Mre.IR8)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(Mre.IR8)$weight, edge.lty=1, main="Mre11ID IR8 Network", sub="FR Tree Layout",
            rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#With ggraph

ggraph(Mre.IR8,layout=Lyt.RT)+#FR[-isolated,])+
  geom_edge_fan(color="gray65", width=1, alpha=0.4)+
  geom_node_point(aes(color=V(Mre.IR8)$Annotation, size=V(Mre.IR8)$Mre.FC.IR8),alpha=0.75)+
  geom_node_text(aes(label=name),color="gray25",size=4)+
  scale_color_manual(limits=V(Mre.IR8)$Annotation, values=V(Mre.IR8)$An.color)+
  scale_size_continuous(range = c(5, 15)) +
  theme_void()+
  title("Mre11 IR8 Network")+
  guides(color=guide_legend(override.aes = list(size=7)))+
  theme(legend.text=element_text(size=12,face = "bold"),
        legend.title = element_text(size=13,face = "bold"))



#Community Detection
#clp <- igraph::cluster_optimal(Mre.IR8)
class(clp)
V(Mre.IR8)$community <- clp$membership
colrs <- adjustcolor( c("gray50", "tomato", "gold", "yellowgreen"), alpha=.6)


plot.igraph(clp,Mre.IR8,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(Mre.IR8)$An.color, vertex.frame.color=V(Mre.IR8)$An.color, vertex.shape="circle",
            vertex.size==V(Mre.IR8)$Mre.FC.IR8, vertex.label=V(Mre.IR8)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(Mre.IR8)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#Distances edges highlight
dist.from.Mre11 <- distances(Mre.IR8, v=V(Mre.IR8)[name=="MRE11A"],
                             to=V(Mre.IR8), weights=E(Mre.IR8)$weight)
dist.from.Mre11<-dist.from.Mre11[,!is.infinite(colSums(dist.from.Mre11))]

dist.from.Mre11
range(dist.from.Mre11)

#Set colors to plot the distances:
oranges <- colorRampPalette(c("dark red", "gold"))
col <- oranges(max(dist.from.Mre11))
col <- col[dist.from.Mre11+1]
plot(Mre.IR8, layout=Lyt.FR[-isolated,],vertex.color=col, vertex.label=dist.from.Mre11,
     vertex.label.color="white",vertex.shape="circle",
     vertex.size=5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(Mre.IR8)$weight, edge.lty=1, main="Mre11 IR8 Full Nuclear Network", sub="FR Layout",
     rescale=T)


#Shortest path from Mre11 to Nucleolin
MreNcl.path <- shortest_paths(Mre.IR8,
                              from = V(Mre.IR8)[name=="MRE11A"],
                              to = V(Mre.IR8)[name=="NCL"],
                              output = "both") # both path nodes and edges
MreNcl.path
# Generate edge color variable to plot the path:
ecol <- rep("gray80", ecount(Mre.IR8))
ecol[unlist(MreNcl.path$epath)] <- "orange"
# Generate edge width variable to plot the path:
ew <- rep(2, ecount(Mre.IR8))
ew[unlist(MreNcl.path$epath)] <- 4
# Generate node color variable to plot the path:
vcol <- rep("gray40", vcount(Mre.IR8))
vcol[unlist(MreNcl.path$vpath)] <- "gold"
plot(Mre.IR8, layout=Lyt.FR[-isolated,], vertex.color=vcol, edge.color=ecol,
     edge.width=ew,vertex.label.color="white",vertex.shape="circle",
     vertex.size=5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(Mre.IR8)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)



#Plot based on Mre11 as incident
inc.edges <- incident(Mre.IR8, V(Mre.IR8)[name=="MRE11A"], mode="all")
# Set colors to plot the selected edges.
ecol <- rep("gray80", ecount(Mre.IR8))
ecol[inc.edges] <- "orange"
vcol <- rep("gray40", vcount(Mre.IR8))
vcol[V(Mre.IR8)$name=="Mre11"] <- "gold"

plot(Mre.IR8, layout=Lyt.FR[-isolated,], vertex.color=vcol, edge.color=ecol, 
     edge.width=ew,vertex.label.color="white",vertex.shape="circle",
     vertex.size=5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(Mre.IR8)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)


#Plot Degree Distribution
Mre.IR8.Deg<-degree_distribution(Mre.IR8, cumulative = T, mode="all")
plot( x=0:max(igraph::degree(Mre.IR8)), y=1-Mre.IR8.Deg, pch=19, cex=1.2, col="orange",
      xlab="Degree", ylab="Cumulative Frequency")

#Network Descriptives
#Number of nodes and edges
vcount(Mre.IR8) #[1] 287
ecount(Mre.IR8)#[1]  2177
#Density (Number of edges in relation to the total number of edges in a full connected graph)
edge_density(Mre.IR8,loops = F) #[1]  0.05304452
#Transitivity (Total number of conected triangles)
transitivity(Mre.IR8,type="global") #[1] 0.6119276
#modularity
mod<-V(Mre.IR8)$Mre.FuzzyCMeans #lets check for fuzzy Cmeans modularity
mod<-mod+1#The non clustered interactors that appears as 0, needs to be transformed to 1, so we add 1 to every cluster
modularity(Mre.IR8,membership = mod)
#S6.FuzzyCMeans modularity [1] -0.003455341
#Ku.FuzzyCMeans modularity [1] 0.00309126
#Mre.FuzzyCMeans modularity [1] -0.003669612

#Assorativity
assortativity_nominal(Mre.IR8,factor(V(Mre.IR8)$Mre.FuzzyCMeans),directed=F)
#Annotation [1] [1] 0.2375432
#S6.FuzzyCMeans [1] -0.005140002
#Ku.FuzzyCMeans [1]  0.005481397
#Mre.FuzzyCMeans [1] -0.00488008

#Diameter (longest geodesic distance), shortest path between two most peripheral nodes
diameter(Mre.IR8,directed = F, weights = E(Mre.IR8)$weights)#[1] 3.896
diam<-get_diameter(Mre.IR8,directed = F,weights = E(Mre.IR8)$weights)
diam #+ 8/287 vertices, named, from 9afb3d6: [1] PELP1  WDR18  GNL2   IFI16  NXF1   NUP88  DDX19A SEPHS1
#Plot Diameter
vcol <- rep("gray50", vcount(Mre.IR8))
vcol[diam] <- "gold"
ecol <- rep("gray90", ecount(Mre.IR8))
ecol[E(Mre.IR8, path=diam)] <- "orange"
# E(net, path=diam) finds edges along a path, here 'diam'
plot(Mre.IR8, vertex.color=vcol, edge.color=ecol, vertex.label.color="black",edge.arrow.mode=0,layout=Lyt.FR[-isolated,],vertex.size=5)

#Node topological descriptives
degree<-igraph::centr_degree(Mre.IR8,mode = "all",normalized = T)
degree<-degree$res
eigen_cent<-igraph::eigen_centrality(Mre.IR8)
eigen_cent<-eigen_cent$vector
clos_cent <- igraph::closeness(Mre.IR8, normalized=T)
betw_cent <- igraph::betweenness(Mre.IR8, normalized=T)
hs <- hub_score(Mre.IR8, weights=E(Mre.IR8)$weight)$vector
as <- authority_score(Mre.IR8, weights=E(Mre.IR8)$weight)$vector
Trans<-transitivity(Mre.IR8,type="local")
Trans
#write results
Mre.IR8.Centrality<-data.frame(NAME=V(Mre.IR8)$name, UniprotID=V(Mre.IR8)$Uniprot, Degree=degree, Eigen= eigen_cent, Closeness= clos_cent, Betweeness= betw_cent, Hubs=hs, Authority=as,Transitivity=Trans)
Mre.IR8.Centrality
write.csv(Mre.IR8.Centrality,"./Mre.IR8.Centrality.csv")

#Cocitation, Number of shared interactions among nodes
Cocitation<-cocitation(Mre.IR8)
Cocitation#This matrix is similar to an adjacency matrix
write.csv(Cocitation,"./Mre.IR8.cocitation.csv")
#Plot with degree as node size
ggraph(Mre.IR8,layout=Lyt.FR[-isolated,])+
  geom_edge_fan(color="gray85", width=1, alpha=0.4)+
  geom_node_point(aes(color=V(Mre.IR8)$Annotation, size=degree),alpha=0.75)+
  geom_node_text(aes(label=name),color="gray25",size=4)+
  scale_color_manual(limits=V(Mre.IR8)$Annotation, values=V(Mre.IR8)$An.color)+
  scale_size_continuous(range = c(5, 15)) +
  theme_void()+
  title("Mre11 IR8 Network")+
  guides(color=guide_legend(override.aes = list(size=7)))+
  theme(legend.text=element_text(size=12,face = "bold"),
        legend.title = element_text(size=13,face = "bold"))


#Modularity Analysis

#Edge_Betweeness Modularity
Mre.IR8.ebmod=edge.betweenness.community(Mre.IR8,weights=E(Mre.IR8)$weight,directed = F) 
Mre.IR8.ebmod# IGRAPH clustering edge betweenness, groups: 61, mod: 0.42

#
plot.igraph(Mre.IR8.ebmod,Mre.IR8,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(Mre.IR8)$An.color, vertex.frame.color=V(Mre.IR8)$An.color, vertex.shape="circle",
            vertex.size=5, vertex.label=V(Mre.IR8)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(Mre.IR8)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)

#
plot(Mre.IR8.ebmod, Mre.IR8,layout=Lyt.FR[-isolated,], alpha=0.6, 
     vertex.color=V(Mre.IR8)$An.color, vertex.frame.color=V(Mre.IR8)$An.color, vertex.shape="circle", 
     vertex.size=5, vertex.label=V(Mre.IR8)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black", 
     edge.width=E(Mre.IR8)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout", 
     rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)
#Dendogram of the Edge_Betweeness
dendPlot(Mre.IR8.ebmod, mode="hclust")

#Louvain modularity
Mre.IR8.Louvain.mod<-cluster_louvain(Mre.IR8,weights=E(Mre.IR8)$weight)
Mre.IR8.Louvain.mod #IGRAPH clustering multi level, groups: 8, mod: 0.52

#Plot the Louvain communities
plot(Mre.IR8.Louvain.mod, Mre.IR8,layout=Lyt.FR[-isolated,], alpha=0.6, 
     vertex.color=V(Mre.IR8)$An.color, vertex.frame.color=V(Mre.IR8)$An.color, vertex.shape="circle", 
     vertex.size=5, vertex.label=V(Mre.IR8)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black", 
     edge.width=E(Mre.IR8)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout", 
     rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#Generate HeatMaps of the adjacency matrixes
#Let's use the neighbor list for generating the Adjaency Matrix
SensorID_FNN #Let's use the edge list that contain the total amount of possible interactors for the three sensors
v<-levels(as.factor(unlist(SensorID_FNN[,1:2]))) #get vertex names from the network list
n<-length(v)#get the number of vertex in the network
e<-matrix(match(as.character(unlist(SensorID_FNN[,1:2])),v),ncol = 2)#Edgelist with row index instead of Names
w<-Mre_IR8$weight #Edge weights, only from the Mre_IR8 network to quit those edges that do not exist at Mre_IR8

# M is our co-expression network adjacency matrix.

# distringuish between i,j and j,i we want to make M symmetric (undirected) by
# considering an undirected network.
Mre.IR8.M <- matrix(0, n, n)                  # set up a co-expression matrix
Mre.IR8.M[e] <- w                             # fill it in with edge weights
Mre.IR8.M <- Mre.IR8.M + t(Mre.IR8.M)                         # make this symmetric
any(is.na(Mre.IR8.M))#check no NAs in the matrix
dimnames(Mre.IR8.M) <- list(v, v)             # label the vertices
# We let A_ij = 1 if there is an interaction. 
#A <- 1*(M > 0) #right now we keep the weights as they are
Mre.IR8.M #The adjacency Matrix not as an igraph object, easier to manipulate, and to subset for smaller Sensor specific networks
quantile(Mre.IR8.M,c(0,0.25,0.5,0.75,1))
#
#       0%   25%   50%   75%  100% 
#     0.000 0.000 0.000 0.000 0.991 

#Creating the HeatMap
# make a heatmap image with colors: white below 0,
# grey from 0.4 to 1 with darker closer to 1
breaks <- c(0, seq(0.4, 0.991, length=11)) # breaks in the color bins
cols <- grey(1-c(0,seq(0.5,1,length=10))) # colors

#Red Blue Colors
library(RColorBrewer)
breaks <- c(seq(-1, 1, length=11))
breaks
cols <- gray(1-c(0,seq(0,1,length=10))) # colors
cols2 <- hcl.colors(10 , palette = "Blues3" , alpha = 1 , rev = TRUE)
cols3 <- colorRamp("red" , space = "rgb" , interpolate = "linear")
HMCol<- brewer.pal(10,name = "RdBu")
HMCol<- rev(HMCol)
HMCol


# color the vertex for CCNB1 blue, its neighbors red, and the others white
SensorIDBind <- which(v == "Mre11") # index for ccnb1
#SensorIDBind<-as.factor(SensorIDBind)
vcols <- rep("white",n)
vcols[SensorIDBind] <- "blue"
#Subsetting Mre.IR8.M by doing a loop asking for look into every row and column the index in SensorIDBind worked
for (i in SensorIDBind) {
  vcols[which(Mre.IR8.M[,i] > 0 | Mre.IR8.M[i,])] <- "red"
}
##Color annotation for functional annotation
MatSet<-subset(Mat,Mat$NAME%in%v)# select only those proteins that are cointained in the node list of our matrixes
match(v,MatSet$NAME) # Check that the indexes of the Annotation Dataframe and the node list are in the same position 
match(rownames(Mre.IR8.M),MatSet$NAME) #Check that the row names of the Adj Matrix and the Annotation dataframe are in the same order
match(colnames(Mre.IR8.M),MatSet$NAME) #Check that the column names of the Adj Matrix and the Annotation dataframe are in the same order
AnCol<-MatSet$An.color#Vector of colors
AnCol<-as.character(AnCol)#transform it into a character vector, so the colors can be read

#Subsetting Mre.IR8.M using the index vector in SensorIDBind do not work
#vcols[which(Mre.IR8.M[,SensorIDBind] > 0 | Mre.IR8.M[SensorIDBind,])] <- "red"

#Generate a dendogram for keeping the proteins in the sample place
#Mre.IR8.dist<-dist(Mre.IR8.M,method = "euclidian")
#Mre.IR8.Dend<-hclust(Mre.IR8.dist, method = "ward.D" )
#plot(Mre.IR8.Dend)

#Cut the dendogram
#cut.height <- 10 # try numbers from 4 to 7
#FNN.hclusts<- cutree(FNN.Dend, cut.height)
#plot(FNN.hclusts)

# now actually make the heat map
par(mar=rep(0,4)) # make plot margins 0
heatmap(Mre.IR8.M, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=cols, breaks=breaks,frame=T, cexRow = 0.5 , cexCol = 0.5, main = "Mre11 IR8 AdjMat", ylab="Total Protein Interactors",xlab="Total Protein Interactors")
#
plot.new()
legend("topleft", c("Cell Cycle","Chromatin Org","DNA Repair","Nuclear Pore","p53","Proteasome","Ribosome Biogenesis","RNA Metabolism","Transcription Reg","UPR"), 
       fill=c("coral1","darkviolet","firebrick1","deepskyblue3","tomato","grey","lightgreen","gold","slateblue1","maroon"),  
       bty="n", inset=0, xpd=T,  border=F, cex =2,y.intersp = 0.7,x.intersp = 0.2,pt.cex = 1)


#Now we substract the MrenonIR8 AdjMatrix to the MreIR6 AdjMatrix
MreIR8.AdjMat.Change<-Mre.IR8.M-Mre.nonIR5.M
range(MreIR8.AdjMat.Change)#Let's check the magnitude of the substracted weights
#[1]-0.813  0.668
any(is.na(MreIR8.AdjMat.Change))
#Generate the breaks and color functions
breaks <- c(seq(-0.5, 0.5, length=6))
breaks
HMCol<- brewer.pal(5,name = "RdBu")
HMCol<- rev(HMCol)
HMCol
[1] "#2166AC" "#4393C3" "#92C5DE" "#D1E5F0" "#F7F7F7" "#FDDBC7" "#F4A582" "#D6604D" "#B2182B"

#Plot the HeatMap
par(mar=rep(0,4)) # make plot margins 0
heatmap(MreIR8.AdjMat.Change, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=HMCol, breaks=breaks,frame=T, cexRow = 0.5 , cexCol = 0.5, main = "Mre11 IR8 AdjMat (minus nonIR5 edges)", ylab="Total Protein Interactors",xlab="Total Protein Interactors")


#Now we substract the MreIR5 AdjMatrix to the MreIR8 AdjMatrix
MreIR8.AdjMat.Change2<-Mre.IR8.M-Mre.IR2.M
range(MreIR8.AdjMat.Change)#Let's check the magnitude of the substracted weights
#[1]-0.813  0.668
any(is.na(MreIR8.AdjMat.Change))
#Generate the breaks and color functions
breaks <- c(seq(-0.5, 0.5, length=10))
breaks
HMCol<- brewer.pal(9,name = "RdBu")
HMCol<- rev(HMCol)
HMCol
[1] "#2166AC" "#4393C3" "#92C5DE" "#D1E5F0" "#F7F7F7" "#FDDBC7" "#F4A582" "#D6604D" "#B2182B"

#Plot the HeatMap
par(mar=rep(0,4)) # make plot margins 0
heatmap(MreIR8.AdjMat.Change2, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=HMCol, breaks=breaks,frame=T, cexRow = 0.5 , cexCol = 0.5, main = "Mre11 IR8 AdjMat (minus IR2 edges)", ylab="Total Protein Interactors",xlab="Total Protein Interactors")
#
#check with the FNN
# now actually make the heat map
par(mar=rep(0,4)) # make plot margins 0
heatmap(FNN.M, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=0:1, frame=T, cexRow = 0.5 , cexCol = 0.5)
#
legend("topleft", c("Neighbors Mre11","Mre11"), fill=c("red","blue"),  
       bty="n", inset=0, xpd=T,  border=F)


####################### Mre11 IR24h Networks
###################################################

#Mre11 IR24 Networks


#Let's keep the Mat as the Node attribute table and load the Neighbor list as the matrixes
Mat# Node attribute matrix for all the networks

#Mre11ID IR24 nodes Network Edge Lists
Mre_IR24<-read.delim("./Mre11_Files/Mre_IR24h.tsv",sep = "\t",col.names = )
str(Mre_IR24)
head(Mre_IR24)

#Generate the igraph network object
Mre.IR24<-graph_from_data_frame(d=Mre_IR24,vertices = Mat, directed = F )#Generate the igraph object with the edge and node annotation attributes
str(vertex_attr(Mre.IR24))#Check node attributes
str(edge_attr(Mre.IR24))#Check for edge attributes
Mre.IR24<-igraph::simplify(Mre.IR24, remove.multiple=TRUE, remove.loops=TRUE) #Remove repetitive edges and self-connecting edges
Mre.IR24
#Small plot
plot(Mre.IR24)

#Generate an adjacency Matrix that will be used after
Mre.IR24_AdjMat<-as_adjacency_matrix(Mre.IR24,type="both",attr = "weight",names=TRUE,sparse = TRUE)
Mre.IR24_AdjMat<-as.matrix(Mre.IR24_AdjMat)
write.csv(Mre.IR24_AdjMat,"./Mre.IR24_AdjMat.csv")

#Generate a vertex size based on degree and count
vsize<-igraph::degree(Mre.IR24)/vcount(Mre.IR24)
#Identify and remove disconnected nodes
isolated<-which(igraph::degree(Mre.IR24)==0)
isolated
length(isolated)# [1] 650
#remove isolated nodes
Mre.IR24<-igraph::delete_vertices(Mre.IR24,isolated)


#Generate layouts for the graph
#Do not run again, to not to change the layouts that I already have
Lyt.Nice<-layout_nicely(Mre.IR24,dim = 3)
Lyt.GOpt<-layout_with_graphopt(Mre.IR24, niter=5000, mass = 30, charge = E(Mre.IR24)$weight)
Lyt.DrL<-layout_with_drl(Mre.IR24, weights = E(Mre.IR24)$weight)
Lyt.FR<-layout.fruchterman.reingold(Mre.IR24, niter=10000, weights=E(Mre.IR24)$weight)#This one looks like the best option
Lyt.KK<-layout_with_kk(Mre.IR24, maxiter=10000, weights=E(Mre.IR24)$weight) #very large separation of peripheral nodes
Lyt.GEM<-layout_with_gem(Mre.IR24, maxiter=10000) #Bad separation
#Lyt.DH<-layout_with_dh(Mre.IR24, maxiter=150, fineiter=max(10, log2(833)))
Lyt.LGL<-layout.lgl(Mre.IR24, maxit=150,maxdelta=833,area=(length(E(Mre.IR24)))^2,coolexp=1.5,repulserad=((length(E(Mre.IR24)))^2)*(length(E(Mre.IR24))))#Do not work


#Star and Tree layouts

V(Mre.IR24)[name=="MRE11A"]
which(V(Mre.IR24)[name=="Mre11"])# 228 #Check what is the rowindex for your desired node
which(V(Mre.IR24)$name=="MRE11A")# 73 #Check what is the rowindex for your desired node
which(V(Mre.IR24)$name=="MRE11A")# 73
which(V(Mre.IR24)$name=="NCL")# 77
Lyt.Star<-layout_as_star(Mre.IR24,center=V(Mre.IR24)[name=="MRE11A"]) #it has a special order argument where we can assign an order to the vertex, hence maybe the OPLS-DA:VIP value order may be used as ordering vector

Lyt.Tree<-layout_as_tree(Mre.IR24,root=c(73))
Lyt.RT<-layout.reingold.tilford(Mre.IR24,root=c(73))# This one looks good, but the indexing do not allow to highlight the correct nodes
Lyt.sugi<-layout_with_sugiyama(Mre.IR24, maxiter=100, weights=E(Mre.IR242)$weight)#This one looks good for undirected tree layouts

#Check what layout looks fine

plot.igraph(Mre.IR24,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(Mre.IR24)$An.color, vertex.frame.color=V(Mre.IR24)$An.color, vertex.shape="circle",
            vertex.size=5, vertex.label=V(Mre.IR24)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(Mre.IR24)$weight, edge.lty=1, main="Mre11ID IR24 Network", sub="FR Tree Layout",
            rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#With ggraph

ggraph(Mre.IR24,layout=Lyt.FR[-isolated,])+
  geom_edge_fan(color="gray65", width=1, alpha=0.4)+
  geom_node_point(aes(color=V(Mre.IR24)$Annotation, size=V(Mre.IR24)$Mre.FC.IR24),alpha=0.75)+
  geom_node_text(aes(label=name),color="gray25",size=4)+
  scale_color_manual(limits=V(Mre.IR24)$Annotation, values=V(Mre.IR24)$An.color)+
  scale_size_continuous(range = c(5, 15)) +
  theme_void()+
  title("Mre11 IR24 Network")+
  guides(color=guide_legend(override.aes = list(size=7)))+
  theme(legend.text=element_text(size=12,face = "bold"),
        legend.title = element_text(size=13,face = "bold"))



#Community Detection
#clp <- igraph::cluster_optimal(Mre.IR24)
class(clp)
V(Mre.IR24)$community <- clp$membership
colrs <- adjustcolor( c("gray50", "tomato", "gold", "yellowgreen"), alpha=.6)


plot.igraph(clp,Mre.IR24,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(Mre.IR24)$An.color, vertex.frame.color=V(Mre.IR24)$An.color, vertex.shape="circle",
            vertex.size==V(Mre.IR24)$Mre.FC.IR24, vertex.label=V(Mre.IR24)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(Mre.IR24)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#Distances edges highlight
dist.from.Mre11 <- distances(Mre.IR24, v=V(Mre.IR24)[name=="MRE11A"],
                             to=V(Mre.IR24), weights=E(Mre.IR24)$weight)
dist.from.Mre11<-dist.from.Mre11[,!is.infinite(colSums(dist.from.Mre11))]

dist.from.Mre11
range(dist.from.Mre11)

#Set colors to plot the distances:
oranges <- colorRampPalette(c("dark red", "gold"))
col <- oranges(max(dist.from.Mre11))
col <- col[dist.from.Mre11+1]
plot(Mre.IR24, layout=Lyt.FR[-isolated,],vertex.color=col, vertex.label=dist.from.Mre11,
     vertex.label.color="white",vertex.shape="circle",
     vertex.size=5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(Mre.IR24)$weight, edge.lty=1, main="Mre11 IR24 Full Nuclear Network", sub="FR Layout",
     rescale=T)


#Shortest path from Mre11 to Nucleolin
MreNcl.path <- shortest_paths(Mre.IR24,
                              from = V(Mre.IR24)[name=="MRE11A"],
                              to = V(Mre.IR24)[name=="NCL"],
                              output = "both") # both path nodes and edges
MreNcl.path
# Generate edge color variable to plot the path:
ecol <- rep("gray80", ecount(Mre.IR24))
ecol[unlist(MreNcl.path$epath)] <- "orange"
# Generate edge width variable to plot the path:
ew <- rep(2, ecount(Mre.IR24))
ew[unlist(MreNcl.path$epath)] <- 4
# Generate node color variable to plot the path:
vcol <- rep("gray40", vcount(Mre.IR24))
vcol[unlist(MreNcl.path$vpath)] <- "gold"
plot(Mre.IR24, layout=Lyt.FR[-isolated,], vertex.color=vcol, edge.color=ecol,
     edge.width=ew,vertex.label.color="white",vertex.shape="circle",
     vertex.size=5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(Mre.IR24)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)



#Plot based on Mre11 as incident
inc.edges <- incident(Mre.IR24, V(Mre.IR24)[name=="MRE11A"], mode="all")
# Set colors to plot the selected edges.
ecol <- rep("gray80", ecount(Mre.IR24))
ecol[inc.edges] <- "orange"
vcol <- rep("gray40", vcount(Mre.IR24))
vcol[V(Mre.IR24)$name=="Mre11"] <- "gold"

plot(Mre.IR24, layout=Lyt.FR[-isolated,], vertex.color=vcol, edge.color=ecol, 
     edge.width=ew,vertex.label.color="white",vertex.shape="circle",
     vertex.size=V(Mre.IR24)$Mre.FC.IR24, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(Mre.IR24)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)


#Plot Degree Distribution
Mre.IR24.Deg<-degree_distribution(Mre.IR24, cumulative = T, mode="all")
plot( x=0:max(igraph::degree(Mre.IR24)), y=1-Mre.IR24.Deg, pch=19, cex=1.2, col="orange",
      xlab="Degree", ylab="Cumulative Frequency")

#Network Descriptives
#Number of nodes and edges
vcount(Mre.IR24) #[1] 183
ecount(Mre.IR24)#[1]  1081
#Density (Number of edges in relation to the total number of edges in a full connected graph)
edge_density(Mre.IR24,loops = F) #[1] 0.06491323
#Transitivity (Total number of conected triangles)
transitivity(Mre.IR24,type="global") #[1] 0.5857033
#modularity
mod<-V(Mre.IR24)$Mre.FuzzyCMeans #lets check for fuzzy Cmeans modularity
mod<-mod+1#The non clustered interactors that appears as 0, needs to be transformed to 1, so we add 1 to every cluster
modularity(Mre.IR24,membership = mod)
#S6.FuzzyCMeans modularity [1]-0.01201007
#Ku.FuzzyCMeans modularity [1] 0.01677876
#Mre.FuzzyCMeans modularity [1] -0.001419695

#Assorativity
assortativity_nominal(Mre.IR24,factor(V(Mre.IR24)$Mre.FuzzyCMeans),directed=F)
#Annotation [1] [1] 0.2580939
#S6.FuzzyCMeans [1] -0.0183629  
#Ku.FuzzyCMeans [1] 0.02987433
#Mre.FuzzyCMeans [1] -0.002330708

#Diameter (longest geodesic distance), shortest path between two most peripheral nodes
diameter(Mre.IR24,directed = F, weights = E(Mre.IR24)$weights)#[1] 4.588
diam<-get_diameter(Mre.IR24,directed = F,weights = E(Mre.IR24)$weights)
diam #+ 9/183 vertices, named, from eb248ab: [1] EXOSC8 MTREX  RPL5   EFTUD2 RPL13  CDC5L  SFPQ   CPSF7  NUDT21
#Plot Diameter
vcol <- rep("gray50", vcount(Mre.IR24))
vcol[diam] <- "gold"
ecol <- rep("gray90", ecount(Mre.IR24))
ecol[E(Mre.IR24, path=diam)] <- "orange"
# E(net, path=diam) finds edges along a path, here 'diam'
plot(Mre.IR24, vertex.color=vcol, edge.color=ecol, vertex.label.color="black",edge.arrow.mode=0,layout=Lyt.FR[-isolated,],vertex.size=5)

#Node topological descriptives
degree<-igraph::centr_degree(Mre.IR24,mode = "all",normalized = T)
degree<-degree$res
eigen_cent<-igraph::eigen_centrality(Mre.IR24)
eigen_cent<-eigen_cent$vector
clos_cent <- igraph::closeness(Mre.IR24, normalized=T)
betw_cent <- igraph::betweenness(Mre.IR24, normalized=T)
hs <- hub_score(Mre.IR24, weights=E(Mre.IR24)$weight)$vector
as <- authority_score(Mre.IR24, weights=E(Mre.IR24)$weight)$vector
Trans<-transitivity(Mre.IR24,type="local")
Trans
#write results
Mre.IR24.Centrality<-data.frame(NAME=V(Mre.IR24)$name, UniprotID=V(Mre.IR24)$Uniprot, Degree=degree, Eigen= eigen_cent, Closeness= clos_cent, Betweeness= betw_cent, Hubs=hs, Authority=as,Transitivity=Trans)
Mre.IR24.Centrality
write.csv(Mre.IR24.Centrality,"./Mre.IR24.Centrality.csv")

#Cocitation, Number of shared interactions among nodes
Cocitation<-cocitation(Mre.IR24)
Cocitation#This matrix is similar to an adjacency matrix
write.csv(Cocitation,"./Mre.IR24.cocitation.csv")
#Plot with degree as node size
ggraph(Mre.IR24,layout=Lyt.FR[-isolated,])+
  geom_edge_fan(color="gray85", width=1, alpha=0.4)+
  geom_node_point(aes(color=V(Mre.IR24)$Annotation, size=degree),alpha=0.75)+
  geom_node_text(aes(label=name),color="gray25",size=4)+
  scale_color_manual(limits=V(Mre.IR24)$Annotation, values=V(Mre.IR24)$An.color)+
  scale_size_continuous(range = c(5, 15)) +
  theme_void()+
  title("Mre11 IR24 Network")+
  guides(color=guide_legend(override.aes = list(size=7)))+
  theme(legend.text=element_text(size=12,face = "bold"),
        legend.title = element_text(size=13,face = "bold"))


#Modularity Analysis

#Edge_Betweeness Modularity
Mre.IR24.ebmod=edge.betweenness.community(Mre.IR24,weights=E(Mre.IR24)$weight,directed = F) 
Mre.IR24.ebmod# IGRAPH clustering edge betweenness, groups: 33, mod: 0.56

#
plot.igraph(Mre.IR24.ebmod,Mre.IR24,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(Mre.IR24)$An.color, vertex.frame.color=V(Mre.IR24)$An.color, vertex.shape="circle",
            vertex.size=5, vertex.label=V(Mre.IR24)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(Mre.IR24)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)

#
plot(Mre.IR24.ebmod, Mre.IR24,layout=Lyt.FR[-isolated,], alpha=0.6, 
     vertex.color=V(Mre.IR24)$An.color, vertex.frame.color=V(Mre.IR24)$An.color, vertex.shape="circle", 
     vertex.size=5, vertex.label=V(Mre.IR24)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black", 
     edge.width=E(Mre.IR24)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout", 
     rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)
#Dendogram of the Edge_Betweeness
dendPlot(Mre.IR24.ebmod, mode="hclust")

#Louvain modularity
Mre.IR24.Louvain.mod<-cluster_louvain(Mre.IR24,weights=E(Mre.IR24)$weight)
Mre.IR24.Louvain.mod #IGRAPH clustering multi level, groups: 9, mod: 0.59

#Plot the Louvain communities
plot(Mre.IR24.Louvain.mod, Mre.IR24,layout=Lyt.FR[-isolated,], alpha=0.6, 
     vertex.color=V(Mre.IR24)$An.color, vertex.frame.color=V(Mre.IR24)$An.color, vertex.shape="circle", 
     vertex.size=5, vertex.label=V(Mre.IR24)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black", 
     edge.width=E(Mre.IR24)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout", 
     rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#Generate HeatMaps of the adjacency matrixes
#Let's use the neighbor list for generating the Adjaency Matrix
SensorID_FNN #Let's use the edge list that contain the total amount of possible interactors for the three sensors
v<-levels(as.factor(unlist(SensorID_FNN[,1:2]))) #get vertex names from the network list
n<-length(v)#get the number of vertex in the network
e<-matrix(match(as.character(unlist(SensorID_FNN[,1:2])),v),ncol = 2)#Edgelist with row index instead of Names
w<-Mre_IR24$weight #Edge weights, only from the Mre_IR24 network to quit those edges that do not exist at Mre_IR24

# M is our co-expression network adjacency matrix.

# distringuish between i,j and j,i we want to make M symmetric (undirected) by
# considering an undirected network.
Mre.IR24.M <- matrix(0, n, n)                  # set up a co-expression matrix
Mre.IR24.M[e] <- w                             # fill it in with edge weights
Mre.IR24.M <- Mre.IR24.M + t(Mre.IR24.M)                         # make this symmetric
any(is.na(Mre.IR24.M))#check no NAs in the matrix
dimnames(Mre.IR24.M) <- list(v, v)             # label the vertices
# We let A_ij = 1 if there is an interaction. 
#A <- 1*(M > 0) #right now we keep the weights as they are
Mre.IR24.M #The adjacency Matrix not as an igraph object, easier to manipulate, and to subset for smaller Sensor specific networks
quantile(Mre.IR24.M,c(0,0.25,0.5,0.75,1))
#
#      0%  25%  50%  75% 100% 
#    0.00 0.00 0.00 0.00 0.99   

#Creating the HeatMap
# make a heatmap image with colors: white below 0,
# grey from 0.4 to 1 with darker closer to 1
breaks <- c(0, seq(0.4, 0.991, length=11)) # breaks in the color bins
cols <- grey(1-c(0,seq(0.5,1,length=10))) # colors

#Red Blue Colors
library(RColorBrewer)
breaks <- c(seq(-1, 1, length=11))
breaks
cols <- gray(1-c(0,seq(0,1,length=10))) # colors
cols2 <- hcl.colors(10 , palette = "Blues3" , alpha = 1 , rev = TRUE)
cols3 <- colorRamp("red" , space = "rgb" , interpolate = "linear")
HMCol<- brewer.pal(10,name = "RdBu")
HMCol<- rev(HMCol)
HMCol


# color the vertex for CCNB1 blue, its neighbors red, and the others white
SensorIDBind <- which(v == "Mre11") # index for ccnb1
#SensorIDBind<-as.factor(SensorIDBind)
vcols <- rep("white",n)
vcols[SensorIDBind] <- "blue"
#Subsetting Mre.IR24.M by doing a loop asking for look into every row and column the index in SensorIDBind worked
for (i in SensorIDBind) {
  vcols[which(Mre.IR24.M[,i] > 0 | Mre.IR24.M[i,])] <- "red"
}
##Color annotation for functional annotation
MatSet<-subset(Mat,Mat$NAME%in%v)# select only those proteins that are cointained in the node list of our matrixes
match(v,MatSet$NAME) # Check that the indexes of the Annotation Dataframe and the node list are in the same position 
match(rownames(Mre.IR24.M),MatSet$NAME) #Check that the row names of the Adj Matrix and the Annotation dataframe are in the same order
match(colnames(Mre.IR24.M),MatSet$NAME) #Check that the column names of the Adj Matrix and the Annotation dataframe are in the same order
AnCol<-MatSet$An.color#Vector of colors
AnCol<-as.character(AnCol)#transform it into a character vector, so the colors can be read

#Subsetting Mre.IR24.M using the index vector in SensorIDBind do not work
#vcols[which(Mre.IR24.M[,SensorIDBind] > 0 | Mre.IR24.M[SensorIDBind,])] <- "red"

#Generate a dendogram for keeping the proteins in the sample place
#Mre.IR24.dist<-dist(Mre.IR24.M,method = "euclidian")
#Mre.IR24.Dend<-hclust(Mre.IR24.dist, method = "ward.D" )
#plot(Mre.IR24.Dend)

#Cut the dendogram
#cut.height <- 10 # try numbers from 4 to 7
#FNN.hclusts<- cutree(FNN.Dend, cut.height)
#plot(FNN.hclusts)

# now actually make the heat map
par(mar=rep(0,4)) # make plot margins 0
heatmap(Mre.IR24.M, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=cols, breaks=breaks,frame=T, cexRow = 0.5 , cexCol = 0.5, main = "Mre11 IR24 AdjMat", ylab="Total Protein Interactors",xlab="Total Protein Interactors")
#
plot.new()
legend("topleft", c("Cell Cycle","Chromatin Org","DNA Repair","Nuclear Pore","p53","Proteasome","Ribosome Biogenesis","RNA Metabolism","Transcription Reg","UPR"), 
       fill=c("coral1","darkviolet","firebrick1","deepskyblue3","tomato","grey","lightgreen","gold","slateblue1","maroon"),  
       bty="n", inset=0, xpd=T,  border=F, cex =2,y.intersp = 0.7,x.intersp = 0.2,pt.cex = 1)


#Now we substract the MrenonIR24 AdjMatrix to the MreIR6 AdjMatrix
MreIR24.AdjMat.Change<-Mre.IR24.M-Mre.nonIR5.M
range(MreIR24.AdjMat.Change)#Let's check the magnitude of the substracted weights
#[1] -0.806  0.805[1] -0.806  0.805
any(is.na(MreIR24.AdjMat.Change))
#Generate the breaks and color functions
breaks <- c(seq(-0.5, 0.5, length=6))
breaks
HMCol<- brewer.pal(5,name = "RdBu")
HMCol<- rev(HMCol)
HMCol
[1] "#2166AC" "#4393C3" "#92C5DE" "#D1E5F0" "#F7F7F7" "#FDDBC7" "#F4A582" "#D6604D" "#B2182B"

#Plot the HeatMap
par(mar=rep(0,4)) # make plot margins 0
heatmap(MreIR24.AdjMat.Change, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=HMCol, breaks=breaks,frame=T, cexRow = 0.5 , cexCol = 0.5, main = "Mre11 IR24 AdjMat (minus nonIR24 edges)", ylab="Total Protein Interactors",xlab="Total Protein Interactors")


#Now we substract the MreIR5 AdjMatrix to the MreIR24 AdjMatrix
MreIR24.AdjMat.Change2<-Mre.IR24.M-Mre.IR8.M
range(MreIR24.AdjMat.Change)#Let's check the magnitude of the substracted weights
#[1]-0.813  0.668
any(is.na(MreIR24.AdjMat.Change))
#Generate the breaks and color functions
breaks <- c(seq(-0.5, 0.5, length=6))
breaks
HMCol<- brewer.pal(5,name = "RdBu")
HMCol<- rev(HMCol)
HMCol
[1] "#2166AC" "#4393C3" "#92C5DE" "#D1E5F0" "#F7F7F7" "#FDDBC7" "#F4A582" "#D6604D" "#B2182B"

#Plot the HeatMap
par(mar=rep(0,4)) # make plot margins 0
heatmap(MreIR24.AdjMat.Change2, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=HMCol, breaks=breaks,frame=T, cexRow = 0.5 , cexCol = 0.5, main = "Mre11 IR24 AdjMat (minus IR8h edges)", ylab="Total Protein Interactors",xlab="Total Protein Interactors")

#
#check with the FNN
# now actually make the heat map
par(mar=rep(0,4)) # make plot margins 0
heatmap(FNN.M, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=0:1, frame=T, cexRow = 0.5 , cexCol = 0.5)
#
legend("topleft", c("Neighbors Mre11","Mre11"), fill=c("red","blue"),  
       bty="n", inset=0, xpd=T,  border=F)


####################### Mre11 nonIR24h Networks
###################################################

#Mre11 nonIR24 Networks


#Let's keep the Mat as the Node attribute table and load the Neighbor list as the matrixes
Mat# Node attribute matrix for all the networks

#Mre11ID nonIR24 nodes Network Edge Lists
Mre_nonIR24<-read.delim("./Mre11_Files/Mre11_nonIR24.tsv",sep = "\t",col.names = )
str(Mre_nonIR24)
head(Mre_nonIR24)

#Generate the igraph network object
Mre.nonIR24<-graph_from_data_frame(d=Mre_nonIR24,vertices = Mat, directed = F )#Generate the igraph object with the edge and node annotation attributes
str(vertex_attr(Mre.nonIR24))#Check node attributes
str(edge_attr(Mre.nonIR24))#Check for edge attributes
Mre.nonIR24<-igraph::simplify(Mre.nonIR24, remove.multiple=TRUE, remove.loops=TRUE) #Remove repetitive edges and self-connecting edges
Mre.nonIR24
#Small plot
plot(Mre.nonIR24)

#Generate an adjacency Matrix that will be used after
Mre.nonIR24_AdjMat<-as_adjacency_matrix(Mre.nonIR24,type="both",attr = "weight",names=TRUE,sparse = TRUE)
Mre.nonIR24_AdjMat<-as.matrix(Mre.nonIR24_AdjMat)
write.csv(Mre.nonIR24_AdjMat,"./Mre.nonIR24_AdjMat.csv")

#Generate a vertex size based on degree and count
vsize<-igraph::degree(Mre.nonIR24)/vcount(Mre.nonIR24)
#Identify and remove disconnected nodes
isolated<-which(igraph::degree(Mre.nonIR24)==0)
isolated
length(isolated)# [1] 754
#remove isolated nodes
Mre.nonIR24<-igraph::delete_vertices(Mre.nonIR24,isolated)


#Generate layouts for the graph
#Do not run again, to not to change the layouts that I already have
Lyt.Nice<-layout_nicely(Mre.nonIR24,dim = 3)
Lyt.GOpt<-layout_with_graphopt(Mre.nonIR24, niter=5000, mass = 30, charge = E(Mre.nonIR24)$weight)
Lyt.DrL<-layout_with_drl(Mre.nonIR24, weights = E(Mre.nonIR24)$weight)
Lyt.FR<-layout.fruchterman.reingold(Mre.nonIR24, niter=10000, weights=E(Mre.nonIR24)$weight)#This one looks like the best option
Lyt.KK<-layout_with_kk(Mre.nonIR24, maxiter=10000, weights=E(Mre.nonIR24)$weight) #very large separation of peripheral nodes
Lyt.GEM<-layout_with_gem(Mre.nonIR24, maxiter=10000) #Bad separation
#Lyt.DH<-layout_with_dh(Mre.nonIR24, maxiter=150, fineiter=max(10, log2(833)))
Lyt.LGL<-layout.lgl(Mre.nonIR24, maxit=150,maxdelta=833,area=(length(E(Mre.nonIR24)))^2,coolexp=1.5,repulserad=((length(E(Mre.nonIR24)))^2)*(length(E(Mre.nonIR24))))#Do not work


#Star and Tree layouts

V(Mre.nonIR24)[name=="MRE11A"] #There is no more MRe11 in this time point
which(V(Mre.nonIR24)[name=="MRE11A"])# 228 #Check what is the rowindex for your desired node
which(V(Mre.nonIR24)$name=="NCL")# 176 #Check what is the rowindex for your desired node
which(V(Mre.nonIR24)$name=="NCL")# 34
which(V(Mre.nonIR24)$name=="MRE11A")# 722
Lyt.Star<-layout_as_star(Mre.nonIR24,center=V(Mre.nonIR24)[name=="NCL"]) #it has a special order argument where we can assign an order to the vertex, hence maybe the OPLS-DA:VIP value order may be used as ordering vector

Lyt.Tree<-layout_as_tree(Mre.nonIR24,root=c(34))
Lyt.RT<-layout.reingold.tilford(Mre.nonIR24,root=c(34))# This one looks good, but the indexing do not allow to highlight the correct nodes
Lyt.sugi<-layout_with_sugiyama(Mre.nonIR24, maxiter=100, weights=E(Mre.nonIR242)$weight)#This one looks good for undirected tree layouts

#Check what layout looks fine

plot.igraph(Mre.nonIR24,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(Mre.nonIR24)$An.color, vertex.frame.color=V(Mre.nonIR24)$An.color, vertex.shape="circle",
            vertex.size=10, vertex.label=V(Mre.nonIR24)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(Mre.nonIR24)$weight, edge.lty=1, main="Mre11ID nonIR24 Network", sub="FR Tree Layout",
            rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#With ggraph

ggraph(Mre.nonIR24,layout=Lyt.FR[-isolated,])+
  geom_edge_fan(color="gray65", width=1, alpha=0.4)+
  geom_node_point(aes(color=V(Mre.nonIR24)$Annotation, size=V(Mre.nonIR24)$Mre.FC.nonIR24),alpha=0.75)+
  geom_node_text(aes(label=name),color="gray25",size=4)+
  scale_color_manual(limits=V(Mre.nonIR24)$Annotation, values=V(Mre.nonIR24)$An.color)+
  scale_size_continuous(range = c(5, 15)) +
  theme_void()+
  title("Mre11 nonIR24 Network")+
  guides(color=guide_legend(override.aes = list(size=7)))+
  theme(legend.text=element_text(size=12,face = "bold"),
        legend.title = element_text(size=13,face = "bold"))



#Community Detection
#clp <- igraph::cluster_optimal(Mre.nonIR24)
class(clp)
V(Mre.nonIR24)$community <- clp$membership
colrs <- adjustcolor( c("gray50", "tomato", "gold", "yellowgreen"), alpha=.6)


plot.igraph(clp,Mre.nonIR24,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(Mre.nonIR24)$An.color, vertex.frame.color=V(Mre.nonIR24)$An.color, vertex.shape="circle",
            vertex.size==V(Mre.nonIR24)$Mre.FC.nonIR24, vertex.label=V(Mre.nonIR24)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(Mre.nonIR24)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#Distances edges highlight
dist.from.Mre11 <- distances(Mre.nonIR24, v=V(Mre.nonIR24)[name=="NCL"],
                             to=V(Mre.nonIR24), weights=E(Mre.nonIR24)$weight)
dist.from.Mre11<-dist.from.Mre11[,!is.infinite(colSums(dist.from.Mre11))]

dist.from.Mre11
range(dist.from.Mre11)

#Set colors to plot the distances:
oranges <- colorRampPalette(c("dark red", "gold"))
col <- oranges(max(dist.from.Mre11))
col <- col[dist.from.Mre11+1]
plot(Mre.nonIR24, layout=Lyt.FR[-isolated,],vertex.color=col, vertex.label=dist.from.Mre11,
     vertex.label.color="white",vertex.shape="circle",
     vertex.size=5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(Mre.nonIR24)$weight, edge.lty=1, main="Mre11 nonIR24 Full Nuclear Network", sub="FR Layout",
     rescale=T)


#Shortest path from Mre11 to Nucleolin
MreNcl.path <- shortest_paths(Mre.nonIR24,
                              from = V(Mre.nonIR24)[name=="XRCC5"],
                              to = V(Mre.nonIR24)[name=="NCL"],
                              output = "both") # both path nodes and edges
MreNcl.path
# Generate edge color variable to plot the path:
ecol <- rep("gray80", ecount(Mre.nonIR24))
ecol[unlist(MreNcl.path$epath)] <- "orange"
# Generate edge width variable to plot the path:
ew <- rep(2, ecount(Mre.nonIR24))
ew[unlist(MreNcl.path$epath)] <- 4
# Generate node color variable to plot the path:
vcol <- rep("gray40", vcount(Mre.nonIR24))
vcol[unlist(MreNcl.path$vpath)] <- "gold"
plot(Mre.nonIR24, layout=Lyt.FR[-isolated,], vertex.color=vcol, edge.color=ecol,
     edge.width=ew,vertex.label.color="white",vertex.shape="circle",
     vertex.size=5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(Mre.nonIR24)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)



#Plot based on Mre11 as incident
inc.edges <- incident(Mre.nonIR24, V(Mre.nonIR24)[name=="NCL"], mode="all")
# Set colors to plot the selected edges.
ecol <- rep("gray80", ecount(Mre.nonIR24))
ecol[inc.edges] <- "orange"
vcol <- rep("gray40", vcount(Mre.nonIR24))
vcol[V(Mre.nonIR24)$name=="Mre11"] <- "gold"

plot(Mre.nonIR24, layout=Lyt.FR[-isolated,], vertex.color=vcol, edge.color=ecol, 
     edge.width=ew,vertex.label.color="white",vertex.shape="circle",
     vertex.size=5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(Mre.nonIR24)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)


#Plot Degree Distribution
Mre.nonIR24.Deg<-degree_distribution(Mre.nonIR24, cumulative = T, mode="all")
plot( x=0:max(igraph::degree(Mre.nonIR24)), y=1-Mre.nonIR24.Deg, pch=19, cex=1.2, col="orange",
      xlab="Degree", ylab="Cumulative Frequency")

#Network Descriptives
#Number of nodes and edges
vcount(Mre.nonIR24) #[1] 79
ecount(Mre.nonIR24)#[1] 185
#Density (Number of edges in relation to the total number of edges in a full connected graph)
edge_density(Mre.nonIR24,loops = F) #[1]0.06004544
#Transitivity (Total number of conected triangles)
transitivity(Mre.nonIR24,type="global") #[1] 0.5466786
#modularity
mod<-V(Mre.nonIR24)$S6.FuzzyCMeans #lets check for fuzzy Cmeans modularity
mod<-mod+1#The non clustered interactors that appears as 0, needs to be transformed to 1, so we add 1 to every cluster
modularity(Mre.nonIR24,membership = mod)
#S6.FuzzyCMeans modularity [1]  0.06037984
#Ku.FuzzyCMeans modularity [1]-0.02515705
#Mre.FuzzyCMeans modularity [1] -0.005273923

#Assorativity
assortativity_nominal(Mre.nonIR24,factor(V(Mre.nonIR24)$Annotation),directed=F)
#Annotation [1] [1] 0.1663444
#S6.FuzzyCMeans [1] 0.09615429  
#Ku.FuzzyCMeans [1] -0.04335566
#Mre.FuzzyCMeans [1]-0.04872452

#Diameter (longest geodesic distance), shortest path between two most peripheral nodes
diameter(Mre.nonIR24,directed = F, weights = E(Mre.nonIR24)$weights)#[1] 4.888
diam<-get_diameter(Mre.nonIR24,directed = F,weights = E(Mre.nonIR24)$weights)
diam #+ 10/79 vertices, named, from 53e27c9: [1] UTP4    EXOSC10 RPL5    NPM1    ILF3    PRKDC   RUVBL2  POLR2E  DDB1    WDR5    
#Plot Diameter
vcol <- rep("gray50", vcount(Mre.nonIR24))
vcol[diam] <- "gold"
ecol <- rep("gray90", ecount(Mre.nonIR24))
ecol[E(Mre.nonIR24, path=diam)] <- "orange"
# E(net, path=diam) finds edges along a path, here 'diam'
plot(Mre.nonIR24, vertex.color=vcol, edge.color=ecol, vertex.label.color="black",edge.arrow.mode=0,layout=Lyt.FR[-isolated,],vertex.size=5)

#Node topological descriptives
degree<-igraph::centr_degree(Mre.nonIR24,mode = "all",normalized = T)
degree<-degree$res
eigen_cent<-igraph::eigen_centrality(Mre.nonIR24)
eigen_cent<-eigen_cent$vector
clos_cent <- igraph::closeness(Mre.nonIR24, normalized=T)
betw_cent <- igraph::betweenness(Mre.nonIR24, normalized=T)
hs <- hub_score(Mre.nonIR24, weights=E(Mre.nonIR24)$weight)$vector
as <- authority_score(Mre.nonIR24, weights=E(Mre.nonIR24)$weight)$vector
Trans<-transitivity(Mre.nonIR24,type="local")
Trans
#write results
Mre.nonIR24.Centrality<-data.frame(NAME=V(Mre.nonIR24)$name, UniprotID=V(Mre.nonIR24)$Uniprot, Degree=degree, Eigen= eigen_cent, Closeness= clos_cent, Betweeness= betw_cent, Hubs=hs, Authority=as,Transitivity=Trans)
Mre.nonIR24.Centrality
write.csv(Mre.nonIR24.Centrality,"./Mre.nonIR24.Centrality.csv")

#Cocitation, Number of shared interactions among nodes
Cocitation<-cocitation(Mre.nonIR24)
Cocitation#This matrix is similar to an adjacency matrix
write.csv(Cocitation,"./Mre.nonIR24.cocitation.csv")
#Plot with degree as node size
ggraph(Mre.nonIR24,layout=Lyt.FR[-isolated,])+
  geom_edge_fan(color="gray85", width=1, alpha=0.4)+
  geom_node_point(aes(color=V(Mre.nonIR24)$Annotation, size=degree),alpha=0.75)+
  geom_node_text(aes(label=name),color="gray25",size=4)+
  scale_color_manual(limits=V(Mre.nonIR24)$Annotation, values=V(Mre.nonIR24)$An.color)+
  scale_size_continuous(range = c(5, 15)) +
  theme_void()+
  title("Mre11 nonIR24 Network")+
  guides(color=guide_legend(override.aes = list(size=7)))+
  theme(legend.text=element_text(size=12,face = "bold"),
        legend.title = element_text(size=13,face = "bold"))


#Modularity Analysis

#Edge_Betweeness Modularity
Mre.nonIR24.ebmod=edge.betweenness.community(Mre.nonIR24,weights=E(Mre.nonIR24)$weight,directed = F) 
Mre.nonIR24.ebmod# IGRAPH clustering edge betweenness, groups: 13, mod: 0.63

#
plot.igraph(Mre.nonIR24.ebmod,Mre.nonIR24,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(Mre.nonIR24)$An.color, vertex.frame.color=V(Mre.nonIR24)$An.color, vertex.shape="circle",
            vertex.size=5, vertex.label=V(Mre.nonIR24)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(Mre.nonIR24)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)

#
plot(Mre.nonIR24.ebmod, Mre.nonIR24,layout=Lyt.FR[-isolated,], alpha=0.6, 
     vertex.color=V(Mre.nonIR24)$An.color, vertex.frame.color=V(Mre.nonIR24)$An.color, vertex.shape="circle", 
     vertex.size=5, vertex.label=V(Mre.nonIR24)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black", 
     edge.width=E(Mre.nonIR24)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout", 
     rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)
#Dendogram of the Edge_Betweeness
dendPlot(Mre.nonIR24.ebmod, mode="hclust")

#Louvain modularity
Mre.nonIR24.Louvain.mod<-cluster_louvain(Mre.nonIR24,weights=E(Mre.nonIR24)$weight)
Mre.nonIR24.Louvain.mod #GRAPH clustering multi level, groups: 9, mod: 0.65

#Plot the Louvain communities
plot(Mre.nonIR24.Louvain.mod, Mre.nonIR24,layout=Lyt.FR[-isolated,], alpha=0.6, 
     vertex.color=V(Mre.nonIR24)$An.color, vertex.frame.color=V(Mre.nonIR24)$An.color, vertex.shape="circle", 
     vertex.size=5, vertex.label=V(Mre.nonIR24)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black", 
     edge.width=E(Mre.nonIR24)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout", 
     rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#Generate HeatMaps of the adjacency matrixes
#Let's use the neighbor list for generating the Adjaency Matrix
SensorID_FNN #Let's use the edge list that contain the total amount of possible interactors for the three sensors
v<-levels(as.factor(unlist(SensorID_FNN[,1:2]))) #get vertex names from the network list
n<-length(v)#get the number of vertex in the network
e<-matrix(match(as.character(unlist(SensorID_FNN[,1:2])),v),ncol = 2)#Edgelist with row index instead of Names
w<-Mre_nonIR24$weight #Edge weights, only from the Mre_nonIR24 network to quit those edges that do not exist at Mre_nonIR24

# M is our co-expression network adjacency matrix.

# distringuish between i,j and j,i we want to make M symmetric (undirected) by
# considering an undirected network.
Mre.nonIR24.M <- matrix(0, n, n)                  # set up a co-expression matrix
Mre.nonIR24.M[e] <- w                             # fill it in with edge weights
Mre.nonIR24.M <- Mre.nonIR24.M + t(Mre.nonIR24.M)                         # make this symmetric
any(is.na(Mre.nonIR24.M))#check no NAs in the matrix
dimnames(Mre.nonIR24.M) <- list(v, v)             # label the vertices
# We let A_ij = 1 if there is an interaction. 
#A <- 1*(M > 0) #right now we keep the weights as they are
Mre.nonIR24.M #The adjacency Matrix not as an igraph object, easier to manipulate, and to subset for smaller Sensor specific networks
quantile(Mre.nonIR24.M,c(0,0.25,0.5,0.75,1))
#
#      0%  25%  50%  75% 100% 
#    0.00 0.00 0.00 0.00 0.99   

#Creating the HeatMap
# make a heatmap image with colors: white below 0,
# grey from 0.4 to 1 with darker closer to 1
breaks <- c(0, seq(0.4, 0.991, length=11)) # breaks in the color bins
cols <- grey(1-c(0,seq(0.5,1,length=10))) # colors

#Red Blue Colors
library(RColorBrewer)
breaks <- c(seq(-1, 1, length=11))
breaks
cols <- gray(1-c(0,seq(0,1,length=10))) # colors
cols2 <- hcl.colors(10 , palette = "Blues3" , alpha = 1 , rev = TRUE)
cols3 <- colorRamp("red" , space = "rgb" , interpolate = "linear")
HMCol<- brewer.pal(10,name = "RdBu")
HMCol<- rev(HMCol)
HMCol


# color the vertex for CCNB1 blue, its neighbors red, and the others white
SensorIDBind <- which(v == "Mre11") # index for ccnb1
#SensorIDBind<-as.factor(SensorIDBind)
vcols <- rep("white",n)
vcols[SensorIDBind] <- "blue"
#Subsetting Mre.nonIR24.M by doing a loop asking for look into every row and column the index in SensorIDBind worked
for (i in SensorIDBind) {
  vcols[which(Mre.nonIR24.M[,i] > 0 | Mre.nonIR24.M[i,])] <- "red"
}
##Color annotation for functional annotation
MatSet<-subset(Mat,Mat$NAME%in%v)# select only those proteins that are cointained in the node list of our matrixes
match(v,MatSet$NAME) # Check that the indexes of the Annotation Dataframe and the node list are in the same position 
match(rownames(Mre.nonIR24.M),MatSet$NAME) #Check that the row names of the Adj Matrix and the Annotation dataframe are in the same order
match(colnames(Mre.nonIR24.M),MatSet$NAME) #Check that the column names of the Adj Matrix and the Annotation dataframe are in the same order
AnCol<-MatSet$An.color#Vector of colors
AnCol<-as.character(AnCol)#transform it into a character vector, so the colors can be read

#Subsetting Mre.nonIR24.M using the index vector in SensorIDBind do not work
#vcols[which(Mre.nonIR24.M[,SensorIDBind] > 0 | Mre.nonIR24.M[SensorIDBind,])] <- "red"

#Generate a dendogram for keeping the proteins in the sample place
#Mre.nonIR24.dist<-dist(Mre.nonIR24.M,method = "euclidian")
#Mre.nonIR24.Dend<-hclust(Mre.nonIR24.dist, method = "ward.D" )
#plot(Mre.nonIR24.Dend)

#Cut the dendogram
#cut.height <- 10 # try numbers from 4 to 7
#FNN.hclusts<- cutree(FNN.Dend, cut.height)
#plot(FNN.hclusts)

# now actually make the heat map
par(mar=rep(0,4)) # make plot margins 0
heatmap(Mre.nonIR24.M, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=cols, breaks=breaks,frame=T, cexRow = 0.5 , cexCol = 0.5, main = "Mre11 nonIR24 AdjMat", ylab="Total Protein Interactors",xlab="Total Protein Interactors")
#
plot.new()
legend("topleft", c("Cell Cycle","Chromatin Org","DNA Repair","Nuclear Pore","p53","Proteasome","Ribosome Biogenesis","RNA Metabolism","Transcription Reg","UPR"), 
       fill=c("coral1","darkviolet","firebrick1","deepskyblue3","tomato","grey","lightgreen","gold","slateblue1","maroon"),  
       bty="n", inset=0, xpd=T,  border=F, cex =2,y.intersp = 0.7,x.intersp = 0.2,pt.cex = 1)


#Now we substract the MrenonnonIR24 AdjMatrix to the MreIR6 AdjMatrix
MrenonIR24.AdjMat.Change<-Mre.nonIR24.M-Mre.nonIR5.M
range(MrenonIR24.AdjMat.Change)#Let's check the magnitude of the substracted weights
#[1] -0.806  0.805[1] -0.806  0.805
any(is.na(MrenonIR24.AdjMat.Change))
#Generate the breaks and color functions
breaks <- c(seq(-0.5, 0.5, length=6))
breaks
HMCol<- brewer.pal(5,name = "RdBu")
HMCol<- rev(HMCol)
HMCol
[1] "#2166AC" "#4393C3" "#92C5DE" "#D1E5F0" "#F7F7F7" "#FDDBC7" "#F4A582" "#D6604D" "#B2182B"
[1] "#0571B0" "#92C5DE" "#F7F7F7" "#F4A582" "#CA0020"
#Plot the HeatMap
par(mar=rep(0,4)) # make plot margins 0
heatmap(MrenonIR24.AdjMat.Change, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=HMCol, breaks=breaks,frame=T, cexRow = 0.5 , cexCol = 0.5, main = "Mre11 nonIR24 AdjMat (minus nonIR5 edges)", ylab="Total Protein Interactors",xlab="Total Protein Interactors")


#Now we substract the MreIR5 AdjMatrix to the MrenonIR24 AdjMatrix
MrenonIR24.AdjMat.Change2<-Mre.nonIR24.M-Mre.IR24.M
range(MrenonIR24.AdjMat.Change)#Let's check the magnitude of the substracted weights
#[1]-0.813  0.668
any(is.na(MrenonIR24.AdjMat.Change))
#Generate the breaks and color functions
breaks <- c(seq(-0.5, 0.5, length=10))
breaks
HMCol<- brewer.pal(9,name = "RdBu")
HMCol<- rev(HMCol)
HMCol
[1] "#2166AC" "#4393C3" "#92C5DE" "#D1E5F0" "#F7F7F7" "#FDDBC7" "#F4A582" "#D6604D" "#B2182B"

#Plot the HeatMap
par(mar=rep(0,4)) # make plot margins 0
heatmap(MrenonIR24.AdjMat.Change2, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=HMCol, breaks=breaks,frame=T, cexRow = 0.5 , cexCol = 0.5, main = "Mre11 IR24 AdjMat (minus IR24 edges)", ylab="Total Protein Interactors",xlab="Total Protein Interactors")

#
#check with the FNN
# now actually make the heat map
par(mar=rep(0,4)) # make plot margins 0
heatmap(FNN.M, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=0:1, frame=T, cexRow = 0.5 , cexCol = 0.5)
#
legend("topleft", c("Neighbors Mre11","Mre11"), fill=c("red","blue"),  
       bty="n", inset=0, xpd=T,  border=F)

#
plot.new()

#################################################
#Mre80ID FNN nodes Network Edge Lists: This network contains all the Mre11ID interctors, in contrast with the FNN 
Mre_FNN<-read.delim("./Mre11_Files/Mre_FNN.tsv",sep = "\t",col.names = )
str(Mre_FNN)
head(Mre_FNN)


#Generate the igraph network object
Mre.FNN<-graph_from_data_frame(d=Mre_FNN,vertices = Mat, directed = F )#Generate the igraph object with the edge and node annotation attributes
str(vertex_attr(Mre.FNN))#Check node attributes
str(edge_attr(Mre.FNN))#Check for edge attributes
Mre.FNN<-igraph::simplify(Mre.FNN, remove.multiple=TRUE, remove.loops=TRUE) #Remove repetitive edges and self-connecting edges
Mre.FNN
#Small plot
plot(Mre.FNN)

#Generate an adjacency Matrix that will be used after
Mre.FNN_AdjMat<-as_adjacency_matrix(Mre.FNN,type="both",attr = "weight",names=TRUE,sparse = TRUE)
Mre.FNN_AdjMat<-as.matrix(Mre.FNN_AdjMat)
write.csv(Mre.FNN_AdjMat,"./Mre.FNN_AdjMat.csv")

#Generate a vertex size based on degree and count
vsize<-igraph::degree(Mre.FNN)/vcount(Mre.FNN)
#Identify and remove disconnected nodes
isolated<-which(igraph::degree(Mre.FNN)==0)
isolated
length(isolated)# [1] 599
#remove isolated nodes
Mre.FNN<-igraph::delete_vertices(Mre.FNN,isolated)


#Generate layouts for the graph
#Do not run again, to not to change the layouts that I already have
Lyt.Nice<-layout_nicely(Mre.FNN,dim = 3)
Lyt.GOpt<-layout_with_graphopt(Mre.FNN, niter=5000, mass = 30, charge = E(Mre.FNN)$weight)
Lyt.DrL<-layout_with_drl(Mre.FNN, weights = E(Mre.FNN)$weight)
Lyt.FR<-layout.fruchterman.reingold(Mre.FNN, niter=10000, weights=E(Mre.FNN)$weight)#This one looks like the best option
Lyt.KK<-layout_with_kk(Mre.FNN, maxiter=10000, weights=E(Mre.FNN)$weight) #very large separation of peripheral nodes
Lyt.GEM<-layout_with_gem(Mre.FNN, maxiter=10000) #Bad separation
#Lyt.DH<-layout_with_dh(Mre.FNN, maxiter=150, fineiter=max(10, log2(833)))
Lyt.LGL<-layout.lgl(Mre.FNN, maxit=150,maxdelta=833,area=(length(E(Mre.FNN)))^2,coolexp=1.5,repulserad=((length(E(Mre.FNN)))^2)*(length(E(Mre.FNN))))#Do not work


#Star and Tree layouts

V(Mre.FNN)[name=="SIRT6"]
which(V(Mre.FNN)[name=="SIRT6"])# 610 #Check what is the rowindex for your desired node
which(V(Mre.FNN)$name=="SIRT6")# 610 #Check what is the rowindex for your desired node
which(V(Mre.FNN)$name=="MRE11A")# 327
which(V(Mre.FNN)$name=="XRCC5")# 722
Lyt.Star<-layout_as_star(Mre.FNN,center=V(Mre.FNN)[name=="SIRT6"]) #it has a special order argument where we can assign an order to the vertex, hence maybe the OPLS-DA:VIP value order may be used as ordering vector

Lyt.Tree<-layout_as_tree(Mre.FNN,root=c(190))
Lyt.RT<-layout.reingold.tilford(Mre.FNN,root=c(190))# This one looks good, but the indexing do not allow to highlight the correct nodes
Lyt.sugi<-layout_with_sugiyama(Mre.FNN, maxiter=100, weights=E(Mre.FNN2)$weight)#This one looks good for undirected tree layouts

#Check what layout looks fine

plot.igraph(Mre.FNN,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(Mre.FNN)$An.color, vertex.frame.color=V(Mre.FNN)$An.color, vertex.shape="circle",
            vertex.size=5, vertex.label=V(Mre.FNN)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(Mre.FNN)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)

#With ggraph

ggraph(Mre.FNN,layout=Lyt.FR[-isolated,])+
  geom_edge_fan(color="gray65", width=1, alpha=0.4)+
  geom_node_point(aes(color=V(Mre.FNN)$Annotation, size=V(Mre.FNN)$Mre.FC.nonIR5),alpha=0.75)+
  geom_node_text(aes(label=name),color="gray25",size=4)+
  scale_color_manual(limits=V(Mre.FNN)$Annotation, values=V(Mre.FNN)$An.color)+
  scale_size_continuous(range = c(5, 15)) +
  theme_void()+
  title("Sirt6 nonIR5 Network")+
  guides(color=guide_legend(override.aes = list(size=7)))+
  theme(legend.text=element_text(size=12,face = "bold"),
        legend.title = element_text(size=13,face = "bold"))



#Community Detection
#clp <- igraph::cluster_optimal(Mre.FNN)
class(clp)
V(Mre.FNN)$community <- clp$membership
colrs <- adjustcolor( c("gray50", "tomato", "gold", "yellowgreen"), alpha=.6)


plot.igraph(clp,Mre.FNN,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(Mre.FNN)$An.color, vertex.frame.color=V(Mre.FNN)$An.color, vertex.shape="circle",
            vertex.size==V(Mre.FNN)$Mre.FC.nonIR5, vertex.label=V(Mre.FNN)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(Mre.FNN)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#Distances edges highlight
dist.from.Sirt6 <- distances(Mre.FNN, v=V(Mre.FNN)[name=="SIRT6"],
                             to=V(Mre.FNN), weights=E(Mre.FNN)$weight)
dist.from.Sirt6<-dist.from.Sirt6[,!is.infinite(colSums(dist.from.Sirt6))]

dist.from.Sirt6
range(dist.from.Sirt6)

#Set colors to plot the distances:
oranges <- colorRampPalette(c("dark red", "gold"))
col <- oranges(max(dist.from.Sirt6))
col <- col[dist.from.Sirt6+1]
plot(Mre.FNN, layout=Lyt.FR[-isolated,],vertex.color=col, vertex.label=dist.from.Sirt6,
     vertex.label.color="white",vertex.shape="circle",
     vertex.size=V(Mre.FNN)$Mre.FC.nonIR5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(Mre.FNN)$weight, edge.lty=1, main="Sirt6 nonIR5 Full Nuclear Network", sub="FR Layout",
     rescale=T)


#Shortest path from Sirt6 to Nucleolin
MreNcl.path <- shortest_paths(Mre.FNN,
                              from = V(Mre.FNN)[name=="SIRT6"],
                              to = V(Mre.FNN)[name=="NCL"],
                              output = "both") # both path nodes and edges
# Generate edge color variable to plot the path:
ecol <- rep("gray80", ecount(Mre.FNN))
ecol[unlist(MreNcl.path$epath)] <- "orange"
# Generate edge width variable to plot the path:
ew <- rep(2, ecount(Mre.FNN))
ew[unlist(MreNcl.path$epath)] <- 4
# Generate node color variable to plot the path:
vcol <- rep("gray40", vcount(Mre.FNN))
vcol[unlist(MreNcl.path$vpath)] <- "gold"
plot(Mre.FNN, layout=Lyt.FR[-isolated,], vertex.color=vcol, edge.color=ecol,
     edge.width=ew,vertex.label.color="white",vertex.shape="circle",
     vertex.size=5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(Mre.FNN)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)



#Plot based on Sirt6 as incident
inc.edges <- incident(Mre.FNN, V(Mre.FNN)[name=="SIRT6"], mode="all")
# Set colors to plot the selected edges.
ecol <- rep("gray80", ecount(Mre.FNN))
ecol[inc.edges] <- "orange"
vcol <- rep("gray40", vcount(Mre.FNN))
vcol[V(Mre.FNN)$name=="SIRT6"] <- "gold"

plot(Mre.FNN, layout=Lyt.FR[-isolated,], vertex.color=vcol, edge.color=ecol, 
     edge.width=ew,vertex.label.color="white",vertex.shape="circle",
     vertex.size=V(Mre.FNN)$Mre.FC.nonIR5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(Mre.FNN)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)


#Plot Degree Distribution
Mre.FNN.Deg<-degree_distribution(Mre.FNN, cumulative = T, mode="all")
plot( x=0:max(igraph::degree(Mre.FNN)), y=1-Mre.FNN.Deg, pch=19, cex=1.2, col="orange",
      xlab="Degree", ylab="Cumulative Frequency")

#Network Descriptives
#Number of nodes and edges
vcount(Mre.FNN) #[1] 493
ecount(Mre.FNN)#[1] 4454
#Density (Number of edges in relation to the total number of edges in a full connected graph)
edge_density(Mre.FNN,loops = F) #[1]0.03672554
#Transitivity (Total number of conected triangles)
transitivity(Mre.FNN,type="global") #[1] 0.5981853
#modularity
mod<-V(Mre.FNN)$S6.FuzzyCMeans #lets check for fuzzy Cmeans modularity
mod<-mod+1#The non clustered interactors that appears as 0, needs to be transformed to 1, so we add 1 to every cluster
modularity(Mre.FNN,membership = mod)
#S6.FuzzyCMeans modularity [1] 0.001156511
#Ku.FuzzyCMeans modularity [1]0.01212958
#Mre.FuzzyCMeans modularity [1]0.008387515

#Assorativity
assortativity_nominal(Mre.FNN,factor(V(Mre.FNN)$Annotation),directed=F)
#Annotation [1] [1]0.2523593
#S6.FuzzyCMeans [1] 0.001760367
#Ku.FuzzyCMeans [1] 0.02227818
#Mre.FuzzyCMeans [1]0.01030464

#Diameter (longest geodesic distance), shortest path between two most peripheral nodes
diameter(Mre.FNN,directed = F, weights = E(Mre.FNN)$weights)#3.87
diam<-get_diameter(Mre.FNN,directed = F,weights = E(Mre.FNN)$weights)
diam #+  8/493 vertices, named, from dda6b83: [1] PELP1  WDR18  RPL5   NCL    PPP1CC CORO1C ACTR2  ACTR3  
#Plot Diameter
vcol <- rep("gray50", vcount(Mre.FNN))
vcol[diam] <- "gold"
ecol <- rep("gray90", ecount(Mre.FNN))
ecol[E(Mre.FNN, path=diam)] <- "orange"
# E(net, path=diam) finds edges along a path, here 'diam'
plot(Mre.FNN, vertex.color=vcol, edge.color=ecol, vertex.label.color="black",edge.arrow.mode=0,layout=Lyt.FR[-isolated,],vertex.size=V(Mre.FNN)$Mre.FC.nonIR5)

#Node topological descriptives
degree<-igraph::centr_degree(Mre.FNN,mode = "all",normalized = T)
degree<-degree$res
eigen_cent<-igraph::eigen_centrality(Mre.FNN)
eigen_cent<-eigen_cent$vector
clos_cent <- igraph::closeness(Mre.FNN, normalized=T)
betw_cent <- igraph::betweenness(Mre.FNN, normalized=T)
hs <- hub_score(Mre.FNN, weights=E(Mre.FNN)$weight)$vector
as <- authority_score(Mre.FNN, weights=E(Mre.FNN)$weight)$vector
Trans<-transitivity(Mre.FNN,type="local")
Trans
#write results
Mre.FNN.Centrality<-data.frame(NAME=V(Mre.FNN)$name, UniprotID=V(Mre.FNN)$Uniprot, Degree=degree, Eigen= eigen_cent, Closeness= clos_cent, Betweeness= betw_cent, Hubs=hs, Authority=as,Transitivity=Trans)
Mre.FNN.Centrality
write.csv(Mre.FNN.Centrality,"./Mre.FNN.Centrality.csv")

#Cocitation, Number of shared interactions among nodes
Cocitation<-cocitation(Mre.FNN)
Cocitation#This matrix is similar to an adjacency matrix
write.csv(Cocitation,"./Mre.FNN.cocitation.csv")
#Plot with degree as node size
ggraph(Mre.FNN,layout=Lyt.FR[-isolated,])+
  geom_edge_fan(color="gray85", width=1, alpha=0.4)+
  geom_node_point(aes(color=V(Mre.FNN)$Annotation, size=degree),alpha=0.75)+
  geom_node_text(aes(label=name),color="gray25",size=4)+
  scale_color_manual(limits=V(Mre.FNN)$Annotation, values=V(Mre.FNN)$An.color)+
  scale_size_continuous(range = c(5, 15)) +
  theme_void()+
  title("Sirt6 nonIR5 Network")+
  guides(color=guide_legend(override.aes = list(size=7)))+
  theme(legend.text=element_text(size=12,face = "bold"),
        legend.title = element_text(size=13,face = "bold"))


#Modularity Analysis

#Edge_Betweeness Modularity
Mre.FNN.ebmod=edge.betweenness.community(Mre.FNN,weights=E(Mre.FNN)$weight,directed = F) 
Mre.FNN.ebmod# groups: 15, mod: 0.48
Mre.FNN.ebmod=edge.betweenness.community(Mre.FNN,weights=E(Mre.FNN)$weight,directed = F) 
Mre.FNN.ebmod #37 groups, good enough mod 0.47
#
plot.igraph(Mre.FNN.ebmod,Mre.FNN,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(Mre.FNN)$An.color, vertex.frame.color=V(Mre.FNN)$An.color, vertex.shape="circle",
            vertex.size=5, vertex.label=V(Mre.FNN)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(Mre.FNN)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)

#
plot(Mre.FNN.ebmod, Mre.FNN,layout=Lyt.FR[-isolated,], alpha=0.6, 
     vertex.color=V(Mre.FNN)$An.color, vertex.frame.color=V(Mre.FNN)$An.color, vertex.shape="circle", 
     vertex.size=5, vertex.label=V(Mre.FNN)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black", 
     edge.width=E(Mre.FNN)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout", 
     rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)
#Dendogram of the Edge_Betweeness
dendPlot(Mre.FNN.ebmod, mode="hclust")

#Louvain modularity
Mre.FNN.Louvain.mod<-cluster_louvain(Mre.FNN,weights=E(Mre.FNN)$weight)
Mre.FNN.Louvain.mod #IGRAPH clustering multi level, groups: 11, mod: 0.43

#Plot the Louvain communities
plot(Mre.FNN.Louvain.mod, Mre.FNN,layout=Lyt.FR[-isolated,], alpha=0.6, 
     vertex.color=V(Mre.FNN)$An.color, vertex.frame.color=V(Mre.FNN)$An.color, vertex.shape="circle", 
     vertex.size=5, vertex.label=V(Mre.FNN)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black", 
     edge.width=E(Mre.FNN)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout", 
     rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#Generate HeatMaps of the adjacency matrixes
#Let's use the neighbor list for generating the Adjaency Matrix
SensorID_FNN #Let's use the edge list that contain the total amount of possible interactors for the three sensors
v<-levels(as.factor(unlist(SensorID_FNN[,1:2]))) #get vertex names from the network list
n<-length(v)#get the number of vertex in the network
e<-matrix(match(as.character(unlist(SensorID_FNN[,1:2])),v),ncol = 2)#Edgelist with row index instead of Names
w<-Mre_FNN$weight #Edge weights, only from the Mre_FNN network to quit those edges that do not exist at Mre_FNN

# M is our co-expression network adjacency matrix.

# distringuish between i,j and j,i we want to make M symmetric (undirected) by
# considering an undirected network.
Mre.FNN.M <- matrix(0, n, n)                  # set up a co-expression matrix
Mre.FNN.M[e] <- w                             # fill it in with edge weights
Mre.FNN.M <- Mre.FNN.M + t(Mre.FNN.M)                         # make this symmetric
any(is.na(Mre.FNN.M))#check no NAs in the matrix
dimnames(Mre.FNN.M) <- list(v, v)             # label the vertices
# We let A_ij = 1 if there is an interaction. 
#A <- 1*(M > 0) #right now we keep the weights as they are
Mre.FNN.M #The adjacency Matrix not as an igraph object, easier to manipulate, and to subset for smaller Sensor specific networks
quantile(Mre.FNN.M,c(0,0.25,0.5,0.75,1))
#
#  0%   25%   50%   75%  100% 
#0.000 0.000 0.000 0.000 0.991 

#Creating the HeatMap
# make a heatmap image with colors: white below 0,
# grey from 0.4 to 1 with darker closer to 1
breaks <- c(0, seq(0.4, 0.991, length=11)) # breaks in the color bins
cols <- grey(1-c(0,seq(0.5,1,length=10))) # colors

#Red Blue Colors
library(RColorBrewer)
breaks <- c(seq(-1, 1, length=10))
breaks
cols <- gray(1-c(0,seq(0,1,length=10))) # colors
cols2 <- hcl.colors(10 , palette = "Blues3" , alpha = 1 , rev = TRUE)
cols3 <- colorRamp("red" , space = "rgb" , interpolate = "linear")
HMCol<- brewer.pal(10,name = "RdBu")
HMCol<- rev(HMCol)
HMCol


# color the vertex for CCNB1 blue, its neighbors red, and the others white
SensorIDBind <- which(v == "SIRT6") # index for ccnb1
#SensorIDBind<-as.factor(SensorIDBind)
vcols <- rep("white",n)
vcols[SensorIDBind] <- "blue"
#Subsetting Mre.FNN.M by doing a loop asking for look into every row and column the index in SensorIDBind worked
for (i in SensorIDBind) {
  vcols[which(Mre.FNN.M[,i] > 0 | Mre.FNN.M[i,])] <- "red"
}
##Color annotation for functional annotation
MatSet<-subset(Mat,Mat$NAME%in%v)# select only those proteins that are cointained in the node list of our matrixes
match(v,MatSet$NAME) # Check that the indexes of the Annotation Dataframe and the node list are in the same position 
match(rownames(Mre.FNN.M),MatSet$NAME) #Check that the row names of the Adj Matrix and the Annotation dataframe are in the same order
match(colnames(Mre.FNN.M),MatSet$NAME) #Check that the column names of the Adj Matrix and the Annotation dataframe are in the same order
AnCol<-MatSet$An.color#Vector of colors
AnCol<-as.character(AnCol)#transform it into a character vector, so the colors can be read

#Subsetting Mre.FNN.M using the index vector in SensorIDBind do not work
#vcols[which(Mre.FNN.M[,SensorIDBind] > 0 | Mre.FNN.M[SensorIDBind,])] <- "red"

#Generate a dendogram for keeping the proteins in the sample place
#Mre.FNN.dist<-dist(Mre.FNN.M,method = "euclidian")
#Mre.FNN.Dend<-hclust(Mre.FNN.dist, method = "ward.D" )
#plot(Mre.FNN.Dend)

#Cut the dendogram
#cut.height <- 10 # try numbers from 4 to 7
#FNN.hclusts<- cutree(FNN.Dend, cut.height)
#plot(FNN.hclusts)

# now actually make the heat map
par(mar=rep(0,4)) # make plot margins 0
heatmap(Mre.FNN.M, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=cols, breaks=breaks,frame=T, cexRow = 0.5 , cexCol = 0.5, main = "Sirt6 nonIR5 AdjMat", ylab="Total Protein Interactors",xlab="Total Protein Interactors")

#
plot.new()
legend("topleft", c("Cell Cycle","Chromatin Org","DNA Repair","Nuclear Pore","p53","Proteasome","Ribosome Biogenesis","RNA Metabolism","Transcription Reg","UPR"), 
       fill=c("coral1","darkviolet","firebrick1","deepskyblue3","tomato","grey","lightgreen","gold","slateblue1","maroon"),  
       bty="n", inset=0, xpd=T,  border=F, cex =2,y.intersp = 0.7,x.intersp = 0.2,pt.cex = 1)

#check with the FNN
# now actually make the heat map
par(mar=rep(0,4)) # make plot margins 0
heatmap(FNN.M, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=0:1, frame=T, cexRow = 0.5 , cexCol = 0.5)
#
legend("topleft", c("Neighbors SIRT6","SIRT6"), fill=c("red","blue"),  
       bty="n", inset=0, xpd=T,  border=F)

#
plot.new()

######################################################Graphics Using the Network and Node Level Metrics
#########################################
#Graphics Using the Network and Node Level Metrics
setwd("C:/InteractomeProteome/Network Analysis/SensorID_NetworkAnalysis")
#NetworkLevel Metrics
Net.Wide.DS<-read.csv("./NetworkMetrics_WideDS.csv",header = TRUE) #Network Level Metrics Wide format
#Class as factor for discrete variables
Net.Wide.DS$Sensor<-as.factor(Net.Wide.DS$Sensor)
Net.Wide.DS$Time.Point<-as.factor(Net.Wide.DS$Time.Point)
Net.Wide.DS$Time<-as.factor(Net.Wide.DS$Time)
str(Net.Wide.DS)
#Generate Ordering
Net.Wide.DS<-  Net.Wide.DS %>%
  mutate(Order2 = case_when(
    (Sensor== "Ku80") ~ 2,
    (Sensor== "Sirt6") ~ 1,
    (Sensor== "Mre11") ~ 3,
  ))

Net.Wide.DS$Time.Point<-reorder(Net.Wide.DS$Time.Point,Net.Wide.DS$Order)
Net.Wide.DS<-Net.Wide.DS[order(Net.Wide.DS$Order),]
Net.Wide.DS$Sensor<-reorder(Net.Wide.DS$Sensor,Net.Wide.DS$Order)
Net.Wide.DS<-Net.Wide.DS[order(Net.Wide.DS$Order2),]
Net.Wide.DS

Net.Long.DS<-read.csv("./NetworkMetrics_LongDS.csv",header = TRUE)#Network level Metrics Long format
#Class as factor for discrete variables
Net.Long.DS$Sensor<-as.factor(Net.Long.DS$Sensor)
Net.Long.DS$Time.Point<-as.factor(Net.Long.DS$Time.Point)
Net.Long.DS$Time<-as.factor(Net.Long.DS$Time)
str(Net.Long.DS)
#Generate Ordering
Net.Long.DS<-  Net.Long.DS %>%
  mutate(Order2 = case_when(
    (Sensor== "Ku80") ~ 3,
    (Sensor== "Sirt6") ~ 2,
    (Sensor== "Mre11") ~ 4,
    (Sensor== "All") ~ 1
  ))

Net.Long.DS$Time.Point<-reorder(Net.Long.DS$Time.Point,Net.Long.DS$Order)
Net.Long.DS<-Net.Long.DS[order(Net.Long.DS$Order),]
Net.Long.DS$Sensor<-reorder(Net.Long.DS$Sensor,Net.Long.DS$Order)
Net.Long.DS<-Net.Long.DS[order(Net.Long.DS$Order2),]
Net.Long.DS



#Plotings

#Plot All Measures (Do not work nice)
ggplot(Net.Long.DS, aes(x = Time.Point, y = Value, color = Sensor, group=Sensor)) + 
  geom_line(alpha=0.01,stat="identity",linewidth=0.01)+
  geom_point(size=0.0001)+
  stat_smooth(method = "loess",geom = "smooth",se=FALSE)+
  facet_grid(Sensor~.~Measure)+
  labs(x="Time Point", y="No. of Interactors")+
  scale_fill_manual(aesthetics = "color",values=c("purple","seagreen1","steelblue1","tomato"))+
  scale_x_discrete("Time Points",labels=c("nonIR5","IR5","IR30","IR2hrs","IR8hrs","IR24h","nonIR24h"))+
  theme_classic()+  
  theme(axis.text=element_text(size=14,face = "bold"),
        axis.title=element_text(size=18,face="bold"),
        legend.text=element_text(size=18,face = "bold"),
        legend.title = element_text(size=18,face = "bold"),
        strip.text.x = element_text(size = 18,face = "bold"),
        strip.text.y = element_text(size = 148))+
  ylim(0,10)


str(Net.Wide.DS)

#Plot No. of Nodes
ggplot(Net.Wide.DS, aes(x = Time.Point, y =No.Nodes , fill = Sensor, group=Sensor)) + 
  geom_bar(alpha=0.6,stat="identity",position = "dodge")+
  labs(x="Time Point", y="No. of Connected Nodes")+
  scale_fill_manual(aesthetics = "fill",values=c("purple","seagreen1","steelblue1","tomato"))+
  scale_x_discrete("Time Points",labels=c("FNN","nonIR5","IR5","IR30","IR2hrs","IR8hrs","IR24h","nonIR24h"))+
  theme_classic()+  
  theme(axis.text=element_text(size=14,face = "bold"),
        axis.title=element_text(size=18,face="bold"),
        legend.text=element_text(size=18,face = "bold"),
        legend.title = element_text(size=18,face = "bold"),
        strip.text.x = element_text(size = 18,face = "bold"),
        strip.text.y = element_text(size = 148))
#

subset(Net.Wide.DS,Net.Wide.DS$Sensor=="Sirt6"|Net.Wide.DS$Sensor=="Ku80"|Net.Wide.DS$Sensor=="Mre11")%>%
  ggplot(aes(x = Time.Point, y =No.Nodes , fill = Sensor, group=Sensor)) + 
  geom_bar(alpha=0.6,stat="identity",position = "dodge")+
  labs(x="Time Point", y="No. of Connected Nodes")+
  scale_fill_manual(aesthetics = "fill",values=c("seagreen1","steelblue1","tomato"))+
  scale_x_discrete("Time Points",labels=c("FNN","nonIR5","IR5","IR30","IR2hrs","IR8hrs","IR24h","nonIR24h"))+
  theme_classic()+  
  theme(axis.text=element_text(size=14,face = "bold"),
        axis.title=element_text(size=18,face="bold"),
        legend.text=element_text(size=18,face = "bold"),
        legend.title = element_text(size=18,face = "bold"),
        strip.text.x = element_text(size = 18,face = "bold"),
        strip.text.y = element_text(size = 148))




#Plot No. of Edges
subset(Net.Wide.DS,Net.Wide.DS$Sensor=="Sirt6"|Net.Wide.DS$Sensor=="Ku80"|Net.Wide.DS$Sensor=="Mre11")%>%
  ggplot(aes(x = Time.Point, y =No.Edges , fill = Sensor, group=Sensor)) + 
  geom_bar(alpha=0.6,stat="identity",position = "stack")+
  labs(x="Time Point", y="No. of Edges")+
  scale_fill_manual(aesthetics = "fill",values=c("seagreen1","steelblue1","tomato"))+
  scale_x_discrete("Time Points",labels=c("FNN","nonIR5","IR5","IR30","IR2hrs","IR8hrs","IR24h","nonIR24h"))+
  theme_classic()+  
  theme(axis.text=element_text(size=14,face = "bold"),
        axis.title=element_text(size=18,face="bold"),
        legend.text=element_text(size=18,face = "bold"),
        legend.title = element_text(size=18,face = "bold"),
        strip.text.x = element_text(size = 18,face = "bold"),
        strip.text.y = element_text(size = 148))
#

subset(Net.Wide.DS,Net.Wide.DS$Sensor=="Sirt6"|Net.Wide.DS$Sensor=="Ku80"|Net.Wide.DS$Sensor=="Mre11")%>%
  ggplot(aes(x = Time.Point, y =No.Edges , fill = Sensor, group=Sensor)) + 
  geom_bar(alpha=0.6,stat="identity",position = "dodge")+
  labs(x="Time Point", y="No. of Edges")+
  scale_fill_manual(aesthetics = "fill",values=c("seagreen1","steelblue1","tomato"))+
  scale_x_discrete("Time Points",labels=c("FNN","nonIR5","IR5","IR30","IR2hrs","IR8hrs","IR24h","nonIR24h"))+
  theme_classic()+  
  theme(axis.text=element_text(size=14,face = "bold"),
        axis.title=element_text(size=18,face="bold"),
        legend.text=element_text(size=18,face = "bold"),
        legend.title = element_text(size=18,face = "bold"),
        strip.text.x = element_text(size = 18,face = "bold"),
        strip.text.y = element_text(size = 148))


#Plot Edge_density
ggplot(Net.Wide.DS, aes(x = Time.Point, y =Edge_density , color = Sensor, group=Sensor)) + 
  geom_path(alpha=0.7,stat="identity",linewidth=1.5)+
  geom_point(size=4)+
  labs(x="Time Point", y="Edge Density")+
  scale_fill_manual(aesthetics = "color",values=c("purple","seagreen1","steelblue1","tomato"))+
  scale_x_discrete("Time Points",labels=c("FNN","nonIR5","IR5","IR30","IR2hrs","IR8hrs","IR24h","nonIR24h"))+
  theme_classic()+  
  theme(axis.text=element_text(size=14,face = "bold"),
        axis.title=element_text(size=18,face="bold"),
        legend.text=element_text(size=18,face = "bold"),
        legend.title = element_text(size=18,face = "bold"),
        strip.text.x = element_text(size = 18,face = "bold"),
        strip.text.y = element_text(size = 148))
#

subset(Net.Wide.DS,Net.Wide.DS$Sensor=="Sirt6"|Net.Wide.DS$Sensor=="Ku80"|Net.Wide.DS$Sensor=="Mre11")%>%
  ggplot(aes(x = Time.Point, y =Edge_density , color = Sensor, group=Sensor)) + 
  geom_path(alpha=0.7,stat="identity",linewidth=1.5)+
  geom_point(size=4)+
  labs(x="Time Point", y="Edge Density")+
  scale_fill_manual(aesthetics = "color",values=c("seagreen1","steelblue1","tomato"))+
  scale_x_discrete("Time Points",labels=c("nonIR5","IR5","IR30","IR2hrs","IR8hrs","IR24h","nonIR24h"))+
  theme_classic()+  
  theme(axis.text=element_text(size=14,face = "bold"),
        axis.title=element_text(size=18,face="bold"),
        legend.text=element_text(size=18,face = "bold"),
        legend.title = element_text(size=18,face = "bold"),
        strip.text.x = element_text(size = 18,face = "bold"),
        strip.text.y = element_text(size = 148))

#Plot Transitivity
ggplot(Net.Wide.DS, aes(x = Time.Point, y =Transitivity , color = Sensor, group=Sensor)) + 
  geom_path(alpha=0.7,stat="identity",linewidth=1.5)+
  geom_point(size=4)+
  labs(x="Time Point", y="Transitivity")+
  scale_fill_manual(aesthetics = "color",values=c("purple","seagreen1","steelblue1","tomato"))+
  scale_x_discrete("Time Points",labels=c("FNN","nonIR5","IR5","IR30","IR2hrs","IR8hrs","IR24h","nonIR24h"))+
  theme_classic()+  
  theme(axis.text=element_text(size=14,face = "bold"),
        axis.title=element_text(size=18,face="bold"),
        legend.text=element_text(size=18,face = "bold"),
        legend.title = element_text(size=18,face = "bold"),
        strip.text.x = element_text(size = 18,face = "bold"),
        strip.text.y = element_text(size = 148))
#

subset(Net.Wide.DS,Net.Wide.DS$Sensor=="Sirt6"|Net.Wide.DS$Sensor=="Ku80"|Net.Wide.DS$Sensor=="Mre11")%>%
  ggplot(aes(x = Time.Point, y =Transitivity , color = Sensor, group=Sensor)) + 
  geom_path(alpha=0.7,stat="identity",linewidth=1.5)+
  geom_point(size=4)+
  labs(x="Time Point", y="Transitivity")+
  scale_fill_manual(aesthetics = "color",values=c("seagreen1","steelblue1","tomato"))+
  scale_x_discrete("Time Points",labels=c("nonIR5","IR5","IR30","IR2hrs","IR8hrs","IR24h","nonIR24h"))+
  theme_classic()+  
  theme(axis.text=element_text(size=14,face = "bold"),
        axis.title=element_text(size=18,face="bold"),
        legend.text=element_text(size=18,face = "bold"),
        legend.title = element_text(size=18,face = "bold"),
        strip.text.x = element_text(size = 18,face = "bold"),
        strip.text.y = element_text(size = 148))

#Plot Diameter
ggplot(Net.Wide.DS, aes(x = Time.Point, y =Diameter , color = Sensor, group=Sensor)) + 
  geom_path(alpha=0.7,stat="identity",linewidth=1.5)+
  geom_point(size=4)+
  labs(x="Time Point", y="Diameter")+
  scale_fill_manual(aesthetics = "color",values=c("purple","seagreen1","steelblue1","tomato"))+
  scale_x_discrete("Time Points",labels=c("FNN","nonIR5","IR5","IR30","IR2hrs","IR8hrs","IR24h","nonIR24h"))+
  theme_classic()+  
  theme(axis.text=element_text(size=14,face = "bold"),
        axis.title=element_text(size=18,face="bold"),
        legend.text=element_text(size=18,face = "bold"),
        legend.title = element_text(size=18,face = "bold"),
        strip.text.x = element_text(size = 18,face = "bold"),
        strip.text.y = element_text(size = 148))
#

subset(Net.Wide.DS,Net.Wide.DS$Sensor=="Sirt6"|Net.Wide.DS$Sensor=="Ku80"|Net.Wide.DS$Sensor=="Mre11")%>%
  ggplot(aes(x = Time.Point, y =Diameter , color = Sensor, group=Sensor)) + 
  geom_path(alpha=0.7,stat="identity",linewidth=1.5)+
  geom_point(size=4)+
  labs(x="Time Point", y="Diameter")+
  scale_fill_manual(aesthetics = "color",values=c("seagreen1","steelblue1","tomato"))+
  scale_x_discrete("Time Points",labels=c("nonIR5","IR5","IR30","IR2hrs","IR8hrs","IR24h","nonIR24h"))+
  theme_classic()+  
  theme(axis.text=element_text(size=14,face = "bold"),
        axis.title=element_text(size=18,face="bold"),
        legend.text=element_text(size=18,face = "bold"),
        legend.title = element_text(size=18,face = "bold"),
        strip.text.x = element_text(size = 18,face = "bold"),
        strip.text.y = element_text(size = 148))

#Plot FuzzyCMeans modularity
ggplot(Net.Wide.DS, aes(x = Time.Point, y =Modularity  , color = Sensor, group=Sensor)) + 
  geom_path(alpha=0.7,stat="identity",linewidth=1.5)+
  geom_point(size=4)+
  labs(x="Time Point", y="FuzzyCMeans modularity")+
  scale_fill_manual(aesthetics = "color",values=c("purple","seagreen1","steelblue1","tomato"))+
  scale_x_discrete("Time Points",labels=c("FNN","nonIR5","IR5","IR30","IR2hrs","IR8hrs","IR24h","nonIR24h"))+
  theme_classic()+  
  theme(axis.text=element_text(size=14,face = "bold"),
        axis.title=element_text(size=18,face="bold"),
        legend.text=element_text(size=18,face = "bold"),
        legend.title = element_text(size=18,face = "bold"),
        strip.text.x = element_text(size = 18,face = "bold"),
        strip.text.y = element_text(size = 148))
#

#Plot S6.FuzzyCMeans modularity
ggplot(Net.Wide.DS, aes(x = Time.Point, y =S6.FuzzyCMeans.modularity  , color = Sensor, group=Sensor)) + 
  geom_path(alpha=0.7,stat="identity",linewidth=1.5)+
  geom_point(size=4)+
  labs(x="Time Point", y="S6 FuzzyCMeans modularity")+
  scale_fill_manual(aesthetics = "color",values=c("purple","seagreen1","steelblue1","tomato"))+
  scale_x_discrete("Time Points",labels=c("FNN","nonIR5","IR5","IR30","IR2hrs","IR8hrs","IR24h","nonIR24h"))+
  theme_classic()+  
  theme(axis.text=element_text(size=14,face = "bold"),
        axis.title=element_text(size=18,face="bold"),
        legend.text=element_text(size=18,face = "bold"),
        legend.title = element_text(size=18,face = "bold"),
        strip.text.x = element_text(size = 18,face = "bold"),
        strip.text.y = element_text(size = 148))
#

subset(Net.Wide.DS,Net.Wide.DS$Sensor=="Sirt6"|Net.Wide.DS$Sensor=="Ku80"|Net.Wide.DS$Sensor=="Mre11")%>%
  ggplot(aes(x = Time.Point, y =S6.FuzzyCMeans.modularity  , color = Sensor, group=Sensor)) + 
  geom_path(alpha=0.7,stat="identity",linewidth=1.5)+
  geom_point(size=4)+
  labs(x="Time Point", y="S6 FuzzyCMeans modularity")+
  scale_fill_manual(aesthetics = "color",values=c("seagreen1","steelblue1","tomato"))+
  scale_x_discrete("Time Points",labels=c("nonIR5","IR5","IR30","IR2hrs","IR8hrs","IR24h","nonIR24h"))+
  theme_classic()+  
  theme(axis.text=element_text(size=14,face = "bold"),
        axis.title=element_text(size=18,face="bold"),
        legend.text=element_text(size=18,face = "bold"),
        legend.title = element_text(size=18,face = "bold"),
        strip.text.x = element_text(size = 18,face = "bold"),
        strip.text.y = element_text(size = 148))

#Plot Ku.FuzzyCMeans modularity
ggplot(Net.Wide.DS, aes(x = Time.Point, y =Ku.FuzzyCMeans.modularity  , color = Sensor, group=Sensor)) + 
  geom_path(alpha=0.7,stat="identity",linewidth=1.5)+
  geom_point(size=4)+
  labs(x="Time Point", y="Ku FuzzyCMeans modularity")+
  scale_fill_manual(aesthetics = "color",values=c("purple","seagreen1","steelblue1","tomato"))+
  scale_x_discrete("Time Points",labels=c("FNN","nonIR5","IR5","IR30","IR2hrs","IR8hrs","IR24h","nonIR24h"))+
  theme_classic()+  
  theme(axis.text=element_text(size=14,face = "bold"),
        axis.title=element_text(size=18,face="bold"),
        legend.text=element_text(size=18,face = "bold"),
        legend.title = element_text(size=18,face = "bold"),
        strip.text.x = element_text(size = 18,face = "bold"),
        strip.text.y = element_text(size = 148))
#

subset(Net.Wide.DS,Net.Wide.DS$Sensor=="Sirt6"|Net.Wide.DS$Sensor=="Ku80"|Net.Wide.DS$Sensor=="Mre11")%>%
  ggplot(aes(x = Time.Point, y =Ku.FuzzyCMeans.modularity  , color = Sensor, group=Sensor)) + 
  geom_path(alpha=0.7,stat="identity",linewidth=1.5)+
  geom_point(size=4)+
  labs(x="Time Point", y="Ku FuzzyCMeans modularity")+
  scale_fill_manual(aesthetics = "color",values=c("seagreen1","steelblue1","tomato"))+
  scale_x_discrete("Time Points",labels=c("nonIR5","IR5","IR30","IR2hrs","IR8hrs","IR24h","nonIR24h"))+
  theme_classic()+  
  theme(axis.text=element_text(size=14,face = "bold"),
        axis.title=element_text(size=18,face="bold"),
        legend.text=element_text(size=18,face = "bold"),
        legend.title = element_text(size=18,face = "bold"),
        strip.text.x = element_text(size = 18,face = "bold"),
        strip.text.y = element_text(size = 148))

#Plot Mre.FuzzyCMeans modularity
ggplot(Net.Wide.DS, aes(x = Time.Point, y =Mre.FuzzyCMeans.modularity  , color = Sensor, group=Sensor)) + 
  geom_path(alpha=0.7,stat="identity",linewidth=1.5)+
  geom_point(size=4)+
  labs(x="Time Point", y="Mre FuzzyCMeans modularity")+
  scale_fill_manual(aesthetics = "color",values=c("purple","seagreen1","steelblue1","tomato"))+
  scale_x_discrete("Time Points",labels=c("FNN","nonIR5","IR5","IR30","IR2hrs","IR8hrs","IR24h","nonIR24h"))+
  theme_classic()+  
  theme(axis.text=element_text(size=14,face = "bold"),
        axis.title=element_text(size=18,face="bold"),
        legend.text=element_text(size=18,face = "bold"),
        legend.title = element_text(size=18,face = "bold"),
        strip.text.x = element_text(size = 18,face = "bold"),
        strip.text.y = element_text(size = 148))
#

subset(Net.Wide.DS,Net.Wide.DS$Sensor=="Sirt6"|Net.Wide.DS$Sensor=="Ku80"|Net.Wide.DS$Sensor=="Mre11")%>%
  ggplot(aes(x = Time.Point, y =Mre.FuzzyCMeans.modularity  , color = Sensor, group=Sensor)) + 
  geom_path(alpha=0.7,stat="identity",linewidth=1.5)+
  geom_point(size=4)+
  labs(x="Time Point", y="Mre FuzzyCMeans modularity")+
  scale_fill_manual(aesthetics = "color",values=c("seagreen1","steelblue1","tomato"))+
  scale_x_discrete("Time Points",labels=c("nonIR5","IR5","IR30","IR2hrs","IR8hrs","IR24h","nonIR24h"))+
  theme_classic()+  
  theme(axis.text=element_text(size=14,face = "bold"),
        axis.title=element_text(size=18,face="bold"),
        legend.text=element_text(size=18,face = "bold"),
        legend.title = element_text(size=18,face = "bold"),
        strip.text.x = element_text(size = 18,face = "bold"),
        strip.text.y = element_text(size = 148))

#Plot Assortativity. Annotation
ggplot(Net.Wide.DS, aes(x = Time.Point, y =Assortativity..Annotation  , color = Sensor, group=Sensor)) + 
  geom_path(alpha=0.7,stat="identity",linewidth=1.5)+
  geom_point(size=4)+
  labs(x="Time Point", y="Assortativity Annotation")+
  scale_fill_manual(aesthetics = "color",values=c("purple","seagreen1","steelblue1","tomato"))+
  scale_x_discrete("Time Points",labels=c("FNN","nonIR5","IR5","IR30","IR2hrs","IR8hrs","IR24h","nonIR24h"))+
  theme_classic()+  
  theme(axis.text=element_text(size=14,face = "bold"),
        axis.title=element_text(size=18,face="bold"),
        legend.text=element_text(size=18,face = "bold"),
        legend.title = element_text(size=18,face = "bold"),
        strip.text.x = element_text(size = 18,face = "bold"),
        strip.text.y = element_text(size = 148))
#

subset(Net.Wide.DS,Net.Wide.DS$Sensor=="Sirt6"|Net.Wide.DS$Sensor=="Ku80"|Net.Wide.DS$Sensor=="Mre11")%>%
  ggplot(aes(x = Time.Point, y =Assortativity..Annotation  , color = Sensor, group=Sensor)) + 
  geom_path(alpha=0.7,stat="identity",linewidth=1.5)+
  geom_point(size=4)+
  labs(x="Time Point", y="Assortativity Annotation")+
  scale_fill_manual(aesthetics = "color",values=c("seagreen1","steelblue1","tomato"))+
  scale_x_discrete("Time Points",labels=c("nonIR5","IR5","IR30","IR2hrs","IR8hrs","IR24h","nonIR24h"))+
  theme_classic()+  
  theme(axis.text=element_text(size=14,face = "bold"),
        axis.title=element_text(size=18,face="bold"),
        legend.text=element_text(size=18,face = "bold"),
        legend.title = element_text(size=18,face = "bold"),
        strip.text.x = element_text(size = 18,face = "bold"),
        strip.text.y = element_text(size = 148))

#Plot Assortativity FuzzyCmeans
ggplot(Net.Wide.DS, aes(x = Time.Point, y =Assortativity  , color = Sensor, group=Sensor)) + 
  geom_path(alpha=0.7,stat="identity",linewidth=1.5)+
  geom_point(size=4)+
  labs(x="Time Point", y="Assortativity FuzzyCmeans")+
  scale_fill_manual(aesthetics = "color",values=c("purple","seagreen1","steelblue1","tomato"))+
  scale_x_discrete("Time Points",labels=c("FNN","nonIR5","IR5","IR30","IR2hrs","IR8hrs","IR24h","nonIR24h"))+
  theme_classic()+  
  theme(axis.text=element_text(size=14,face = "bold"),
        axis.title=element_text(size=18,face="bold"),
        legend.text=element_text(size=18,face = "bold"),
        legend.title = element_text(size=18,face = "bold"),
        strip.text.x = element_text(size = 18,face = "bold"),
        strip.text.y = element_text(size = 148))


#Plot Assortativity..S6FuzzyCmeans
ggplot(Net.Wide.DS, aes(x = Time.Point, y =Assortativity..S6FuzzyCmeans  , color = Sensor, group=Sensor)) + 
  geom_path(alpha=0.7,stat="identity",linewidth=1.5)+
  geom_point(size=4)+
  labs(x="Time Point", y="Assortativity S6FuzzyCmeans")+
  scale_fill_manual(aesthetics = "color",values=c("purple","seagreen1","steelblue1","tomato"))+
  scale_x_discrete("Time Points",labels=c("FNN","nonIR5","IR5","IR30","IR2hrs","IR8hrs","IR24h","nonIR24h"))+
  theme_classic()+  
  theme(axis.text=element_text(size=14,face = "bold"),
        axis.title=element_text(size=18,face="bold"),
        legend.text=element_text(size=18,face = "bold"),
        legend.title = element_text(size=18,face = "bold"),
        strip.text.x = element_text(size = 18,face = "bold"),
        strip.text.y = element_text(size = 148))
#

subset(Net.Wide.DS,Net.Wide.DS$Sensor=="Sirt6"|Net.Wide.DS$Sensor=="Ku80"|Net.Wide.DS$Sensor=="Mre11")%>%
  ggplot(aes(x = Time.Point, y =Assortativity..S6FuzzyCmeans  , color = Sensor, group=Sensor)) + 
  geom_path(alpha=0.7,stat="identity",linewidth=1.5)+
  geom_point(size=4)+
  labs(x="Time Point", y="Assortativity S6FuzzyCmeansn")+
  scale_fill_manual(aesthetics = "color",values=c("seagreen1","steelblue1","tomato"))+
  scale_x_discrete("Time Points",labels=c("nonIR5","IR5","IR30","IR2hrs","IR8hrs","IR24h","nonIR24h"))+
  theme_classic()+  
  theme(axis.text=element_text(size=14,face = "bold"),
        axis.title=element_text(size=18,face="bold"),
        legend.text=element_text(size=18,face = "bold"),
        legend.title = element_text(size=18,face = "bold"),
        strip.text.x = element_text(size = 18,face = "bold"),
        strip.text.y = element_text(size = 148))

#Plot Assortativity..KuFuzzyCmeans
ggplot(Net.Wide.DS, aes(x = Time.Point, y =Assortativity..KuFuzzyCmeans  , color = Sensor, group=Sensor)) + 
  geom_path(alpha=0.7,stat="identity",linewidth=1.5)+
  geom_point(size=4)+
  labs(x="Time Point", y="Assortativity KuFuzzyCmeans")+
  scale_fill_manual(aesthetics = "color",values=c("purple","seagreen1","steelblue1","tomato"))+
  scale_x_discrete("Time Points",labels=c("FNN","nonIR5","IR5","IR30","IR2hrs","IR8hrs","IR24h","nonIR24h"))+
  theme_classic()+  
  theme(axis.text=element_text(size=14,face = "bold"),
        axis.title=element_text(size=18,face="bold"),
        legend.text=element_text(size=18,face = "bold"),
        legend.title = element_text(size=18,face = "bold"),
        strip.text.x = element_text(size = 18,face = "bold"),
        strip.text.y = element_text(size = 148))
#

subset(Net.Wide.DS,Net.Wide.DS$Sensor=="Sirt6"|Net.Wide.DS$Sensor=="Ku80"|Net.Wide.DS$Sensor=="Mre11")%>%
  ggplot(aes(x = Time.Point, y =Assortativity..KuFuzzyCmeans  , color = Sensor, group=Sensor)) + 
  geom_path(alpha=0.7,stat="identity",linewidth=1.5)+
  geom_point(size=4)+
  labs(x="Time Point", y="Assortativity..KuFuzzyCmeans")+
  scale_fill_manual(aesthetics = "color",values=c("seagreen1","steelblue1","tomato"))+
  scale_x_discrete("Time Points",labels=c("nonIR5","IR5","IR30","IR2hrs","IR8hrs","IR24h","nonIR24h"))+
  theme_classic()+  
  theme(axis.text=element_text(size=14,face = "bold"),
        axis.title=element_text(size=18,face="bold"),
        legend.text=element_text(size=18,face = "bold"),
        legend.title = element_text(size=18,face = "bold"),
        strip.text.x = element_text(size = 18,face = "bold"),
        strip.text.y = element_text(size = 148))

#Plot Assortativity..MreFuzzyCmeans
ggplot(Net.Wide.DS, aes(x = Time.Point, y =Assortativity..MreFuzzyCmeans  , color = Sensor, group=Sensor)) + 
  geom_path(alpha=0.7,stat="identity",linewidth=1.5)+
  geom_point(size=4)+
  labs(x="Time Point", y="Assortativity MreFuzzyCmeans")+
  scale_fill_manual(aesthetics = "color",values=c("purple","seagreen1","steelblue1","tomato"))+
  scale_x_discrete("Time Points",labels=c("FNN","nonIR5","IR5","IR30","IR2hrs","IR8hrs","IR24h","nonIR24h"))+
  theme_classic()+  
  theme(axis.text=element_text(size=14,face = "bold"),
        axis.title=element_text(size=18,face="bold"),
        legend.text=element_text(size=18,face = "bold"),
        legend.title = element_text(size=18,face = "bold"),
        strip.text.x = element_text(size = 18,face = "bold"),
        strip.text.y = element_text(size = 148))
#

subset(Net.Wide.DS,Net.Wide.DS$Sensor=="Sirt6"|Net.Wide.DS$Sensor=="Ku80"|Net.Wide.DS$Sensor=="Mre11")%>%
  ggplot(aes(x = Time.Point, y =Assortativity..MreFuzzyCmeans  , color = Sensor, group=Sensor)) + 
  geom_path(alpha=0.7,stat="identity",linewidth=1.5)+
  geom_point(size=4)+
  labs(x="Time Point", y="Assortativity..MreFuzzyCmeans")+
  scale_fill_manual(aesthetics = "color",values=c("seagreen1","steelblue1","tomato"))+
  scale_x_discrete("Time Points",labels=c("nonIR5","IR5","IR30","IR2hrs","IR8hrs","IR24h","nonIR24h"))+
  theme_classic()+  
  theme(axis.text=element_text(size=14,face = "bold"),
        axis.title=element_text(size=18,face="bold"),
        legend.text=element_text(size=18,face = "bold"),
        legend.title = element_text(size=18,face = "bold"),
        strip.text.x = element_text(size = 18,face = "bold"),
        strip.text.y = element_text(size = 148))
#################
########################################################Node Metrics Level HeatMap
########################
#Node Level Metrics HeatMaps
library(ComplexHeatmap)
library(circlize)
library(dplyr)
#Read the composed matrix
NodeMat<-read.csv("./NodeMetrics_Matrix.csv",row.names = "NAME")
str(NodeMat)
row.names(NodeMat)
#Annotation Matrix
DF<-NodeMat[,2:13]
DF$Annotation<-as.factor(DF$Annotation)
GO_Anno<-as.matrix(subset(DF,select = Annotation))
ha_anno<-rowAnnotation(Annotation=as.factor(GO_Anno))#,col=Col_Anno)#annotation_legend_param=)#,labels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22"))
ha_anno

#Column Annotations
colAnno<-read.csv("./ColAnnoDF.csv",row.names = 1)
Col_Anno<-HeatmapAnnotation(df=colAnno, col = list(Treatment=c("IR"="red","nonIR"="blue"),
                                                   Sensor=c("Ku80"="seagreen1","Mre11"="steelblue1","Sirt6"="tomato")))  

#Generate the Metrics specific matrices
EigenMat<-NodeMat[,56:76]
head(EigenMat)
ClosenessMat<-NodeMat[,77:97]
head(ClosenessMat)
BetweenessMat<-NodeMat[,98:118]
head(BetweenessMat)
HubsMat<-NodeMat[,119:139]
head(HubsMat)
AuthorityMat<-NodeMat[,140:160]
head(AuthorityMat)
TransitivityMat<-NodeMat[,161:181]
head(TransitivityMat)

#Degree HeatMap
DegreeMat<-NodeMat[,34:54]
row_sums <- rowSums(DegreeMat)
row_sums
DegreeMat_filtered <- DegreeMat[row_sums != 0, ]
DegreeMat_filtered<-as.matrix(DegreeMat_filtered)
head(DegreeMat)
DegreeMat<-as.matrix(DegreeMat)
DegreeMat<-t(scale(t(DegreeMat),center = TRUE,scale = TRUE))
str(DegreeMat)
any(is.na(DegreeMat))
is.na(DegreeMat_filtered)%>%table()
is.infinite(DegreeMat_filtered)%>%table()
dim(DegreeMat)
row.names(DegreeMat)
mean(DegreeMat)#[1] [1] -3.747166e-18
median(DegreeMat)#[1] -0.3267613
range(DegreeMat) #1] -0.4287193  8.5510674
quantile(DegreeMat,c(0,0.25,0.5,0.75,1))
#              0%        25%        50%        75%       100% 
#-0.4287193 -0.3874102 -0.3267613 -0.2418331  8.5510674 

breaks <- c(seq(min(DegreeMat_filtered),max(DegreeMat_filtered),length=11))
#Color function for HeatMaps
col_fun1<-colorRamp2(c(min(DegreeMat_filtered),max(DegreeMat_filtered)),c("#F7F7F7","#B2182B"))
col_fun1(seq(breaks)) 

#Generate a dendogram for keeping the proteins in the sample place
DD<-dist(DegreeMat_filtered,method = "euclidean")
#DD<-hclust(DD, method = "average" )
plot(DD)

#Annotation Matrix
DF<-NodeMat[,1:12]
DF<-subset(DF,rownames(DF)%in%rownames(DegreeMat_filtered))
DF$Annotation<-as.factor(DF$Annotation)
GO_Anno<-as.matrix(subset(DF,select = Annotation))
ha_anno<-rowAnnotation(Annotation=as.factor(GO_Anno))#,col=Col_Anno)#annotation_legend_param=)#,labels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22"))
ha_anno

require(ComplexHeatmap)
#Plot heat map

Heatmap(DegreeMat_filtered, name = "Node Degree" ,  col = col_fun1 , cluster_rows = T , 
        cluster_columns = FALSE , clustering_distance_rows = DD , clustering_distance_columns = "pearson" ,
        clustering_method_rows = "average" , clustering_method_columns = "average" , column_names_centered = TRUE , 
        column_names_gp = gpar(fontsize=15) , row_names_gp = gpar(fontsize=3) ,top_annotation = Col_Anno, show_column_names = T,
        row_gap = unit(1, "mm") , 
        left_annotation = ha_anno)

#EigenCentrality HeatMap
EigenMat<-NodeMat[,55:75]
head(EigenMat)
EigenMat<-as.matrix(EigenMat)
#EigenMat<-scale(EigenMat,center = TRUE,scale = TRUE)
str(EigenMat)
any(is.na(EigenMat))
is.na(EigenMat)%>%table()
is.infinite(EigenMat)%>%table()
dim(EigenMat)
#
row_sums <- rowSums(EigenMat)
row_sums
EigenMat_filtered <- EigenMat[row_sums != 0, ]

mean(EigenMat)#[1] 1] 0.03962607
median(EigenMat)#[1] [1] 0
range(EigenMat) #1] [1] 0.000000 1.093504
quantile(EigenMat,c(0,0.25,0.5,0.75,1))
#              0%      25%      50%      75%     100% 
#0.000000 0.000000 0.000000 0.000000 1.093504 
breaks <- c(seq(min(EigenMat_filtered),max(EigenMat_filtered),length=11))
#Color function for HeatMaps
col_fun1<-colorRamp2(c(min(EigenMat_filtered),max(EigenMat_filtered)),c("#F7F7F7","#B2182B"))
col_fun1(seq(breaks)) 

#Generate a dendogram for keeping the proteins in the sample place
DD<-dist(EigenMat_filtered,method = "euclidean")
#DD<-hclust(DD, method = "average" )
plot(DD)
# Row annotations
DF<-NodeMat[,1:12]
DF<-subset(DF,rownames(DF)%in%rownames(EigenMat_filtered))
DF$Annotation<-as.factor(DF$Annotation)
#
GO_Anno<-as.matrix(subset(DF,select = Annotation))
ha_anno<-rowAnnotation(Annotation=as.factor(GO_Anno))#,col=Col_Anno)#annotation_legend_param=)#,labels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22"))
ha_anno


#Plot the Heat map
Heatmap(EigenMat_filtered, name = "Eigen Cent" ,  col = col_fun1 , cluster_rows = T , 
        cluster_columns = FALSE , clustering_distance_rows = DD , clustering_distance_columns = "pearson" ,
        clustering_method_rows = "average" , clustering_method_columns = "average" , column_names_centered = TRUE , 
        column_names_gp = gpar(fontsize=15) , row_names_gp = gpar(fontsize=3) ,top_annotation = Col_Anno, show_column_names = T,
        show_row_names = TRUE,
        left_annotation = ha_anno)

#ClosenessCentrality HeatMap
ClosenessMat<-NodeMat[,76:96]
head(ClosenessMat)
ClosenessMat<-as.matrix(ClosenessMat)
#ClosenessMat<-scale(ClosenessMat,center = TRUE,scale = TRUE)
str(ClosenessMat)
any(is.na(ClosenessMat))
is.na(ClosenessMat)%>%table()
is.infinite(ClosenessMat)%>%table()
dim(ClosenessMat)

#
row_sums <- rowSums(ClosenessMat)
row_sums
ClosenessMat_filtered <- ClosenessMat[row_sums != 0, ]

mean(ClosenessMat)#[1]  0.1566645
median(ClosenessMat)#[1] 0
range(ClosenessMat) #1] 0.000000 2.375297
quantile(ClosenessMat,c(0,0.25,0.5,0.75,1))
#      0%      25%      50%      75%     100% 
#   0.000000 0.000000 0.000000 0.000000 2.375297 

breaks <- c(seq(min(ClosenessMat_filtered ),max(ClosenessMat_filtered ),length=11))
#Color function for HeatMaps
col_fun1<-colorRamp2(c(min(ClosenessMat_filtered ),max(ClosenessMat_filtered )),c("#F7F7F7","#B2182B"))
col_fun1(seq(breaks)) 

#Generate a dendogram for keeping the proteins in the sample place
DD<-dist(ClosenessMat_filtered ,method = "euclidean")
#DD<-hclust(DD, method = "average" )
plot(DD)

#Row annotations
DF<-NodeMat[,1:12]
DF<-subset(DF,rownames(DF)%in%rownames(ClosenessMat_filtered))

DF$Annotation<-as.factor(DF$Annotation)
GO_Anno<-as.matrix(subset(DF,select = Annotation))
ha_anno<-rowAnnotation(Annotation=as.factor(GO_Anno))#,col=Col_Anno)#annotation_legend_param=)#,labels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22"))
ha_anno

#Plot the Heat map
Heatmap(ClosenessMat_filtered , name = "Closeness" ,  col = col_fun1 , cluster_rows = T , 
        cluster_columns = FALSE , clustering_distance_rows = DD , clustering_distance_columns = "pearson" ,
        clustering_method_rows = "average" , clustering_method_columns = "average" , column_names_centered = TRUE , 
        column_names_gp = gpar(fontsize=15) , row_names_gp = gpar(fontsize=3) ,top_annotation = Col_Anno, show_column_names = T,
        show_row_names = TRUE,
        left_annotation = ha_anno)

##BetweenessCentrality HeatMap
BetweenessMat<-NodeMat[,97:117]
head(BetweenessMat)
BetweenessMat<-as.matrix(BetweenessMat)
#BetweenessMat<-scale(BetweenessMat,center = TRUE,scale = TRUE)
str(BetweenessMat)
any(is.na(BetweenessMat))
is.na(BetweenessMat)%>%table()
is.infinite(BetweenessMat)%>%table()
dim(BetweenessMat)

row_sums <- rowSums(BetweenessMat)
row_sums
BetweenessMat_filtered <- BetweenessMat[row_sums != 0, ]


mean(BetweenessMat)#[1] 0.002341562
median(BetweenessMat)#[1]0
range(BetweenessMat) #1] 0.0000000 0.4771895
quantile(BetweenessMat,c(0,0.25,0.5,0.75,1))
#      0%       25%       50%       75%      100% 
#  0.0000000 0.0000000 0.0000000 0.0000000 0.4771895 

breaks <- c(seq(min(BetweenessMat_filtered),0.002341562,length=11))
#Color function for HeatMaps
col_fun1<-colorRamp2(c(min(BetweenessMat),0.002341562),c("#F7F7F7","#B2182B"))
col_fun1(seq(breaks)) 

#Generate a dendogram for keeping the proteins in the sample place
DD<-dist(BetweenessMat_filtered,method = "euclidean")
#DD<-hclust(DD, method = "average" )
plot(DD)

#Row annotations
DF<-NodeMat[,1:12]
DF<-subset(DF,rownames(DF)%in%rownames(BetweenessMat_filtered))
DF$Annotation<-as.factor(DF$Annotation)
GO_Anno<-as.matrix(subset(DF,select = Annotation))
ha_anno<-rowAnnotation(Annotation=as.factor(GO_Anno))#,col=Col_Anno)#annotation_legend_param=)#,labels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22"))
ha_anno

#Plot the heatmap
Heatmap(BetweenessMat_filtered, name = "Betweeness Cent" ,  col = col_fun1 , cluster_rows = T , 
        cluster_columns = FALSE , clustering_distance_rows = DD , clustering_distance_columns = "pearson" ,
        clustering_method_rows = "average" , clustering_method_columns = "average" , column_names_centered = TRUE , 
        column_names_gp = gpar(fontsize=15) , row_names_gp = gpar(fontsize=3) ,top_annotation = Col_Anno, show_column_names = T,
        show_row_names = TRUE,
        left_annotation = ha_anno)

##HubCentrality HeatMap
HubsMat<-NodeMat[,118:138]
head(HubsMat)
HubsMat<-as.matrix(HubsMat)
#HubsMat<-t(scale(t(HubsMat),center = TRUE,scale = TRUE))
str(HubsMat)
any(is.na(HubsMat))
is.na(HubsMat)%>%table()
is.infinite(HubsMat)%>%table()
dim(HubsMat)

row_sums <- rowSums(HubsMat)
row_sums
HubsMat_filtered <- HubsMat[row_sums != 0, ]

mean(HubsMat)#[1]  0.03170686
median(HubsMat)#[1] 0
range(HubsMat) #1] 0 1
quantile(HubsMat,c(0,0.25,0.5,0.75,1))

breaks <- c(seq(min(HubsMat_filtered),max(HubsMat_filtered),length=11))
#Color function for HeatMaps
col_fun1<-colorRamp2(c(min(HubsMat_filtered),max(HubsMat_filtered)),c("#F7F7F7","#B2182B"))
col_fun1(seq(breaks)) 

#Generate a dendogram for keeping the proteins in the sample place
DD<-dist(HubsMat_filtered,method = "euclidean")
#DD<-hclust(DD, method = "average" )
plot(DD)

#Row annotations
DF<-NodeMat[,1:12]
DF<-subset(DF,rownames(DF)%in%rownames(HubsMat_filtered))
DF$Annotation<-as.factor(DF$Annotation)
GO_Anno<-as.matrix(subset(DF,select = Annotation))
ha_anno<-rowAnnotation(Annotation=as.factor(GO_Anno))#,col=Col_Anno)#annotation_legend_param=)#,labels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22"))
ha_anno

#Plot the Heatmap

Heatmap(HubsMat_filtered, name = "Hub Degree" ,  col = col_fun1 , cluster_rows = T , 
        cluster_columns = FALSE , clustering_distance_rows = DD , clustering_distance_columns = "pearson" ,
        clustering_method_rows = "average" , clustering_method_columns = "average" , column_names_centered = TRUE , 
        column_names_gp = gpar(fontsize=15) , row_names_gp = gpar(fontsize=3) ,top_annotation = Col_Anno, show_column_names = T,
        show_row_names = TRUE,
        left_annotation = ha_anno)

##AuthorityCentrality HeatMap
AuthorityMat<-NodeMat[,139:159]
head(AuthorityMat)
AuthorityMat<-as.matrix(AuthorityMat)
#AuthorityMat<-t(scale(t(AuthorityMat),center = TRUE,scale = TRUE))
str(AuthorityMat)
any(is.na(AuthorityMat))
is.na(AuthorityMat)%>%table()
is.infinite(AuthorityMat)%>%table()
dim(AuthorityMat)

mean(AuthorityMat)#[1] 0.03170686
median(AuthorityMat)#[1]0
range(AuthorityMat) #1] 0 1
quantile(AuthorityMat,c(0,0.25,0.5,0.75,1))
#      0%  25%  50%  75% 100% 
#     0    0    0    0    1 

row_sums <- rowSums(AuthorityMat)
row_sums
AuthorityMat_filtered <- AuthorityMat[row_sums != 0, ]

breaks <- c(seq(min(AuthorityMat_filtered),max(AuthorityMat_filtered),length=11))
#Color function for HeatMaps
col_fun1<-colorRamp2(c(min(AuthorityMat_filtered),max(AuthorityMat_filtered)),c("#F7F7F7","#B2182B"))
col_fun1(seq(breaks)) 

#Generate a dendogram for keeping the proteins in the sample place
DD<-dist(AuthorityMat_filtered,method = "euclidean")
#DD<-hclust(DD, method = "average" )
plot(DD)

#Row annotations
DF<-NodeMat[,1:12]
DF<-subset(DF,rownames(DF)%in%rownames(AuthorityMat_filtered))
DF$Annotation<-as.factor(DF$Annotation)
GO_Anno<-as.matrix(subset(DF,select = Annotation))
ha_anno<-rowAnnotation(Annotation=as.factor(GO_Anno))#,col=Col_Anno)#annotation_legend_param=)#,labels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22"))
ha_anno

#Plot the heat map
Heatmap(AuthorityMat_filtered, name = "Authority degree" ,  col = col_fun1 , cluster_rows = T , 
        cluster_columns = FALSE , clustering_distance_rows = DD , clustering_distance_columns = "pearson" ,
        clustering_method_rows = "average" , clustering_method_columns = "average" , column_names_centered = TRUE , 
        column_names_gp = gpar(fontsize=15) , row_names_gp = gpar(fontsize=3) ,top_annotation = Col_Anno, show_column_names = T,
        show_row_names = TRUE,
        left_annotation = ha_anno)

##Transitivity Centrality HeatMap
TransitivityMat<-NodeMat[,160:180]
head(TransitivityMat)
TransitivityMat<-as.matrix(TransitivityMat)
#TransitivityMat<-t(scale(t(TransitivityMat),center = TRUE,scale = TRUE))
str(TransitivityMat)
any(is.na(TransitivityMat))
is.na(TransitivityMat)%>%table()
is.infinite(TransitivityMat)%>%table()
dim(TransitivityMat)

row_sums <- rowSums(TransitivityMat)
row_sums
TransitivityMat_filtered <- TransitivityMat[row_sums != 0, ]

mean(TransitivityMat)#[1]  0.09529262
median(TransitivityMat)#[1] 0
range(TransitivityMat) #1] 0 1
quantile(TransitivityMat,c(0,0.25,0.5,0.75,1))
#      0%  25%  50%  75% 100% 
#      0    0    0    0    1 

breaks <- c(seq(min(TransitivityMat_filtered),max(TransitivityMat_filtered),length=11))
#Color function for HeatMaps
col_fun1<-colorRamp2(c(min(TransitivityMat_filtered),max(TransitivityMat_filtered)),c("#F7F7F7","#B2182B"))
col_fun1(seq(breaks)) 

#Generate a dendogram for keeping the proteins in the sample place
DD<-dist(TransitivityMat_filtered,method = "euclidean")
#DD<-hclust(DD, method = "average" )
plot(DD)

#Row annotations
DF<-NodeMat[,1:12]
DF<-subset(DF,rownames(DF)%in%rownames(TransitivityMat_filtered))
DF$Annotation<-as.factor(DF$Annotation)
GO_Anno<-as.matrix(subset(DF,select = Annotation))
ha_anno<-rowAnnotation(Annotation=as.factor(GO_Anno))#,col=Col_Anno)#annotation_legend_param=)#,labels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22"))
ha_anno

#Plot the Heatmap
Heatmap(TransitivityMat_filtered, name = "Transitivity" ,  col = col_fun1 , cluster_rows = T , 
        cluster_columns = FALSE , clustering_distance_rows = DD , clustering_distance_columns = "pearson" ,
        clustering_method_rows = "average" , clustering_method_columns = "average" , column_names_centered = TRUE , 
        column_names_gp = gpar(fontsize=15) , row_names_gp = gpar(fontsize=3) ,top_annotation = Col_Anno, show_column_names = T,
        show_row_names = TRUE,
        left_annotation = ha_anno)


