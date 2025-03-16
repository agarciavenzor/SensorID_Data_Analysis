
### This R Script contains the Network analysis using the igraph package of the SensorID interactomes treated with the shRNA targeting Nucleolin (shNcl), 
#and the control shRNA, shScramble.

#Intall the required package
install.packages(c("igraph", "bipartite", "asnipe", "assortnet", "ggplot2", "ggmap", "rnetcarto"
                   , "ecodist", "igraphdata", "statnet", "RColorBrewer", "tidyverse"))
#Load Packages
require(igraph)
require(igraphdata)
require(RColorBrewer)
require(tidyverse)
require(ggplot2)
require(dplyr)
require(Matrix)
require(MASS)
require(ggraph)
require(ComplexHeatmap)
require(forcats)

##Select Working Dir
setwd("C:/./shNcl_Network_Analysis")

# Load the Node metadata with a matrix format, it will contain the functional annotations of the proteins, and the Fold Change values versus the negative controls.
#Load the Matrix
Mat<-read.csv("./SensorID_shNclExp_Matrix.csv")
#Convert into categorical variables the columns with annotations, or other information regarding the specific T-Test performed
Mat$DNA_Repair<-as.factor(Mat$DNA_Repair)
Mat$Annotation<-as.factor(Mat$Annotation)
Mat$An.color<-as.factor(Mat$An.color)
Mat$An.number<-as.factor(Mat$An.number)
Mat$Sensor<-as.factor(Mat$Sensor)
Mat$Sirt6_shScr.Treatment<-as.factor(Mat$Sirt6_shScr.Treatment)
Mat$Sirt6_shNcl.Treatment<-as.factor(Mat$Sirt6_shNcl.Treatment)
Mat$Sirt6_DiffBound.Treatment<-as.factor(Mat$Sirt6_DiffBound.Treatment)
Mat$Ku80ID_shScr.Treatment<-as.factor(Mat$Ku80ID_shScr.Treatment)
Mat$Ku80ID_shNcl.Treatment<-as.factor(Mat$Ku80ID_shNcl.Treatment)
Mat$Ku80ID_DiffBound.Treatment<-as.factor(Mat$Ku80ID_DiffBound.Treatment)
Mat$Mre11_shScr.Treatment<-as.factor(Mat$Mre11_shScr.Treatment)
Mat$Mre11_shNcl.Treatment<-as.factor(Mat$Mre11_shNcl.Treatment)
Mat$Mre11_DiffBound.Treatment<-as.factor(Mat$Mre11_DiffBound.Treatment)
Mat$Sirt6_shScr.T.Test<-as.factor(Mat$Sirt6_shScr.T.Test)
Mat$Sirt6_shNcl.T.Test<-as.factor(Mat$Sirt6_shNcl.T.Test)
Mat$Sirt6_DiffBound.T.Test<-as.factor(Mat$Sirt6_DiffBound.T.Test)
Mat$Ku80ID_shScr.T.Test<-as.factor(Mat$Ku80ID_shScr.T.Test)
Mat$Ku80ID_shNcl.T.Test<-as.factor(Mat$Ku80ID_shNcl.T.Test)
Mat$Ku80ID_DiffBound.T.Test<-as.factor(Mat$Ku80ID_DiffBound.T.Test)
Mat$Mre11_shScr.T.Test<-as.factor(Mat$Mre11_shScr.T.Test)
Mat$Mre11_shNcl.T.Test<-as.factor(Mat$Mre11_shNcl.T.Test)
Mat$Mre11_DiffBound.T.Test<-as.factor(Mat$Mre11_DiffBound.T.Test)
Mat$Sirt6ID_Presence<-as.factor(Mat$Sirt6ID_Presence)
Mat$Ku80_Prescence<-as.factor(Mat$Ku80_Prescence)
Mat$Mre11_Prescence<-as.factor(Mat$Mre11_Prescence)
str(Mat)


### Similar to the procedure conducted in the analysis of the SensorID dynamic interactomes, here we will also start with the Full nuclear network.
##This will allow us to generate a single layout for all the networks, and similar dendograms for the position of the interactors in the networks and heatmaps.

#Full Network lists and files
#Node files
Mat#This may be the best node attribute list that we may have
###
#SensorID All nodes Network Edge Lists
SensorID_FNN<-read.delim("./SensorID_shNclExp_CompleteNetwork_EdgeList.tsv",sep = "\t",col.names = )
str(SensorID_FNN)
head(SensorID_FNN)

#Generate the igraph network object
FNN<-graph_from_data_frame(d=SensorID_FNN,vertices = Mat, directed = F )#Generate the igraph object with the edge and node annotation attributes
str(vertex_attr(FNN))#Check node attributes
str(edge_attr(FNN))#Check for edge attributes
FNN<-igraph::simplify(FNN, remove.multiple=TRUE, remove.loops=TRUE,edge.attr.comb = "mean") #Remove repetitive edges and self-connecting edges
FNN
#Small plot
plot(FNN)
#Generate an adjacency Matrix that will be used after
FNN_AdjMat<-as_adjacency_matrix(FNN,type="both",attr = "weight1",names=TRUE,sparse = TRUE)
FNN_AdjMat

#Generate a vertex size based on degree and count
vsize<-igraph::degree(FNN)/vcount(FNN)
#Identify and remove disconnected nodes
isolated<-which(igraph::degree(FNN)==0)
isolated
length(isolated)# 95 isolated nodes
#remove isolated nodes
FNN<-igraph::delete_vertices(FNN,isolated)

#Generate layouts for the graph
Lyt.Nice<-layout_nicely(FNN,dim = 3)
Lyt.GOpt<-layout_with_graphopt(FNN, niter=5000, mass = 30, charge = E(FNN)$weight)
Lyt.DrL<-layout_with_drl(FNN, weights = E(FNN)$weight)
LytFR<-layout_with_fr(FNN, niter=10000, weights=E(FNN)$weight)#This one looks like the best option
Lyt.KK<-layout_with_kk(FNN, maxiter=10000, weights=E(FNN)$weight) #very large separation of peripheral nodes
Lyt.GEM<-layout_with_gem(FNN, maxiter=10000) #Bad separation
Lyt.LGL<-layout_with_lgl(FNN, maxit=150,maxdelta=833,area=(length(E(FNN)))^2,coolexp=1.5,repulserad=((length(E(FNN)))^2)*(length(E(FNN))))#Do not work

#Star and Tree layouts these should be used when the specific SensorID networks are being analyzed
V(FNN)[name=="SIRT6"]
which(V(FNN)[name=="SIRT6"])# 610 #Check what is the rowindex for your desired node
which(V(FNN)$name=="SIRT6")# 960 #Check what is the rowindex for your desired node
which(V(FNN)$name=="MRE11A")# 130
which(V(FNN)$name=="XRCC5")# 141
Lyt.Star<-layout_as_star(FNN,center=V(FNN)[name=="SIRT6"]) #it has a special order argument where we can assign an order to the vertex, hence maybe the OPLS-DA:VIP value order may be used as ordering vector
Lyt.Tree<-layout_as_tree(FNN,root=c(960)) #For Sirt6ID networks
Lyt.RT<-layout.reingold.tilford(FNN,root=c(960))# For Sirt6ID networks
Lyt.sugi<-layout_with_sugiyama(FNN, maxiter=100, weights=E(FNN)$weight)#This one looks good for undirected tree layouts

#Check what layout looks fine
plot.igraph(FNN,layout=Lyt.Tree[-isolated,], alpha=0.6,
            vertex.color=V(FNN)$An.color, vertex.frame.color=V(FNN)$An.color, vertex.shape="circle",
            vertex.size=vsize*100, vertex.label=V(FNN)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(FNN)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)

#With ggraph
ggraph(FNN,layout=Lyt.FR)+
  geom_edge_fan(color="gray85", aes(width=(E(FNN)$weight1)/3), alpha=0.4)+
  geom_node_point(aes(color=V(FNN)$Annotation), size=10,alpha=0.75)+
  geom_node_text(aes(label=name),color="gray25",size=4)+
  scale_color_manual(limits=V(FNN)$Annotation, values=V(FNN)$An.color)+
  theme_void()

#Community Detection by the clique algorithm
clp <- igraph::cluster_optimal(FNN)
class(clp)
V(FNN)$community <- clp$membership
colrs <- adjustcolor( c("gray50", "tomato", "gold", "yellowgreen"), alpha=.6)

#Plot the communities detected previously 
plot.igraph(clp,FNN,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(FNN)$An.color, vertex.frame.color=V(FNN)$An.color, vertex.shape="circle",
            vertex.size=5, vertex.label=V(FNN)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(FNN)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


# Determine how distant is each node in the network to a desired node, in this case Sirt6
dist.from.Sirt6 <- distances(FNN, v=V(FNN)[name=="SIRT6"],
                             to=V(FNN), weights=E(FNN)$weight)
dist.from.Sirt6<-dist.from.Sirt6[,!is.infinite(colSums(dist.from.Sirt6))]
dist.from.Sirt6

#Set colors to plot the distances:
oranges <- colorRampPalette(c("dark red", "gold"))
col <- oranges(max(dist.from.Sirt6))
col <- col[dist.from.Sirt6+1]
#Plot the distance of each node in the network to Sirt6ID
plot(FNN, layout=Lyt.Tree[-isolated,],vertex.color=col, vertex.label=dist.from.Sirt6,
     vertex.label.color="white",vertex.shape="circle",
     vertex.size=5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(FNN)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)

#Shortest path from Sirt6 to Nucleolin
S6Ncl.path <- shortest_paths(FNN,
                             from = V(FNN)[name=="SIRT6"],
                             to = V(FNN)[name=="NCL"],
                             output = "both") # both path nodes and edges
# Generate edge color variable to plot the path:
ecol <- rep("gray80", ecount(FNN))
ecol[unlist(S6Ncl.path$epath)] <- "orange"
# Generate edge width variable to plot the path:
ew <- rep(2, ecount(FNN))
ew[unlist(S6Ncl.path$epath)] <- 4
# Generate node color variable to plot the path:
vcol <- rep("gray40", vcount(FNN))
vcol[unlist(S6Ncl.path$vpath)] <- "gold"

#Plot the path from Sirt6 to Nucleolin
plot(FNN,layout=Lyt.Tree[-isolated,], vertex.color=vcol, edge.color=ecol,
     edge.width=ew,vertex.label.color="white",vertex.shape="circle",
     vertex.size=5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(FNN)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)

#Identify Sirt6 incident edges and neighbor nodes, this can be applied to any node of the network.
inc.edges <- incident(FNN, V(FNN)[name=="SIRT6"], mode="all")
# Set colors to plot the selected edges.
ecol <- rep("gray80", ecount(FNN))
ecol[inc.edges] <- "orange"
vcol <- V(FNN)$An.color
vcol[V(FNN)$name=="SIRT6"] <- "gold"

#Plot the Sirt6 incident edges and the neighbor nodes
plot(FNN,layout=Lyt.Tree[-isolated,], vertex.color=vcol, edge.color=ecol, layout=Lyt.FR[-isolated,],
     edge.width=ew,vertex.label.color="white",vertex.shape="circle",
     vertex.size=5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(FNN)$weight, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)

#Network-level metrics
#Number of nodes and edges
vcount(FNN) #1129
ecount(FNN)#13828
#Density (Number of edges in relation to the total number of edges in a full connected graph)
edge_density(FNN,loops = F) #0.02171632
#Transitivity (Total number of conected triangles)
transitivity(FNN,type="global") #[1]0.511229
#modularity
vertex_attr_names(FNN)
#Modularity Measurements
#Annotation
mod<-V(FNN)$An.number #lets check for functional annotation
mod<-as.numeric(mod)
mod<-mod+1#The non clustered interactors that appears as 0, needs to be transformed to 1, so we add 1 to every cluster
modularity(FNN,membership = mod)
#Annotation modularity [1] 0.2792875
#Nuclear compartment
mod<-V(FNN)$CellCompartment #
mod<-as.factor(mod)
mod<-as.numeric(mod)
mod
mod<-mod+1#The non clustered interactors that appears as 0, needs to be transformed to 1, so we add 1 to every cluster
modularity(FNN,membership = mod)
#NucLoc modularity [1] [1] 0.07863669
#DNA_Repair
mod<-V(FNN)$DNA_Repair #lets check for DNA repair pathway
mod<-as.factor(mod)
mod<-as.numeric(mod)
mod
mod<-mod+1#The non clustered interactors that appears as 0, needs to be transformed to 1, so we add 1 to every cluster
modularity(FNN,membership = mod)
#DNA repair modularity [1] [1] 0.02559662

##Assorativity Measurements
#Assorativity
assortativity_nominal(FNN,factor(V(FNN)$Annotation),directed=F)
#Annotation [1] 0.3457159
assortativity_nominal(FNN,factor(V(FNN)$CellCompartment),directed=F)
#CellCompartment [1] 0.1330889
assortativity_nominal(FNN,factor(V(FNN)$DNA_Repair),directed=F)
#DNA_Repair [1] 0.04582501

#Diameter (longest geodesic distance), shortest path between two most peripheral nodes
diameter(FNN,directed = F, weights = E(FNN)$weight1)#[1] 4.989
diam<-get_diameter(FNN,directed = F,weights = E(FNN)$weight1)
diam #+ 8/1129 vertices, named, from 7a6d679:
#[1] TFB2M  POLRMT TFB1M  DAP3   RPS5   EIF4A1 NFYB   NFYC  
#Color variables for the diameter of the networks
vcol <- rep("gray50", vcount(FNN))
vcol[diam] <- "gold"
ecol <- rep("gray90", ecount(FNN))
ecol[E(FNN, path=diam)] <- "orange"

#Plot the diameter on the network
plot(FNN, layout=Lyt.Tree[-isolated,], vertex.color=vcol, edge.color=ecol, vertex.label.color="black",edge.arrow.mode=0)

#Node-level metrics descriptives
degree<-igraph::centr_degree(FNN,mode = "all",normalized = T)
degree<-degree$res
degree
eigen_cent<-igraph::eigen_centrality(FNN)
eigen_cent<-eigen_cent$vector
clos_cent <- igraph::closeness(FNN, normalized=T)
betw_cent <- igraph::betweenness(FNN, normalized=T)
hits <- hits_scores(FNN, weights=E(FNN)$weight)
hs <-hits$hub
authority<-hits$authority
Trans<-transitivity(FNN,type="local")
Trans
#write results
FNN.Centrality<-data.frame(NAME=V(FNN)$name, UniprotID=V(FNN)$UniprotID, Degree=degree, Eigen= eigen_cent, Closeness= clos_cent, Betweeness= betw_cent, Hubs=hs, Authority=authority,Transitivity=Trans)
FNN.Centrality

write.csv(FNN.Centrality,"./FNN.Centrality.csv")

#Cocitation, Number of shared interactions among nodes
Cocitation<-cocitation(FNN)
Cocitation#This matrix is similar to an adjacency matrix

#Plot with degree as node size
A1<-ggraph(FNN,layout=Lyt.Tree)+
  geom_edge_fan(color="gray85", aes(width=(E(FNN)$weight1)/3), alpha=0.4)+
  geom_node_point(aes(color=V(FNN)$Annotation, size=10+degree),alpha=0.75)+
  geom_node_text(aes(label=name),color="gray25",size=4)+
  scale_color_manual(limits=V(FNN)$Annotation, values=V(FNN)$An.color)+
  theme_void()

A1

#Igraph build-in modularity algortithms, they do not calculate the Connectivity and participation metrics.
#Edge_Betweness Modularity

FNN.ebmod=cluster_edge_betweenness(FNN,weights=E(FNN)$weight1,directed = F) 
FNN.ebmod #IGRAPH clustering edge betweenness, groups: 50, mod: 0.53
# Plot the Edge-betweness modules
plot(FNN.ebmod, FNN,layout=Lyt.Tree, alpha=0.6, 
     vertex.color=V(FNN)$An.color, vertex.frame.color=V(FNN)$An.color, vertex.shape="circle", 
     vertex.size=5, vertex.label=V(FNN)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black", 
     edge.width=E(FNN)$weight1, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout", 
     rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#Louvain modularity
FNN.Louvain.mod<-cluster_louvain(FNN,weights=E(FNN)$weight1)
FNN.Louvain.mod #The result shows a total of groups: 18, mod: 0.6

#Plot the Louvain modules
plot(FNN.Louvain.mod, FNN,layout=Lyt.Tree, alpha=0.6, 
     vertex.color=V(FNN)$An.color, vertex.frame.color=V(FNN)$An.color, vertex.shape="circle", 
     vertex.size=5, vertex.label=V(FNN)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black", 
     edge.width=E(FNN)$weight1, edge.lty=1, main="SensorID Full Nuclear Network", sub="FR Tree Layout", 
     rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#Generate HeatMaps of the adjacency matrixes
#Let's use the neighbor list for generating the Adjaency Matrix
SensorID_FNN
v<-levels(as.factor(unlist(SensorID_FNN[,1:2]))) #get vertex names from the network list
n<-length(v)#get the number of vertex in the network
e<-matrix(match(as.character(unlist(SensorID_FNN[,1:2])),v),ncol = 2)#Edgelist with row index instead of Names
w<-SensorID_FNN$weight1 #Edge weights

# M is our co-expression network adjacency matrix.

# distringuish between i,j and j,i we want to make M symmetric (undirected) by
# considering an undirected network.
FNN.M <- matrix(0, n, n)                  # set up a co-expression matrix
FNN.M[e] <- w                             # fill it in with edge weights
FNN.M <- FNN.M + t(FNN.M)                         # make this symmetric
dimnames(FNN.M) <- list(v, v)             # label the vertices
# We let A_ij = 1 if there is an interaction. 
#A <- 1*(M > 0) #right now we keep the weights as they are
FNN.M #The adjacency Matrix not as an igraph object, easier to manipulate, and to subset for smaller Sensor specific networks
quantile(FNN.M,c(0,0.25,0.5,0.75,1))
#
#0%   25%   50%   75%  100% 
#0.000 0.000 0.000 0.000 0.991 

#Creating the HeatMap
# make a heatmap image with colors: white below 0,
# grey from 0.4 to 1 with darker closer to 1
breaks <- c(0, seq(0.2, 1, length=11)) # breaks in the color bins
cols <- grey(1-c(0,seq(0.5,1,length=10))) # colors

# color the vertex for SensorID nodes in blue, its neighbors red, and the others white
SensorIDBind <- which(v == "SIRT6"|v == "MRE11A"|v == "XRCC5"|v=="NCL") # index for Sirt6ID
#SensorIDBind<-as.factor(SensorIDBind)
vcols <- rep("white",n)
vcols[SensorIDBind] <- "blue"
#Subsetting FNN.M by doing a loop asking for look into every row and column the index in SensorIDBind worked
for (i in SensorIDBind) {
  vcols[which(FNN.M[,i] > 0 | FNN.M[i,])] <- "red"
}

#Generate a dendogram for keeping the proteins in the sample place
FNN.dist<-dist(FNN.M,method = "binary")
FNN.Dend<-hclust(FNN.dist, method = "ward.D2" )
plot(FNN.Dend)

#Cut the dendogram
cut.height <- 10 # try numbers from 4 to 7
FNN.hclusts<- cutree(FNN.Dend, cut.height)
plot(FNN.hclusts)

# now actually make the heat map
par(mar=rep(0,4)) # make plot margins 0
heatmap(FNN.M, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=vcols, RowSideColors=vcols,  col=0:1, frame=T)+
  legend("topleft", c("Neighbors(SIRT6, MRE11, KU80)","SIRT6, MRE11, XRCC5"), fill=c("red","blue"),  
         bty="n", inset=0, xpd=T,  border=F)
##Color annotation for functional annotation
MatSet<-subset(Mat,Mat$NAME%in%v)# select only those proteins that are cointained in the node list of our matrixes
match(v,MatSet$NAME) # Check that the indexes of the Annotation Dataframe and the node list are in the same position 
match(rownames(FNN.M),MatSet$NAME) #Check that the row names of the Adj Matrix and the Annotation dataframe are in the same order
match(colnames(FNN.M),MatSet$NAME) #Check that the column names of the Adj Matrix and the Annotation dataframe are in the same order
AnCol<-MatSet$An.color#Vector of colors
AnCol<-as.character(AnCol)#transform it into a character vector, so the colors can be read


#
plot.new()
legend("topleft", c("Cell Cycle","Chromatin Org","DNA Repair","Nuclear Pore","p53","Proteasome","Ribosome Biogenesis","RNA Metabolism","Transcription Reg","UPR","Calcium Signaling"), 
       fill=c("coral1","darkviolet","firebrick1","deepskyblue3","tomato","grey","lightgreen","gold","slateblue1","maroon","sienna2"),
       bty="n", inset=0, xpd=T,  border=F, cex =2,y.intersp = 0.7,x.intersp = 0.2,pt.cex = 1)

#check with the FNN
# now actually make the heat map
par(mar=rep(0,4)) # make plot margins 0
heatmap(FNN.M, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=0:1, frame=T, cexRow = 0.5 , cexCol = 0.5)

legend("topleft", c("Cell Cycle","Chromatin Org","DNA Repair","Nuclear Pore","p53","Proteasome","Ribosome Biogenesis","RNA Metabolism","Transcription Reg","UPR","Calcium Signaling"), 
       fill=c("coral1","darkviolet","firebrick1","deepskyblue3","tomato","grey","lightgreen","gold","slateblue1","maroon","sienna2"), bty="n", inset=0, xpd=T,  border=F, cex =2,y.intersp = 0.7,x.intersp = 0.2,pt.cex = 1)



#
dev.off()
plot.new()

###Sirt6 Networks of the shNcl experiment
###Sirt6 shScr Network: The Sirt6ID IR2h interactome after transduction with the control shRNA, shScr.

#S6ID All nodes Network Edge Lists
S6ID_S6shScr<-read.delim("./Sirt6ID_shScrlvsCneg_Network_EdgeList.tsv",sep = "\t",col.names = )
str(S6ID_S6shScr)
head(S6ID_S6shScr)

#Generate the igraph network object
S6shScr<-graph_from_data_frame(d=S6ID_S6shScr,vertices = Mat, directed = F )#Generate the igraph object with the edge and node annotation attributes
str(vertex_attr(S6shScr))#Check node attributes
str(edge_attr(S6shScr))#Check for edge attributes
S6shScr<-igraph::simplify(S6shScr, remove.multiple=TRUE, remove.loops=TRUE,edge.attr.comb = "first") #Remove repetitive edges and self-connecting edges
S6shScr

#Small plot
plot(S6shScr)
#Generate an adjacency Matrix that will be used after
S6shScr_AdjMat<-as_adjacency_matrix(S6shScr,type="both",attr = "weight1",names=TRUE,sparse = TRUE)
S6shScr_AdjMat
#Generate a vertex size based on degree and count
vsize<-igraph::degree(S6shScr)/vcount(S6shScr)
#Identify and remove disconnected nodes
isolated<-which(igraph::degree(S6shScr)==0)
isolated
length(isolated)# 484
#remove isolated nodes
S6shScr<-igraph::delete_vertices(S6shScr,isolated)

#Generate layouts for the graph
#Lyt.Nice<-layout_nicely(S6shScr,dim = 3)
#Lyt.GOpt<-layout_with_graphopt(S6shScr, niter=5000, mass = 30, charge = E(S6shScr)$weight)
#Lyt.DrL<-layout_with_drl(S6shScr, weights = E(S6shScr)$weight)
#LytFR<-layout_with_fr(S6shScr, niter=10000, weights=E(S6shScr)$weight)#This one looks like the best option
#Lyt.KK<-layout_with_kk(S6shScr, maxiter=10000, weights=E(S6shScr)$weight) #very large separation of peripheral nodes
#Lyt.GEM<-layout_with_gem(S6shScr, maxiter=10000) #Bad separation
#Lyt.LGL<-layout_with_lgl(S6shScr, maxit=150,maxdelta=833,area=(length(E(S6shScr)))^2,coolexp=1.5,repulserad=((length(E(S6shScr)))^2)*(length(E(S6shScr))))#Do not work


#Star and Tree layouts
V(S6shScr)[name=="SIRT6"]
which(V(S6shScr)[name=="SIRT6"])# 610 #Check what is the rowindex for your desired node
which(V(S6shScr)$name=="SIRT6")# 665 #Check what is the rowindex for your desired node
which(V(S6shScr)$name=="MRE11A")# 104
which(V(S6shScr)$name=="XRCC5")# 114
Lyt.Star<-layout_as_star(S6shScr,center=V(S6shScr)[name=="SIRT6"]) #it has a special order argument where we can assign an order to the vertex, hence maybe the OPLS-DA:VIP value order may be used as ordering vector

Lyt.Tree<-layout_as_tree(S6shScr,root=c(665))
Lyt.RT<-layout.reingold.tilford(S6shScr,root=c(665))# This one looks good, but the indexing do not allow to highlight the correct nodes
Lyt.sugi<-layout_with_sugiyama(S6shScr, maxiter=100, weights=E(S6shScr)$weight)#This one looks good for undirected tree layouts

#Check what layout looks fine

plot.igraph(S6shScr,layout=Lyt.Tree, alpha=0.6,
            vertex.color=V(S6shScr)$An.color, vertex.frame.color=V(S6shScr)$An.color, vertex.shape="circle",
            vertex.size=vsize*100, vertex.label=V(S6shScr)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(S6shScr)$weight1, edge.lty=1, main="S6ID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)

#With ggraph

ggraph(S6shScr,layout=Lyt.Tree)+
  geom_edge_fan(color="gray85", aes(width=(E(S6shScr)$weight1)/3), alpha=0.4)+
  geom_node_point(aes(color=V(S6shScr)$Annotation), size=10,alpha=0.75)+
  geom_node_text(aes(label=name),color="gray25",size=4)+
  scale_color_manual(limits=V(S6shScr)$Annotation, values=V(S6shScr)$An.color)+
  theme_void()

#Community Detection
clp <- igraph::cluster_optimal(S6shScr)
class(clp)
V(S6shScr)$community <- clp$membership
colrs <- adjustcolor( c("gray50", "tomato", "gold", "yellowgreen"), alpha=.6)

# Plot the clique detected communities
plot.igraph(clp,S6shScr,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(S6shScr)$An.color, vertex.frame.color=V(S6shScr)$An.color, vertex.shape="circle",
            vertex.size=5, vertex.label=V(S6shScr)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(S6shScr)$weight, edge.lty=1, main="S6ID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#Calculate the Distances of each node in the network to Sirt6
dist.from.Sirt6 <- distances(S6shScr, v=V(S6shScr)[name=="SIRT6"],
                             to=V(S6shScr), weights=E(S6shScr)$weight)
dist.from.Sirt6<-dist.from.Sirt6[,!is.infinite(colSums(dist.from.Sirt6))]

dist.from.Sirt6
range(dist.from.Sirt6)

#Set colors to plot the distances:
oranges <- colorRampPalette(c("dark red", "gold"))
col <- oranges(max(dist.from.Sirt6))
col <- col[dist.from.Sirt6+1]
plot(S6shScr, layout=Lyt.Tree[-isolated,],vertex.color=col, vertex.label=dist.from.Sirt6,
     vertex.label.color="white",vertex.shape="circle",
     vertex.size=5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(S6shScr)$weight, edge.lty=1, main="S6ID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)


#Shortest path from Sirt6 to Nucleolin
S6Ncl.path <- shortest_paths(S6shScr,
                             from = V(S6shScr)[name=="SIRT6"],
                             to = V(S6shScr)[name=="NCL"],
                             output = "both") # both path nodes and edges
# Generate edge color variable to plot the path:
ecol <- rep("gray80", ecount(S6shScr))
ecol[unlist(S6Ncl.path$epath)] <- "orange"
# Generate edge width variable to plot the path:
ew <- rep(2, ecount(S6shScr))
ew[unlist(S6Ncl.path$epath)] <- 4
# Generate node color variable to plot the path:
vcol <- rep("gray40", vcount(S6shScr))
vcol[unlist(S6Ncl.path$vpath)] <- "gold"
#Plot the Distance between Sirt6 and Nucleolin in the network
plot(S6shScr,layout=Lyt.Tree, vertex.color=vcol, edge.color=ecol,
     edge.width=ew,vertex.label.color="white",vertex.shape="circle",
     vertex.size=5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(S6shScr)$weight1, edge.lty=1, main="S6ID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)

#Plot the Sirt6 incident edges and neighbot nodes
inc.edges <- incident(S6shScr, V(S6shScr)[name=="SIRT6"], mode="all")
inc.edges
# Set colors to plot the selected edges.
ecol <- rep("gray80", ecount(S6shScr))
ecol[inc.edges] <- "orange"
vcol <- V(S6shScr)$An.color
vcol[V(S6shScr)$name=="SIRT6"] <- "gold"

#Plot the Sirt6 incident edges and neighbot nodes
plot(S6shScr,layout=Lyt.Tree, vertex.color=vcol, edge.color=ecol, layout=Lyt.FR[-isolated,],
     edge.width=ew,vertex.label.color="white",vertex.shape="circle",
     vertex.size=5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(S6shScr)$weight1, edge.lty=1, main="S6ID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)

#Plot Degree Distribution
S6shScr.Deg<-degree_distribution(S6shScr, cumulative = T, mode="all")
plot( x=0:max(igraph::degree(S6shScr)), y=1-S6shScr.Deg, pch=19, cex=1.2, col="orange",
      xlab="Degree", ylab="Cumulative Frequency")

#Network-level Descriptive metrics.
#Number of nodes and edges
vcount(S6shScr) #772
ecount(S6shScr)#7953
#Density (Number of edges in relation to the total number of edges in a full connected graph)
edge_density(S6shScr,loops = F) #0.02672325
#Transitivity (Total number of conected triangles)
transitivity(S6shScr,type="global") #[1]0.5343984
#modularity
vertex_attr_names(S6shScr)
#Modularity Measurements
#Annotation
mod<-V(S6shScr)$An.number #lets check for Annotation
mod<-as.numeric(mod)
mod<-mod+1#The non clustered interactors that appears as 0, needs to be transformed to 1, so we add 1 to every cluster
modularity(S6shScr,membership = mod)
#Annotation modularity [1] 0.2630792
#Nuclear compartment
mod<-V(S6shScr)$CellCompartment #lets check for Nuclear compartment modularity
mod<-as.factor(mod)
mod<-as.numeric(mod)
mod
mod<-mod+1#The non clustered interactors that appears as 0, needs to be transformed to 1, so we add 1 to every cluster
modularity(S6shScr,membership = mod)
#NucLoc modularity [1] [1] 0.08559185
#DNA_Repair
mod<-V(S6shScr)$DNA_Repair #lets check for DNA_Repair pathway modularity
mod<-as.factor(mod)
mod<-as.numeric(mod)
mod
mod<-mod+1#The non clustered interactors that appears as 0, needs to be transformed to 1, so we add 1 to every cluster
modularity(S6shScr,membership = mod)
#NucLoc modularity [1] [1] 0.02786796

##Assorativity Measurements
#Assorativity
assortativity_nominal(S6shScr,factor(V(S6shScr)$Annotation),directed=F)
#Annotation [1] 0.3302574
assortativity_nominal(S6shScr,factor(V(S6shScr)$CellCompartment),directed=F)
#CellCompartment [1] 0.1346954
assortativity_nominal(S6shScr,factor(V(S6shScr)$DNA_Repair),directed=F)
#DNA_Repair [1] 0.05219048


#Diameter (longest geodesic distance), shortest path between two most peripheral nodes
diameter(S6shScr,directed = F, weights = E(S6shScr)$weight1)#[1] 5.876
diam<-get_diameter(S6shScr,directed = F,weights = E(S6shScr)$weight1)
diam #+ 11/772 vertices, named, from 3333c65:
#[1] RAP1B   RAB5C   RAP1A   CDC42   RACGAP1 AURKB   RPS16   DAP3    TFB1M   POLRMT  TFB2M  
vcol <- rep("gray50", vcount(S6shScr))
vcol[diam] <- "gold"
ecol <- rep("gray90", ecount(S6shScr))
ecol[E(S6shScr, path=diam)] <- "orange"
#Plot the Diameter in the network
plot(S6shScr, layout=Lyt.Tree, vertex.color=vcol, edge.color=ecol, vertex.label.color="black",edge.arrow.mode=0)

#Node-level descriptive metrics.
degree<-igraph::centr_degree(S6shScr,mode = "all",normalized = T)
degree<-degree$res
degree
eigen_cent<-igraph::eigen_centrality(S6shScr)
eigen_cent<-eigen_cent$vector
clos_cent <- igraph::closeness(S6shScr, normalized=T)
betw_cent <- igraph::betweenness(S6shScr, normalized=T)
hits <- hits_scores(S6shScr, weights=E(S6shScr)$weight)
hs <-hits$hub
authority<-hits$authority
Trans<-transitivity(S6shScr,type="local")
Trans
#write results
S6shScr.Centrality<-data.frame(NAME=V(S6shScr)$name, UniprotID=V(S6shScr)$UniprotID, Degree=degree, Eigen= eigen_cent, Closeness= clos_cent, Betweeness= betw_cent, Hubs=hs, Authority=authority,Transitivity=Trans)
S6shScr.Centrality

write.csv(S6shScr.Centrality,"./S6shScr.Centrality.csv")

#Cocitation, Number of shared interactions among nodes
Cocitation<-cocitation(S6shScr)
Cocitation#This matrix is similar to an adjacency matrix

#Plot with degree as node size
A1<-ggraph(S6shScr,layout=Lyt.Tree)+
  geom_edge_fan(color="gray85", aes(width=(E(S6shScr)$weight1)/3), alpha=0.4)+
  geom_node_point(aes(color=V(S6shScr)$Annotation, size=10+degree),alpha=0.75)+
  geom_node_text(aes(label=name),color="gray25",size=4)+
  scale_color_manual(limits=V(S6shScr)$Annotation, values=V(S6shScr)$An.color)+
  theme_void()
A1
#Modularity Analysis
#Edge_Betweeness Modularity
S6shScr.ebmod=cluster_edge_betweenness(S6shScr,weights=E(S6shScr)$weight1,directed = F) 
S6shScr.ebmod #IGRAPH clustering edge betweenness, groups: 93, mod: 0.58
# Plot the Edge-betweness modules
plot(S6shScr.ebmod, S6shScr,layout=Lyt.Tree, alpha=0.6, 
     vertex.color=V(S6shScr)$An.color, vertex.frame.color=V(S6shScr)$An.color, vertex.shape="circle", 
     vertex.size=5, vertex.label=V(S6shScr)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black", 
     edge.width=E(S6shScr)$weight1, edge.lty=1, main="S6ID Full Nuclear Network", sub="FR Tree Layout", 
     rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)

#Dendogram of the Edge_Betweeness
plot_dendrogram(S6shScr.ebmod, mode="hclust")

#Louvain modularity
S6shScr.Louvain.mod<-cluster_louvain(S6shScr,weights=E(S6shScr)$weight1)
S6shScr.Louvain.mod #The result shows a total of IGRAPH clustering multi level, groups: 19, mod: 0.61

#Plot the Louvain modules
plot(S6shScr.Louvain.mod, S6shScr,layout=Lyt.Tree, alpha=0.6, 
     vertex.color=V(S6shScr)$An.color, vertex.frame.color=V(S6shScr)$An.color, vertex.shape="circle", 
     vertex.size=5, vertex.label=V(S6shScr)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black", 
     edge.width=E(S6shScr)$weight1, edge.lty=1, main="S6ID Full Nuclear Network", sub="FR Tree Layout", 
     rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)

#Generate HeatMaps of the adjacency matrixes
#Let's use the neighbor list for generating the Adjaency Matrix
SensorID_FNN
v<-levels(as.factor(unlist(SensorID_FNN[,1:2]))) #get vertex names from the network list
n<-length(v)#get the number of vertex in the network
e<-matrix(match(as.character(unlist(SensorID_FNN[,1:2])),v),ncol = 2)#Edgelist with row index instead of Names
w<-S6ID_S6shScr$weight1 #Edge weights

# M is our co-expression network adjacency matrix.

# distringuish between i,j and j,i we want to make M symmetric (undirected) by
# considering an undirected network.
S6shScr.M <- matrix(0, n, n)                  # set up a co-expression matrix
S6shScr.M[e] <- w                             # fill it in with edge weights
S6shScr.M <- S6shScr.M + t(S6shScr.M)                         # make this symmetric
dimnames(S6shScr.M) <- list(v, v)             # label the vertices
# We let A_ij = 1 if there is an interaction. 
#A <- 1*(M > 0) #right now we keep the weights as they are
S6shScr.M #The adjacency Matrix not as an igraph object, easier to manipulate, and to subset for smaller Sensor specific networks
quantile(S6shScr.M,c(0,0.25,0.5,0.75,1))
#
#0%   25%   50%   75%  100% 
#0.000 0.000 0.000 0.000 0.991 

#Creating the HeatMap
# make a heatmap image with colors: white below 0,
# grey from 0.4 to 1 with darker closer to 1
breaks <- c(0, seq(0.05, 1, length=11)) # breaks in the color bins
cols <- grey(1-c(0,seq(0.5,1,length=10))) # colors

# color the vertex for Sirt6 blue, its neighbors red, and the others white
S6IDBind <- which(v == "SIRT6") # index for Sirt6
#S6IDBind<-as.factor(S6IDBind)
vcols <- rep("white",n)
vcols[S6IDBind] <- "blue"
#Subsetting S6shScr.M by doing a loop asking for look into every row and column the index in S6IDBind worked
for (i in S6IDBind) {
  vcols[which(S6shScr.M[,i] > 0 | S6shScr.M[i,])] <- "red"
}
#Subsetting S6shScr.M using the index vector in S6IDBind do not work
#vcols[which(S6shScr.M[,S6IDBind] > 0 | S6shScr.M[S6IDBind,])] <- "red"

#Generate a dendogram for keeping the proteins in the sample place
#S6shScr.dist<-dist(S6shScr.M,method = "euclidian")
#S6shScr.Dend<-hclust(S6shScr.dist, method = "ward.D2" )
#plot(S6shScr.Dend)

#Cut the dendogram
#cut.height <- 10 # try numbers from 4 to 7
#S6shScr.hclusts<- cutree(S6shScr.Dend, cut.height)
#plot(S6shScr.hclusts)

# now actually make the heat map
par(mar=rep(0,4)) # make plot margins 0
heatmap(S6shScr.M, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=vcols, RowSideColors=vcols,  col=0:1, frame=T)+
  legend("topleft", c("Neighbors(SIRT6, MRE11, KU80)","SIRT6, MRE11, XRCC5"), fill=c("red","blue"),  
         bty="n", inset=0, xpd=T,  border=F)
##Color annotation for functional annotation
MatSet<-subset(Mat,Mat$NAME%in%v)# select only those proteins that are cointained in the node list of our matrixes
match(v,MatSet$NAME) # Check that the indexes of the Annotation Dataframe and the node list are in the same position 
match(rownames(S6shScr.M),MatSet$NAME) #Check that the row names of the Adj Matrix and the Annotation dataframe are in the same order
match(colnames(S6shScr.M),MatSet$NAME) #Check that the column names of the Adj Matrix and the Annotation dataframe are in the same order
AnCol<-MatSet$An.color#Vector of colors
AnCol<-as.character(AnCol)#transform it into a character vector, so the colors can be read


#
plot.new()
legend("topleft", c("Cell Cycle","Chromatin Org","DNA Repair","Nuclear Pore","p53","Proteasome","Ribosome Biogenesis","RNA Metabolism","Transcription Reg","UPR","Calcium Signaling"), 
       fill=c("coral1","darkviolet","firebrick1","deepskyblue3","tomato","grey","lightgreen","gold","slateblue1","maroon","sienna2"),
       bty="n", inset=0, xpd=T,  border=F, cex =2,y.intersp = 0.7,x.intersp = 0.2,pt.cex = 1)

#check with the S6shScr
# now actually make the heat map
par(mar=rep(0,4)) # make plot margins 0
heatmap(S6shScr.M, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=0:1, frame=T, cexRow = 0.5 , cexCol = 0.5)

legend("topleft", c("Cell Cycle","Chromatin Org","DNA Repair","Nuclear Pore","p53","Proteasome","Ribosome Biogenesis","RNA Metabolism","Transcription Reg","UPR","Calcium Signaling"), 
       fill=c("coral1","darkviolet","firebrick1","deepskyblue3","tomato","grey","lightgreen","gold","slateblue1","maroon","sienna2"), bty="n", inset=0, xpd=T,  border=F, cex =2,y.intersp = 0.7,x.intersp = 0.2,pt.cex = 1)



###Sirt6ID shNcl Networks: Sirt6ID IR2h interactome after treatment with the shRNA targeting Nucleolin: shNcl
#S6ID All nodes Network Edge Lists
S6ID_S6shNcl<-read.delim("./Sirt6ID_shNclvsCneg_Network_EdgeList.tsv",sep = "\t",col.names = )
str(S6ID_S6shNcl)
head(S6ID_S6shNcl)

#Generate the igraph network object
S6shNcl<-graph_from_data_frame(d=S6ID_S6shNcl,vertices = Mat, directed = F )#Generate the igraph object with the edge and node annotation attributes
str(vertex_attr(S6shNcl))#Check node attributes
str(edge_attr(S6shNcl))#Check for edge attributes
S6shNcl<-igraph::simplify(S6shNcl, remove.multiple=TRUE, remove.loops=TRUE,edge.attr.comb = "first") #Remove repetitive edges and self-connecting edges
S6shNcl
vcount(S6shNcl)

#Small plot
plot(S6shNcl)
#Generate an adjacency Matrix that will be used after
S6shNcl_AdjMat<-as_adjacency_matrix(S6shNcl,type="both",attr = "weight1",names=TRUE,sparse = TRUE)
S6shNcl_AdjMat
#Generate a vertex size based on degree and count
vsize<-igraph::degree(S6shNcl)/vcount(S6shNcl)
#Identify and remove disconnected nodes
isolated<-which(igraph::degree(S6shNcl)==0)
isolated
length(isolated)# 637
#remove isolated nodes
S6shNcl<-igraph::delete_vertices(S6shNcl,isolated)

igraph::vcount(S6shNcl)
#Generate layouts for the graph
#Lyt.Nice<-layout_nicely(S6shNcl,dim = 3)
#Lyt.GOpt<-layout_with_graphopt(S6shNcl, niter=5000, mass = 30, charge = E(S6shNcl)$weight)
#Lyt.DrL<-layout_with_drl(S6shNcl, weights = E(S6shNcl)$weight)
#LytFR<-layout_with_fr(S6shNcl, niter=10000, weights=E(S6shNcl)$weight)#This one looks like the best option
#Lyt.KK<-layout_with_kk(S6shNcl, maxiter=10000, weights=E(S6shNcl)$weight) #very large separation of peripheral nodes
#Lyt.GEM<-layout_with_gem(S6shNcl, maxiter=10000) #Bad separation
#Lyt.DH<-layout_with_dh(S6shNcl, maxiter=150, fineiter=max(10, log2(833)))
#Lyt.LGL<-layout_with_lgl(S6shNcl, maxit=150,maxdelta=833,area=(length(E(S6shNcl)))^2,coolexp=1.5,repulserad=((length(E(S6shNcl)))^2)*(length(E(S6shNcl))))#Do not work

#Star and Tree layouts
V(S6shNcl)[name=="SIRT6"]
which(V(S6shNcl)[name=="SIRT6"])# 558 #Check what is the rowindex for your desired node
which(V(S6shNcl)$name=="SIRT6")# 558 #Check what is the rowindex for your desired node
which(V(S6shNcl)$name=="MRE11A")# 115
which(V(S6shNcl)$name=="XRCC5")# 124
Lyt.Star<-layout_as_star(S6shNcl,center=V(S6shNcl)[name=="SIRT6"]) #it has a special order argument where we can assign an order to the vertex, hence maybe the OPLS-DA:VIP value order may be used as ordering vector
Lyt.Tree<-layout_as_tree(S6shNcl,root=c(558))
Lyt.RT<-layout.reingold.tilford(S6shNcl,root=c(558))# This one looks good, but the indexing do not allow to highlight the correct nodes
Lyt.sugi<-layout_with_sugiyama(S6shNcl, maxiter=100, weights=E(S6shNcl)$weight)#This one looks good for undirected tree layouts

#Check what layout looks fine

plot.igraph(S6shNcl,layout=Lyt.Tree, alpha=0.6,
            vertex.color=V(S6shNcl)$An.color, vertex.frame.color=V(S6shNcl)$An.color, vertex.shape="circle",
            vertex.size=vsize*100, vertex.label=V(S6shNcl)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(S6shNcl)$weight1, edge.lty=1, main="S6ID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)

#With ggraph
ggraph(S6shNcl,layout=Lyt.Tree)+
  geom_edge_fan(color="gray85", aes(width=(E(S6shNcl)$weight1)/3), alpha=0.4)+
  geom_node_point(aes(color=V(S6shNcl)$Annotation), size=10,alpha=0.75)+
  geom_node_text(aes(label=name),color="gray25",size=4)+
  scale_color_manual(limits=V(S6shNcl)$Annotation, values=V(S6shNcl)$An.color)+
  theme_void()

#igraph clique Community Detection
clp <- igraph::cluster_optimal(S6shNcl)
class(clp)
V(S6shNcl)$community <- clp$membership
colrs <- adjustcolor( c("gray50", "tomato", "gold", "yellowgreen"), alpha=.6)

#Plot the clique detected communities
plot.igraph(clp,S6shNcl,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(S6shNcl)$An.color, vertex.frame.color=V(S6shNcl)$An.color, vertex.shape="circle",
            vertex.size=5, vertex.label=V(S6shNcl)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(S6shNcl)$weight, edge.lty=1, main="S6ID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)

#Calculate the distances from every node in the network to Sirt6 
dist.from.Sirt6 <- distances(S6shNcl, v=V(S6shNcl)[name=="SIRT6"],
                             to=V(S6shNcl), weights=E(S6shNcl)$weight1)
dist.from.Sirt6<-dist.from.Sirt6[,!is.infinite(colSums(dist.from.Sirt6))]

dist.from.Sirt6
range(dist.from.Sirt6)

#Set colors to plot the distances:
oranges <- colorRampPalette(c("dark red", "gold"))
col <- oranges(max(dist.from.Sirt6))
col <- col[dist.from.Sirt6+1]
plot(S6shNcl, layout=Lyt.Tree,vertex.color=col, vertex.label=dist.from.Sirt6,
     vertex.label.color="white",vertex.shape="circle",
     vertex.size=5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(S6shNcl)$weight1, edge.lty=1, main="S6ID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)


#Shortest path from Sirt6 to Nucleolin
S6Ncl.path <- shortest_paths(S6shNcl,
                             from = V(S6shNcl)[name=="SIRT6"],
                             to = V(S6shNcl)[name=="NCL"],
                             output = "both") # both path nodes and edges
# Generate edge color variable to plot the path:
ecol <- rep("gray80", ecount(S6shNcl))
ecol[unlist(S6Ncl.path$epath)] <- "orange"
# Generate edge width variable to plot the path:
ew <- rep(2, ecount(S6shNcl))
ew[unlist(S6Ncl.path$epath)] <- 4
# Generate node color variable to plot the path:
vcol <- rep("gray40", vcount(S6shNcl))
vcol[unlist(S6Ncl.path$vpath)] <- "gold"
#Plot the path from Sirt6 to Nucleolin in the network
plot(S6shNcl,layout=Lyt.Tree, vertex.color=vcol, edge.color=ecol,
     edge.width=ew,vertex.label.color="white",vertex.shape="circle",
     vertex.size=5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(S6shNcl)$weight1, edge.lty=1, main="S6ID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)

# Identify the Sirt6 incident edges and neighbor nodes
inc.edges <- incident(S6shNcl, V(S6shNcl)[name=="SIRT6"], mode="all")
inc.edges
# Set colors to plot the selected edges.
ecol <- rep("gray80", ecount(S6shNcl))
ecol[inc.edges] <- "orange"
vcol <- V(S6shNcl)$An.color
vcol[V(S6shNcl)$name=="SIRT6"] <- "gold"
#Plot the Sirt6 incident edges and neighbor nodes
plot(S6shNcl,layout=Lyt.Tree, vertex.color=vcol, edge.color=ecol, layout=Lyt.FR[-isolated,],
     edge.width=ew,vertex.label.color="white",vertex.shape="circle",
     vertex.size=5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(S6shNcl)$weight1, edge.lty=1, main="S6ID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)

#Plot Degree Distribution
S6shNcl.Deg<-degree_distribution(S6shNcl, cumulative = T, mode="all")
plot( x=0:max(igraph::degree(S6shNcl)), y=1-S6shNcl.Deg, pch=19, cex=1.2, col="orange",
      xlab="Degree", ylab="Cumulative Frequency")

#Network-level metrics
#Number of nodes and edges
vcount(S6shNcl) #619
ecount(S6shNcl)#5021
#Density (Number of edges in relation to the total number of edges in a full connected graph)
edge_density(S6shNcl,loops = F) #0.02625071
#Transitivity (Total number of conected triangles)
transitivity(S6shNcl,type="global") #[1]0.4961447
#modularity
vertex_attr_names(S6shNcl)
#Modularity Measurements
#Annotation
mod<-V(S6shNcl)$An.number #lets check for Annotation
mod<-as.numeric(mod)
mod<-mod+1#The non clustered interactors that appears as 0, needs to be transformed to 1, so we add 1 to every cluster
modularity(S6shNcl,membership = mod)
#Annotation modularity [1] 0.2630792
#Nuclear compartment
mod<-V(S6shNcl)$CellCompartment #lets check for Nuclear compartment modularity
mod<-as.factor(mod)
mod<-as.numeric(mod)
mod
mod<-mod+1#The non clustered interactors that appears as 0, needs to be transformed to 1, so we add 1 to every cluster
modularity(S6shNcl,membership = mod)
#NucLoc modularity [1] [1] 0.07009921
#DNA_Repair
mod<-V(S6shNcl)$DNA_Repair #lets check for DNA_Repair modularity
mod<-as.factor(mod)
mod<-as.numeric(mod)
mod
mod<-mod+1#The non clustered interactors that appears as 0, needs to be transformed to 1, so we add 1 to every cluster
modularity(S6shNcl,membership = mod)
#DNA_Repair modularity [1] [1] 0.02742711


##Assorativity Measurements
#Assorativity
assortativity_nominal(S6shNcl,factor(V(S6shNcl)$Annotation),directed=F)
#Annotation [1] 0.356348
assortativity_nominal(S6shNcl,factor(V(S6shNcl)$CellCompartment),directed=F)
#CellCompartment [1] 0.1357395
assortativity_nominal(S6shNcl,factor(V(S6shNcl)$DNA_Repair),directed=F)
#DNA_Repair [1]0.04571016


#Diameter (longest geodesic distance), shortest path between two most peripheral nodes
diameter(S6shNcl,directed = F, weights = E(S6shNcl)$weight1)#[1] 4.835
diam<-get_diameter(S6shNcl,directed = F,weights = E(S6shNcl)$weight1)
diam #12/619 vertices, named, from 3f82947:
#[1] RAP1B  RAB5C  RAP1A  CDC42  EPHA2  EGFR   XRCC5  YBX1   MRPS27 GRSF1  NDUFS3 UQCRC2  
vcol <- rep("gray50", vcount(S6shNcl))
vcol[diam] <- "gold"
ecol <- rep("gray90", ecount(S6shNcl))
ecol[E(S6shNcl, path=diam)] <- "orange"
#Plot the diameter of the network
plot(S6shNcl, layout=Lyt.Tree, vertex.color=vcol, edge.color=ecol, vertex.label.color="black",edge.arrow.mode=0)

#Node-level metrics
degree<-igraph::centr_degree(S6shNcl,mode = "all",normalized = T)
degree<-degree$res
degree
eigen_cent<-igraph::eigen_centrality(S6shNcl)
eigen_cent<-eigen_cent$vector
clos_cent <- igraph::closeness(S6shNcl, normalized=T)
betw_cent <- igraph::betweenness(S6shNcl, normalized=T)
hits <- hits_scores(S6shNcl, weights=E(S6shNcl)$weight)
hs <-hits$hub
authority<-hits$authority
Trans<-transitivity(S6shNcl,type="local")
Trans
#write results
S6shNcl.Centrality<-data.frame(NAME=V(S6shNcl)$name, UniprotID=V(S6shNcl)$UniprotID, Degree=degree, Eigen= eigen_cent, Closeness= clos_cent, Betweeness= betw_cent, Hubs=hs, Authority=authority,Transitivity=Trans)
S6shNcl.Centrality

write.csv(S6shNcl.Centrality,"./S6shNcl.Centrality.csv")

#Cocitation, Number of shared interactions among nodes
Cocitation<-cocitation(S6shNcl)
Cocitation#This matrix is similar to an adjacency matrix

#Plot with degree as node size
A1<-ggraph(S6shNcl,layout=Lyt.Tree)+
  geom_edge_fan(color="gray85", aes(width=(E(S6shNcl)$weight1)/3), alpha=0.4)+
  geom_node_point(aes(color=V(S6shNcl)$Annotation, size=10+degree),alpha=0.75)+
  geom_node_text(aes(label=name),color="gray25",size=4)+
  scale_color_manual(limits=V(S6shNcl)$Annotation, values=V(S6shNcl)$An.color)+
  theme_void()

A1
#Modularity Analysis

#Edge_Betweeness Modularity

S6shNcl.ebmod=cluster_edge_betweenness(S6shNcl,weights=E(S6shNcl)$weight1,directed = F) 
S6shNcl.ebmod #IGRAPH clustering edge betweenness, groups: 37, mod: 0.58

#Plot the Edge betweness modularity
plot(S6shNcl.ebmod, S6shNcl,layout=Lyt.Tree, alpha=0.6, 
     vertex.color=V(S6shNcl)$An.color, vertex.frame.color=V(S6shNcl)$An.color, vertex.shape="circle", 
     vertex.size=5, vertex.label=V(S6shNcl)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black", 
     edge.width=E(S6shNcl)$weight1, edge.lty=1, main="S6ID Full Nuclear Network", sub="FR Tree Layout", 
     rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)
#Dendogram of the Edge_Betweeness
plot_dendrogram(S6shNcl.ebmod, mode="hclust")

#Louvain modularity
S6shNcl.Louvain.mod<-cluster_louvain(S6shNcl,weights=E(S6shNcl)$weight1)
S6shNcl.Louvain.mod #The result shows a total of IGRAPH clustering multi level, groups: 11, mod: 0.62

#Plot the Louvain modules
plot(S6shNcl.Louvain.mod, S6shNcl,layout=Lyt.Tree, alpha=0.6, 
     vertex.color=V(S6shNcl)$An.color, vertex.frame.color=V(S6shNcl)$An.color, vertex.shape="circle", 
     vertex.size=5, vertex.label=V(S6shNcl)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black", 
     edge.width=E(S6shNcl)$weight1, edge.lty=1, main="S6ID Full Nuclear Network", sub="FR Tree Layout", 
     rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#Generate HeatMaps of the adjacency matrixes
#Let's use the neighbor list for generating the Adjaency Matrix
SensorID_FNN
v<-levels(as.factor(unlist(SensorID_FNN[,1:2]))) #get vertex names from the network list
n<-length(v)#get the number of vertex in the network
e<-matrix(match(as.character(unlist(SensorID_FNN[,1:2])),v),ncol = 2)#Edgelist with row index instead of Names
w<-S6ID_S6shNcl$weight1 #Edge weights

# M is our co-expression network adjacency matrix.

# distringuish between i,j and j,i we want to make M symmetric (undirected) by
# considering an undirected network.
S6shNcl.M <- matrix(0, n, n)                  # set up a co-expression matrix
S6shNcl.M[e] <- w                             # fill it in with edge weights
S6shNcl.M <- S6shNcl.M + t(S6shNcl.M)                         # make this symmetric
dimnames(S6shNcl.M) <- list(v, v)             # label the vertices
# We let A_ij = 1 if there is an interaction. 
#A <- 1*(M > 0) #right now we keep the weights as they are
S6shNcl.M #The adjacency Matrix not as an igraph object, easier to manipulate, and to subset for smaller Sensor specific networks
quantile(S6shNcl.M,c(0,0.25,0.5,0.75,1))
#
#0%   25%   50%   75%  100% 
#0.000 0.000 0.000 0.000 0.991 

#Creating the HeatMap
# make a heatmap image with colors: white below 0,
# grey from 0.4 to 1 with darker closer to 1
breaks <- c(0, seq(0.05, 1, length=11)) # breaks in the color bins
cols <- grey(1-c(0,seq(0.5,1,length=10))) # colors

# color the vertex for CCNB1 blue, its neighbors red, and the others white
S6IDBind <- which(v == "SIRT6") # index for ccnb1
#S6IDBind<-as.factor(S6IDBind)
vcols <- rep("white",n)
vcols[S6IDBind] <- "blue"
#Subsetting S6shNcl.M by doing a loop asking for look into every row and column the index in S6IDBind worked
for (i in S6IDBind) {
  vcols[which(S6shNcl.M[,i] > 0 | S6shNcl.M[i,])] <- "red"
}
#Subsetting S6shNcl.M using the index vector in S6IDBind do not work
#vcols[which(S6shNcl.M[,S6IDBind] > 0 | S6shNcl.M[S6IDBind,])] <- "red"

#Generate a dendogram for keeping the proteins in the sample place
#S6shNcl.dist<-dist(S6shNcl.M,method = "euclidian")
#S6shNcl.Dend<-hclust(S6shNcl.dist, method = "ward.D2" )
#plot(S6shNcl.Dend)
#?hclust #methods= "ward.D", "ward.D2" this one works nice, "single", "complete", "average","mcquitty","median" (= WPGMC) or "centroid" (= UPGMC).

#Cut the dendogram
#cut.height <- 10 # try numbers from 4 to 7
#S6shNcl.hclusts<- cutree(S6shNcl.Dend, cut.height)
#plot(S6shNcl.hclusts)

# now actually make the heat map
par(mar=rep(0,4)) # make plot margins 0
heatmap(S6shNcl.M, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=vcols, RowSideColors=vcols,  col=0:1, frame=T)+
  legend("topleft", c("Neighbors(SIRT6, MRE11, KU80)","SIRT6, MRE11, XRCC5"), fill=c("red","blue"),  
         bty="n", inset=0, xpd=T,  border=F)
##Color annotation for functional annotation
MatSet<-subset(Mat,Mat$NAME%in%v)# select only those proteins that are cointained in the node list of our matrixes
match(v,MatSet$NAME) # Check that the indexes of the Annotation Dataframe and the node list are in the same position 
match(rownames(S6shNcl.M),MatSet$NAME) #Check that the row names of the Adj Matrix and the Annotation dataframe are in the same order
match(colnames(S6shNcl.M),MatSet$NAME) #Check that the column names of the Adj Matrix and the Annotation dataframe are in the same order
AnCol<-MatSet$An.color#Vector of colors
AnCol<-as.character(AnCol)#transform it into a character vector, so the colors can be read


#
plot.new()
legend("topleft", c("Cell Cycle","Chromatin Org","DNA Repair","Nuclear Pore","p53","Proteasome","Ribosome Biogenesis","RNA Metabolism","Transcription Reg","UPR","Calcium Signaling"), 
       fill=c("coral1","darkviolet","firebrick1","deepskyblue3","tomato","grey","lightgreen","gold","slateblue1","maroon","sienna2"),
       bty="n", inset=0, xpd=T,  border=F, cex =2,y.intersp = 0.7,x.intersp = 0.2,pt.cex = 1)

#check with the S6shNcl
# now actually make the heat map
par(mar=rep(0,4)) # make plot margins 0
heatmap(S6shNcl.M, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=0:1, frame=T, cexRow = 0.5 , cexCol = 0.5)

legend("topleft", c("Cell Cycle","Chromatin Org","DNA Repair","Nuclear Pore","p53","Proteasome","Ribosome Biogenesis","RNA Metabolism","Transcription Reg","UPR","Calcium Signaling"), 
       fill=c("coral1","darkviolet","firebrick1","deepskyblue3","tomato","grey","lightgreen","gold","slateblue1","maroon","sienna2"), bty="n", inset=0, xpd=T,  border=F, cex =2,y.intersp = 0.7,x.intersp = 0.2,pt.cex = 1)


#Now we substract the S6shScr AdjMatrix to the S6shNcl AdjMatrix
#lets add 1 to all the weights of the shNcl matrix
S6shNcl.M.Plus1<-S6shNcl.M+1
range(S6shNcl.M.Plus1)
S6shNcl.M.Plus1.Change<-S6shNcl.M.Plus1-S6shScr.M
S6shNcl.M.Change<-S6shNcl.M-S6shScr.M
range(S6shNcl.M.Change)#Let's check the magnitude of the substracted weights
#[1] -0.716  0.705
any(is.na(S6shNcl.M.Change))
#Generate the breaks and color functions
breaks <- c(seq(-0.5, 0.5, length=10))
breaks
HMCol<- brewer.pal(9,name = "RdBu")
HMCol<- rev(HMCol)
HMCol
[1] "#2166AC" "#4393C3" "#92C5DE" "#D1E5F0" "#F7F7F7" "#FDDBC7" "#F4A582" "#D6604D" "#B2182B"

#Plot the HeatMap
par(mar=rep(0,4)) # make plot margins 0
heatmap(S6shNcl.M.Change, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=HMCol, breaks=breaks,frame=T, cexRow = 0.5 , cexCol = 0.5, main = "Sirt6 shNcl (minus shScr edges)", ylab="Total Protein Interactors",xlab="Total Protein Interactors")


#
dev.off()
plot.new()



####Ku80 Networks
##
## Ku80ID shScr Networks: The Ku80ID shScr network is the Ku80ID IR2h after treatment with the shRNA control shScr

#KuID All nodes Network Edge Lists
KuID_KushScr<-read.delim("./Ku80ID_shScrlvsCneg_Network_EdgeList.tsv",sep = "\t",col.names = )
str(KuID_KushScr)
head(KuID_KushScr)

#Generate the igraph network object
KushScr<-graph_from_data_frame(d=KuID_KushScr,vertices = Mat, directed = F )#Generate the igraph object with the edge and node annotation attributes
str(vertex_attr(KushScr))#Check node attributes
str(edge_attr(KushScr))#Check for edge attributes
KushScr<-igraph::simplify(KushScr, remove.multiple=TRUE, remove.loops=TRUE,edge.attr.comb = "first") #Remove repetitive edges and self-connecting edges
KushScr
vcount(KushScr)

#Small plot
plot(KushScr)
#Generate an adjacency Matrix that will be used after
KushScr_AdjMat<-as_adjacency_matrix(KushScr,type="both",attr = "weight1",names=TRUE,sparse = TRUE)
KushScr_AdjMat
#Generate a vertex size based on degree and count
vsize<-igraph::degree(KushScr)/vcount(KushScr)
#Identify and remove disconnected nodes
isolated<-which(igraph::degree(KushScr)==0)
isolated
length(isolated)# [1] 463
#remove isolated nodes
KushScr<-igraph::delete_vertices(KushScr,isolated)

vcount(KushScr)
#Generate layouts for the graph
#Lyt.Nice<-layout_nicely(KushScr,dim = 3)
#Lyt.GOpt<-layout_with_graphopt(KushScr, niter=5000, mass = 30, charge = E(KushScr)$weight)
#Lyt.DrL<-layout_with_drl(KushScr, weights = E(KushScr)$weight)
#LytFR<-layout_with_fr(KushScr, niter=10000, weights=E(KushScr)$weight)#This one looks like the best option
#Lyt.KK<-layout_with_kk(KushScr, maxiter=10000, weights=E(KushScr)$weight) #very large separation of peripheral nodes
#Lyt.GEM<-layout_with_gem(KushScr, maxiter=10000) #Bad separation
#Lyt.DH<-layout_with_dh(KushScr, maxiter=150, fineiter=max(10, log2(833)))
#Lyt.LGL<-layout_with_lgl(KushScr, maxit=150,maxdelta=833,area=(length(E(KushScr)))^2,coolexp=1.5,repulserad=((length(E(KushScr)))^2)*(length(E(KushScr))))#Do not work


#Star and Tree layouts

V(KushScr)[name=="Ku80"]
which(V(KushScr)[name=="SIRT6"])# 610 #Check what is the rowindex for your desired node
which(V(KushScr)$name=="SIRT6")# 665 #Check what is the rowindex for your desired node
which(V(KushScr)$name=="MRE11A")# 104
which(V(KushScr)$name=="XRCC5")# 127
Lyt.Star<-layout_as_star(KushScr,center=V(KushScr)[name=="Ku80"]) #it has a special order argument where we can assign an order to the vertex, hence maybe the OPLS-DA:VIP value order may be used as ordering vector

Lyt.Tree<-layout_as_tree(KushScr,root=c(127))
Lyt.RT<-layout.reingold.tilford(KushScr,root=c(127))# This one looks good, but the indexing do not allow to highlight the correct nodes
Lyt.sugi<-layout_with_sugiyama(KushScr, maxiter=100, weights=E(KushScr)$weight)#This one looks good for undirected tree layouts

#Check what layout looks fine

plot.igraph(KushScr,layout=Lyt.Tree, alpha=0.6,
            vertex.color=V(KushScr)$An.color, vertex.frame.color=V(KushScr)$An.color, vertex.shape="circle",
            vertex.size=vsize*100, vertex.label=V(KushScr)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(KushScr)$weight1, edge.lty=1, main="KuID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)

#With ggraph

ggraph(KushScr,layout=Lyt.Tree)+
  geom_edge_fan(color="gray85", aes(width=(E(KushScr)$weight1)/3), alpha=0.4)+
  geom_node_point(aes(color=V(KushScr)$Annotation), size=10,alpha=0.75)+
  geom_node_text(aes(label=name),color="gray25",size=4)+
  scale_color_manual(limits=V(KushScr)$Annotation, values=V(KushScr)$An.color)+
  theme_void()

#Community Detection
clp <- igraph::cluster_optimal(KushScr)
class(clp)
V(KushScr)$community <- clp$membership
colrs <- adjustcolor( c("gray50", "tomato", "gold", "yellowgreen"), alpha=.6)


plot.igraph(clp,KushScr,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(KushScr)$An.color, vertex.frame.color=V(KushScr)$An.color, vertex.shape="circle",
            vertex.size=5, vertex.label=V(KushScr)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(KushScr)$weight, edge.lty=1, main="KuID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#Calculate the distance from every node to Ku80
dist.from.Ku80 <- distances(KushScr, v=V(KushScr)[name=="XRCC5"],
                            to=V(KushScr), weights=E(KushScr)$weight1)
dist.from.Ku80<-dist.from.Ku80[,!is.infinite(colSums(dist.from.Ku80))]

dist.from.Ku80
range(dist.from.Ku80)

#Set colors to plot the distances:
oranges <- colorRampPalette(c("dark red", "gold"))
col <- oranges(max(dist.from.Ku80))
col <- col[dist.from.Ku80+1]
#Plot the distances in the network
plot(KushScr, layout=Lyt.Tree,vertex.color=col, vertex.label=dist.from.Ku80,
     vertex.label.color="white",vertex.shape="circle",
     vertex.size=5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(KushScr)$weight1, edge.lty=1, main="KuID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)


#Shortest path from Ku80 to Nucleolin
KuNcl.path <- shortest_paths(KushScr,
                             from = V(KushScr)[name=="XRCC5"],
                             to = V(KushScr)[name=="NCL"],
                             output = "both") # both path nodes and edges
# Generate edge color variable to plot the path:
ecol <- rep("gray80", ecount(KushScr))
ecol[unlist(KuNcl.path$epath)] <- "orange"
# Generate edge width variable to plot the path:
ew <- rep(2, ecount(KushScr))
ew[unlist(KuNcl.path$epath)] <- 4
# Generate node color variable to plot the path:
vcol <- rep("gray40", vcount(KushScr))
vcol[unlist(KuNcl.path$vpath)] <- "gold"
#Plot the path from Ku80 to Nucleolin
plot(KushScr,layout=Lyt.Tree, vertex.color=vcol, edge.color=ecol,
     edge.width=ew,vertex.label.color="white",vertex.shape="circle",
     vertex.size=5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(KushScr)$weight1, edge.lty=1, main="KuID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)

#Identify the Ku80 incident edges and neighbor nodes
inc.edges <- incident(KushScr, V(KushScr)[name=="XRCC5"], mode="all")
inc.edges
# Set colors to plot the selected edges.
ecol <- rep("gray80", ecount(KushScr))
ecol[inc.edges] <- "orange"
vcol <- V(KushScr)$An.color
vcol[V(KushScr)$name=="XRCC5"] <- "gold"

# Plot the Ku80 incident edges and neighbor nodes
plot(KushScr,layout=Lyt.Tree, vertex.color=vcol, edge.color=ecol, layout=Lyt.FR[-isolated,],
     edge.width=ew,vertex.label.color="white",vertex.shape="circle",
     vertex.size=5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(KushScr)$weight1, edge.lty=1, main="KuID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)

#Plot Degree Distribution
KushScr.Deg<-degree_distribution(KushScr, cumulative = T, mode="all")
plot( x=0:max(igraph::degree(KushScr)), y=1-KushScr.Deg, pch=19, cex=1.2, col="orange",
      xlab="Degree", ylab="Cumulative Frequency")

#Network-level Descriptive metrics
#Number of nodes and edges
vcount(KushScr) #793
ecount(KushScr)#9309
#Density (Number of edges in relation to the total number of edges in a full connected graph)
edge_density(KushScr,loops = F) #0.02964385
#Transitivity (Total number of conected triangles)
transitivity(KushScr,type="global") #0.5731592
#modularity
vertex_attr_names(KushScr)
#Modularity Measurements
#Annotation
mod<-V(KushScr)$An.number #lets check for Annotation modularity
mod<-as.numeric(mod)
mod<-mod+1#The non clustered interactors that appears as 0, needs to be transformed to 1, so we add 1 to every cluster
modularity(KushScr,membership = mod)
#Annotation modularity [1] 0.2314044
#Nuclear compartment
mod<-V(KushScr)$CellCompartment #lets check for Nuclear compartment modularity
mod<-as.factor(mod)
mod<-as.numeric(mod)
mod
mod<-mod+1#The non clustered interactors that appears as 0, needs to be transformed to 1, so we add 1 to every cluster
modularity(KushScr,membership = mod)
#NucLoc modularity [1] [1] 0.0818481
#DNA_Repair
mod<-V(KushScr)$DNA_Repair #lets check for DNA_Repair modularity
mod<-as.factor(mod)
mod<-as.numeric(mod)
mod
mod<-mod+1#The non clustered interactors that appears as 0, needs to be transformed to 1, so we add 1 to every cluster
modularity(KushScr,membership = mod)
#NucLoc modularity [1] [1] 0.02854773


##Assorativity Measurements
#Assorativity
assortativity_nominal(KushScr,factor(V(KushScr)$Annotation),directed=F)
#Annotation [1] 0.3059801
assortativity_nominal(KushScr,factor(V(KushScr)$CellCompartment),directed=F)
#CellCompartment [1]  0.1308041
assortativity_nominal(KushScr,factor(V(KushScr)$DNA_Repair),directed=F)
#DNA_Repair [1] 0.05079576


#Diameter (longest geodesic distance), shortest path between two most peripheral nodes
diameter(KushScr,directed = F, weights = E(KushScr)$weight1)#[1] 4.998
diam<-get_diameter(KushScr,directed = F,weights = E(KushScr)$weight1)
diam #+ + 10/793 vertices, named, from 3bdaf56:
#[1] ARPC4  ACTR2  CORO1C PPP1CC NCL    RPL13A DAP3   TFB1M  POLRMT TFB2M 
vcol <- rep("gray50", vcount(KushScr))
vcol[diam] <- "gold"
ecol <- rep("gray90", ecount(KushScr))
ecol[E(KushScr, path=diam)] <- "orange"

# Plot the diameter of the network
plot(KushScr, layout=Lyt.Tree, vertex.color=vcol, edge.color=ecol, vertex.label.color="black",edge.arrow.mode=0)

#Node-level metrics
degree<-igraph::centr_degree(KushScr,mode = "all",normalized = T)
degree<-degree$res
degree
eigen_cent<-igraph::eigen_centrality(KushScr)
eigen_cent<-eigen_cent$vector
clos_cent <- igraph::closeness(KushScr, normalized=T)
betw_cent <- igraph::betweenness(KushScr, normalized=T)
hits <- hits_scores(KushScr, weights=E(KushScr)$weight)
hs <-hits$hub
authority<-hits$authority
Trans<-transitivity(KushScr,type="local")
Trans
#write results
KushScr.Centrality<-data.frame(NAME=V(KushScr)$name, UniprotID=V(KushScr)$UniprotID, Degree=degree, Eigen= eigen_cent, Closeness= clos_cent, Betweeness= betw_cent, Hubs=hs, Authority=authority,Transitivity=Trans)
KushScr.Centrality

write.csv(KushScr.Centrality,"./KushScr.Centrality.csv")

#Cocitation, Number of shared interactions among nodes
Cocitation<-cocitation(KushScr)
Cocitation#This matrix is similar to an adjacency matrix

#Plot with degree as node size
A1<-ggraph(KushScr,layout=Lyt.Tree)+
  geom_edge_fan(color="gray85", aes(width=(E(KushScr)$weight1)/3), alpha=0.4)+
  geom_node_point(aes(color=V(KushScr)$Annotation, size=10+degree),alpha=0.75)+
  geom_node_text(aes(label=name),color="gray25",size=4)+
  scale_color_manual(limits=V(KushScr)$Annotation, values=V(KushScr)$An.color)+
  theme_void()

A1
#Modularity Analysis

#Edge_Betweeness Modularity

KushScr.ebmod=cluster_edge_betweenness(KushScr,weights=E(KushScr)$weight1,directed = F) 
KushScr.ebmod #IGRAPH clustering edge betweenness, groups: 121, mod: 0.47

# Plot the Edge_Betweeness modules
plot(KushScr.ebmod, KushScr,layout=Lyt.Tree, alpha=0.6, 
     vertex.color=V(KushScr)$An.color, vertex.frame.color=V(KushScr)$An.color, vertex.shape="circle", 
     vertex.size=5, vertex.label=V(KushScr)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black", 
     edge.width=E(KushScr)$weight1, edge.lty=1, main="KuID Full Nuclear Network", sub="FR Tree Layout", 
     rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)
#Dendogram of the Edge_Betweeness
plot_dendrogram(KushScr.ebmod, mode="hclust")

#Louvain modularity
KushScr.Louvain.mod<-cluster_louvain(KushScr,weights=E(KushScr)$weight1)
KushScr.Louvain.mod #The result shows a total of IGRAPH clustering multi level, groups: 19, mod: 0.55

#Plot the Louvain modules
plot(KushScr.Louvain.mod, KushScr,layout=Lyt.Tree, alpha=0.6, 
     vertex.color=V(KushScr)$An.color, vertex.frame.color=V(KushScr)$An.color, vertex.shape="circle", 
     vertex.size=5, vertex.label=V(KushScr)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black", 
     edge.width=E(KushScr)$weight1, edge.lty=1, main="KuID Full Nuclear Network", sub="FR Tree Layout", 
     rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#Generate HeatMaps of the adjacency matrixes
#Let's use the neighbor list for generating the Adjaency Matrix
SensorID_FNN
v<-levels(as.factor(unlist(SensorID_FNN[,1:2]))) #get vertex names from the network list
n<-length(v)#get the number of vertex in the network
e<-matrix(match(as.character(unlist(SensorID_FNN[,1:2])),v),ncol = 2)#Edgelist with row index instead of Names
w<-KuID_KushScr$weight1 #Edge weights

# M is our co-expression network adjacency matrix.

# distringuish between i,j and j,i we want to make M symmetric (undirected) by
# considering an undirected network.
KushScr.M <- matrix(0, n, n)                  # set up a co-expression matrix
KushScr.M[e] <- w                             # fill it in with edge weights
KushScr.M <- KushScr.M + t(KushScr.M)                         # make this symmetric
dimnames(KushScr.M) <- list(v, v)             # label the vertices
# We let A_ij = 1 if there is an interaction. 
#A <- 1*(M > 0) #right now we keep the weights as they are
KushScr.M #The adjacency Matrix not as an igraph object, easier to manipulate, and to subset for smaller Sensor specific networks
quantile(KushScr.M,c(0,0.25,0.5,0.75,1))
#
#0%   25%   50%   75%  100% 
#0.000 0.000 0.000 0.000 0.991 

#Creating the HeatMap
# make a heatmap image with colors: white below 0,
# grey from 0.4 to 1 with darker closer to 1
breaks <- c(0, seq(0.05, 1, length=11)) # breaks in the color bins
cols <- grey(1-c(0,seq(0.1,1,length=10))) # colors

# color the vertex for CCNB1 blue, its neighbors red, and the others white
KuIDBind <- which(v == "XRCC5") # index for ccnb1
#KuIDBind<-as.factor(KuIDBind)
vcols <- rep("white",n)
vcols[KuIDBind] <- "blue"
#Subsetting KushScr.M by doing a loop asking for look into every row and column the index in KuIDBind worked
for (i in KuIDBind) {
  vcols[which(KushScr.M[,i] > 0 | KushScr.M[i,])] <- "red"
}
#Subsetting KushScr.M using the index vector in KuIDBind do not work
#vcols[which(KushScr.M[,KuIDBind] > 0 | KushScr.M[KuIDBind,])] <- "red"

#Generate a dendogram for keeping the proteins in the sample place
#KushScr.dist<-dist(KushScr.M,method = "euclidian")
#KushScr.Dend<-hclust(KushScr.dist, method = "ward.D2" )
#plot(KushScr.Dend)
#?hclust #methods= "ward.D", "ward.D2" this one works nice, "single", "complete", "average","mcquitty","median" (= WPGMC) or "centroid" (= UPGMC).

#Cut the dendogram
#cut.height <- 10 # try numbers from 4 to 7
#KushScr.hclusts<- cutree(KushScr.Dend, cut.height)
#plot(KushScr.hclusts)

# now actually make the heat map
par(mar=rep(0,4)) # make plot margins 0
heatmap(KushScr.M, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=vcols, RowSideColors=vcols,  col=0:1, frame=T)+
  legend("topleft", c("Neighbors(Ku80, MRE11, KU80)","Ku80, MRE11, XRCC5"), fill=c("red","blue"),  
         bty="n", inset=0, xpd=T,  border=F)
##Color annotation for functional annotation
MatSet<-subset(Mat,Mat$NAME%in%v)# select only those proteins that are cointained in the node list of our matrixes
match(v,MatSet$NAME) # Check that the indexes of the Annotation Dataframe and the node list are in the same position 
match(rownames(KushScr.M),MatSet$NAME) #Check that the row names of the Adj Matrix and the Annotation dataframe are in the same order
match(colnames(KushScr.M),MatSet$NAME) #Check that the column names of the Adj Matrix and the Annotation dataframe are in the same order
AnCol<-MatSet$An.color#Vector of colors
AnCol<-as.character(AnCol)#transform it into a character vector, so the colors can be read


#
plot.new()
legend("topleft", c("Cell Cycle","Chromatin Org","DNA Repair","Nuclear Pore","p53","Proteasome","Ribosome Biogenesis","RNA Metabolism","Transcription Reg","UPR","Calcium Signaling"), 
       fill=c("coral1","darkviolet","firebrick1","deepskyblue3","tomato","grey","lightgreen","gold","slateblue1","maroon","sienna2"),
       bty="n", inset=0, xpd=T,  border=F, cex =2,y.intersp = 0.7,x.intersp = 0.2,pt.cex = 1)

#check with the KushScr
# now actually make the heat map
par(mar=rep(0,4)) # make plot margins 0
heatmap(KushScr.M, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=0:1, frame=T, cexRow = 0.5 , cexCol = 0.5)

legend("topleft", c("Cell Cycle","Chromatin Org","DNA Repair","Nuclear Pore","p53","Proteasome","Ribosome Biogenesis","RNA Metabolism","Transcription Reg","UPR","Calcium Signaling"), 
       fill=c("coral1","darkviolet","firebrick1","deepskyblue3","tomato","grey","lightgreen","gold","slateblue1","maroon","sienna2"), bty="n", inset=0, xpd=T,  border=F, cex =2,y.intersp = 0.7,x.intersp = 0.2,pt.cex = 1)





###############
###Ku80 shNcl Networks: Ku80ID IR2h interactome after treatment with the shRNA targeting Nucleolin (shNcl)

#KuID All nodes Network Edge Lists
KuID_KushNcl<-read.delim("./Ku80ID_shNclvsCneg_Network_EdgeList.tsv",sep = "\t",col.names = )
str(KuID_KushNcl)
head(KuID_KushNcl)

#Generate the igraph network object
KushNcl<-graph_from_data_frame(d=KuID_KushNcl,vertices = Mat, directed = F )#Generate the igraph object with the edge and node annotation attributes
str(vertex_attr(KushNcl))#Check node attributes
str(edge_attr(KushNcl))#Check for edge attributes
KushNcl<-igraph::simplify(KushNcl, remove.multiple=TRUE, remove.loops=TRUE,edge.attr.comb = "first") #Remove repetitive edges and self-connecting edges
KushNcl
vcount(KushNcl)

#Small plot
plot(KushNcl)
#Generate an adjacency Matrix that will be used after
KushNcl_AdjMat<-as_adjacency_matrix(KushNcl,type="both",attr = "weight1",names=TRUE,sparse = TRUE)
KushNcl_AdjMat
#Generate a vertex size based on degree and count
vsize<-igraph::degree(KushNcl)/vcount(KushNcl)
#Identify and remove disconnected nodes
isolated<-which(igraph::degree(KushNcl)==0)
isolated
length(isolated)# 923
#remove isolated nodes
KushNcl<-igraph::delete_vertices(KushNcl,isolated)

igraph::vcount(KushNcl)

#Generate layouts for the graph
#Lyt.Nice<-layout_nicely(KushNcl,dim = 3)
#Lyt.GOpt<-layout_with_graphopt(KushNcl, niter=5000, mass = 30, charge = E(KushNcl)$weight)
#Lyt.DrL<-layout_with_drl(KushNcl, weights = E(KushNcl)$weight)
#LytFR<-layout_with_fr(KushNcl, niter=10000, weights=E(KushNcl)$weight)#This one looks like the best option
#Lyt.KK<-layout_with_kk(KushNcl, maxiter=10000, weights=E(KushNcl)$weight) #very large separation of peripheral nodes
#Lyt.GEM<-layout_with_gem(KushNcl, maxiter=10000) #Bad separation
#Lyt.DH<-layout_with_dh(KushNcl, maxiter=150, fineiter=max(10, log2(833)))
#Lyt.LGL<-layout_with_lgl(KushNcl, maxit=150,maxdelta=833,area=(length(E(KushNcl)))^2,coolexp=1.5,repulserad=((length(E(KushNcl)))^2)*(length(E(KushNcl))))#Do not work


#Star and Tree layouts

V(KushNcl)[name=="XRCC5"]
which(V(KushNcl)[name=="XRCC5"])# 558 #Check what is the rowindex for your desired node
which(V(KushNcl)$name=="SIRT6")# 558 #Check what is the rowindex for your desired node
which(V(KushNcl)$name=="MRE11A")# 115
which(V(KushNcl)$name=="XRCC5")# 83
Lyt.Star<-layout_as_star(KushNcl,center=V(KushNcl)[name=="Ku80"]) #it has a special order argument where we can assign an order to the vertex, hence maybe the OPLS-DA:VIP value order may be used as ordering vector

Lyt.Tree<-layout_as_tree(KushNcl,root=c(83))
Lyt.RT<-layout.reingold.tilford(KushNcl,root=c(83))# This one looks good, but the indexing do not allow to highlight the correct nodes
Lyt.sugi<-layout_with_sugiyama(KushNcl, maxiter=100, weights=E(KushNcl)$weight)#This one looks good for undirected tree layouts

#Check what layout looks fine

plot.igraph(KushNcl,layout=Lyt.Tree, alpha=0.6,
            vertex.color=V(KushNcl)$An.color, vertex.frame.color=V(KushNcl)$An.color, vertex.shape="circle",
            vertex.size=vsize*100, vertex.label=V(KushNcl)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(KushNcl)$weight1, edge.lty=1, main="KuID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)
#
legend( x="bottomright", levels(Mat$Annotation), pch=21,
        col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)

#With ggraph

ggraph(KushNcl,layout=Lyt.Tree)+
  geom_edge_fan(color="gray85", aes(width=(E(KushNcl)$weight1)/3), alpha=0.4)+
  geom_node_point(aes(color=V(KushNcl)$Annotation), size=10,alpha=0.75)+
  geom_node_text(aes(label=name),color="gray25",size=4)+
  scale_color_manual(limits=V(KushNcl)$Annotation, values=V(KushNcl)$An.color)+
  theme_void()

#Community Detection
clp <- igraph::cluster_optimal(KushNcl)
class(clp)
V(KushNcl)$community <- clp$membership
colrs <- adjustcolor( c("gray50", "tomato", "gold", "yellowgreen"), alpha=.6)

#Plot the clique detected communities
plot.igraph(clp,KushNcl,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(KushNcl)$An.color, vertex.frame.color=V(KushNcl)$An.color, vertex.shape="circle",
            vertex.size=5, vertex.label=V(KushNcl)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(KushNcl)$weight, edge.lty=1, main="KuID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#Calculate the distance from Ku80 to any other node in the network
dist.from.Ku80 <- distances(KushNcl, v=V(KushNcl)[name=="XRCC5"],
                            to=V(KushNcl), weights=E(KushNcl)$weight1)
dist.from.Ku80<-dist.from.Ku80[,!is.infinite(colSums(dist.from.Ku80))]

dist.from.Ku80
range(dist.from.Ku80)

#Set colors to plot the distances:
oranges <- colorRampPalette(c("dark red", "gold"))
col <- oranges(max(dist.from.Ku80))
col <- col[dist.from.Ku80+1]
#Plot the distances from Ku80 in the network
plot(KushNcl, layout=Lyt.Tree,vertex.color=col, vertex.label=dist.from.Ku80,
     vertex.label.color="white",vertex.shape="circle",
     vertex.size=5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(KushNcl)$weight1, edge.lty=1, main="KuID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)


#Shortest path from Ku80 to Nucleolin
KuNcl.path <- shortest_paths(KushNcl,
                             from = V(KushNcl)[name=="XRCC5"],
                             to = V(KushNcl)[name=="NCL"],
                             output = "both") # both path nodes and edges
# Generate edge color variable to plot the path:
ecol <- rep("gray80", ecount(KushNcl))
ecol[unlist(KuNcl.path$epath)] <- "orange"
# Generate edge width variable to plot the path:
ew <- rep(2, ecount(KushNcl))
ew[unlist(KuNcl.path$epath)] <- 4
# Generate node color variable to plot the path:
vcol <- rep("gray40", vcount(KushNcl))
vcol[unlist(KuNcl.path$vpath)] <- "gold"
#Plot the path between Ku80 and Nucleolin in the network
plot(KushNcl,layout=Lyt.Tree, vertex.color=vcol, edge.color=ecol,
     edge.width=ew,vertex.label.color="white",vertex.shape="circle",
     vertex.size=5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(KushNcl)$weight1, edge.lty=1, main="KuID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)

#Identify the Ku80 incident edges and neighbor nodes
inc.edges <- incident(KushNcl, V(KushNcl)[name=="XRCC5"], mode="all")
inc.edges
# Set colors to plot the selected edges.
ecol <- rep("gray80", ecount(KushNcl))
ecol[inc.edges] <- "orange"
vcol <- V(KushNcl)$An.color
vcol[V(KushNcl)$name=="XRCC5"] <- "gold"
#Plot the Ku80 incident edges and neighbor nodes
plot(KushNcl,layout=Lyt.Tree, vertex.color=vcol, edge.color=ecol, layout=Lyt.FR[-isolated,],
     edge.width=ew,vertex.label.color="white",vertex.shape="circle",
     vertex.size=5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(KushNcl)$weight1, edge.lty=1, main="KuID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)

#Plot Degree Distribution
KushNcl.Deg<-degree_distribution(KushNcl, cumulative = T, mode="all")
plot( x=0:max(igraph::degree(KushNcl)), y=1-KushNcl.Deg, pch=19, cex=1.2, col="orange",
      xlab="Degree", ylab="Cumulative Frequency")

#Network=level Descriptive metrics
#Number of nodes and edges
vcount(KushNcl) #333
ecount(KushNcl)#2169
#Density (Number of edges in relation to the total number of edges in a full connected graph)
edge_density(KushNcl,loops = F) #0.03923803
#Transitivity (Total number of conected triangles)
transitivity(KushNcl,type="global") #[1]0.60805
#modularity
vertex_attr_names(KushNcl)
#Modularity Measurements
#Annotation
mod<-V(KushNcl)$An.number #lets check for Annotation modularity
mod<-as.numeric(mod)
mod<-mod+1#The non clustered interactors that appears as 0, needs to be transformed to 1, so we add 1 to every cluster
modularity(KushNcl,membership = mod)
#Annotation modularity [1] 0.2469959
#Nuclear compartment
mod<-V(KushNcl)$CellCompartment #lets check for Nuclear compartment modularity
mod<-as.factor(mod)
mod<-as.numeric(mod)
mod
mod<-mod+1#The non clustered interactors that appears as 0, needs to be transformed to 1, so we add 1 to every cluster
modularity(KushNcl,membership = mod)
#NucLoc modularity [1] [1] 0.04141173
#DNA_Repair
mod<-V(KushNcl)$DNA_Repair #lets check for DNA_Repair modularity
mod<-as.factor(mod)
mod<-as.numeric(mod)
mod
mod<-mod+1#The non clustered interactors that appears as 0, needs to be transformed to 1, so we add 1 to every cluster
modularity(KushNcl,membership = mod)
#DNA_Repair modularity [1] [1] 0.03024746


##Assorativity Measurements
#Assorativity
assortativity_nominal(KushNcl,factor(V(KushNcl)$Annotation),directed=F)
#Annotation [1]  0.3271191
assortativity_nominal(KushNcl,factor(V(KushNcl)$CellCompartment),directed=F)
#CellCompartment [1]  0.08459237
assortativity_nominal(KushNcl,factor(V(KushNcl)$DNA_Repair),directed=F)
#DNA_Repair [1]0.04960411


#Diameter (longest geodesic distance), shortest path between two most peripheral nodes
diameter(KushNcl,directed = F, weights = E(KushNcl)$weight1)#[1] 4.172
diam<-get_diameter(KushNcl,directed = F,weights = E(KushNcl)$weight1)
diam #+ 9/333 vertices, named, from 1eb80a0:
#[1] ARPC4  ACTR2  CFL1   EGFR   XRCC5  PHB    MT-CO2 UQCRC2 PTGES2

vcol <- rep("gray50", vcount(KushNcl))
vcol[diam] <- "gold"
ecol <- rep("gray90", ecount(KushNcl))
ecol[E(KushNcl, path=diam)] <- "orange"
# E(net, path=diam) finds edges along a path, here 'diam'
plot(KushNcl, layout=Lyt.Tree, vertex.color=vcol, edge.color=ecol, vertex.label.color="black",edge.arrow.mode=0)

#Node-level metrics
degree<-igraph::centr_degree(KushNcl,mode = "all",normalized = T)
degree<-degree$res
degree
eigen_cent<-igraph::eigen_centrality(KushNcl)
eigen_cent<-eigen_cent$vector
clos_cent <- igraph::closeness(KushNcl, normalized=T)
betw_cent <- igraph::betweenness(KushNcl, normalized=T)
hits <- hits_scores(KushNcl, weights=E(KushNcl)$weight)
hs <-hits$hub
authority<-hits$authority
Trans<-transitivity(KushNcl,type="local")
Trans
#write results
KushNcl.Centrality<-data.frame(NAME=V(KushNcl)$name, UniprotID=V(KushNcl)$UniprotID, Degree=degree, Eigen= eigen_cent, Closeness= clos_cent, Betweeness= betw_cent, Hubs=hs, Authority=authority,Transitivity=Trans)
KushNcl.Centrality

write.csv(KushNcl.Centrality,"./KushNcl.Centrality.csv")

#Cocitation, Number of shared interactions among nodes
Cocitation<-cocitation(KushNcl)
Cocitation#This matrix is similar to an adjacency matrix

#Plot with degree as node size
A1<-ggraph(KushNcl,layout=Lyt.Tree)+
  geom_edge_fan(color="gray85", aes(width=(E(KushNcl)$weight1)/3), alpha=0.4)+
  geom_node_point(aes(color=V(KushNcl)$Annotation, size=10+degree),alpha=0.75)+
  geom_node_text(aes(label=name),color="gray25",size=4)+
  scale_color_manual(limits=V(KushNcl)$Annotation, values=V(KushNcl)$An.color)+
  theme_void()

A1
#Modularity Analysis

#Edge_Betweeness Modularity

KushNcl.ebmod=cluster_edge_betweenness(KushNcl,weights=E(KushNcl)$weight1,directed = F) 
KushNcl.ebmod #IGRAPH clustering edge betweenness, groups: 47, mod: 0.44
#

vertex_attr_names(KushNcl)
# Plot the Edge_Betweeness Modules
plot(KushNcl.ebmod, KushNcl,layout=Lyt.Tree, alpha=0.6, 
     vertex.color=V(KushNcl)$An.color, vertex.frame.color=V(KushNcl)$An.color, vertex.shape="circle", 
     vertex.size=5, vertex.label=V(KushNcl)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black", 
     edge.width=E(KushNcl)$weight1, edge.lty=1, main="KuID Full Nuclear Network", sub="FR Tree Layout", 
     rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)
#Dendogram of the Edge_Betweeness
plot_dendrogram(KushNcl.ebmod, mode="hclust")

#Louvain modularity
KushNcl.Louvain.mod<-cluster_louvain(KushNcl,weights=E(KushNcl)$weight1)
KushNcl.Louvain.mod #The result shows a total of IGRAPH clustering multi level,groups: 14, mod: 0.5

#Plot the Louvain Modules
plot(KushNcl.Louvain.mod, KushNcl,layout=Lyt.Tree, alpha=0.6, 
     vertex.color=V(KushNcl)$An.color, vertex.frame.color=V(KushNcl)$An.color, vertex.shape="circle", 
     vertex.size=5, vertex.label=V(KushNcl)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black", 
     edge.width=E(KushNcl)$weight1, edge.lty=1, main="KuID Full Nuclear Network", sub="FR Tree Layout", 
     rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#Generate HeatMaps of the adjacency matrixes
#Let's use the neighbor list for generating the Adjaency Matrix
SensorID_FNN
v<-levels(as.factor(unlist(SensorID_FNN[,1:2]))) #get vertex names from the network list
n<-length(v)#get the number of vertex in the network
e<-matrix(match(as.character(unlist(SensorID_FNN[,1:2])),v),ncol = 2)#Edgelist with row index instead of Names
w<-KuID_KushNcl$weight1 #Edge weights

# M is our co-expression network adjacency matrix.

# distringuish between i,j and j,i we want to make M symmetric (undirected) by
# considering an undirected network.
KushNcl.M <- matrix(0, n, n)                  # set up a co-expression matrix
KushNcl.M[e] <- w                             # fill it in with edge weights
KushNcl.M <- KushNcl.M + t(KushNcl.M)                         # make this symmetric
dimnames(KushNcl.M) <- list(v, v)             # label the vertices
# We let A_ij = 1 if there is an interaction. 
#A <- 1*(M > 0) #right now we keep the weights as they are
KushNcl.M #The adjacency Matrix not as an igraph object, easier to manipulate, and to subset for smaller Sensor specific networks
quantile(KushNcl.M,c(0,0.25,0.5,0.75,1))
#
#0%   25%   50%   75%  100% 
#0.000 0.000 0.000 0.000 0.991 

#Creating the HeatMap
# make a heatmap image with colors: white below 0,
# grey from 0.4 to 1 with darker closer to 1
breaks <- c(0, seq(0.05, 1, length=11)) # breaks in the color bins
cols <- grey(1-c(0,seq(0.5,1,length=10))) # colors

# color the vertex for CCNB1 blue, its neighbors red, and the others white
KuIDBind <- which(v == "XRCC5") # index for ccnb1
#KuIDBind<-as.factor(KuIDBind)
vcols <- rep("white",n)
vcols[KuIDBind] <- "blue"
#Subsetting KushNcl.M by doing a loop asking for look into every row and column the index in KuIDBind worked
for (i in KuIDBind) {
  vcols[which(KushNcl.M[,i] > 0 | KushNcl.M[i,])] <- "red"
}
#Subsetting KushNcl.M using the index vector in KuIDBind do not work
#vcols[which(KushNcl.M[,KuIDBind] > 0 | KushNcl.M[KuIDBind,])] <- "red"

#Generate a dendogram for keeping the proteins in the sample place
#KushNcl.dist<-dist(KushNcl.M,method = "euclidian")
#KushNcl.Dend<-hclust(KushNcl.dist, method = "ward.D2" )
#plot(KushNcl.Dend)
#?hclust #methods= "ward.D", "ward.D2" this one works nice, "single", "complete", "average","mcquitty","median" (= WPGMC) or "centroid" (= UPGMC).

#Cut the dendogram
#cut.height <- 10 # try numbers from 4 to 7
#KushNcl.hclusts<- cutree(KushNcl.Dend, cut.height)
#plot(KushNcl.hclusts)

# now actually make the heat map
par(mar=rep(0,4)) # make plot margins 0
heatmap(KushNcl.M, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=vcols, RowSideColors=vcols,  col=0:1, frame=T)+
  legend("topleft", c("Neighbors(Ku80, MRE11, KU80)","Ku80, MRE11, XRCC5"), fill=c("red","blue"),  
         bty="n", inset=0, xpd=T,  border=F)
##Color annotation for functional annotation
MatSet<-subset(Mat,Mat$NAME%in%v)# select only those proteins that are cointained in the node list of our matrixes
match(v,MatSet$NAME) # Check that the indexes of the Annotation Dataframe and the node list are in the same position 
match(rownames(KushNcl.M),MatSet$NAME) #Check that the row names of the Adj Matrix and the Annotation dataframe are in the same order
match(colnames(KushNcl.M),MatSet$NAME) #Check that the column names of the Adj Matrix and the Annotation dataframe are in the same order
AnCol<-MatSet$An.color#Vector of colors
AnCol<-as.character(AnCol)#transform it into a character vector, so the colors can be read


#
plot.new()
legend("topleft", c("Cell Cycle","Chromatin Org","DNA Repair","Nuclear Pore","p53","Proteasome","Ribosome Biogenesis","RNA Metabolism","Transcription Reg","UPR","Calcium Signaling"), 
       fill=c("coral1","darkviolet","firebrick1","deepskyblue3","tomato","grey","lightgreen","gold","slateblue1","maroon","sienna2"),
       bty="n", inset=0, xpd=T,  border=F, cex =2,y.intersp = 0.7,x.intersp = 0.2,pt.cex = 1)

#check with the KushNcl
# now actually make the heat map
par(mar=rep(0,4)) # make plot margins 0
heatmap(KushNcl.M, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=0:1, frame=T, cexRow = 0.5 , cexCol = 0.5)

legend("topleft", c("Cell Cycle","Chromatin Org","DNA Repair","Nuclear Pore","p53","Proteasome","Ribosome Biogenesis","RNA Metabolism","Transcription Reg","UPR","Calcium Signaling"), 
       fill=c("coral1","darkviolet","firebrick1","deepskyblue3","tomato","grey","lightgreen","gold","slateblue1","maroon","sienna2"), bty="n", inset=0, xpd=T,  border=F, cex =2,y.intersp = 0.7,x.intersp = 0.2,pt.cex = 1)


#Now we substract the KushScr AdjMatrix to the KushNcl AdjMatrix
#lets add 1 to all the weights of the shNcl matrix
KushNcl.M.Plus1<-KushNcl.M+1
range(KushNcl.M.Plus1)
KushNcl.M.Plus1.Change<-KushNcl.M.Plus1-KushScr.M
KushNcl.M.Change<-KushNcl.M-KushScr.M
range(KushNcl.M.Change)#Let's check the magnitude of the substracted weights
#[1] -0.716  0.705
any(is.na(KushNcl.M.Change))
#Generate the breaks and color functions
breaks <- c(seq(-0.5, 0.5, length=10))
breaks
HMCol<- brewer.pal(9,name = "RdBu")
HMCol<- rev(HMCol)
HMCol
[1] "#2166AC" "#4393C3" "#92C5DE" "#D1E5F0" "#F7F7F7" "#FDDBC7" "#F4A582" "#D6604D" "#B2182B"

#Plot the HeatMap
par(mar=rep(0,4)) # make plot margins 0
heatmap(KushNcl.M.Change, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=HMCol, breaks=breaks,frame=T, cexRow = 0.5 , cexCol = 0.5, main = "Ku80 shNcl (minus shScr edges)", ylab="Total Protein Interactors",xlab="Total Protein Interactors")


#
dev.off()
plot.new()


#############################################  Mre11 Networks

### Mre11 shScr Networks: The Mre11ID IR2h interactome after treatment with the control shRNA (shScr)

#MreID All nodes Network Edge Lists
MreID_MreshScr<-read.delim("./Mre11ID_shScrlvsCneg_Network_EsgeList.tsv",sep = "\t",col.names = )
str(MreID_MreshScr)
head(MreID_MreshScr)

#Generate the igraph network object
MreshScr<-graph_from_data_frame(d=MreID_MreshScr,vertices = Mat, directed = F )#Generate the igraph object with the edge and node annotation attributes
str(vertex_attr(MreshScr))#Check node attributes
str(edge_attr(MreshScr))#Check for edge attributes
MreshScr<-igraph::simplify(MreshScr, remove.multiple=TRUE, remove.loops=TRUE,edge.attr.comb = "first") #Remove repetitive edges and self-connecting edges
MreshScr
vcount(MreshScr)

#Small plot
plot(MreshScr)
#Generate an adjacency Matrix that will be used after
MreshScr_AdjMat<-as_adjacency_matrix(MreshScr,type="both",attr = "weight1",names=TRUE,sparse = TRUE)
MreshScr_AdjMat
#Generate a vertex size based on degree and count
vsize<-igraph::degree(MreshScr)/vcount(MreshScr)
#Identify and remove disconnected nodes
isolated<-which(igraph::degree(MreshScr)==0)
isolated
length(isolated)# [1] 463
#remove isolated nodes
MreshScr<-igraph::delete_vertices(MreshScr,isolated)

vcount(MreshScr)
#Generate layouts for the graph
#Lyt.Nice<-layout_nicely(MreshScr,dim = 3)
#Lyt.GOpt<-layout_with_graphopt(MreshScr, niter=5000, mass = 30, charge = E(MreshScr)$weight)
#Lyt.DrL<-layout_with_drl(MreshScr, weights = E(MreshScr)$weight)
#LytFR<-layout_with_fr(MreshScr, niter=10000, weights=E(MreshScr)$weight)#This one looks like the best option
#Lyt.KK<-layout_with_kk(MreshScr, maxiter=10000, weights=E(MreshScr)$weight) #very large separation of peripheral nodes
#Lyt.GEM<-layout_with_gem(MreshScr, maxiter=10000) #Bad separation
#Lyt.DH<-layout_with_dh(MreshScr, maxiter=150, fineiter=max(10, log2(833)))
#Lyt.LGL<-layout_with_lgl(MreshScr, maxit=150,maxdelta=833,area=(length(E(MreshScr)))^2,coolexp=1.5,repulserad=((length(E(MreshScr)))^2)*(length(E(MreshScr))))#Do not work


#Star and Tree layouts
V(MreshScr)[name=="Mre11A"]
which(V(MreshScr)[name=="SIRT6"])# 610 #Check what is the rowindex for your desired node
which(V(MreshScr)$name=="SIRT6")# 665 #Check what is the rowindex for your desired node
which(V(MreshScr)$name=="MRE11A")# 80
which(V(MreshScr)$name=="XRCC5")# 89
Lyt.Star<-layout_as_star(MreshScr,center=V(MreshScr)[name=="MRE11A"]) #it has a special order argument where we can assign an order to the vertex, hence maybe the OPLS-DA:VIP value order may be used as ordering vector
Lyt.Tree<-layout_as_tree(MreshScr,root=c(80))
Lyt.RT<-layout.reingold.tilford(MreshScr,root=c(80))# This one looks good, but the indexing do not allow to highlight the correct nodes
Lyt.sugi<-layout_with_sugiyama(MreshScr, maxiter=100, weights=E(MreshScr)$weight1)#This one looks good for undirected tree layouts

#Check what layout looks fine

plot.igraph(MreshScr,layout=Lyt.Tree, alpha=0.6,
            vertex.color=V(MreshScr)$An.color, vertex.frame.color=V(MreshScr)$An.color, vertex.shape="circle",
            vertex.size=vsize*100, vertex.label=V(MreshScr)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(MreshScr)$weight1, edge.lty=1, main="MreID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)

legend( x="bottomright", levels(Mat$Annotation), pch=21,
        col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)

#With ggraph

ggraph(MreshScr,layout=Lyt.Tree)+
  geom_edge_fan(color="gray85", aes(width=(E(MreshScr)$weight1)/3), alpha=0.4)+
  geom_node_point(aes(color=V(MreshScr)$Annotation), size=10,alpha=0.75)+
  geom_node_text(aes(label=name),color="gray25",size=4)+
  scale_color_manual(limits=V(MreshScr)$Annotation, values=V(MreshScr)$An.color)+
  theme_void()

#Community Detection by the clique algoright in igraph
clp <- igraph::cluster_optimal(MreshScr)
class(clp)
V(MreshScr)$community <- clp$membership
colrs <- adjustcolor( c("gray50", "tomato", "gold", "yellowgreen"), alpha=.6)

#Plot the clique communities
plot.igraph(clp,MreshScr,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(MreshScr)$An.color, vertex.frame.color=V(MreshScr)$An.color, vertex.shape="circle",
            vertex.size=5, vertex.label=V(MreshScr)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(MreshScr)$weight, edge.lty=1, main="MreID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#Calculate the distances from Ku80 to every node in the network
dist.from.Mre11 <- distances(MreshScr, v=V(MreshScr)[name=="MRE11A"],
                             to=V(MreshScr), weights=E(MreshScr)$weight1)
dist.from.Mre11<-dist.from.Mre11[,!is.infinite(colSums(dist.from.Mre11))]

dist.from.Mre11
range(dist.from.Mre11)

# Set colors to plot the distances:
oranges <- colorRampPalette(c("dark red", "gold"))
col <- oranges(max(dist.from.Mre11))
col <- col[dist.from.Mre11+1]
# Plot the distances into the network
plot(MreshScr, layout=Lyt.Tree,vertex.color=col, vertex.label=dist.from.Mre11,
     vertex.label.color="white",vertex.shape="circle",
     vertex.size=5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(MreshScr)$weight1, edge.lty=1, main="MreID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)


# Shortest path from Mre11 to Nucleolin
MreNcl.path <- shortest_paths(MreshScr,
                              from = V(MreshScr)[name=="MRE11A"],
                              to = V(MreshScr)[name=="NCL"],
                              output = "both") # both path nodes and edges

# Generate edge color variable to plot the path:
ecol <- rep("gray80", ecount(MreshScr))
ecol[unlist(MreNcl.path$epath)] <- "orange"
# Generate edge width variable to plot the path:
ew <- rep(2, ecount(MreshScr))
ew[unlist(MreNcl.path$epath)] <- 4
# Generate node color variable to plot the path:
vcol <- rep("gray40", vcount(MreshScr))
vcol[unlist(MreNcl.path$vpath)] <- "gold"

#Plot the path between Mre11 and Ncl
plot(MreshScr,layout=Lyt.Tree, vertex.color=vcol, edge.color=ecol,
     edge.width=ew,vertex.label.color="white",vertex.shape="circle",
     vertex.size=5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(MreshScr)$weight1, edge.lty=1, main="MreID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)

#Identify the Mre11 incident edges and neighbor nodes
inc.edges <- incident(MreshScr, V(MreshScr)[name=="MRE11A"], mode="all")
inc.edges
# Set colors to plot the selected edges.
ecol <- rep("gray80", ecount(MreshScr))
ecol[inc.edges] <- "orange"
vcol <- V(MreshScr)$An.color
vcol[V(MreshScr)$name=="MRE11A"] <- "gold"

#Plot the Mre11 incident edges and neighbor nodes in the network
plot(MreshScr,layout=Lyt.Tree, vertex.color=vcol, edge.color=ecol, layout=Lyt.FR[-isolated,],
     edge.width=ew,vertex.label.color="white",vertex.shape="circle",
     vertex.size=5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(MreshScr)$weight1, edge.lty=1, main="MreID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)

#Plot Degree Distribution
MreshScr.Deg<-degree_distribution(MreshScr, cumulative = T, mode="all")
plot( x=0:max(igraph::degree(MreshScr)), y=1-MreshScr.Deg, pch=19, cex=1.2, col="orange",
      xlab="Degree", ylab="Cumulative Frequency")

#Network-level Descriptive metrics
#Number of nodes and edges
vcount(MreshScr) #325
ecount(MreshScr)#2405
#Density (Number of edges in relation to the total number of edges in a full connected graph)
edge_density(MreshScr,loops = F) #00.04567901
#Transitivity (Total number of conected triangles)
transitivity(MreshScr,type="global") #0.6381993
#modularity
vertex_attr_names(MreshScr)
#Modularity Measurements
#Annotation
mod<-V(MreshScr)$An.number #lets check for Annotation modularity
mod<-as.numeric(mod)
mod<-mod+1#The non clustered interactors that appears as 0, needs to be transformed to 1, so we add 1 to every cluster
modularity(MreshScr,membership = mod)
#Annotation modularity [1] 0.2471295
#Nuclear compartment
mod<-V(MreshScr)$CellCompartment #lets check for Nuclear compartment modularity
mod<-as.factor(mod)
mod<-as.numeric(mod)
mod
mod<-mod+1#The non clustered interactors that appears as 0, needs to be transformed to 1, so we add 1 to every cluster
modularity(MreshScr,membership = mod)
#NucLoc modularity [1] [1] 0.06386081
#DNA_Repair
mod<-V(MreshScr)$DNA_Repair #lets check for DNA_Repair pathway modularity
mod<-as.factor(mod)
mod<-as.numeric(mod)
mod
mod<-mod+1#The non clustered interactors that appears as 0, needs to be transformed to 1, so we add 1 to every cluster
modularity(MreshScr,membership = mod)
#NucLoc modularity [1] [1] 0.04835223


##Assorativity Measurements
#Assorativity
assortativity_nominal(MreshScr,factor(V(MreshScr)$Annotation),directed=F)
#Annotation [1]  0.3194816
assortativity_nominal(MreshScr,factor(V(MreshScr)$CellCompartment),directed=F)
#CellCompartment [1]   0.1047945
assortativity_nominal(MreshScr,factor(V(MreshScr)$DNA_Repair),directed=F)
#DNA_Repair [1] 0.08193346


#Diameter (longest geodesic distance), shortest path between two most peripheral nodes
diameter(MreshScr,directed = F, weights = E(MreshScr)$weight1)#[1] 3.631
diam<-get_diameter(MreshScr,directed = F,weights = E(MreshScr)$weight1)
diam #+ 8/325 vertices, named, from c19dd0f:
#[1] SLC7A5 SLC3A2 EGFR   NCL    RPL13A MRPL21 DDX28  YARS2 
vcol <- rep("gray50", vcount(MreshScr))
vcol[diam] <- "gold"
ecol <- rep("gray90", ecount(MreshScr))
ecol[E(MreshScr, path=diam)] <- "orange"

#Plot the diameter into the Network
plot(MreshScr, layout=Lyt.Tree, vertex.color=vcol, edge.color=ecol, vertex.label.color="black",edge.arrow.mode=0)

#Node-level metrics
degree<-igraph::centr_degree(MreshScr,mode = "all",normalized = T)
degree<-degree$res
degree
eigen_cent<-igraph::eigen_centrality(MreshScr)
eigen_cent<-eigen_cent$vector
clos_cent <- igraph::closeness(MreshScr, normalized=T)
betw_cent <- igraph::betweenness(MreshScr, normalized=T)
hits <- hits_scores(MreshScr, weights=E(MreshScr)$weight)
hs <-hits$hub
authority<-hits$authority
Trans<-transitivity(MreshScr,type="local")
Trans
#write results
MreshScr.Centrality<-data.frame(NAME=V(MreshScr)$name, UniprotID=V(MreshScr)$UniprotID, Degree=degree, Eigen= eigen_cent, Closeness= clos_cent, Betweeness= betw_cent, Hubs=hs, Authority=authority,Transitivity=Trans)
MreshScr.Centrality

write.csv(MreshScr.Centrality,"./MreshScr.Centrality.csv")

#Cocitation, Number of shared interactions among nodes
Cocitation<-cocitation(MreshScr)
Cocitation#This matrix is similar to an adjacency matrix

#Plot with degree as node size
A1<-ggraph(MreshScr,layout=Lyt.Tree)+
  geom_edge_fan(color="gray85", aes(width=(E(MreshScr)$weight1)/3), alpha=0.4)+
  geom_node_point(aes(color=V(MreshScr)$Annotation, size=10+degree),alpha=0.75)+
  geom_node_text(aes(label=name),color="gray25",size=4)+
  scale_color_manual(limits=V(MreshScr)$Annotation, values=V(MreshScr)$An.color)+
  theme_void()

A1
#Modularity Analysis
#Edge_Betweeness Modularity

MreshScr.ebmod=cluster_edge_betweenness(MreshScr,weights=E(MreshScr)$weight1,directed = F) 
MreshScr.ebmod #IGRAPH clustering edge betweenness, groups: 36, mod: 0.46
#

vertex_attr_names(MreshScr)
# Plot the Edge_Betweeness Modules
plot(MreshScr.ebmod, MreshScr,layout=Lyt.Tree, alpha=0.6, 
     vertex.color=V(MreshScr)$An.color, vertex.frame.color=V(MreshScr)$An.color, vertex.shape="circle", 
     vertex.size=5, vertex.label=V(MreshScr)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black", 
     edge.width=E(MreshScr)$weight1, edge.lty=1, main="MreID Full Nuclear Network", sub="FR Tree Layout", 
     rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)
#Dendogram of the Edge_Betweeness
plot_dendrogram(MreshScr.ebmod, mode="hclust")

#Louvain modularity
MreshScr.Louvain.mod<-cluster_louvain(MreshScr,weights=E(MreshScr)$weight1)
MreshScr.Louvain.mod #The result shows a total of IIGRAPH clustering multi level, groups: 10, mod: 0.51

#Plot the Louvain Modules
plot(MreshScr.Louvain.mod, MreshScr,layout=Lyt.Tree, alpha=0.6, 
     vertex.color=V(MreshScr)$An.color, vertex.frame.color=V(MreshScr)$An.color, vertex.shape="circle", 
     vertex.size=5, vertex.label=V(MreshScr)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black", 
     edge.width=E(MreshScr)$weight1, edge.lty=1, main="MreID Full Nuclear Network", sub="FR Tree Layout", 
     rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#Generate HeatMaps of the adjacency matrixes
#Let's use the neighbor list for generating the Adjaency Matrix
SensorID_FNN
v<-levels(as.factor(unlist(SensorID_FNN[,1:2]))) #get vertex names from the network list
n<-length(v)#get the number of vertex in the network
e<-matrix(match(as.character(unlist(SensorID_FNN[,1:2])),v),ncol = 2)#Edgelist with row index instead of Names
w<-MreID_MreshScr$weight1 #Edge weights

# M is our co-expression network adjacency matrix.

# distringuish between i,j and j,i we want to make M symmetric (undirected) by
# considering an undirected network.
MreshScr.M <- matrix(0, n, n)                  # set up a co-expression matrix
MreshScr.M[e] <- w                             # fill it in with edge weights
MreshScr.M <- MreshScr.M + t(MreshScr.M)                         # make this symmetric
dimnames(MreshScr.M) <- list(v, v)             # label the vertices
# We let A_ij = 1 if there is an interaction. 
#A <- 1*(M > 0) #right now we keep the weights as they are
MreshScr.M #The adjacency Matrix not as an igraph object, easier to manipulate, and to subset for smaller Sensor specific networks
quantile(MreshScr.M,c(0,0.25,0.5,0.75,1))
#
#0%   25%   50%   75%  100% 
#0.000 0.000 0.000 0.000 0.991 

#Creating the HeatMap
# make a heatmap image with colors: white below 0,
# grey from 0.4 to 1 with darker closer to 1
breaks <- c(0, seq(0.05, 1, length=11)) # breaks in the color bins
cols <- grey(1-c(0,seq(0.1,1,length=10))) # colors

# color the vertex for CCNB1 blue, its neighbors red, and the others white
MreIDBind <- which(v == "MRE11A") # index for ccnb1
#MreIDBind<-as.factor(MreIDBind)
vcols <- rep("white",n)
vcols[MreIDBind] <- "blue"
#Subsetting MreshScr.M by doing a loop asking for look into every row and column the index in MreIDBind worked
for (i in MreIDBind) {
  vcols[which(MreshScr.M[,i] > 0 | MreshScr.M[i,])] <- "red"
}
#Subsetting MreshScr.M using the index vector in MreIDBind do not work
#vcols[which(MreshScr.M[,MreIDBind] > 0 | MreshScr.M[MreIDBind,])] <- "red"

#Generate a dendogram for keeping the proteins in the sample place
#MreshScr.dist<-dist(MreshScr.M,method = "euclidian")
#MreshScr.Dend<-hclust(MreshScr.dist, method = "ward.D2" )
#plot(MreshScr.Dend)
#?hclust #methods= "ward.D", "ward.D2" this one works nice, "single", "complete", "average","mcquitty","median" (= WPGMC) or "centroid" (= UPGMC).

#Cut the dendogram
#cut.height <- 10 # try numbers from 4 to 7
#MreshScr.hclusts<- cutree(MreshScr.Dend, cut.height)
#plot(MreshScr.hclusts)

# now actually make the heat map
par(mar=rep(0,4)) # make plot margins 0
heatmap(MreshScr.M, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=vcols, RowSideColors=vcols,  col=0:1, frame=T)+
  legend("topleft", c("Neighbors(Mre11, MRE11, Mre11)","Mre11, MRE11, MRE11A"), fill=c("red","blue"),  
         bty="n", inset=0, xpd=T,  border=F)
##Color annotation for functional annotation
MatSet<-subset(Mat,Mat$NAME%in%v)# select only those proteins that are cointained in the node list of our matrixes
match(v,MatSet$NAME) # Check that the indexes of the Annotation Dataframe and the node list are in the same position 
match(rownames(MreshScr.M),MatSet$NAME) #Check that the row names of the Adj Matrix and the Annotation dataframe are in the same order
match(colnames(MreshScr.M),MatSet$NAME) #Check that the column names of the Adj Matrix and the Annotation dataframe are in the same order
AnCol<-MatSet$An.color#Vector of colors
AnCol<-as.character(AnCol)#transform it into a character vector, so the colors can be read


#
plot.new()
legend("topleft", c("Cell Cycle","Chromatin Org","DNA Repair","Nuclear Pore","p53","Proteasome","Ribosome Biogenesis","RNA Metabolism","Transcription Reg","UPR","Calcium Signaling"), 
       fill=c("coral1","darkviolet","firebrick1","deepskyblue3","tomato","grey","lightgreen","gold","slateblue1","maroon","sienna2"),
       bty="n", inset=0, xpd=T,  border=F, cex =2,y.intersp = 0.7,x.intersp = 0.2,pt.cex = 1)

#check with the MreshScr
# now actually make the heat map
par(mar=rep(0,4)) # make plot margins 0
heatmap(MreshScr.M, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=0:1, frame=T, cexRow = 0.5 , cexCol = 0.5)

legend("topleft", c("Cell Cycle","Chromatin Org","DNA Repair","Nuclear Pore","p53","Proteasome","Ribosome Biogenesis","RNA Metabolism","Transcription Reg","UPR","Calcium Signaling"), 
       fill=c("coral1","darkviolet","firebrick1","deepskyblue3","tomato","grey","lightgreen","gold","slateblue1","maroon","sienna2"), bty="n", inset=0, xpd=T,  border=F, cex =2,y.intersp = 0.7,x.intersp = 0.2,pt.cex = 1)






###Mre11 shNcl Networks: Mre11ID IR2h interactome after treatment with the shRNA targeting Nucleolin (shNcl)

#MreID All nodes Network Edge Lists
MreID_MreshNcl<-read.delim("./Mre11ID_shNclvsCneg_Network_EdgeList.tsv",sep = "\t",col.names = )
str(MreID_MreshNcl)
head(MreID_MreshNcl)

#Generate the igraph network object
MreshNcl<-graph_from_data_frame(d=MreID_MreshNcl,vertices = Mat, directed = F )#Generate the igraph object with the edge and node annotation attributes
str(vertex_attr(MreshNcl))#Check node attributes
str(edge_attr(MreshNcl))#Check for edge attributes
MreshNcl<-igraph::simplify(MreshNcl, remove.multiple=TRUE, remove.loops=TRUE,edge.attr.comb = "first") #Remove repetitive edges and self-connecting edges
MreshNcl
vcount(MreshNcl)

#Small plot
plot(MreshNcl)
#Generate an adjacency Matrix that will be used after
MreshNcl_AdjMat<-as_adjacency_matrix(MreshNcl,type="both",attr = "weight1",names=TRUE,sparse = TRUE)
MreshNcl_AdjMat
#Generate a vertex size based on degree and count
vsize<-igraph::degree(MreshNcl)/vcount(MreshNcl)
#Identify and remove disconnected nodes
isolated<-which(igraph::degree(MreshNcl)==0)
isolated
length(isolated)# 923
#remove isolated nodes
MreshNcl<-igraph::delete_vertices(MreshNcl,isolated)

igraph::vcount(MreshNcl)
#Generate layouts for the graph
#Lyt.Nice<-layout_nicely(MreshNcl,dim = 3)
#Lyt.GOpt<-layout_with_graphopt(MreshNcl, niter=5000, mass = 30, charge = E(MreshNcl)$weight)
#Lyt.DrL<-layout_with_drl(MreshNcl, weights = E(MreshNcl)$weight)
#LytFR<-layout_with_fr(MreshNcl, niter=10000, weights=E(MreshNcl)$weight)#This one looks like the best option
#Lyt.KK<-layout_with_kk(MreshNcl, maxiter=10000, weights=E(MreshNcl)$weight) #very large separation of peripheral nodes
#Lyt.GEM<-layout_with_gem(MreshNcl, maxiter=10000) #Bad separation
#Lyt.DH<-layout_with_dh(MreshNcl, maxiter=150, fineiter=max(10, log2(833)))
#Lyt.LGL<-layout_with_lgl(MreshNcl, maxit=150,maxdelta=833,area=(length(E(MreshNcl)))^2,coolexp=1.5,repulserad=((length(E(MreshNcl)))^2)*(length(E(MreshNcl))))#Do not work

#Star and Tree layouts
V(MreshNcl)[name=="MRE11A"]
which(V(MreshNcl)[name=="MRE11A"])# 558 #Check what is the rowindex for your desired node
which(V(MreshNcl)$name=="SIRT6")# 558 #Check what is the rowindex for your desired node
which(V(MreshNcl)$name=="MRE11A")# 115
which(V(MreshNcl)$name=="XRCC5")# 125
Lyt.Star<-layout_as_star(MreshNcl,center=V(MreshNcl)[name=="Mre11"]) #it has a special order argument where we can assign an order to the vertex, hence maybe the OPLS-DA:VIP value order may be used as ordering vector
Lyt.Tree<-layout_as_tree(MreshNcl,root=c(115))
Lyt.RT<-layout.reingold.tilford(MreshNcl,root=c(115))# This one looks good, but the indexing do not allow to highlight the correct nodes
Lyt.sugi<-layout_with_sugiyama(MreshNcl, maxiter=100, weights=E(MreshNcl)$weight1)#This one looks good for undirected tree layouts

#Check what layout looks fine

plot.igraph(MreshNcl,layout=Lyt.Tree, alpha=0.6,
            vertex.color=V(MreshNcl)$An.color, vertex.frame.color=V(MreshNcl)$An.color, vertex.shape="circle",
            vertex.size=vsize*100, vertex.label=V(MreshNcl)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(MreshNcl)$weight1, edge.lty=1, main="MreID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)
#
legend( x="bottomright", levels(Mat$Annotation), pch=21,
        col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)

#With ggraph

ggraph(MreshNcl,layout=Lyt.Tree)+
  geom_edge_fan(color="gray85", aes(width=(E(MreshNcl)$weight1)/3), alpha=0.4)+
  geom_node_point(aes(color=V(MreshNcl)$Annotation), size=10,alpha=0.75)+
  geom_node_text(aes(label=name),color="gray25",size=4)+
  scale_color_manual(limits=V(MreshNcl)$Annotation, values=V(MreshNcl)$An.color)+
  theme_void()

#igraph clique Community Detection
clp <- igraph::cluster_optimal(MreshNcl)
class(clp)
V(MreshNcl)$community <- clp$membership
colrs <- adjustcolor( c("gray50", "tomato", "gold", "yellowgreen"), alpha=.6)
#Plot the clique detected communities
plot.igraph(clp,MreshNcl,layout=Lyt.FR[-isolated,], alpha=0.6,
            vertex.color=V(MreshNcl)$An.color, vertex.frame.color=V(MreshNcl)$An.color, vertex.shape="circle",
            vertex.size=5, vertex.label=V(MreshNcl)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
            edge.width=E(MreshNcl)$weight, edge.lty=1, main="MreID Full Nuclear Network", sub="FR Tree Layout",
            rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#Calculate the distances between Mre11 and every node in the network
dist.from.Mre11 <- distances(MreshNcl, v=V(MreshNcl)[name=="MRE11A"],
                             to=V(MreshNcl), weights=E(MreshNcl)$weight1)
dist.from.Mre11<-dist.from.Mre11[,!is.infinite(colSums(dist.from.Mre11))]
dist.from.Mre11
range(dist.from.Mre11)

#Set colors to plot the distances:
oranges <- colorRampPalette(c("dark red", "gold"))
col <- oranges(max(dist.from.Mre11))
col <- col[dist.from.Mre11+1]
#Plot the distances into the network
plot(MreshNcl, layout=Lyt.Tree,vertex.color=col, vertex.label=dist.from.Mre11,
     vertex.label.color="white",vertex.shape="circle",
     vertex.size=5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(MreshNcl)$weight1, edge.lty=1, main="MreID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)


#Identify the shortest path from Mre11 to Nucleolin
MreNcl.path <- shortest_paths(MreshNcl,
                              from = V(MreshNcl)[name=="MRE11A"],
                              to = V(MreshNcl)[name=="NCL"],
                              output = "both") # both path nodes and edges
# Generate edge color variable to plot the path:
ecol <- rep("gray80", ecount(MreshNcl))
ecol[unlist(MreNcl.path$epath)] <- "orange"
# Generate edge width variable to plot the path:
ew <- rep(2, ecount(MreshNcl))
ew[unlist(MreNcl.path$epath)] <- 4
# Generate node color variable to plot the path:
vcol <- rep("gray40", vcount(MreshNcl))
vcol[unlist(MreNcl.path$vpath)] <- "gold"
#Plot the shortest path from Mre11 to Nucleolin
plot(MreshNcl,layout=Lyt.Tree, vertex.color=vcol, edge.color=ecol,
     edge.width=ew,vertex.label.color="white",vertex.shape="circle",
     vertex.size=5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(MreshNcl)$weight1, edge.lty=1, main="MreID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)

#Identify the Mre11 incident edges and neighbor nodes
inc.edges <- incident(MreshNcl, V(MreshNcl)[name=="MRE11A"], mode="all")
inc.edges
# Set colors to plot the selected edges.
ecol <- rep("gray80", ecount(MreshNcl))
ecol[inc.edges] <- "orange"
vcol <- V(MreshNcl)$An.color
vcol[V(MreshNcl)$name=="MRE11A"] <- "gold"
#Plot the Mre11 incident edges and neighbor nodes
plot(MreshNcl,layout=Lyt.Tree, vertex.color=vcol, edge.color=ecol, layout=Lyt.FR[-isolated,],
     edge.width=ew,vertex.label.color="white",vertex.shape="circle",
     vertex.size=5, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black",
     edge.width=E(MreshNcl)$weight1, edge.lty=1, main="MreID Full Nuclear Network", sub="FR Tree Layout",
     rescale=T)

#Plot Degree Distribution
MreshNcl.Deg<-degree_distribution(MreshNcl, cumulative = T, mode="all")
plot( x=0:max(igraph::degree(MreshNcl)), y=1-MreshNcl.Deg, pch=19, cex=1.2, col="orange",
      xlab="Degree", ylab="Cumulative Frequency")

#Network-level Descriptive metrics
#Number of nodes and edges
vcount(MreshNcl) #453
ecount(MreshNcl)#3567
#Density (Number of edges in relation to the total number of edges in a full connected graph)
edge_density(MreshNcl,loops = F) #0.03484147
#Transitivity (Total number of conected triangles)
transitivity(MreshNcl,type="global") #[1]0.6244292
#modularity
vertex_attr_names(MreshNcl)
#Modularity Measurements
#Annotation
mod<-V(MreshNcl)$An.number #lets check for Annotation modularity
mod<-as.numeric(mod)
mod<-mod+1#The non clustered interactors that appears as 0, needs to be transformed to 1, so we add 1 to every cluster
modularity(MreshNcl,membership = mod)
#Annotation modularity [1] 0.2548033
#Nuclear compartment
mod<-V(MreshNcl)$CellCompartment #lets check for Nuclear compartment modularity
mod<-as.factor(mod)
mod<-as.numeric(mod)
mod
mod<-mod+1#The non clustered interactors that appears as 0, needs to be transformed to 1, so we add 1 to every cluster
modularity(MreshNcl,membership = mod)
#NucLoc modularity [1] [1] 00.06548353
#DNA_Repair
mod<-V(MreshNcl)$DNA_Repair #lets check for DNA_Repair pathway modularity
mod<-as.factor(mod)
mod<-as.numeric(mod)
mod
mod<-mod+1#The non clustered interactors that appears as 0, needs to be transformed to 1, so we add 1 to every cluster
modularity(MreshNcl,membership = mod)
#DNA_Repair modularity [1] [1] 0.0243871


##Assorativity Measurements
#Assorativity
assortativity_nominal(MreshNcl,factor(V(MreshNcl)$Annotation),directed=F)
#Annotation [1]  0.342344
assortativity_nominal(MreshNcl,factor(V(MreshNcl)$CellCompartment),directed=F)
#CellCompartment [1]   0.116465
assortativity_nominal(MreshNcl,factor(V(MreshNcl)$DNA_Repair),directed=F)
#DNA_Repair [1]0.04483984


#Diameter (longest geodesic distance), shortest path between two most peripheral nodes
diameter(MreshNcl,directed = F, weights = E(MreshNcl)$weight1)#[1] 5.205
diam<-get_diameter(MreshNcl,directed = F,weights = E(MreshNcl)$weight1)
diam #+ 11/453 vertices, named, from 7d84b45:
#[1] PTGES2 UQCRC2 NDUFS3 GRSF1  MRPS27 RBM39  RPL27  PPP1CC WDR82  CPSF3  CSTF3 

vcol <- rep("gray50", vcount(MreshNcl))
vcol[diam] <- "gold"
ecol <- rep("gray90", ecount(MreshNcl))
ecol[E(MreshNcl, path=diam)] <- "orange"
#Plot the diameter into the network
plot(MreshNcl, layout=Lyt.Tree, vertex.color=vcol, edge.color=ecol, vertex.label.color="black",edge.arrow.mode=0)

#Node-level metrics
degree<-igraph::centr_degree(MreshNcl,mode = "all",normalized = T)
degree<-degree$res
degree
eigen_cent<-igraph::eigen_centrality(MreshNcl)
eigen_cent<-eigen_cent$vector
clos_cent <- igraph::closeness(MreshNcl, normalized=T)
betw_cent <- igraph::betweenness(MreshNcl, normalized=T)
hits <- hits_scores(MreshNcl, weights=E(MreshNcl)$weight)
hs <-hits$hub
authority<-hits$authority
Trans<-transitivity(MreshNcl,type="local")
Trans
#write results
MreshNcl.Centrality<-data.frame(NAME=V(MreshNcl)$name, UniprotID=V(MreshNcl)$UniprotID, Degree=degree, Eigen= eigen_cent, Closeness= clos_cent, Betweeness= betw_cent, Hubs=hs, Authority=authority,Transitivity=Trans)
MreshNcl.Centrality

write.csv(MreshNcl.Centrality,"./MreshNcl.Centrality.csv")

#Cocitation, Number of shared interactions among nodes
Cocitation<-cocitation(MreshNcl)
Cocitation#This matrix is similar to an adjacency matrix

#Plot with degree as node size
A1<-ggraph(MreshNcl,layout=Lyt.Tree)+
  geom_edge_fan(color="gray85", aes(width=(E(MreshNcl)$weight1)/3), alpha=0.4)+
  geom_node_point(aes(color=V(MreshNcl)$Annotation, size=10+degree),alpha=0.75)+
  geom_node_text(aes(label=name),color="gray25",size=4)+
  scale_color_manual(limits=V(MreshNcl)$Annotation, values=V(MreshNcl)$An.color)+
  theme_void()

A1
#Modularity Analysis
#Edge_Betweeness Modularity
MreshNcl.ebmod=cluster_edge_betweenness(MreshNcl,weights=E(MreshNcl)$weight1,directed = F) 
MreshNcl.ebmod #IGRAPH clustering edge betweenness, groups: 48, mod: 0.43

# Plot the Edge_Betweeness modules
plot(MreshNcl.ebmod, MreshNcl,layout=Lyt.Tree, alpha=0.6, 
     vertex.color=V(MreshNcl)$An.color, vertex.frame.color=V(MreshNcl)$An.color, vertex.shape="circle", 
     vertex.size=5, vertex.label=V(MreshNcl)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black", 
     edge.width=E(MreshNcl)$weight1, edge.lty=1, main="MreID Full Nuclear Network", sub="FR Tree Layout", 
     rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)
#Dendogram of the Edge_Betweeness
plot_dendrogram(MreshNcl.ebmod, mode="hclust")

#Louvain modularity
MreshNcl.Louvain.mod<-cluster_louvain(MreshNcl,weights=E(MreshNcl)$weight1)
MreshNcl.Louvain.mod #The result shows a total of IGRAPH clustering multi level,groups: 16, mod: 0.48

#Plot the Louvain modules
plot(MreshNcl.Louvain.mod, MreshNcl,layout=Lyt.Tree, alpha=0.6, 
     vertex.color=V(MreshNcl)$An.color, vertex.frame.color=V(MreshNcl)$An.color, vertex.shape="circle", 
     vertex.size=5, vertex.label=V(MreshNcl)$name, edge.color="gray50",vertex.label.cex=1, vertex.label.font=2, vertex.label.color="black", 
     edge.width=E(MreshNcl)$weight1, edge.lty=1, main="MreID Full Nuclear Network", sub="FR Tree Layout", 
     rescale=T)+
  legend( x="bottomright", levels(Mat$Annotation), pch=21,
          col=levels(Mat$An.color), pt.bg=levels(Mat$An.color), pt.cex=2, cex=.8, bty="n", ncol=1)


#Generate HeatMaps of the adjacency matrixes
#Let's use the neighbor list for generating the Adjaency Matrix
SensorID_FNN
v<-levels(as.factor(unlist(SensorID_FNN[,1:2]))) #get vertex names from the network list
n<-length(v)#get the number of vertex in the network
e<-matrix(match(as.character(unlist(SensorID_FNN[,1:2])),v),ncol = 2)#Edgelist with row index instead of Names
w<-MreID_MreshNcl$weight1 #Edge weights

# M is our co-expression network adjacency matrix.

# distringuish between i,j and j,i we want to make M symmetric (undirected) by
# considering an undirected network.
MreshNcl.M <- matrix(0, n, n)                  # set up a co-expression matrix
MreshNcl.M[e] <- w                             # fill it in with edge weights
MreshNcl.M <- MreshNcl.M + t(MreshNcl.M)                         # make this symmetric
dimnames(MreshNcl.M) <- list(v, v)             # label the vertices
# We let A_ij = 1 if there is an interaction. 
#A <- 1*(M > 0) #right now we keep the weights as they are
MreshNcl.M #The adjacency Matrix not as an igraph object, easier to manipulate, and to subset for smaller Sensor specific networks
quantile(MreshNcl.M,c(0,0.25,0.5,0.75,1))
#
#0%   25%   50%   75%  100% 
#0.000 0.000 0.000 0.000 0.991 

#Creating the HeatMap
# make a heatmap image with colors: white below 0,
# grey from 0.4 to 1 with darker closer to 1
breaks <- c(0, seq(0.05, 1, length=11)) # breaks in the color bins
cols <- grey(1-c(0,seq(0.5,1,length=10))) # colors

# color the vertex for CCNB1 blue, its neighbors red, and the others white
MreIDBind <- which(v == "MRE11A") # index for ccnb1
#MreIDBind<-as.factor(MreIDBind)
vcols <- rep("white",n)
vcols[MreIDBind] <- "blue"
#Subsetting MreshNcl.M by doing a loop asking for look into every row and column the index in MreIDBind worked
for (i in MreIDBind) {
  vcols[which(MreshNcl.M[,i] > 0 | MreshNcl.M[i,])] <- "red"
}
#Subsetting MreshNcl.M using the index vector in MreIDBind do not work
#vcols[which(MreshNcl.M[,MreIDBind] > 0 | MreshNcl.M[MreIDBind,])] <- "red"

#Generate a dendogram for keeping the proteins in the sample place
#MreshNcl.dist<-dist(MreshNcl.M,method = "euclidian")
#MreshNcl.Dend<-hclust(MreshNcl.dist, method = "ward.D2" )
#plot(MreshNcl.Dend)
#?hclust #methods= "ward.D", "ward.D2" this one works nice, "single", "complete", "average","mcquitty","median" (= WPGMC) or "centroid" (= UPGMC).

#Cut the dendogram
#cut.height <- 10 # try numbers from 4 to 7
#MreshNcl.hclusts<- cutree(MreshNcl.Dend, cut.height)
#plot(MreshNcl.hclusts)

# now actually make the heat map
par(mar=rep(0,4)) # make plot margins 0
heatmap(MreshNcl.M, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=vcols, RowSideColors=vcols,  col=0:1, frame=T)+
  legend("topleft", c("Neighbors(Mre11, MRE11, Mre11)","Mre11, MRE11, MRE11A"), fill=c("red","blue"),  
         bty="n", inset=0, xpd=T,  border=F)
##Color annotation for functional annotation
MatSet<-subset(Mat,Mat$NAME%in%v)# select only those proteins that are cointained in the node list of our matrixes
match(v,MatSet$NAME) # Check that the indexes of the Annotation Dataframe and the node list are in the same position 
match(rownames(MreshNcl.M),MatSet$NAME) #Check that the row names of the Adj Matrix and the Annotation dataframe are in the same order
match(colnames(MreshNcl.M),MatSet$NAME) #Check that the column names of the Adj Matrix and the Annotation dataframe are in the same order
AnCol<-MatSet$An.color#Vector of colors
AnCol<-as.character(AnCol)#transform it into a character vector, so the colors can be read


#
plot.new()
legend("topleft", c("Cell Cycle","Chromatin Org","DNA Repair","Nuclear Pore","p53","Proteasome","Ribosome Biogenesis","RNA Metabolism","Transcription Reg","UPR","Calcium Signaling"), 
       fill=c("coral1","darkviolet","firebrick1","deepskyblue3","tomato","grey","lightgreen","gold","slateblue1","maroon","sienna2"),
       bty="n", inset=0, xpd=T,  border=F, cex =2,y.intersp = 0.7,x.intersp = 0.2,pt.cex = 1)

#check with the MreshNcl
# now actually make the heat map
par(mar=rep(0,4)) # make plot margins 0
heatmap(MreshNcl.M, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=0:1, frame=T, cexRow = 0.5 , cexCol = 0.5)

legend("topleft", c("Cell Cycle","Chromatin Org","DNA Repair","Nuclear Pore","p53","Proteasome","Ribosome Biogenesis","RNA Metabolism","Transcription Reg","UPR","Calcium Signaling"), 
       fill=c("coral1","darkviolet","firebrick1","deepskyblue3","tomato","grey","lightgreen","gold","slateblue1","maroon","sienna2"), bty="n", inset=0, xpd=T,  border=F, cex =2,y.intersp = 0.7,x.intersp = 0.2,pt.cex = 1)


#Now we substract the MreshScr AdjMatrix to the MreshNcl AdjMatrix
#lets add 1 to all the weights of the shNcl matrix
MreshNcl.M.Plus1<-MreshNcl.M+1
range(MreshNcl.M.Plus1)
MreshNcl.M.Plus1.Change<-MreshNcl.M.Plus1-MreshScr.M
MreshNcl.M.Change<-MreshNcl.M-MreshScr.M
range(MreshNcl.M.Change)#Let's check the magnitude of the substracted weights
#[1] -0.716  0.705
any(is.na(MreshNcl.M.Change))
#Generate the breaks and color functions
breaks <- c(seq(-0.5, 0.5, length=10))
breaks
HMCol<- brewer.pal(9,name = "RdBu")
HMCol<- rev(HMCol)
HMCol
[1] "#2166AC" "#4393C3" "#92C5DE" "#D1E5F0" "#F7F7F7" "#FDDBC7" "#F4A582" "#D6604D" "#B2182B"

#Plot the HeatMap
par(mar=rep(0,4)) # make plot margins 0
heatmap(MreshNcl.M.Change, Rowv = as.dendrogram(FNN.Dend), symm=TRUE, 
        ColSideColors=AnCol, RowSideColors=AnCol,  col=HMCol, breaks=breaks,frame=T, cexRow = 0.5 , cexCol = 0.5, main = "Mre11 shNcl (minus shScr edges)", ylab="Total Protein Interactors",xlab="Total Protein Interactors")
#
dev.off()
plot.new()

#############################################################################################################################################
###### From the Igraph analysis, the resulting tables were merged in excel and usted to plot different network metrics
####################################### Network metrics plots################################################################################
#Set the working directory
setwd("C:/InteractomeProteome/LimmaAnalysis/shNcl/SensorID_All_RawData/Networks/igraph_Analysis")
#open the Network metrix data frame
SensorID_NetMet<-read.csv("./Results_Tables/Network_Metrics.csv")
SensorID_NetMet$Sensor<-as.factor(SensorID_NetMet$Sensor)
SensorID_NetMet$Treatment<-as.factor(SensorID_NetMet$Treatment)
str(SensorID_NetMet)

#Generate an ordered table
#Generate order vectors
SensorID_NetMet<-  SensorID_NetMet %>%
  mutate(Order.Sensor = case_when(
    (Sensor=="All")~1,
    (Sensor== "Ku80ID") ~ 3,
    (Sensor== "Sirt6ID") ~ 2,
    (Sensor== "Mre11ID") ~ 4
  ))
#
SensorID_NetMet<-  SensorID_NetMet %>%
  mutate(Order.Treatment = case_when(
    (Treatment=="All")~1,
    (Treatment== "shScramble") ~ 2,
    (Treatment== "shNcl") ~ 3,
  ))

#
SensorID_NetMet<-  SensorID_NetMet %>%
  mutate(Order.Network = case_when(
    (Network=="FNN")~1,
    (Network=="Sirt6ID_shScr")~2,
    (Network== "Sirt6ID_shNcl") ~ 3,
    (Network=="Ku80ID_shScr")~4,
    (Network== "Ku80ID_shNcl") ~ 5,
    (Network=="Mre11ID_shScr")~6,
    (Network=="Mre11ID_shNcl")~7
  ))

unique(SensorID_NetMet$Network)


SensorID_NetMet<-SensorID_NetMet[order(SensorID_NetMet$Order.Sensor),]
SensorID_NetMet$Sensor<-reorder(SensorID_NetMet$Sensor,SensorID_NetMet$Order.Sensor)
SensorID_NetMet<-SensorID_NetMet[order(SensorID_NetMet$Order.Treatment),]
SensorID_NetMet$Treatment<-reorder(SensorID_NetMet$Treatment,SensorID_NetMet$Order.Treatment)
#SensorID_NetMet<-SensorID_NetMet[order(SensorID_NetMet$Order.Annotation),]
#SensorID_NetMet$Annotation<-reorder(SensorID_NetMet$Annotation,SensorID_NetMet$Order.Annotation)
str(SensorID_NetMet)

SensorID_NetMet<-as.data.frame(SensorID_NetMet)

#Plot Connected Components 
SensorID_NetMet%>%
  mutate(Network=fct_reorder(Network,Order.Network))%>%
  ggplot(aes(x = Network, y =Connected.Components  , color = Sensor, fill=Sensor, group=Network)) + 
  geom_bar(stat = "identity", alpha=0.6,)+
  geom_text(aes(label=Connected.Components), position = position_stack(vjust = 0.5), size = 3, 
            colour = "gray10")+
  labs(x="Treatment", y="No of Connected Components")+
  scale_fill_manual(aesthetics = "color",values=c("purple","tomato","seagreen1","steelblue1"))+
  scale_fill_manual(aesthetics = "fill",values=c("purple","tomato","seagreen1","steelblue1"))+
  #scale_x_discrete("Time Points",labels=c("nonIR5","IR5","IR30","IR2hrs","IR8hrs","IR24h","nonIR24h"))+
  theme_bw()+  
  theme(axis.text=element_text(size=14,face = "bold"),
        axis.title=element_text(size=18,face="bold"),
        legend.text=element_text(size=18,face = "bold"),
        legend.title = element_text(size=18,face = "bold"),
        strip.text.x = element_text(size = 12,face = "bold"),
        strip.text.y = element_text(size = 12))

#Plot No.Connected.Interactors 
SensorID_NetMet%>%
  mutate(Network=fct_reorder(Network,Order.Network))%>%
  ggplot(aes(x = Network, y =No.Connected.Interactors  , color = Sensor, fill=Sensor, group=Network)) + 
  geom_bar(stat = "identity", alpha=0.6,)+
  geom_text(aes(label=No.Connected.Interactors), position = position_stack(vjust = 0.5), size = 3, 
            colour = "gray10")+
  labs(x="Treatment", y="No of Connected Interactors")+
  scale_fill_manual(aesthetics = "color",values=c("purple","tomato","seagreen1","steelblue1"))+
  scale_fill_manual(aesthetics = "fill",values=c("purple","tomato","seagreen1","steelblue1"))+
  #scale_x_discrete("Time Points",labels=c("nonIR5","IR5","IR30","IR2hrs","IR8hrs","IR24h","nonIR24h"))+
  theme_bw()+  
  theme(axis.text=element_text(size=14,face = "bold"),
        axis.title=element_text(size=18,face="bold"),
        legend.text=element_text(size=18,face = "bold"),
        legend.title = element_text(size=18,face = "bold"),
        strip.text.x = element_text(size = 12,face = "bold"),
        strip.text.y = element_text(size = 12))

#Plot No.Edges 
SensorID_NetMet%>%
  mutate(Network=fct_reorder(Network,Order.Network))%>%
  ggplot(aes(x = Network, y =No.Edges  , color = Sensor, fill=Sensor, group=Network)) + 
  geom_bar(stat = "identity", alpha=0.6,)+
  geom_text(aes(label=No.Edges), position = position_stack(vjust = 0.5), size = 3, 
            colour = "gray10")+
  labs(x="Treatment", y="No of Edges")+
  scale_fill_manual(aesthetics = "color",values=c("purple","tomato","seagreen1","steelblue1"))+
  scale_fill_manual(aesthetics = "fill",values=c("purple","tomato","seagreen1","steelblue1"))+
  #scale_x_discrete("Time Points",labels=c("nonIR5","IR5","IR30","IR2hrs","IR8hrs","IR24h","nonIR24h"))+
  theme_bw()+  
  theme(axis.text=element_text(size=14,face = "bold"),
        axis.title=element_text(size=18,face="bold"),
        legend.text=element_text(size=18,face = "bold"),
        legend.title = element_text(size=18,face = "bold"),
        strip.text.x = element_text(size = 12,face = "bold"),
        strip.text.y = element_text(size = 12))

#Plot Edge.Density 
SensorID_NetMet%>%
  mutate(Network=fct_reorder(Network,Order.Network))%>%
  ggplot(aes(x = Network, y =Edge.Density  , color = Sensor, fill=Sensor, group=Network,
             label=sprintf("%0.4f", round(Edge.Density, digits = 4)))) + 
  geom_bar(stat = "identity", alpha=0.6,)+
  geom_text(position = position_stack(vjust = 0.5), size = 3, 
            colour = "gray10")+
  labs(x="Treatment", y="Edge Density")+
  scale_fill_manual(aesthetics = "color",values=c("purple","tomato","seagreen1","steelblue1"))+
  scale_fill_manual(aesthetics = "fill",values=c("purple","tomato","seagreen1","steelblue1"))+
  #scale_x_discrete("Time Points",labels=c("nonIR5","IR5","IR30","IR2hrs","IR8hrs","IR24h","nonIR24h"))+
  theme_bw()+  
  theme(axis.text=element_text(size=14,face = "bold"),
        axis.title=element_text(size=18,face="bold"),
        legend.text=element_text(size=18,face = "bold"),
        legend.title = element_text(size=18,face = "bold"),
        strip.text.x = element_text(size = 12,face = "bold"),
        strip.text.y = element_text(size = 12))

#Plot Network Diameter
SensorID_NetMet%>%
  mutate(Network=fct_reorder(Network,Order.Network))%>%
  ggplot(aes(x = Network, y =Diameter  , color = Sensor, fill=Sensor, group=Network,
             label=sprintf("%0.4f", round(Diameter, digits = 4)))) + 
  geom_bar(stat = "identity", alpha=0.6,)+
  geom_text(position = position_stack(vjust = 0.5), size = 3, 
            colour = "gray10")+
  labs(x="Treatment", y="Network Diameter")+
  scale_fill_manual(aesthetics = "color",values=c("purple","tomato","seagreen1","steelblue1"))+
  scale_fill_manual(aesthetics = "fill",values=c("purple","tomato","seagreen1","steelblue1"))+
  #scale_x_discrete("Time Points",labels=c("nonIR5","IR5","IR30","IR2hrs","IR8hrs","IR24h","nonIR24h"))+
  theme_bw()+  
  theme(axis.text=element_text(size=14,face = "bold"),
        axis.title=element_text(size=18,face="bold"),
        legend.text=element_text(size=18,face = "bold"),
        legend.title = element_text(size=18,face = "bold"),
        strip.text.x = element_text(size = 12,face = "bold"),
        strip.text.y = element_text(size = 12))


#Plot Centralization 
SensorID_NetMet%>%
  mutate(Network=fct_reorder(Network,Order.Network))%>%
  ggplot(aes(x = Network, y =Centralization  , color = Sensor, fill=Sensor, group=Network,
             label=sprintf("%0.4f", round(Centralization, digits = 4)))) + 
  geom_bar(stat = "identity", alpha=0.6,)+
  geom_text(position = position_stack(vjust = 0.5), size = 3, 
            colour = "gray10")+
  labs(x="Treatment", y="Network Centralization")+
  scale_fill_manual(aesthetics = "color",values=c("purple","tomato","seagreen1","steelblue1"))+
  scale_fill_manual(aesthetics = "fill",values=c("purple","tomato","seagreen1","steelblue1"))+
  #scale_x_discrete("Time Points",labels=c("nonIR5","IR5","IR30","IR2hrs","IR8hrs","IR24h","nonIR24h"))+
  theme_bw()+  
  theme(axis.text=element_text(size=14,face = "bold"),
        axis.title=element_text(size=18,face="bold"),
        legend.text=element_text(size=18,face = "bold"),
        legend.title = element_text(size=18,face = "bold"),
        strip.text.x = element_text(size = 12,face = "bold"),
        strip.text.y = element_text(size = 12))


#Plot Mean.Neighbors 
SensorID_NetMet%>%
  mutate(Network=fct_reorder(Network,Order.Network))%>%
  ggplot(aes(x = Network, y =Mean.Neighbors  , color = Sensor, fill=Sensor, group=Network,
             label=sprintf("%0.4f", round(Mean.Neighbors, digits = 4)))) + 
  geom_bar(stat = "identity", alpha=0.6,)+
  geom_text(position = position_stack(vjust = 0.5), size = 3, 
            colour = "gray10")+
  labs(x="Treatment", y="Average No of Neighbors ")+
  scale_fill_manual(aesthetics = "color",values=c("purple","tomato","seagreen1","steelblue1"))+
  scale_fill_manual(aesthetics = "fill",values=c("purple","tomato","seagreen1","steelblue1"))+
  #scale_x_discrete("Time Points",labels=c("nonIR5","IR5","IR30","IR2hrs","IR8hrs","IR24h","nonIR24h"))+
  theme_bw()+  
  theme(axis.text=element_text(size=14,face = "bold"),
        axis.title=element_text(size=18,face="bold"),
        legend.text=element_text(size=18,face = "bold"),
        legend.title = element_text(size=18,face = "bold"),
        strip.text.x = element_text(size = 12,face = "bold"),
        strip.text.y = element_text(size = 12))

#Plot ClusteringCoeff
SensorID_NetMet%>%
  mutate(Network=fct_reorder(Network,Order.Network))%>%
  ggplot(aes(x = Network, y =ClusteringCoeff  , color = Sensor, fill=Sensor, group=Network,
             label=sprintf("%0.4f", round(ClusteringCoeff, digits = 4)))) + 
  geom_bar(stat = "identity", alpha=0.6,)+
  geom_text(position = position_stack(vjust = 0.5), size = 3, 
            colour = "gray10")+
  labs(x="Treatment", y="Network Clsutering Coefficient")+
  scale_fill_manual(aesthetics = "color",values=c("purple","tomato","seagreen1","steelblue1"))+
  scale_fill_manual(aesthetics = "fill",values=c("purple","tomato","seagreen1","steelblue1"))+
  #scale_x_discrete("Time Points",labels=c("nonIR5","IR5","IR30","IR2hrs","IR8hrs","IR24h","nonIR24h"))+
  theme_bw()+  
  theme(axis.text=element_text(size=14,face = "bold"),
        axis.title=element_text(size=18,face="bold"),
        legend.text=element_text(size=18,face = "bold"),
        legend.title = element_text(size=18,face = "bold"),
        strip.text.x = element_text(size = 12,face = "bold"),
        strip.text.y = element_text(size = 12))


#Plot Transitivity

SensorID_NetMet%>%
  mutate(Network=fct_reorder(Network,Order.Network))%>%
  ggplot(aes(x = Network, y =Transitivity  , color = Sensor, fill=Sensor, group=Network,
             label=sprintf("%0.4f", round(Transitivity, digits = 4)))) + 
  geom_bar(stat = "identity", alpha=0.6,)+
  geom_text(position = position_stack(vjust = 0.5), size = 3, 
            colour = "gray10")+
  labs(x="Treatment", y="Network Transitivity")+
  scale_fill_manual(aesthetics = "color",values=c("purple","tomato","seagreen1","steelblue1"))+
  scale_fill_manual(aesthetics = "fill",values=c("purple","tomato","seagreen1","steelblue1"))+
  #scale_x_discrete("Time Points",labels=c("nonIR5","IR5","IR30","IR2hrs","IR8hrs","IR24h","nonIR24h"))+
  theme_bw()+  
  theme(axis.text=element_text(size=14,face = "bold"),
        axis.title=element_text(size=18,face="bold"),
        legend.text=element_text(size=18,face = "bold"),
        legend.title = element_text(size=18,face = "bold"),
        strip.text.x = element_text(size = 12,face = "bold"),
        strip.text.y = element_text(size = 12))

#Plot Annotation modularity

SensorID_NetMet%>%
  mutate(Network=fct_reorder(Network,Order.Network))%>%
  ggplot(aes(x = Network, y =Anno.Modularity  , color = Sensor, fill=Sensor, group=Network,
             label=sprintf("%0.4f", round(Anno.Modularity, digits = 4)))) + 
  geom_bar(stat = "identity", alpha=0.6,)+
  geom_text(position = position_stack(vjust = 0.5), size = 3, 
            colour = "gray10")+
  labs(x="Treatment", y="Network Annotation-dependent Modularity")+
  scale_fill_manual(aesthetics = "color",values=c("purple","tomato","seagreen1","steelblue1"))+
  scale_fill_manual(aesthetics = "fill",values=c("purple","tomato","seagreen1","steelblue1"))+
  #scale_x_discrete("Time Points",labels=c("nonIR5","IR5","IR30","IR2hrs","IR8hrs","IR24h","nonIR24h"))+
  theme_bw()+  
  theme(axis.text=element_text(size=14,face = "bold"),
        axis.title=element_text(size=18,face="bold"),
        legend.text=element_text(size=18,face = "bold"),
        legend.title = element_text(size=18,face = "bold"),
        strip.text.x = element_text(size = 12,face = "bold"),
        strip.text.y = element_text(size = 12))


#Plot Nuclear localization modularity

SensorID_NetMet%>%
  mutate(Network=fct_reorder(Network,Order.Network))%>%
  ggplot(aes(x = Network, y =NucLoc.Modularity  , color = Sensor, fill=Sensor, group=Network,
             label=sprintf("%0.4f", round(NucLoc.Modularity, digits = 4)))) + 
  geom_bar(stat = "identity", alpha=0.6,)+
  geom_text(position = position_stack(vjust = 0.5), size = 3, 
            colour = "gray10")+
  labs(x="Treatment", y="Network Nuclear Localization-dependent Modularity")+
  scale_fill_manual(aesthetics = "color",values=c("purple","tomato","seagreen1","steelblue1"))+
  scale_fill_manual(aesthetics = "fill",values=c("purple","tomato","seagreen1","steelblue1"))+
  #scale_x_discrete("Time Points",labels=c("nonIR5","IR5","IR30","IR2hrs","IR8hrs","IR24h","nonIR24h"))+
  theme_bw()+  
  theme(axis.text=element_text(size=14,face = "bold"),
        axis.title=element_text(size=18,face="bold"),
        legend.text=element_text(size=18,face = "bold"),
        legend.title = element_text(size=18,face = "bold"),
        strip.text.x = element_text(size = 12,face = "bold"),
        strip.text.y = element_text(size = 12))

#Plot DNA repair modularity

SensorID_NetMet%>%
  mutate(Network=fct_reorder(Network,Order.Network))%>%
  ggplot(aes(x = Network, y =DNARep.Modularity  , color = Sensor, fill=Sensor, group=Network,
             label=sprintf("%0.4f", round(DNARep.Modularity, digits = 4)))) + 
  geom_bar(stat = "identity", alpha=0.6,)+
  geom_text(position = position_stack(vjust = 0.5), size = 3, 
            colour = "gray10")+
  labs(x="Treatment", y="Network Nuclear DNA repair proteins-dependent Modularity")+
  scale_fill_manual(aesthetics = "color",values=c("purple","tomato","seagreen1","steelblue1"))+
  scale_fill_manual(aesthetics = "fill",values=c("purple","tomato","seagreen1","steelblue1"))+
  #scale_x_discrete("Time Points",labels=c("nonIR5","IR5","IR30","IR2hrs","IR8hrs","IR24h","nonIR24h"))+
  theme_bw()+  
  theme(axis.text=element_text(size=14,face = "bold"),
        axis.title=element_text(size=18,face="bold"),
        legend.text=element_text(size=18,face = "bold"),
        legend.title = element_text(size=18,face = "bold"),
        strip.text.x = element_text(size = 12,face = "bold"),
        strip.text.y = element_text(size = 12))

#Plot Annotation assortativity

SensorID_NetMet%>%
  mutate(Network=fct_reorder(Network,Order.Network))%>%
  ggplot(aes(x = Network, y =Anno.Assortativity  , color = Sensor, fill=Sensor, group=Network,
             label=sprintf("%0.4f", round(Anno.Assortativity, digits = 4)))) + 
  geom_bar(stat = "identity", alpha=0.6,)+
  geom_text(position = position_stack(vjust = 0.5), size = 3, 
            colour = "gray10")+
  labs(x="Treatment", y="Network Annotation Assorativity")+
  scale_fill_manual(aesthetics = "color",values=c("purple","tomato","seagreen1","steelblue1"))+
  scale_fill_manual(aesthetics = "fill",values=c("purple","tomato","seagreen1","steelblue1"))+
  #scale_x_discrete("Time Points",labels=c("nonIR5","IR5","IR30","IR2hrs","IR8hrs","IR24h","nonIR24h"))+
  theme_bw()+  
  theme(axis.text=element_text(size=14,face = "bold"),
        axis.title=element_text(size=18,face="bold"),
        legend.text=element_text(size=18,face = "bold"),
        legend.title = element_text(size=18,face = "bold"),
        strip.text.x = element_text(size = 12,face = "bold"),
        strip.text.y = element_text(size = 12))


#Plot Nuclear Localization assortativity

SensorID_NetMet%>%
  mutate(Network=fct_reorder(Network,Order.Network))%>%
  ggplot(aes(x = Network, y =NucLoc.Assortativity  , color = Sensor, fill=Sensor, group=Network,
             label=sprintf("%0.4f", round(NucLoc.Assortativity, digits = 4)))) + 
  geom_bar(stat = "identity", alpha=0.6,)+
  geom_text(position = position_stack(vjust = 0.5), size = 3, 
            colour = "gray10")+
  labs(x="Treatment", y="Network Nuclear Localization Assorativity")+
  scale_fill_manual(aesthetics = "color",values=c("purple","tomato","seagreen1","steelblue1"))+
  scale_fill_manual(aesthetics = "fill",values=c("purple","tomato","seagreen1","steelblue1"))+
  #scale_x_discrete("Time Points",labels=c("nonIR5","IR5","IR30","IR2hrs","IR8hrs","IR24h","nonIR24h"))+
  theme_bw()+  
  theme(axis.text=element_text(size=14,face = "bold"),
        axis.title=element_text(size=18,face="bold"),
        legend.text=element_text(size=18,face = "bold"),
        legend.title = element_text(size=18,face = "bold"),
        strip.text.x = element_text(size = 12,face = "bold"),
        strip.text.y = element_text(size = 12))

#Plot DNARep.Assortativity

SensorID_NetMet%>%
  mutate(Network=fct_reorder(Network,Order.Network))%>%
  ggplot(aes(x = Network, y =DNARep.Assortativity  , color = Sensor, fill=Sensor, group=Network,
             label=sprintf("%0.4f", round(DNARep.Assortativity, digits = 4)))) + 
  geom_bar(stat = "identity", alpha=0.6,)+
  geom_text(position = position_stack(vjust = 0.5), size = 3, 
            colour = "gray10")+
  labs(x="Treatment", y="Network DNA repair proteins Assorativity")+
  scale_fill_manual(aesthetics = "color",values=c("purple","tomato","seagreen1","steelblue1"))+
  scale_fill_manual(aesthetics = "fill",values=c("purple","tomato","seagreen1","steelblue1"))+
  #scale_x_discrete("Time Points",labels=c("nonIR5","IR5","IR30","IR2hrs","IR8hrs","IR24h","nonIR24h"))+
  theme_bw()+  
  theme(axis.text=element_text(size=14,face = "bold"),
        axis.title=element_text(size=18,face="bold"),
        legend.text=element_text(size=18,face = "bold"),
        legend.title = element_text(size=18,face = "bold"),
        strip.text.x = element_text(size = 12,face = "bold"),
        strip.text.y = element_text(size = 12))

################################################################################################
