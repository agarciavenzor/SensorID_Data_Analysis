
### This R Script has the objective of using RNetcarto simulated annealing algorithm to identify the module composition of the SensorID networks,
# morevoer, the connectivity (within-module degree) and participation (other-module degree) coefficients to assign specific roles to each node.


#NOTE: The results from the modularity analysis were merged with previous results to generate a node data frame that includes all the results generated 
#From the analysis, this dataframe is the "SensorID_Node_DataFrame.csv". note that the node metadata file already contains the module, connectiviy and 
#Participation values from the first analysis we did.

###
#Install rnetcarto package, and other required packages.
install.packages(c("rnetcarto","igraph","ggplot2"))
#Load required packages
require(igraph)
require(rnetcarto)
require(ggplot2)

#Set the working directory, as the folder that contain the SensorID networks edgelists.
setwd("../SensorID_rNetcarto_Modularity")

#### Create a color palette to generate contrast between the modules
custom_palette <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')
functions_colours <- rainbow(10)

###### Generate an Edge weight to highlight those edges that crosses its own module ######
edge.weights <- function(community, network, weight.within = .001, weight.between = .01) {
  bridges <- crossing(communities = community, graph = network)
  weights <- ifelse(test = bridges, yes = weight.between, no = weight.within)
  return(weights) 
}

###Sirt6ID networks
### Full Sirt6ID Nuclear Network
#### read the data
metadata <- read.csv('./SensorID_Node_DataFrame.csv')#Load the Node metadata table
metadata<-subset(metadata,metadata$Sensor=="Sirt6")#Isolate the node information corresponding to Sirt6ID networks
df <- read.table('./Sirt6_FullNetwork.tsv', header=T) #Load the edge list table
g <- graph_from_data_frame(df, directed = FALSE) #generate an igraph objetct
g <- simplify(g, remove.multiple=TRUE, remove.loops=TRUE, edge.attr.comb = "first") #remove duplicated edges, and self-loops
edge.attributes(g) #Check no NAs in the edge attributes
adj_mat_g <- as_adjacency_matrix(g, attr = 'weight', sparse = FALSE) #convert the igraph object into an adjacency matrix
mod <- netcarto(adj_mat_g) #generate the modules using the netcarto algorithm
#Now generate a function to run the netcarto algorithm 20 times, to improve the modularity results from the first run
for(i in 1:20){
  temp_mod <- netcarto(adj_mat_g)
  if(temp_mod[[2]] > mod[[2]])
    mod <- temp_mod
}

# integrate the module results from netcarto to the igraph network 
com_mod <- make_clusters(
  g,
  membership = mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1,
  algorithm = 'netcarto',
  merges = NULL,
  modularity = mod[[2]]
)
#Generate a color palette for the modules
V(g)$module <- mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1
V(g)$color <- custom_palette[V(g)$module]
#generate an color palette for the functional annotation of each node
mod[[1]]$annotation <- as.factor(metadata[match(mod[[1]]$name, metadata$NAME),]$Annotation)
V(g)$color_border <- functions_colours[as.numeric(mod[[1]][match(V(g)$name, mod[[1]]$name),]$annotation)]

## Generate a weight column that includes the detected modules as an edge attribute
cl_for_layout <- com_mod
E(g)$weight <- edge.weights(cl_for_layout, g)

##### ATTENTION : THE CODE BELOW HAS TO BE RUN AT LEAST ONCE - UNCOMMENT ON FIRST RUN

# # create layout - run this line several times until you find a layout that you like
coord_graphopt <- layout_with_graphopt(g, mass = 20, charge = E(g)$weight)
igraph::plot.igraph(g,layout=coord_graphopt)
#After generating a good layout this can be saved for ploting all the other Sirt6ID networks 
## and we store the layout coordinates to plot the other networks in the same way
coords_nodes <- coord_graphopt
row.names(coords_nodes) <- V(g)$name
coords_nodes
idxs <- match(V(g)$name, row.names(coords_nodes))
coord_graphopt <- coords_nodes[idxs,]
vsize <- igraph::degree(g)/vcount(g)

#Plot the networks with the module as the fill color, and the annotation as the color of the circle line
par(mar=c(0,0,0,0))
pdf(paste0('Sirt6_protein-protein-network-whole.pdf'), width=30, height=30)
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('M = ', format(mod[[2]], digits=2)))
#
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.shape='circle', vertex.frame.color=V(g)$color_border, vertex.frame.width=4, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('whole network M = ', format(mod[[2]], digits = 2)), edge.arrow.size=.2)

legend('bottomleft', legend=levels(mod[[1]]$annotation), pt.lwd=4, cex=2, col=functions_colours, pch=21, pt.bg='white', title='function')

#recover the module results into a new table
cur_mod <- mod[[1]]
cur_mod$name <- as.character(cur_mod$name)
cur_mod$module <- cur_mod$module + 1
cur_mod$color <- custom_palette[cur_mod$module]
cur_mod$treatment <- 'whole'
cur_mod$degree <- degree(g, mode='all')[match(cur_mod$name, names(degree(g, mode='all')))]

#Plot the  Participation coefficient versus the degree 
par(mar=c(4,4,4,4))
ggplot(cur_mod, aes(degree, participation, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('degree-part-whole.pdf'), width = 12, height = 8)

# Plot the connectivity (within-module degree) versus the degree
ggplot(cur_mod, aes(degree, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('degree-conn-whole.pdf'), width = 12, height = 8)

#Plot the connectivity versus the participation
ggplot(cur_mod, aes(participation, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('part-conn-whole.pdf'), width = 12, height = 8)

#Print the results in a new table to merge them with the Node data frames
mod_overall <- cur_mod
write.csv(mod_overall,'./Sirt6_FullNetwork_overall_mod.csv')

########### Work with each one of the time-point networks
##########
###nonIR5
metadata <- read.csv('./SensorID_Node_DataFrame.csv')#Load the Node metadata table
metadata<-subset(metadata,metadata$Sensor=="Sirt6")#Isolate the node information corresponding to Sirt6ID networks
df <- read.table('./S6_nonIR5_Edges.tsv', header=T) #Load the edge list table
g <- graph_from_data_frame(df, directed = FALSE) #generate an igraph objetct
g <- simplify(g, remove.multiple=TRUE, remove.loops=TRUE, edge.attr.comb = "first") #remove duplicated edges, and self-loops
edge.attributes(g) #Check no NAs in the edge attributes
adj_mat_g <- as_adjacency_matrix(g, attr = 'weight', sparse = FALSE) #convert the igraph object into an adjacency matrix
mod <- netcarto(adj_mat_g) #generate the modules using the netcarto algorithm
#Now generate a function to run the netcarto algorithm 20 times, to improve the modularity results from the first run
for(i in 1:20){
  temp_mod <- netcarto(adj_mat_g)
  if(temp_mod[[2]] > mod[[2]])
    mod <- temp_mod
}

# integrate the module results from netcarto to the igraph network 
com_mod <- make_clusters(
  g,
  membership = mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1,
  algorithm = 'netcarto',
  merges = NULL,
  modularity = mod[[2]]
)
#Generate a color palette for the modules
V(g)$module <- mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1
V(g)$color <- custom_palette[V(g)$module]
#generate an color palette for the functional annotation of each node
mod[[1]]$annotation <- as.factor(metadata[match(mod[[1]]$name, metadata$NAME),]$Annotation)
V(g)$color_border <- functions_colours[as.numeric(mod[[1]][match(V(g)$name, mod[[1]]$name),]$annotation)]

## Generate a weight column that includes the detected modules as an edge attribute
cl_for_layout <- com_mod
E(g)$weight <- edge.weights(cl_for_layout, g)

##### ATTENTION : THE CODE BELOW HAS TO BE RUN AT LEAST ONCE - UNCOMMENT ON FIRST RUN

# # create layout - run this line several times until you find a layout that you like
coord_graphopt <- layout_with_graphopt(g, mass = 20, charge = E(g)$weight)
igraph::plot.igraph(g,layout=coord_graphopt)
#After generating a good layout this can be saved for ploting all the other Sirt6ID networks 
## and we store the layout coordinates to plot the other networks in the same way
coords_nodes <- coord_graphopt
row.names(coords_nodes) <- V(g)$name
coords_nodes
idxs <- match(V(g)$name, row.names(coords_nodes))
coord_graphopt <- coords_nodes[idxs,]
vsize <- igraph::degree(g)/vcount(g)

#Plot the networks with the module as the fill color, and the annotation as the color of the circle line
par(mar=c(0,0,0,0))
pdf(paste0('Sirt6_protein-protein-network-whole.pdf'), width=30, height=30)
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('M = ', format(mod[[2]], digits=2)))
#
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.shape='circle', vertex.frame.color=V(g)$color_border, vertex.frame.width=4, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('whole network M = ', format(mod[[2]], digits = 2)), edge.arrow.size=.2)

legend('bottomleft', legend=levels(mod[[1]]$annotation), pt.lwd=4, cex=2, col=functions_colours, pch=21, pt.bg='white', title='function')

#recover the module results into a new table
cur_mod <- mod[[1]]
cur_mod$name <- as.character(cur_mod$name)
cur_mod$module <- cur_mod$module + 1
cur_mod$color <- custom_palette[cur_mod$module]
cur_mod$treatment <- 'nonIR5'
cur_mod$degree <- degree(g, mode='all')[match(cur_mod$name, names(degree(g, mode='all')))]

#Plot the  Participation coefficient versus the degree 
par(mar=c(4,4,4,4))
ggplot(cur_mod, aes(degree, participation, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('degree-part-nonIR5.pdf'), width = 12, height = 8)

# Plot the connectivity (within-module degree) versus the degree
ggplot(cur_mod, aes(degree, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('degree-conn-nonIR5.pdf'), width = 12, height = 8)

#Plot the connectivity versus the participation
ggplot(cur_mod, aes(participation, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('part-conn-nonIR5.pdf'), width = 12, height = 8)

#Print the results in a new table to merge them with the Node data frames
mod_overall <- cur_mod
write.csv(mod_overall,'./Sirt6_nonIR5_mod.csv')

###IR5
metadata <- read.csv('./SensorID_Node_DataFrame.csv')#Load the Node metadata table
metadata<-subset(metadata,metadata$Sensor=="Sirt6")#Isolate the node information corresponding to Sirt6ID networks
df <- read.table('./S6_IR5_Edges.tsv', header=T) #Load the edge list table
g <- graph_from_data_frame(df, directed = FALSE) #generate an igraph objetct
g <- simplify(g, remove.multiple=TRUE, remove.loops=TRUE, edge.attr.comb = "first") #remove duplicated edges, and self-loops
edge.attributes(g) #Check no NAs in the edge attributes
adj_mat_g <- as_adjacency_matrix(g, attr = 'weight', sparse = FALSE) #convert the igraph object into an adjacency matrix
mod <- netcarto(adj_mat_g) #generate the modules using the netcarto algorithm
#Now generate a function to run the netcarto algorithm 20 times, to improve the modularity results from the first run
for(i in 1:20){
  temp_mod <- netcarto(adj_mat_g)
  if(temp_mod[[2]] > mod[[2]])
    mod <- temp_mod
}

# integrate the module results from netcarto to the igraph network 
com_mod <- make_clusters(
  g,
  membership = mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1,
  algorithm = 'netcarto',
  merges = NULL,
  modularity = mod[[2]]
)
#Generate a color palette for the modules
V(g)$module <- mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1
V(g)$color <- custom_palette[V(g)$module]
#generate an color palette for the functional annotation of each node
mod[[1]]$annotation <- as.factor(metadata[match(mod[[1]]$name, metadata$NAME),]$Annotation)
V(g)$color_border <- functions_colours[as.numeric(mod[[1]][match(V(g)$name, mod[[1]]$name),]$annotation)]

## Generate a weight column that includes the detected modules as an edge attribute
cl_for_layout <- com_mod
E(g)$weight <- edge.weights(cl_for_layout, g)

##### ATTENTION : THE CODE BELOW HAS TO BE RUN AT LEAST ONCE - UNCOMMENT ON FIRST RUN

# # create layout - run this line several times until you find a layout that you like
coord_graphopt <- layout_with_graphopt(g, mass = 20, charge = E(g)$weight)
igraph::plot.igraph(g,layout=coord_graphopt)
#After generating a good layout this can be saved for ploting all the other Sirt6ID networks 
## and we store the layout coordinates to plot the other networks in the same way
coords_nodes <- coord_graphopt
row.names(coords_nodes) <- V(g)$name
coords_nodes
idxs <- match(V(g)$name, row.names(coords_nodes))
coord_graphopt <- coords_nodes[idxs,]
vsize <- igraph::degree(g)/vcount(g)

#Plot the networks with the module as the fill color, and the annotation as the color of the circle line
par(mar=c(0,0,0,0))
pdf(paste0('Sirt6_protein-protein-network-whole.pdf'), width=30, height=30)
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('M = ', format(mod[[2]], digits=2)))
#
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.shape='circle', vertex.frame.color=V(g)$color_border, vertex.frame.width=4, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('whole network M = ', format(mod[[2]], digits = 2)), edge.arrow.size=.2)

legend('bottomleft', legend=levels(mod[[1]]$annotation), pt.lwd=4, cex=2, col=functions_colours, pch=21, pt.bg='white', title='function')

#recover the module results into a new table
cur_mod <- mod[[1]]
cur_mod$name <- as.character(cur_mod$name)
cur_mod$module <- cur_mod$module + 1
cur_mod$color <- custom_palette[cur_mod$module]
cur_mod$treatment <- 'IR5'
cur_mod$degree <- degree(g, mode='all')[match(cur_mod$name, names(degree(g, mode='all')))]

#Plot the  Participation coefficient versus the degree 
par(mar=c(4,4,4,4))
ggplot(cur_mod, aes(degree, participation, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('degree-part-IR5.pdf'), width = 12, height = 8)

# Plot the connectivity (within-module degree) versus the degree
ggplot(cur_mod, aes(degree, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('degree-conn-IR5.pdf'), width = 12, height = 8)

#Plot the connectivity versus the participation
ggplot(cur_mod, aes(participation, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('part-conn-IR5.pdf'), width = 12, height = 8)

#Print the results in a new table to merge them with the Node data frames
mod_overall <- cur_mod
write.csv(mod_overall,'./Sirt6_IR5_mod.csv')

###IR30
metadata <- read.csv('./SensorID_Node_DataFrame.csv')#Load the Node metadata table
metadata<-subset(metadata,metadata$Sensor=="Sirt6")#Isolate the node information corresponding to Sirt6ID networks
df <- read.table('./S6_IR30_Edges.tsv', header=T) #Load the edge list table
g <- graph_from_data_frame(df, directed = FALSE) #generate an igraph objetct
g <- simplify(g, remove.multiple=TRUE, remove.loops=TRUE, edge.attr.comb = "first") #remove duplicated edges, and self-loops
edge.attributes(g) #Check no NAs in the edge attributes
adj_mat_g <- as_adjacency_matrix(g, attr = 'weight', sparse = FALSE) #convert the igraph object into an adjacency matrix
mod <- netcarto(adj_mat_g) #generate the modules using the netcarto algorithm
#Now generate a function to run the netcarto algorithm 20 times, to improve the modularity results from the first run
for(i in 1:20){
  temp_mod <- netcarto(adj_mat_g)
  if(temp_mod[[2]] > mod[[2]])
    mod <- temp_mod
}

# integrate the module results from netcarto to the igraph network 
com_mod <- make_clusters(
  g,
  membership = mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1,
  algorithm = 'netcarto',
  merges = NULL,
  modularity = mod[[2]]
)
#Generate a color palette for the modules
V(g)$module <- mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1
V(g)$color <- custom_palette[V(g)$module]
#generate an color palette for the functional annotation of each node
mod[[1]]$annotation <- as.factor(metadata[match(mod[[1]]$name, metadata$NAME),]$Annotation)
V(g)$color_border <- functions_colours[as.numeric(mod[[1]][match(V(g)$name, mod[[1]]$name),]$annotation)]

## Generate a weight column that includes the detected modules as an edge attribute
cl_for_layout <- com_mod
E(g)$weight <- edge.weights(cl_for_layout, g)

##### ATTENTION : THE CODE BELOW HAS TO BE RUN AT LEAST ONCE - UNCOMMENT ON FIRST RUN

# # create layout - run this line several times until you find a layout that you like
coord_graphopt <- layout_with_graphopt(g, mass = 20, charge = E(g)$weight)
igraph::plot.igraph(g,layout=coord_graphopt)
#After generating a good layout this can be saved for ploting all the other Sirt6ID networks 
## and we store the layout coordinates to plot the other networks in the same way
coords_nodes <- coord_graphopt
row.names(coords_nodes) <- V(g)$name
coords_nodes
idxs <- match(V(g)$name, row.names(coords_nodes))
coord_graphopt <- coords_nodes[idxs,]
vsize <- igraph::degree(g)/vcount(g)

#Plot the networks with the module as the fill color, and the annotation as the color of the circle line
par(mar=c(0,0,0,0))
pdf(paste0('Sirt6_protein-protein-network-whole.pdf'), width=30, height=30)
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('M = ', format(mod[[2]], digits=2)))
#
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.shape='circle', vertex.frame.color=V(g)$color_border, vertex.frame.width=4, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('whole network M = ', format(mod[[2]], digits = 2)), edge.arrow.size=.2)

legend('bottomleft', legend=levels(mod[[1]]$annotation), pt.lwd=4, cex=2, col=functions_colours, pch=21, pt.bg='white', title='function')

#recover the module results into a new table
cur_mod <- mod[[1]]
cur_mod$name <- as.character(cur_mod$name)
cur_mod$module <- cur_mod$module + 1
cur_mod$color <- custom_palette[cur_mod$module]
cur_mod$treatment <- 'IR30'
cur_mod$degree <- degree(g, mode='all')[match(cur_mod$name, names(degree(g, mode='all')))]

#Plot the  Participation coefficient versus the degree 
par(mar=c(4,4,4,4))
ggplot(cur_mod, aes(degree, participation, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('degree-part-IR30.pdf'), width = 12, height = 8)

# Plot the connectivity (within-module degree) versus the degree
ggplot(cur_mod, aes(degree, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('degree-conn-IR30.pdf'), width = 12, height = 8)

#Plot the connectivity versus the participation
ggplot(cur_mod, aes(participation, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('part-conn-IR30.pdf'), width = 12, height = 8)

#Print the results in a new table to merge them with the Node data frames
mod_overall <- cur_mod
write.csv(mod_overall,'./Sirt6_IR30_mod.csv')

###IR2h
metadata <- read.csv('./SensorID_Node_DataFrame.csv')#Load the Node metadata table
metadata<-subset(metadata,metadata$Sensor=="Sirt6")#Isolate the node information corresponding to Sirt6ID networks
df <- read.table('./S6_IR2h_Edges.tsv', header=T) #Load the edge list table
g <- graph_from_data_frame(df, directed = FALSE) #generate an igraph objetct
g <- simplify(g, remove.multiple=TRUE, remove.loops=TRUE, edge.attr.comb = "first") #remove duplicated edges, and self-loops
edge.attributes(g) #Check no NAs in the edge attributes
adj_mat_g <- as_adjacency_matrix(g, attr = 'weight', sparse = FALSE) #convert the igraph object into an adjacency matrix
mod <- netcarto(adj_mat_g) #generate the modules using the netcarto algorithm
#Now generate a function to run the netcarto algorithm 20 times, to improve the modularity results from the first run
for(i in 1:20){
  temp_mod <- netcarto(adj_mat_g)
  if(temp_mod[[2]] > mod[[2]])
    mod <- temp_mod
}

# integrate the module results from netcarto to the igraph network 
com_mod <- make_clusters(
  g,
  membership = mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1,
  algorithm = 'netcarto',
  merges = NULL,
  modularity = mod[[2]]
)
#Generate a color palette for the modules
V(g)$module <- mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1
V(g)$color <- custom_palette[V(g)$module]
#generate an color palette for the functional annotation of each node
mod[[1]]$annotation <- as.factor(metadata[match(mod[[1]]$name, metadata$NAME),]$Annotation)
V(g)$color_border <- functions_colours[as.numeric(mod[[1]][match(V(g)$name, mod[[1]]$name),]$annotation)]

## Generate a weight column that includes the detected modules as an edge attribute
cl_for_layout <- com_mod
E(g)$weight <- edge.weights(cl_for_layout, g)

##### ATTENTION : THE CODE BELOW HAS TO BE RUN AT LEAST ONCE - UNCOMMENT ON FIRST RUN

# # create layout - run this line several times until you find a layout that you like
coord_graphopt <- layout_with_graphopt(g, mass = 20, charge = E(g)$weight)
igraph::plot.igraph(g,layout=coord_graphopt)
#After generating a good layout this can be saved for ploting all the other Sirt6ID networks 
## and we store the layout coordinates to plot the other networks in the same way
coords_nodes <- coord_graphopt
row.names(coords_nodes) <- V(g)$name
coords_nodes
idxs <- match(V(g)$name, row.names(coords_nodes))
coord_graphopt <- coords_nodes[idxs,]
vsize <- igraph::degree(g)/vcount(g)

#Plot the networks with the module as the fill color, and the annotation as the color of the circle line
par(mar=c(0,0,0,0))
pdf(paste0('Sirt6_protein-protein-network-whole.pdf'), width=30, height=30)
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('M = ', format(mod[[2]], digits=2)))
#
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.shape='circle', vertex.frame.color=V(g)$color_border, vertex.frame.width=4, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('whole network M = ', format(mod[[2]], digits = 2)), edge.arrow.size=.2)

legend('bottomleft', legend=levels(mod[[1]]$annotation), pt.lwd=4, cex=2, col=functions_colours, pch=21, pt.bg='white', title='function')

#recover the module results into a new table
cur_mod <- mod[[1]]
cur_mod$name <- as.character(cur_mod$name)
cur_mod$module <- cur_mod$module + 1
cur_mod$color <- custom_palette[cur_mod$module]
cur_mod$treatment <- 'IR2h'
cur_mod$degree <- degree(g, mode='all')[match(cur_mod$name, names(degree(g, mode='all')))]

#Plot the  Participation coefficient versus the degree 
par(mar=c(4,4,4,4))
ggplot(cur_mod, aes(degree, participation, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('degree-part-IR2h.pdf'), width = 12, height = 8)

# Plot the connectivity (within-module degree) versus the degree
ggplot(cur_mod, aes(degree, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('degree-conn-IR2h.pdf'), width = 12, height = 8)

#Plot the connectivity versus the participation
ggplot(cur_mod, aes(participation, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('part-conn-IR2h.pdf'), width = 12, height = 8)

#Print the results in a new table to merge them with the Node data frames
mod_overall <- cur_mod
write.csv(mod_overall,'./Sirt6_IR2h_mod.csv')

###IR8h
metadata <- read.csv('./SensorID_Node_DataFrame.csv')#Load the Node metadata table
metadata<-subset(metadata,metadata$Sensor=="Sirt6")#Isolate the node information corresponding to Sirt6ID networks
df <- read.table('./S6_IR8h_Edges.tsv', header=T) #Load the edge list table
g <- graph_from_data_frame(df, directed = FALSE) #generate an igraph objetct
g <- simplify(g, remove.multiple=TRUE, remove.loops=TRUE, edge.attr.comb = "first") #remove duplicated edges, and self-loops
edge.attributes(g) #Check no NAs in the edge attributes
adj_mat_g <- as_adjacency_matrix(g, attr = 'weight', sparse = FALSE) #convert the igraph object into an adjacency matrix
mod <- netcarto(adj_mat_g) #generate the modules using the netcarto algorithm
#Now generate a function to run the netcarto algorithm 20 times, to improve the modularity results from the first run
for(i in 1:20){
  temp_mod <- netcarto(adj_mat_g)
  if(temp_mod[[2]] > mod[[2]])
    mod <- temp_mod
}

# integrate the module results from netcarto to the igraph network 
com_mod <- make_clusters(
  g,
  membership = mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1,
  algorithm = 'netcarto',
  merges = NULL,
  modularity = mod[[2]]
)
#Generate a color palette for the modules
V(g)$module <- mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1
V(g)$color <- custom_palette[V(g)$module]
#generate an color palette for the functional annotation of each node
mod[[1]]$annotation <- as.factor(metadata[match(mod[[1]]$name, metadata$NAME),]$Annotation)
V(g)$color_border <- functions_colours[as.numeric(mod[[1]][match(V(g)$name, mod[[1]]$name),]$annotation)]

## Generate a weight column that includes the detected modules as an edge attribute
cl_for_layout <- com_mod
E(g)$weight <- edge.weights(cl_for_layout, g)

##### ATTENTION : THE CODE BELOW HAS TO BE RUN AT LEAST ONCE - UNCOMMENT ON FIRST RUN

# # create layout - run this line several times until you find a layout that you like
coord_graphopt <- layout_with_graphopt(g, mass = 20, charge = E(g)$weight)
igraph::plot.igraph(g,layout=coord_graphopt)
#After generating a good layout this can be saved for ploting all the other Sirt6ID networks 
## and we store the layout coordinates to plot the other networks in the same way
coords_nodes <- coord_graphopt
row.names(coords_nodes) <- V(g)$name
coords_nodes
idxs <- match(V(g)$name, row.names(coords_nodes))
coord_graphopt <- coords_nodes[idxs,]
vsize <- igraph::degree(g)/vcount(g)

#Plot the networks with the module as the fill color, and the annotation as the color of the circle line
par(mar=c(0,0,0,0))
pdf(paste0('Sirt6_protein-protein-network-whole.pdf'), width=30, height=30)
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('M = ', format(mod[[2]], digits=2)))
#
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.shape='circle', vertex.frame.color=V(g)$color_border, vertex.frame.width=4, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('whole network M = ', format(mod[[2]], digits = 2)), edge.arrow.size=.2)

legend('bottomleft', legend=levels(mod[[1]]$annotation), pt.lwd=4, cex=2, col=functions_colours, pch=21, pt.bg='white', title='function')

#recover the module results into a new table
cur_mod <- mod[[1]]
cur_mod$name <- as.character(cur_mod$name)
cur_mod$module <- cur_mod$module + 1
cur_mod$color <- custom_palette[cur_mod$module]
cur_mod$treatment <- 'IR8h'
cur_mod$degree <- degree(g, mode='all')[match(cur_mod$name, names(degree(g, mode='all')))]

#Plot the  Participation coefficient versus the degree 
par(mar=c(4,4,4,4))
ggplot(cur_mod, aes(degree, participation, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('degree-part-IR8h.pdf'), width = 12, height = 8)

# Plot the connectivity (within-module degree) versus the degree
ggplot(cur_mod, aes(degree, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('degree-conn-IR8h.pdf'), width = 12, height = 8)

#Plot the connectivity versus the participation
ggplot(cur_mod, aes(participation, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('part-conn-IR8h.pdf'), width = 12, height = 8)

#Print the results in a new table to merge them with the Node data frames
mod_overall <- cur_mod
write.csv(mod_overall,'./Sirt6_IR8h_mod.csv')

###IR24h
metadata <- read.csv('./SensorID_Node_DataFrame.csv')#Load the Node metadata table
metadata<-subset(metadata,metadata$Sensor=="Sirt6")#Isolate the node information corresponding to Sirt6ID networks
df <- read.table('./S6_IR24h_Edges.tsv', header=T) #Load the edge list table
g <- graph_from_data_frame(df, directed = FALSE) #generate an igraph objetct
g <- simplify(g, remove.multiple=TRUE, remove.loops=TRUE, edge.attr.comb = "first") #remove duplicated edges, and self-loops
edge.attributes(g) #Check no NAs in the edge attributes
adj_mat_g <- as_adjacency_matrix(g, attr = 'weight', sparse = FALSE) #convert the igraph object into an adjacency matrix
mod <- netcarto(adj_mat_g) #generate the modules using the netcarto algorithm
#Now generate a function to run the netcarto algorithm 20 times, to improve the modularity results from the first run
for(i in 1:20){
  temp_mod <- netcarto(adj_mat_g)
  if(temp_mod[[2]] > mod[[2]])
    mod <- temp_mod
}

# integrate the module results from netcarto to the igraph network 
com_mod <- make_clusters(
  g,
  membership = mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1,
  algorithm = 'netcarto',
  merges = NULL,
  modularity = mod[[2]]
)
#Generate a color palette for the modules
V(g)$module <- mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1
V(g)$color <- custom_palette[V(g)$module]
#generate an color palette for the functional annotation of each node
mod[[1]]$annotation <- as.factor(metadata[match(mod[[1]]$name, metadata$NAME),]$Annotation)
V(g)$color_border <- functions_colours[as.numeric(mod[[1]][match(V(g)$name, mod[[1]]$name),]$annotation)]

## Generate a weight column that includes the detected modules as an edge attribute
cl_for_layout <- com_mod
E(g)$weight <- edge.weights(cl_for_layout, g)

##### ATTENTION : THE CODE BELOW HAS TO BE RUN AT LEAST ONCE - UNCOMMENT ON FIRST RUN

# # create layout - run this line several times until you find a layout that you like
coord_graphopt <- layout_with_graphopt(g, mass = 20, charge = E(g)$weight)
igraph::plot.igraph(g,layout=coord_graphopt)
#After generating a good layout this can be saved for ploting all the other Sirt6ID networks 
## and we store the layout coordinates to plot the other networks in the same way
coords_nodes <- coord_graphopt
row.names(coords_nodes) <- V(g)$name
coords_nodes
idxs <- match(V(g)$name, row.names(coords_nodes))
coord_graphopt <- coords_nodes[idxs,]
vsize <- igraph::degree(g)/vcount(g)

#Plot the networks with the module as the fill color, and the annotation as the color of the circle line
par(mar=c(0,0,0,0))
pdf(paste0('Sirt6_protein-protein-network-whole.pdf'), width=30, height=30)
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('M = ', format(mod[[2]], digits=2)))
#
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.shape='circle', vertex.frame.color=V(g)$color_border, vertex.frame.width=4, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('whole network M = ', format(mod[[2]], digits = 2)), edge.arrow.size=.2)

legend('bottomleft', legend=levels(mod[[1]]$annotation), pt.lwd=4, cex=2, col=functions_colours, pch=21, pt.bg='white', title='function')

#recover the module results into a new table
cur_mod <- mod[[1]]
cur_mod$name <- as.character(cur_mod$name)
cur_mod$module <- cur_mod$module + 1
cur_mod$color <- custom_palette[cur_mod$module]
cur_mod$treatment <- 'IR24h'
cur_mod$degree <- degree(g, mode='all')[match(cur_mod$name, names(degree(g, mode='all')))]

#Plot the  Participation coefficient versus the degree 
par(mar=c(4,4,4,4))
ggplot(cur_mod, aes(degree, participation, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('degree-part-IR24h.pdf'), width = 12, height = 8)

# Plot the connectivity (within-module degree) versus the degree
ggplot(cur_mod, aes(degree, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('degree-conn-IR24h.pdf'), width = 12, height = 8)

#Plot the connectivity versus the participation
ggplot(cur_mod, aes(participation, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('part-conn-IR24h.pdf'), width = 12, height = 8)

#Print the results in a new table to merge them with the Node data frames
mod_overall <- cur_mod
write.csv(mod_overall,'./Sirt6_IR24h_mod.csv')

###nonIR24
metadata <- read.csv('./SensorID_Node_DataFrame.csv')#Load the Node metadata table
metadata<-subset(metadata,metadata$Sensor=="Sirt6")#Isolate the node information corresponding to Sirt6ID networks
df <- read.table('./S6_nonIR24_Edges.tsv', header=T) #Load the edge list table
g <- graph_from_data_frame(df, directed = FALSE) #generate an igraph objetct
g <- simplify(g, remove.multiple=TRUE, remove.loops=TRUE, edge.attr.comb = "first") #remove duplicated edges, and self-loops
edge.attributes(g) #Check no NAs in the edge attributes
adj_mat_g <- as_adjacency_matrix(g, attr = 'weight', sparse = FALSE) #convert the igraph object into an adjacency matrix
mod <- netcarto(adj_mat_g) #generate the modules using the netcarto algorithm
#Now generate a function to run the netcarto algorithm 20 times, to improve the modularity results from the first run
for(i in 1:20){
  temp_mod <- netcarto(adj_mat_g)
  if(temp_mod[[2]] > mod[[2]])
    mod <- temp_mod
}

# integrate the module results from netcarto to the igraph network 
com_mod <- make_clusters(
  g,
  membership = mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1,
  algorithm = 'netcarto',
  merges = NULL,
  modularity = mod[[2]]
)
#Generate a color palette for the modules
V(g)$module <- mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1
V(g)$color <- custom_palette[V(g)$module]
#generate an color palette for the functional annotation of each node
mod[[1]]$annotation <- as.factor(metadata[match(mod[[1]]$name, metadata$NAME),]$Annotation)
V(g)$color_border <- functions_colours[as.numeric(mod[[1]][match(V(g)$name, mod[[1]]$name),]$annotation)]

## Generate a weight column that includes the detected modules as an edge attribute
cl_for_layout <- com_mod
E(g)$weight <- edge.weights(cl_for_layout, g)

##### ATTENTION : THE CODE BELOW HAS TO BE RUN AT LEAST ONCE - UNCOMMENT ON FIRST RUN

# # create layout - run this line several times until you find a layout that you like
coord_graphopt <- layout_with_graphopt(g, mass = 20, charge = E(g)$weight)
igraph::plot.igraph(g,layout=coord_graphopt)
#After generating a good layout this can be saved for ploting all the other Sirt6ID networks 
## and we store the layout coordinates to plot the other networks in the same way
coords_nodes <- coord_graphopt
row.names(coords_nodes) <- V(g)$name
coords_nodes
idxs <- match(V(g)$name, row.names(coords_nodes))
coord_graphopt <- coords_nodes[idxs,]
vsize <- igraph::degree(g)/vcount(g)

#Plot the networks with the module as the fill color, and the annotation as the color of the circle line
par(mar=c(0,0,0,0))
pdf(paste0('Sirt6_protein-protein-network-whole.pdf'), width=30, height=30)
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('M = ', format(mod[[2]], digits=2)))
#
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.shape='circle', vertex.frame.color=V(g)$color_border, vertex.frame.width=4, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('whole network M = ', format(mod[[2]], digits = 2)), edge.arrow.size=.2)

legend('bottomleft', legend=levels(mod[[1]]$annotation), pt.lwd=4, cex=2, col=functions_colours, pch=21, pt.bg='white', title='function')

#recover the module results into a new table
cur_mod <- mod[[1]]
cur_mod$name <- as.character(cur_mod$name)
cur_mod$module <- cur_mod$module + 1
cur_mod$color <- custom_palette[cur_mod$module]
cur_mod$treatment <- 'nonIR24'
cur_mod$degree <- degree(g, mode='all')[match(cur_mod$name, names(degree(g, mode='all')))]

#Plot the  Participation coefficient versus the degree 
par(mar=c(4,4,4,4))
ggplot(cur_mod, aes(degree, participation, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('degree-part-nonIR24.pdf'), width = 12, height = 8)

# Plot the connectivity (within-module degree) versus the degree
ggplot(cur_mod, aes(degree, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('degree-conn-nonIR24.pdf'), width = 12, height = 8)

#Plot the connectivity versus the participation
ggplot(cur_mod, aes(participation, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('part-conn-nonIR24.pdf'), width = 12, height = 8)

#Print the results in a new table to merge them with the Node data frames
mod_overall <- cur_mod
write.csv(mod_overall,'./Sirt6_nonIR24_mod.csv')

#### NOTE: The results of the modularity analysis were manually reviewed, and merged with previous results to generate a single table that resumes all the analysis done. 
##This table is the SensorID_Node_DataFrame.csv. Once this table was generated, the connectivity versus participation plots, from Figure 5 can be easily done in ggplot.

###Ku80ID networks
### Full Ku80ID Nuclear Network
#### read the data
metadata <- read.csv('./SensorID_Node_DataFrame.csv')#Load the Node metadata table
metadata<-subset(metadata,metadata$Sensor=="Ku80")#Isolate the node information corresponding to Ku80ID networks
df <- read.table('./Ku80_FullNetwork.tsv', header=T) #Load the edge list table
g <- graph_from_data_frame(df, directed = FALSE) #generate an igraph objetct
g <- simplify(g, remove.multiple=TRUE, remove.loops=TRUE, edge.attr.comb = "first") #remove duplicated edges, and self-loops
edge.attributes(g) #Check no NAs in the edge attributes
adj_mat_g <- as_adjacency_matrix(g, attr = 'weight', sparse = FALSE) #convert the igraph object into an adjacency matrix
mod <- netcarto(adj_mat_g) #generate the modules using the netcarto algorithm
#Now generate a function to run the netcarto algorithm 20 times, to improve the modularity results from the first run
for(i in 1:20){
  temp_mod <- netcarto(adj_mat_g)
  if(temp_mod[[2]] > mod[[2]])
    mod <- temp_mod
}

# integrate the module results from netcarto to the igraph network 
com_mod <- make_clusters(
  g,
  membership = mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1,
  algorithm = 'netcarto',
  merges = NULL,
  modularity = mod[[2]]
)
#Generate a color palette for the modules
V(g)$module <- mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1
V(g)$color <- custom_palette[V(g)$module]
#generate an color palette for the functional annotation of each node
mod[[1]]$annotation <- as.factor(metadata[match(mod[[1]]$name, metadata$NAME),]$Annotation)
V(g)$color_border <- functions_colours[as.numeric(mod[[1]][match(V(g)$name, mod[[1]]$name),]$annotation)]

## Generate a weight column that includes the detected modules as an edge attribute
cl_for_layout <- com_mod
E(g)$weight <- edge.weights(cl_for_layout, g)

##### ATTENTION : THE CODE BELOW HAS TO BE RUN AT LEAST ONCE - UNCOMMENT ON FIRST RUN

# # create layout - run this line several times until you find a layout that you like
coord_graphopt <- layout_with_graphopt(g, mass = 20, charge = E(g)$weight)
igraph::plot.igraph(g,layout=coord_graphopt)
#After generating a good layout this can be saved for ploting all the other Ku80ID networks 
## and we store the layout coordinates to plot the other networks in the same way
coords_nodes <- coord_graphopt
row.names(coords_nodes) <- V(g)$name
coords_nodes
idxs <- match(V(g)$name, row.names(coords_nodes))
coord_graphopt <- coords_nodes[idxs,]
vsize <- igraph::degree(g)/vcount(g)

#Plot the networks with the module as the fill color, and the annotation as the color of the circle line
par(mar=c(0,0,0,0))
pdf(paste0('Ku80_protein-protein-network-whole.pdf'), width=30, height=30)
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('M = ', format(mod[[2]], digits=2)))
#
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.shape='circle', vertex.frame.color=V(g)$color_border, vertex.frame.width=4, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('whole network M = ', format(mod[[2]], digits = 2)), edge.arrow.size=.2)

legend('bottomleft', legend=levels(mod[[1]]$annotation), pt.lwd=4, cex=2, col=functions_colours, pch=21, pt.bg='white', title='function')

#recover the module results into a new table
cur_mod <- mod[[1]]
cur_mod$name <- as.character(cur_mod$name)
cur_mod$module <- cur_mod$module + 1
cur_mod$color <- custom_palette[cur_mod$module]
cur_mod$treatment <- 'whole'
cur_mod$degree <- degree(g, mode='all')[match(cur_mod$name, names(degree(g, mode='all')))]

#Plot the  Participation coefficient versus the degree 
par(mar=c(4,4,4,4))
ggplot(cur_mod, aes(degree, participation, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('degree-part-whole.pdf'), width = 12, height = 8)

# Plot the connectivity (within-module degree) versus the degree
ggplot(cur_mod, aes(degree, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('degree-conn-whole.pdf'), width = 12, height = 8)

#Plot the connectivity versus the participation
ggplot(cur_mod, aes(participation, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('part-conn-whole.pdf'), width = 12, height = 8)

#Print the results in a new table to merge them with the Node data frames
mod_overall <- cur_mod
write.csv(mod_overall,'./Ku80_FullNetwork_overall_mod.csv')

########### Work with each one of the time-point networks
##########
###nonIR5
metadata <- read.csv('./SensorID_Node_DataFrame.csv')#Load the Node metadata table
metadata<-subset(metadata,metadata$Sensor=="Ku80")#Isolate the node information corresponding to Ku80ID networks
df <- read.table('./Ku_nonIR5_Edges.tsv', header=T) #Load the edge list table
g <- graph_from_data_frame(df, directed = FALSE) #generate an igraph objetct
g <- simplify(g, remove.multiple=TRUE, remove.loops=TRUE, edge.attr.comb = "first") #remove duplicated edges, and self-loops
edge.attributes(g) #Check no NAs in the edge attributes
adj_mat_g <- as_adjacency_matrix(g, attr = 'weight', sparse = FALSE) #convert the igraph object into an adjacency matrix
mod <- netcarto(adj_mat_g) #generate the modules using the netcarto algorithm
#Now generate a function to run the netcarto algorithm 20 times, to improve the modularity results from the first run
for(i in 1:20){
  temp_mod <- netcarto(adj_mat_g)
  if(temp_mod[[2]] > mod[[2]])
    mod <- temp_mod
}

# integrate the module results from netcarto to the igraph network 
com_mod <- make_clusters(
  g,
  membership = mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1,
  algorithm = 'netcarto',
  merges = NULL,
  modularity = mod[[2]]
)
#Generate a color palette for the modules
V(g)$module <- mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1
V(g)$color <- custom_palette[V(g)$module]
#generate an color palette for the functional annotation of each node
mod[[1]]$annotation <- as.factor(metadata[match(mod[[1]]$name, metadata$NAME),]$Annotation)
V(g)$color_border <- functions_colours[as.numeric(mod[[1]][match(V(g)$name, mod[[1]]$name),]$annotation)]

## Generate a weight column that includes the detected modules as an edge attribute
cl_for_layout <- com_mod
E(g)$weight <- edge.weights(cl_for_layout, g)

##### ATTENTION : THE CODE BELOW HAS TO BE RUN AT LEAST ONCE - UNCOMMENT ON FIRST RUN

# # create layout - run this line several times until you find a layout that you like
coord_graphopt <- layout_with_graphopt(g, mass = 20, charge = E(g)$weight)
igraph::plot.igraph(g,layout=coord_graphopt)
#After generating a good layout this can be saved for ploting all the other Ku80ID networks 
## and we store the layout coordinates to plot the other networks in the same way
coords_nodes <- coord_graphopt
row.names(coords_nodes) <- V(g)$name
coords_nodes
idxs <- match(V(g)$name, row.names(coords_nodes))
coord_graphopt <- coords_nodes[idxs,]
vsize <- igraph::degree(g)/vcount(g)

#Plot the networks with the module as the fill color, and the annotation as the color of the circle line
par(mar=c(0,0,0,0))
pdf(paste0('Ku80_protein-protein-network-whole.pdf'), width=30, height=30)
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('M = ', format(mod[[2]], digits=2)))
#
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.shape='circle', vertex.frame.color=V(g)$color_border, vertex.frame.width=4, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('whole network M = ', format(mod[[2]], digits = 2)), edge.arrow.size=.2)

legend('bottomleft', legend=levels(mod[[1]]$annotation), pt.lwd=4, cex=2, col=functions_colours, pch=21, pt.bg='white', title='function')

#recover the module results into a new table
cur_mod <- mod[[1]]
cur_mod$name <- as.character(cur_mod$name)
cur_mod$module <- cur_mod$module + 1
cur_mod$color <- custom_palette[cur_mod$module]
cur_mod$treatment <- 'nonIR5'
cur_mod$degree <- degree(g, mode='all')[match(cur_mod$name, names(degree(g, mode='all')))]

#Plot the  Participation coefficient versus the degree 
par(mar=c(4,4,4,4))
ggplot(cur_mod, aes(degree, participation, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('degree-part-nonIR5.pdf'), width = 12, height = 8)

# Plot the connectivity (within-module degree) versus the degree
ggplot(cur_mod, aes(degree, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('degree-conn-nonIR5.pdf'), width = 12, height = 8)

#Plot the connectivity versus the participation
ggplot(cur_mod, aes(participation, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('part-conn-nonIR5.pdf'), width = 12, height = 8)

#Print the results in a new table to merge them with the Node data frames
mod_overall <- cur_mod
write.csv(mod_overall,'./Ku80_nonIR5_mod.csv')

###IR5
metadata <- read.csv('./SensorID_Node_DataFrame.csv')#Load the Node metadata table
metadata<-subset(metadata,metadata$Sensor=="Ku80")#Isolate the node information corresponding to Ku80ID networks
df <- read.table('./Ku_IR5_Edges.tsv', header=T) #Load the edge list table
g <- graph_from_data_frame(df, directed = FALSE) #generate an igraph objetct
g <- simplify(g, remove.multiple=TRUE, remove.loops=TRUE, edge.attr.comb = "first") #remove duplicated edges, and self-loops
edge.attributes(g) #Check no NAs in the edge attributes
adj_mat_g <- as_adjacency_matrix(g, attr = 'weight', sparse = FALSE) #convert the igraph object into an adjacency matrix
mod <- netcarto(adj_mat_g) #generate the modules using the netcarto algorithm
#Now generate a function to run the netcarto algorithm 20 times, to improve the modularity results from the first run
for(i in 1:20){
  temp_mod <- netcarto(adj_mat_g)
  if(temp_mod[[2]] > mod[[2]])
    mod <- temp_mod
}

# integrate the module results from netcarto to the igraph network 
com_mod <- make_clusters(
  g,
  membership = mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1,
  algorithm = 'netcarto',
  merges = NULL,
  modularity = mod[[2]]
)
#Generate a color palette for the modules
V(g)$module <- mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1
V(g)$color <- custom_palette[V(g)$module]
#generate an color palette for the functional annotation of each node
mod[[1]]$annotation <- as.factor(metadata[match(mod[[1]]$name, metadata$NAME),]$Annotation)
V(g)$color_border <- functions_colours[as.numeric(mod[[1]][match(V(g)$name, mod[[1]]$name),]$annotation)]

## Generate a weight column that includes the detected modules as an edge attribute
cl_for_layout <- com_mod
E(g)$weight <- edge.weights(cl_for_layout, g)

##### ATTENTION : THE CODE BELOW HAS TO BE RUN AT LEAST ONCE - UNCOMMENT ON FIRST RUN

# # create layout - run this line several times until you find a layout that you like
coord_graphopt <- layout_with_graphopt(g, mass = 20, charge = E(g)$weight)
igraph::plot.igraph(g,layout=coord_graphopt)
#After generating a good layout this can be saved for ploting all the other Ku80ID networks 
## and we store the layout coordinates to plot the other networks in the same way
coords_nodes <- coord_graphopt
row.names(coords_nodes) <- V(g)$name
coords_nodes
idxs <- match(V(g)$name, row.names(coords_nodes))
coord_graphopt <- coords_nodes[idxs,]
vsize <- igraph::degree(g)/vcount(g)

#Plot the networks with the module as the fill color, and the annotation as the color of the circle line
par(mar=c(0,0,0,0))
pdf(paste0('Ku80_protein-protein-network-whole.pdf'), width=30, height=30)
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('M = ', format(mod[[2]], digits=2)))
#
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.shape='circle', vertex.frame.color=V(g)$color_border, vertex.frame.width=4, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('whole network M = ', format(mod[[2]], digits = 2)), edge.arrow.size=.2)

legend('bottomleft', legend=levels(mod[[1]]$annotation), pt.lwd=4, cex=2, col=functions_colours, pch=21, pt.bg='white', title='function')

#recover the module results into a new table
cur_mod <- mod[[1]]
cur_mod$name <- as.character(cur_mod$name)
cur_mod$module <- cur_mod$module + 1
cur_mod$color <- custom_palette[cur_mod$module]
cur_mod$treatment <- 'IR5'
cur_mod$degree <- degree(g, mode='all')[match(cur_mod$name, names(degree(g, mode='all')))]

#Plot the  Participation coefficient versus the degree 
par(mar=c(4,4,4,4))
ggplot(cur_mod, aes(degree, participation, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('degree-part-IR5.pdf'), width = 12, height = 8)

# Plot the connectivity (within-module degree) versus the degree
ggplot(cur_mod, aes(degree, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('degree-conn-IR5.pdf'), width = 12, height = 8)

#Plot the connectivity versus the participation
ggplot(cur_mod, aes(participation, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('part-conn-IR5.pdf'), width = 12, height = 8)

#Print the results in a new table to merge them with the Node data frames
mod_overall <- cur_mod
write.csv(mod_overall,'./Ku80_IR5_mod.csv')

###IR30
metadata <- read.csv('./SensorID_Node_DataFrame.csv')#Load the Node metadata table
metadata<-subset(metadata,metadata$Sensor=="Ku80")#Isolate the node information corresponding to Ku80ID networks
df <- read.table('./Ku_IR30_Edges.tsv', header=T) #Load the edge list table
g <- graph_from_data_frame(df, directed = FALSE) #generate an igraph objetct
g <- simplify(g, remove.multiple=TRUE, remove.loops=TRUE, edge.attr.comb = "first") #remove duplicated edges, and self-loops
edge.attributes(g) #Check no NAs in the edge attributes
adj_mat_g <- as_adjacency_matrix(g, attr = 'weight', sparse = FALSE) #convert the igraph object into an adjacency matrix
mod <- netcarto(adj_mat_g) #generate the modules using the netcarto algorithm
#Now generate a function to run the netcarto algorithm 20 times, to improve the modularity results from the first run
for(i in 1:20){
  temp_mod <- netcarto(adj_mat_g)
  if(temp_mod[[2]] > mod[[2]])
    mod <- temp_mod
}

# integrate the module results from netcarto to the igraph network 
com_mod <- make_clusters(
  g,
  membership = mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1,
  algorithm = 'netcarto',
  merges = NULL,
  modularity = mod[[2]]
)
#Generate a color palette for the modules
V(g)$module <- mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1
V(g)$color <- custom_palette[V(g)$module]
#generate an color palette for the functional annotation of each node
mod[[1]]$annotation <- as.factor(metadata[match(mod[[1]]$name, metadata$NAME),]$Annotation)
V(g)$color_border <- functions_colours[as.numeric(mod[[1]][match(V(g)$name, mod[[1]]$name),]$annotation)]

## Generate a weight column that includes the detected modules as an edge attribute
cl_for_layout <- com_mod
E(g)$weight <- edge.weights(cl_for_layout, g)

##### ATTENTION : THE CODE BELOW HAS TO BE RUN AT LEAST ONCE - UNCOMMENT ON FIRST RUN

# # create layout - run this line several times until you find a layout that you like
coord_graphopt <- layout_with_graphopt(g, mass = 20, charge = E(g)$weight)
igraph::plot.igraph(g,layout=coord_graphopt)
#After generating a good layout this can be saved for ploting all the other Ku80ID networks 
## and we store the layout coordinates to plot the other networks in the same way
coords_nodes <- coord_graphopt
row.names(coords_nodes) <- V(g)$name
coords_nodes
idxs <- match(V(g)$name, row.names(coords_nodes))
coord_graphopt <- coords_nodes[idxs,]
vsize <- igraph::degree(g)/vcount(g)

#Plot the networks with the module as the fill color, and the annotation as the color of the circle line
par(mar=c(0,0,0,0))
pdf(paste0('Ku80_protein-protein-network-whole.pdf'), width=30, height=30)
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('M = ', format(mod[[2]], digits=2)))
#
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.shape='circle', vertex.frame.color=V(g)$color_border, vertex.frame.width=4, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('whole network M = ', format(mod[[2]], digits = 2)), edge.arrow.size=.2)

legend('bottomleft', legend=levels(mod[[1]]$annotation), pt.lwd=4, cex=2, col=functions_colours, pch=21, pt.bg='white', title='function')

#recover the module results into a new table
cur_mod <- mod[[1]]
cur_mod$name <- as.character(cur_mod$name)
cur_mod$module <- cur_mod$module + 1
cur_mod$color <- custom_palette[cur_mod$module]
cur_mod$treatment <- 'IR30'
cur_mod$degree <- degree(g, mode='all')[match(cur_mod$name, names(degree(g, mode='all')))]

#Plot the  Participation coefficient versus the degree 
par(mar=c(4,4,4,4))
ggplot(cur_mod, aes(degree, participation, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('degree-part-IR30.pdf'), width = 12, height = 8)

# Plot the connectivity (within-module degree) versus the degree
ggplot(cur_mod, aes(degree, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('degree-conn-IR30.pdf'), width = 12, height = 8)

#Plot the connectivity versus the participation
ggplot(cur_mod, aes(participation, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('part-conn-IR30.pdf'), width = 12, height = 8)

#Print the results in a new table to merge them with the Node data frames
mod_overall <- cur_mod
write.csv(mod_overall,'./Ku80_IR30_mod.csv')

###IR2h
metadata <- read.csv('./SensorID_Node_DataFrame.csv')#Load the Node metadata table
metadata<-subset(metadata,metadata$Sensor=="Ku80")#Isolate the node information corresponding to Ku80ID networks
df <- read.table('./Ku_IR2h_Edges.tsv', header=T) #Load the edge list table
g <- graph_from_data_frame(df, directed = FALSE) #generate an igraph objetct
g <- simplify(g, remove.multiple=TRUE, remove.loops=TRUE, edge.attr.comb = "first") #remove duplicated edges, and self-loops
edge.attributes(g) #Check no NAs in the edge attributes
adj_mat_g <- as_adjacency_matrix(g, attr = 'weight', sparse = FALSE) #convert the igraph object into an adjacency matrix
mod <- netcarto(adj_mat_g) #generate the modules using the netcarto algorithm
#Now generate a function to run the netcarto algorithm 20 times, to improve the modularity results from the first run
for(i in 1:20){
  temp_mod <- netcarto(adj_mat_g)
  if(temp_mod[[2]] > mod[[2]])
    mod <- temp_mod
}

# integrate the module results from netcarto to the igraph network 
com_mod <- make_clusters(
  g,
  membership = mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1,
  algorithm = 'netcarto',
  merges = NULL,
  modularity = mod[[2]]
)
#Generate a color palette for the modules
V(g)$module <- mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1
V(g)$color <- custom_palette[V(g)$module]
#generate an color palette for the functional annotation of each node
mod[[1]]$annotation <- as.factor(metadata[match(mod[[1]]$name, metadata$NAME),]$Annotation)
V(g)$color_border <- functions_colours[as.numeric(mod[[1]][match(V(g)$name, mod[[1]]$name),]$annotation)]

## Generate a weight column that includes the detected modules as an edge attribute
cl_for_layout <- com_mod
E(g)$weight <- edge.weights(cl_for_layout, g)

##### ATTENTION : THE CODE BELOW HAS TO BE RUN AT LEAST ONCE - UNCOMMENT ON FIRST RUN

# # create layout - run this line several times until you find a layout that you like
coord_graphopt <- layout_with_graphopt(g, mass = 20, charge = E(g)$weight)
igraph::plot.igraph(g,layout=coord_graphopt)
#After generating a good layout this can be saved for ploting all the other Ku80ID networks 
## and we store the layout coordinates to plot the other networks in the same way
coords_nodes <- coord_graphopt
row.names(coords_nodes) <- V(g)$name
coords_nodes
idxs <- match(V(g)$name, row.names(coords_nodes))
coord_graphopt <- coords_nodes[idxs,]
vsize <- igraph::degree(g)/vcount(g)

#Plot the networks with the module as the fill color, and the annotation as the color of the circle line
par(mar=c(0,0,0,0))
pdf(paste0('Ku80_protein-protein-network-whole.pdf'), width=30, height=30)
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('M = ', format(mod[[2]], digits=2)))
#
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.shape='circle', vertex.frame.color=V(g)$color_border, vertex.frame.width=4, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('whole network M = ', format(mod[[2]], digits = 2)), edge.arrow.size=.2)

legend('bottomleft', legend=levels(mod[[1]]$annotation), pt.lwd=4, cex=2, col=functions_colours, pch=21, pt.bg='white', title='function')

#recover the module results into a new table
cur_mod <- mod[[1]]
cur_mod$name <- as.character(cur_mod$name)
cur_mod$module <- cur_mod$module + 1
cur_mod$color <- custom_palette[cur_mod$module]
cur_mod$treatment <- 'IR2h'
cur_mod$degree <- degree(g, mode='all')[match(cur_mod$name, names(degree(g, mode='all')))]

#Plot the  Participation coefficient versus the degree 
par(mar=c(4,4,4,4))
ggplot(cur_mod, aes(degree, participation, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('degree-part-IR2h.pdf'), width = 12, height = 8)

# Plot the connectivity (within-module degree) versus the degree
ggplot(cur_mod, aes(degree, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('degree-conn-IR2h.pdf'), width = 12, height = 8)

#Plot the connectivity versus the participation
ggplot(cur_mod, aes(participation, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('part-conn-IR2h.pdf'), width = 12, height = 8)

#Print the results in a new table to merge them with the Node data frames
mod_overall <- cur_mod
write.csv(mod_overall,'./Ku80_IR2h_mod.csv')

###IR8h
metadata <- read.csv('./SensorID_Node_DataFrame.csv')#Load the Node metadata table
metadata<-subset(metadata,metadata$Sensor=="Ku80")#Isolate the node information corresponding to Ku80ID networks
df <- read.table('./Ku_IR8h_Edges.tsv', header=T) #Load the edge list table
g <- graph_from_data_frame(df, directed = FALSE) #generate an igraph objetct
g <- simplify(g, remove.multiple=TRUE, remove.loops=TRUE, edge.attr.comb = "first") #remove duplicated edges, and self-loops
edge.attributes(g) #Check no NAs in the edge attributes
adj_mat_g <- as_adjacency_matrix(g, attr = 'weight', sparse = FALSE) #convert the igraph object into an adjacency matrix
mod <- netcarto(adj_mat_g) #generate the modules using the netcarto algorithm
#Now generate a function to run the netcarto algorithm 20 times, to improve the modularity results from the first run
for(i in 1:20){
  temp_mod <- netcarto(adj_mat_g)
  if(temp_mod[[2]] > mod[[2]])
    mod <- temp_mod
}

# integrate the module results from netcarto to the igraph network 
com_mod <- make_clusters(
  g,
  membership = mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1,
  algorithm = 'netcarto',
  merges = NULL,
  modularity = mod[[2]]
)
#Generate a color palette for the modules
V(g)$module <- mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1
V(g)$color <- custom_palette[V(g)$module]
#generate an color palette for the functional annotation of each node
mod[[1]]$annotation <- as.factor(metadata[match(mod[[1]]$name, metadata$NAME),]$Annotation)
V(g)$color_border <- functions_colours[as.numeric(mod[[1]][match(V(g)$name, mod[[1]]$name),]$annotation)]

## Generate a weight column that includes the detected modules as an edge attribute
cl_for_layout <- com_mod
E(g)$weight <- edge.weights(cl_for_layout, g)

##### ATTENTION : THE CODE BELOW HAS TO BE RUN AT LEAST ONCE - UNCOMMENT ON FIRST RUN

# # create layout - run this line several times until you find a layout that you like
coord_graphopt <- layout_with_graphopt(g, mass = 20, charge = E(g)$weight)
igraph::plot.igraph(g,layout=coord_graphopt)
#After generating a good layout this can be saved for ploting all the other Ku80ID networks 
## and we store the layout coordinates to plot the other networks in the same way
coords_nodes <- coord_graphopt
row.names(coords_nodes) <- V(g)$name
coords_nodes
idxs <- match(V(g)$name, row.names(coords_nodes))
coord_graphopt <- coords_nodes[idxs,]
vsize <- igraph::degree(g)/vcount(g)

#Plot the networks with the module as the fill color, and the annotation as the color of the circle line
par(mar=c(0,0,0,0))
pdf(paste0('Ku80_protein-protein-network-whole.pdf'), width=30, height=30)
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('M = ', format(mod[[2]], digits=2)))
#
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.shape='circle', vertex.frame.color=V(g)$color_border, vertex.frame.width=4, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('whole network M = ', format(mod[[2]], digits = 2)), edge.arrow.size=.2)

legend('bottomleft', legend=levels(mod[[1]]$annotation), pt.lwd=4, cex=2, col=functions_colours, pch=21, pt.bg='white', title='function')

#recover the module results into a new table
cur_mod <- mod[[1]]
cur_mod$name <- as.character(cur_mod$name)
cur_mod$module <- cur_mod$module + 1
cur_mod$color <- custom_palette[cur_mod$module]
cur_mod$treatment <- 'IR8h'
cur_mod$degree <- degree(g, mode='all')[match(cur_mod$name, names(degree(g, mode='all')))]

#Plot the  Participation coefficient versus the degree 
par(mar=c(4,4,4,4))
ggplot(cur_mod, aes(degree, participation, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('degree-part-IR8h.pdf'), width = 12, height = 8)

# Plot the connectivity (within-module degree) versus the degree
ggplot(cur_mod, aes(degree, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('degree-conn-IR8h.pdf'), width = 12, height = 8)

#Plot the connectivity versus the participation
ggplot(cur_mod, aes(participation, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('part-conn-IR8h.pdf'), width = 12, height = 8)

#Print the results in a new table to merge them with the Node data frames
mod_overall <- cur_mod
write.csv(mod_overall,'./Ku80_IR8h_mod.csv')

###IR24h
metadata <- read.csv('./SensorID_Node_DataFrame.csv')#Load the Node metadata table
metadata<-subset(metadata,metadata$Sensor=="Ku80")#Isolate the node information corresponding to Ku80ID networks
df <- read.table('./Ku_IR24h_Edges.tsv', header=T) #Load the edge list table
g <- graph_from_data_frame(df, directed = FALSE) #generate an igraph objetct
g <- simplify(g, remove.multiple=TRUE, remove.loops=TRUE, edge.attr.comb = "first") #remove duplicated edges, and self-loops
edge.attributes(g) #Check no NAs in the edge attributes
adj_mat_g <- as_adjacency_matrix(g, attr = 'weight', sparse = FALSE) #convert the igraph object into an adjacency matrix
mod <- netcarto(adj_mat_g) #generate the modules using the netcarto algorithm
#Now generate a function to run the netcarto algorithm 20 times, to improve the modularity results from the first run
for(i in 1:20){
  temp_mod <- netcarto(adj_mat_g)
  if(temp_mod[[2]] > mod[[2]])
    mod <- temp_mod
}

# integrate the module results from netcarto to the igraph network 
com_mod <- make_clusters(
  g,
  membership = mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1,
  algorithm = 'netcarto',
  merges = NULL,
  modularity = mod[[2]]
)
#Generate a color palette for the modules
V(g)$module <- mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1
V(g)$color <- custom_palette[V(g)$module]
#generate an color palette for the functional annotation of each node
mod[[1]]$annotation <- as.factor(metadata[match(mod[[1]]$name, metadata$NAME),]$Annotation)
V(g)$color_border <- functions_colours[as.numeric(mod[[1]][match(V(g)$name, mod[[1]]$name),]$annotation)]

## Generate a weight column that includes the detected modules as an edge attribute
cl_for_layout <- com_mod
E(g)$weight <- edge.weights(cl_for_layout, g)

##### ATTENTION : THE CODE BELOW HAS TO BE RUN AT LEAST ONCE - UNCOMMENT ON FIRST RUN

# # create layout - run this line several times until you find a layout that you like
coord_graphopt <- layout_with_graphopt(g, mass = 20, charge = E(g)$weight)
igraph::plot.igraph(g,layout=coord_graphopt)
#After generating a good layout this can be saved for ploting all the other Ku80ID networks 
## and we store the layout coordinates to plot the other networks in the same way
coords_nodes <- coord_graphopt
row.names(coords_nodes) <- V(g)$name
coords_nodes
idxs <- match(V(g)$name, row.names(coords_nodes))
coord_graphopt <- coords_nodes[idxs,]
vsize <- igraph::degree(g)/vcount(g)

#Plot the networks with the module as the fill color, and the annotation as the color of the circle line
par(mar=c(0,0,0,0))
pdf(paste0('Ku80_protein-protein-network-whole.pdf'), width=30, height=30)
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('M = ', format(mod[[2]], digits=2)))
#
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.shape='circle', vertex.frame.color=V(g)$color_border, vertex.frame.width=4, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('whole network M = ', format(mod[[2]], digits = 2)), edge.arrow.size=.2)

legend('bottomleft', legend=levels(mod[[1]]$annotation), pt.lwd=4, cex=2, col=functions_colours, pch=21, pt.bg='white', title='function')

#recover the module results into a new table
cur_mod <- mod[[1]]
cur_mod$name <- as.character(cur_mod$name)
cur_mod$module <- cur_mod$module + 1
cur_mod$color <- custom_palette[cur_mod$module]
cur_mod$treatment <- 'IR24h'
cur_mod$degree <- degree(g, mode='all')[match(cur_mod$name, names(degree(g, mode='all')))]

#Plot the  Participation coefficient versus the degree 
par(mar=c(4,4,4,4))
ggplot(cur_mod, aes(degree, participation, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('degree-part-IR24h.pdf'), width = 12, height = 8)

# Plot the connectivity (within-module degree) versus the degree
ggplot(cur_mod, aes(degree, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('degree-conn-IR24h.pdf'), width = 12, height = 8)

#Plot the connectivity versus the participation
ggplot(cur_mod, aes(participation, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('part-conn-IR24h.pdf'), width = 12, height = 8)

#Print the results in a new table to merge them with the Node data frames
mod_overall <- cur_mod
write.csv(mod_overall,'./Ku80_IR24h_mod.csv')

###nonIR24
metadata <- read.csv('./SensorID_Node_DataFrame.csv')#Load the Node metadata table
metadata<-subset(metadata,metadata$Sensor=="Ku80")#Isolate the node information corresponding to Ku80ID networks
df <- read.table('./Ku_nonIR24_Edges.tsv', header=T) #Load the edge list table
g <- graph_from_data_frame(df, directed = FALSE) #generate an igraph objetct
g <- simplify(g, remove.multiple=TRUE, remove.loops=TRUE, edge.attr.comb = "first") #remove duplicated edges, and self-loops
edge.attributes(g) #Check no NAs in the edge attributes
adj_mat_g <- as_adjacency_matrix(g, attr = 'weight', sparse = FALSE) #convert the igraph object into an adjacency matrix
mod <- netcarto(adj_mat_g) #generate the modules using the netcarto algorithm
#Now generate a function to run the netcarto algorithm 20 times, to improve the modularity results from the first run
for(i in 1:20){
  temp_mod <- netcarto(adj_mat_g)
  if(temp_mod[[2]] > mod[[2]])
    mod <- temp_mod
}

# integrate the module results from netcarto to the igraph network 
com_mod <- make_clusters(
  g,
  membership = mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1,
  algorithm = 'netcarto',
  merges = NULL,
  modularity = mod[[2]]
)
#Generate a color palette for the modules
V(g)$module <- mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1
V(g)$color <- custom_palette[V(g)$module]
#generate an color palette for the functional annotation of each node
mod[[1]]$annotation <- as.factor(metadata[match(mod[[1]]$name, metadata$NAME),]$Annotation)
V(g)$color_border <- functions_colours[as.numeric(mod[[1]][match(V(g)$name, mod[[1]]$name),]$annotation)]

## Generate a weight column that includes the detected modules as an edge attribute
cl_for_layout <- com_mod
E(g)$weight <- edge.weights(cl_for_layout, g)

##### ATTENTION : THE CODE BELOW HAS TO BE RUN AT LEAST ONCE - UNCOMMENT ON FIRST RUN

# # create layout - run this line several times until you find a layout that you like
coord_graphopt <- layout_with_graphopt(g, mass = 20, charge = E(g)$weight)
igraph::plot.igraph(g,layout=coord_graphopt)
#After generating a good layout this can be saved for ploting all the other Ku80ID networks 
## and we store the layout coordinates to plot the other networks in the same way
coords_nodes <- coord_graphopt
row.names(coords_nodes) <- V(g)$name
coords_nodes
idxs <- match(V(g)$name, row.names(coords_nodes))
coord_graphopt <- coords_nodes[idxs,]
vsize <- igraph::degree(g)/vcount(g)

#Plot the networks with the module as the fill color, and the annotation as the color of the circle line
par(mar=c(0,0,0,0))
pdf(paste0('Ku80_protein-protein-network-whole.pdf'), width=30, height=30)
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('M = ', format(mod[[2]], digits=2)))
#
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.shape='circle', vertex.frame.color=V(g)$color_border, vertex.frame.width=4, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('whole network M = ', format(mod[[2]], digits = 2)), edge.arrow.size=.2)

legend('bottomleft', legend=levels(mod[[1]]$annotation), pt.lwd=4, cex=2, col=functions_colours, pch=21, pt.bg='white', title='function')

#recover the module results into a new table
cur_mod <- mod[[1]]
cur_mod$name <- as.character(cur_mod$name)
cur_mod$module <- cur_mod$module + 1
cur_mod$color <- custom_palette[cur_mod$module]
cur_mod$treatment <- 'nonIR24'
cur_mod$degree <- degree(g, mode='all')[match(cur_mod$name, names(degree(g, mode='all')))]

#Plot the  Participation coefficient versus the degree 
par(mar=c(4,4,4,4))
ggplot(cur_mod, aes(degree, participation, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('degree-part-nonIR24.pdf'), width = 12, height = 8)

# Plot the connectivity (within-module degree) versus the degree
ggplot(cur_mod, aes(degree, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('degree-conn-nonIR24.pdf'), width = 12, height = 8)

#Plot the connectivity versus the participation
ggplot(cur_mod, aes(participation, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('part-conn-nonIR24.pdf'), width = 12, height = 8)

#Print the results in a new table to merge them with the Node data frames
mod_overall <- cur_mod
write.csv(mod_overall,'./Ku80_nonIR24_mod.csv')

###Mre11ID networks
### Full Mre11ID Nuclear Network
#### read the data
metadata <- read.csv('./SensorID_Node_DataFrame.csv')#Load the Node metadata table
metadata<-subset(metadata,metadata$Sensor=="Mre11")#Isolate the node information corresponding to Mre11ID networks
df <- read.table('./Mre11_FullNetwork.tsv', header=T) #Load the edge list table
g <- graph_from_data_frame(df, directed = FALSE) #generate an igraph objetct
g <- simplify(g, remove.multiple=TRUE, remove.loops=TRUE, edge.attr.comb = "first") #remove duplicated edges, and self-loops
edge.attributes(g) #Check no NAs in the edge attributes
adj_mat_g <- as_adjacency_matrix(g, attr = 'weight', sparse = FALSE) #convert the igraph object into an adjacency matrix
mod <- netcarto(adj_mat_g) #generate the modules using the netcarto algorithm
#Now generate a function to run the netcarto algorithm 20 times, to improve the modularity results from the first run
for(i in 1:20){
  temp_mod <- netcarto(adj_mat_g)
  if(temp_mod[[2]] > mod[[2]])
    mod <- temp_mod
}

# integrate the module results from netcarto to the igraph network 
com_mod <- make_clusters(
  g,
  membership = mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1,
  algorithm = 'netcarto',
  merges = NULL,
  modularity = mod[[2]]
)
#Generate a color palette for the modules
V(g)$module <- mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1
V(g)$color <- custom_palette[V(g)$module]
#generate an color palette for the functional annotation of each node
mod[[1]]$annotation <- as.factor(metadata[match(mod[[1]]$name, metadata$NAME),]$Annotation)
V(g)$color_border <- functions_colours[as.numeric(mod[[1]][match(V(g)$name, mod[[1]]$name),]$annotation)]

## Generate a weight column that includes the detected modules as an edge attribute
cl_for_layout <- com_mod
E(g)$weight <- edge.weights(cl_for_layout, g)

##### ATTENTION : THE CODE BELOW HAS TO BE RUN AT LEAST ONCE - UNCOMMENT ON FIRST RUN

# # create layout - run this line several times until you find a layout that you like
coord_graphopt <- layout_with_graphopt(g, mass = 20, charge = E(g)$weight)
igraph::plot.igraph(g,layout=coord_graphopt)
#After generating a good layout this can be saved for ploting all the other Mre11ID networks 
## and we store the layout coordinates to plot the other networks in the same way
coords_nodes <- coord_graphopt
row.names(coords_nodes) <- V(g)$name
coords_nodes
idxs <- match(V(g)$name, row.names(coords_nodes))
coord_graphopt <- coords_nodes[idxs,]
vsize <- igraph::degree(g)/vcount(g)

#Plot the networks with the module as the fill color, and the annotation as the color of the circle line
par(mar=c(0,0,0,0))
pdf(paste0('Mre11_protein-protein-network-whole.pdf'), width=30, height=30)
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('M = ', format(mod[[2]], digits=2)))
#
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.shape='circle', vertex.frame.color=V(g)$color_border, vertex.frame.width=4, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('whole network M = ', format(mod[[2]], digits = 2)), edge.arrow.size=.2)

legend('bottomleft', legend=levels(mod[[1]]$annotation), pt.lwd=4, cex=2, col=functions_colours, pch=21, pt.bg='white', title='function')

#recover the module results into a new table
cur_mod <- mod[[1]]
cur_mod$name <- as.character(cur_mod$name)
cur_mod$module <- cur_mod$module + 1
cur_mod$color <- custom_palette[cur_mod$module]
cur_mod$treatment <- 'whole'
cur_mod$degree <- degree(g, mode='all')[match(cur_mod$name, names(degree(g, mode='all')))]

#Plot the  Participation coefficient versus the degree 
par(mar=c(4,4,4,4))
ggplot(cur_mod, aes(degree, participation, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('degree-part-whole.pdf'), width = 12, height = 8)

# Plot the connectivity (within-module degree) versus the degree
ggplot(cur_mod, aes(degree, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('degree-conn-whole.pdf'), width = 12, height = 8)

#Plot the connectivity versus the participation
ggplot(cur_mod, aes(participation, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('part-conn-whole.pdf'), width = 12, height = 8)

#Print the results in a new table to merge them with the Node data frames
mod_overall <- cur_mod
write.csv(mod_overall,'./Mre11_FullNetwork_overall_mod.csv')

########### Work with each one of the time-point networks
##########
###nonIR5
metadata <- read.csv('./SensorID_Node_DataFrame.csv')#Load the Node metadata table
metadata<-subset(metadata,metadata$Sensor=="Mre11")#Isolate the node information corresponding to Mre11ID networks
df <- read.table('./Mre_nonIR5_Edges.tsv', header=T) #Load the edge list table
g <- graph_from_data_frame(df, directed = FALSE) #generate an igraph objetct
g <- simplify(g, remove.multiple=TRUE, remove.loops=TRUE, edge.attr.comb = "first") #remove duplicated edges, and self-loops
edge.attributes(g) #Check no NAs in the edge attributes
adj_mat_g <- as_adjacency_matrix(g, attr = 'weight', sparse = FALSE) #convert the igraph object into an adjacency matrix
mod <- netcarto(adj_mat_g) #generate the modules using the netcarto algorithm
#Now generate a function to run the netcarto algorithm 20 times, to improve the modularity results from the first run
for(i in 1:20){
  temp_mod <- netcarto(adj_mat_g)
  if(temp_mod[[2]] > mod[[2]])
    mod <- temp_mod
}

# integrate the module results from netcarto to the igraph network 
com_mod <- make_clusters(
  g,
  membership = mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1,
  algorithm = 'netcarto',
  merges = NULL,
  modularity = mod[[2]]
)
#Generate a color palette for the modules
V(g)$module <- mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1
V(g)$color <- custom_palette[V(g)$module]
#generate an color palette for the functional annotation of each node
mod[[1]]$annotation <- as.factor(metadata[match(mod[[1]]$name, metadata$NAME),]$Annotation)
V(g)$color_border <- functions_colours[as.numeric(mod[[1]][match(V(g)$name, mod[[1]]$name),]$annotation)]

## Generate a weight column that includes the detected modules as an edge attribute
cl_for_layout <- com_mod
E(g)$weight <- edge.weights(cl_for_layout, g)

##### ATTENTION : THE CODE BELOW HAS TO BE RUN AT LEAST ONCE - UNCOMMENT ON FIRST RUN

# # create layout - run this line several times until you find a layout that you like
coord_graphopt <- layout_with_graphopt(g, mass = 20, charge = E(g)$weight)
igraph::plot.igraph(g,layout=coord_graphopt)
#After generating a good layout this can be saved for ploting all the other Mre11ID networks 
## and we store the layout coordinates to plot the other networks in the same way
coords_nodes <- coord_graphopt
row.names(coords_nodes) <- V(g)$name
coords_nodes
idxs <- match(V(g)$name, row.names(coords_nodes))
coord_graphopt <- coords_nodes[idxs,]
vsize <- igraph::degree(g)/vcount(g)

#Plot the networks with the module as the fill color, and the annotation as the color of the circle line
par(mar=c(0,0,0,0))
pdf(paste0('Mre11_protein-protein-network-whole.pdf'), width=30, height=30)
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('M = ', format(mod[[2]], digits=2)))
#
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.shape='circle', vertex.frame.color=V(g)$color_border, vertex.frame.width=4, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('whole network M = ', format(mod[[2]], digits = 2)), edge.arrow.size=.2)

legend('bottomleft', legend=levels(mod[[1]]$annotation), pt.lwd=4, cex=2, col=functions_colours, pch=21, pt.bg='white', title='function')

#recover the module results into a new table
cur_mod <- mod[[1]]
cur_mod$name <- as.character(cur_mod$name)
cur_mod$module <- cur_mod$module + 1
cur_mod$color <- custom_palette[cur_mod$module]
cur_mod$treatment <- 'nonIR5'
cur_mod$degree <- degree(g, mode='all')[match(cur_mod$name, names(degree(g, mode='all')))]

#Plot the  Participation coefficient versus the degree 
par(mar=c(4,4,4,4))
ggplot(cur_mod, aes(degree, participation, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('degree-part-nonIR5.pdf'), width = 12, height = 8)

# Plot the connectivity (within-module degree) versus the degree
ggplot(cur_mod, aes(degree, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('degree-conn-nonIR5.pdf'), width = 12, height = 8)

#Plot the connectivity versus the participation
ggplot(cur_mod, aes(participation, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('part-conn-nonIR5.pdf'), width = 12, height = 8)

#Print the results in a new table to merge them with the Node data frames
mod_overall <- cur_mod
write.csv(mod_overall,'./Mre11_nonIR5_mod.csv')

###IR5
metadata <- read.csv('./SensorID_Node_DataFrame.csv')#Load the Node metadata table
metadata<-subset(metadata,metadata$Sensor=="Mre11")#Isolate the node information corresponding to Mre11ID networks
df <- read.table('./Mre_IR5_Edges.tsv', header=T) #Load the edge list table
g <- graph_from_data_frame(df, directed = FALSE) #generate an igraph objetct
g <- simplify(g, remove.multiple=TRUE, remove.loops=TRUE, edge.attr.comb = "first") #remove duplicated edges, and self-loops
edge.attributes(g) #Check no NAs in the edge attributes
adj_mat_g <- as_adjacency_matrix(g, attr = 'weight', sparse = FALSE) #convert the igraph object into an adjacency matrix
mod <- netcarto(adj_mat_g) #generate the modules using the netcarto algorithm
#Now generate a function to run the netcarto algorithm 20 times, to improve the modularity results from the first run
for(i in 1:20){
  temp_mod <- netcarto(adj_mat_g)
  if(temp_mod[[2]] > mod[[2]])
    mod <- temp_mod
}

# integrate the module results from netcarto to the igraph network 
com_mod <- make_clusters(
  g,
  membership = mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1,
  algorithm = 'netcarto',
  merges = NULL,
  modularity = mod[[2]]
)
#Generate a color palette for the modules
V(g)$module <- mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1
V(g)$color <- custom_palette[V(g)$module]
#generate an color palette for the functional annotation of each node
mod[[1]]$annotation <- as.factor(metadata[match(mod[[1]]$name, metadata$NAME),]$Annotation)
V(g)$color_border <- functions_colours[as.numeric(mod[[1]][match(V(g)$name, mod[[1]]$name),]$annotation)]

## Generate a weight column that includes the detected modules as an edge attribute
cl_for_layout <- com_mod
E(g)$weight <- edge.weights(cl_for_layout, g)

##### ATTENTION : THE CODE BELOW HAS TO BE RUN AT LEAST ONCE - UNCOMMENT ON FIRST RUN

# # create layout - run this line several times until you find a layout that you like
coord_graphopt <- layout_with_graphopt(g, mass = 20, charge = E(g)$weight)
igraph::plot.igraph(g,layout=coord_graphopt)
#After generating a good layout this can be saved for ploting all the other Mre11ID networks 
## and we store the layout coordinates to plot the other networks in the same way
coords_nodes <- coord_graphopt
row.names(coords_nodes) <- V(g)$name
coords_nodes
idxs <- match(V(g)$name, row.names(coords_nodes))
coord_graphopt <- coords_nodes[idxs,]
vsize <- igraph::degree(g)/vcount(g)

#Plot the networks with the module as the fill color, and the annotation as the color of the circle line
par(mar=c(0,0,0,0))
pdf(paste0('Mre11_protein-protein-network-whole.pdf'), width=30, height=30)
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('M = ', format(mod[[2]], digits=2)))
#
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.shape='circle', vertex.frame.color=V(g)$color_border, vertex.frame.width=4, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('whole network M = ', format(mod[[2]], digits = 2)), edge.arrow.size=.2)

legend('bottomleft', legend=levels(mod[[1]]$annotation), pt.lwd=4, cex=2, col=functions_colours, pch=21, pt.bg='white', title='function')

#recover the module results into a new table
cur_mod <- mod[[1]]
cur_mod$name <- as.character(cur_mod$name)
cur_mod$module <- cur_mod$module + 1
cur_mod$color <- custom_palette[cur_mod$module]
cur_mod$treatment <- 'IR5'
cur_mod$degree <- degree(g, mode='all')[match(cur_mod$name, names(degree(g, mode='all')))]

#Plot the  Participation coefficient versus the degree 
par(mar=c(4,4,4,4))
ggplot(cur_mod, aes(degree, participation, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('degree-part-IR5.pdf'), width = 12, height = 8)

# Plot the connectivity (within-module degree) versus the degree
ggplot(cur_mod, aes(degree, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('degree-conn-IR5.pdf'), width = 12, height = 8)

#Plot the connectivity versus the participation
ggplot(cur_mod, aes(participation, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('part-conn-IR5.pdf'), width = 12, height = 8)

#Print the results in a new table to merge them with the Node data frames
mod_overall <- cur_mod
write.csv(mod_overall,'./Mre11_IR5_mod.csv')

###IR30
metadata <- read.csv('./SensorID_Node_DataFrame.csv')#Load the Node metadata table
metadata<-subset(metadata,metadata$Sensor=="Mre11")#Isolate the node information corresponding to Mre11ID networks
df <- read.table('./Mre_IR30_Edges.tsv', header=T) #Load the edge list table
g <- graph_from_data_frame(df, directed = FALSE) #generate an igraph objetct
g <- simplify(g, remove.multiple=TRUE, remove.loops=TRUE, edge.attr.comb = "first") #remove duplicated edges, and self-loops
edge.attributes(g) #Check no NAs in the edge attributes
adj_mat_g <- as_adjacency_matrix(g, attr = 'weight', sparse = FALSE) #convert the igraph object into an adjacency matrix
mod <- netcarto(adj_mat_g) #generate the modules using the netcarto algorithm
#Now generate a function to run the netcarto algorithm 20 times, to improve the modularity results from the first run
for(i in 1:20){
  temp_mod <- netcarto(adj_mat_g)
  if(temp_mod[[2]] > mod[[2]])
    mod <- temp_mod
}

# integrate the module results from netcarto to the igraph network 
com_mod <- make_clusters(
  g,
  membership = mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1,
  algorithm = 'netcarto',
  merges = NULL,
  modularity = mod[[2]]
)
#Generate a color palette for the modules
V(g)$module <- mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1
V(g)$color <- custom_palette[V(g)$module]
#generate an color palette for the functional annotation of each node
mod[[1]]$annotation <- as.factor(metadata[match(mod[[1]]$name, metadata$NAME),]$Annotation)
V(g)$color_border <- functions_colours[as.numeric(mod[[1]][match(V(g)$name, mod[[1]]$name),]$annotation)]

## Generate a weight column that includes the detected modules as an edge attribute
cl_for_layout <- com_mod
E(g)$weight <- edge.weights(cl_for_layout, g)

##### ATTENTION : THE CODE BELOW HAS TO BE RUN AT LEAST ONCE - UNCOMMENT ON FIRST RUN

# # create layout - run this line several times until you find a layout that you like
coord_graphopt <- layout_with_graphopt(g, mass = 20, charge = E(g)$weight)
igraph::plot.igraph(g,layout=coord_graphopt)
#After generating a good layout this can be saved for ploting all the other Mre11ID networks 
## and we store the layout coordinates to plot the other networks in the same way
coords_nodes <- coord_graphopt
row.names(coords_nodes) <- V(g)$name
coords_nodes
idxs <- match(V(g)$name, row.names(coords_nodes))
coord_graphopt <- coords_nodes[idxs,]
vsize <- igraph::degree(g)/vcount(g)

#Plot the networks with the module as the fill color, and the annotation as the color of the circle line
par(mar=c(0,0,0,0))
pdf(paste0('Mre11_protein-protein-network-whole.pdf'), width=30, height=30)
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('M = ', format(mod[[2]], digits=2)))
#
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.shape='circle', vertex.frame.color=V(g)$color_border, vertex.frame.width=4, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('whole network M = ', format(mod[[2]], digits = 2)), edge.arrow.size=.2)

legend('bottomleft', legend=levels(mod[[1]]$annotation), pt.lwd=4, cex=2, col=functions_colours, pch=21, pt.bg='white', title='function')

#recover the module results into a new table
cur_mod <- mod[[1]]
cur_mod$name <- as.character(cur_mod$name)
cur_mod$module <- cur_mod$module + 1
cur_mod$color <- custom_palette[cur_mod$module]
cur_mod$treatment <- 'IR30'
cur_mod$degree <- degree(g, mode='all')[match(cur_mod$name, names(degree(g, mode='all')))]

#Plot the  Participation coefficient versus the degree 
par(mar=c(4,4,4,4))
ggplot(cur_mod, aes(degree, participation, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('degree-part-IR30.pdf'), width = 12, height = 8)

# Plot the connectivity (within-module degree) versus the degree
ggplot(cur_mod, aes(degree, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('degree-conn-IR30.pdf'), width = 12, height = 8)

#Plot the connectivity versus the participation
ggplot(cur_mod, aes(participation, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('part-conn-IR30.pdf'), width = 12, height = 8)

#Print the results in a new table to merge them with the Node data frames
mod_overall <- cur_mod
write.csv(mod_overall,'./Mre11_IR30_mod.csv')

###IR2h
metadata <- read.csv('./SensorID_Node_DataFrame.csv')#Load the Node metadata table
metadata<-subset(metadata,metadata$Sensor=="Mre11")#Isolate the node information corresponding to Mre11ID networks
df <- read.table('./Mre_IR2h_Edges.tsv', header=T) #Load the edge list table
g <- graph_from_data_frame(df, directed = FALSE) #generate an igraph objetct
g <- simplify(g, remove.multiple=TRUE, remove.loops=TRUE, edge.attr.comb = "first") #remove duplicated edges, and self-loops
edge.attributes(g) #Check no NAs in the edge attributes
adj_mat_g <- as_adjacency_matrix(g, attr = 'weight', sparse = FALSE) #convert the igraph object into an adjacency matrix
mod <- netcarto(adj_mat_g) #generate the modules using the netcarto algorithm
#Now generate a function to run the netcarto algorithm 20 times, to improve the modularity results from the first run
for(i in 1:20){
  temp_mod <- netcarto(adj_mat_g)
  if(temp_mod[[2]] > mod[[2]])
    mod <- temp_mod
}

# integrate the module results from netcarto to the igraph network 
com_mod <- make_clusters(
  g,
  membership = mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1,
  algorithm = 'netcarto',
  merges = NULL,
  modularity = mod[[2]]
)
#Generate a color palette for the modules
V(g)$module <- mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1
V(g)$color <- custom_palette[V(g)$module]
#generate an color palette for the functional annotation of each node
mod[[1]]$annotation <- as.factor(metadata[match(mod[[1]]$name, metadata$NAME),]$Annotation)
V(g)$color_border <- functions_colours[as.numeric(mod[[1]][match(V(g)$name, mod[[1]]$name),]$annotation)]

## Generate a weight column that includes the detected modules as an edge attribute
cl_for_layout <- com_mod
E(g)$weight <- edge.weights(cl_for_layout, g)

##### ATTENTION : THE CODE BELOW HAS TO BE RUN AT LEAST ONCE - UNCOMMENT ON FIRST RUN

# # create layout - run this line several times until you find a layout that you like
coord_graphopt <- layout_with_graphopt(g, mass = 20, charge = E(g)$weight)
igraph::plot.igraph(g,layout=coord_graphopt)
#After generating a good layout this can be saved for ploting all the other Mre11ID networks 
## and we store the layout coordinates to plot the other networks in the same way
coords_nodes <- coord_graphopt
row.names(coords_nodes) <- V(g)$name
coords_nodes
idxs <- match(V(g)$name, row.names(coords_nodes))
coord_graphopt <- coords_nodes[idxs,]
vsize <- igraph::degree(g)/vcount(g)

#Plot the networks with the module as the fill color, and the annotation as the color of the circle line
par(mar=c(0,0,0,0))
pdf(paste0('Mre11_protein-protein-network-whole.pdf'), width=30, height=30)
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('M = ', format(mod[[2]], digits=2)))
#
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.shape='circle', vertex.frame.color=V(g)$color_border, vertex.frame.width=4, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('whole network M = ', format(mod[[2]], digits = 2)), edge.arrow.size=.2)

legend('bottomleft', legend=levels(mod[[1]]$annotation), pt.lwd=4, cex=2, col=functions_colours, pch=21, pt.bg='white', title='function')

#recover the module results into a new table
cur_mod <- mod[[1]]
cur_mod$name <- as.character(cur_mod$name)
cur_mod$module <- cur_mod$module + 1
cur_mod$color <- custom_palette[cur_mod$module]
cur_mod$treatment <- 'IR2h'
cur_mod$degree <- degree(g, mode='all')[match(cur_mod$name, names(degree(g, mode='all')))]

#Plot the  Participation coefficient versus the degree 
par(mar=c(4,4,4,4))
ggplot(cur_mod, aes(degree, participation, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('degree-part-IR2h.pdf'), width = 12, height = 8)

# Plot the connectivity (within-module degree) versus the degree
ggplot(cur_mod, aes(degree, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('degree-conn-IR2h.pdf'), width = 12, height = 8)

#Plot the connectivity versus the participation
ggplot(cur_mod, aes(participation, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('part-conn-IR2h.pdf'), width = 12, height = 8)

#Print the results in a new table to merge them with the Node data frames
mod_overall <- cur_mod
write.csv(mod_overall,'./Mre11_IR2h_mod.csv')

###IR8h
metadata <- read.csv('./SensorID_Node_DataFrame.csv')#Load the Node metadata table
metadata<-subset(metadata,metadata$Sensor=="Mre11")#Isolate the node information corresponding to Mre11ID networks
df <- read.table('./Mre_IR8h_Edges.tsv', header=T) #Load the edge list table
g <- graph_from_data_frame(df, directed = FALSE) #generate an igraph objetct
g <- simplify(g, remove.multiple=TRUE, remove.loops=TRUE, edge.attr.comb = "first") #remove duplicated edges, and self-loops
edge.attributes(g) #Check no NAs in the edge attributes
adj_mat_g <- as_adjacency_matrix(g, attr = 'weight', sparse = FALSE) #convert the igraph object into an adjacency matrix
mod <- netcarto(adj_mat_g) #generate the modules using the netcarto algorithm
#Now generate a function to run the netcarto algorithm 20 times, to improve the modularity results from the first run
for(i in 1:20){
  temp_mod <- netcarto(adj_mat_g)
  if(temp_mod[[2]] > mod[[2]])
    mod <- temp_mod
}

# integrate the module results from netcarto to the igraph network 
com_mod <- make_clusters(
  g,
  membership = mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1,
  algorithm = 'netcarto',
  merges = NULL,
  modularity = mod[[2]]
)
#Generate a color palette for the modules
V(g)$module <- mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1
V(g)$color <- custom_palette[V(g)$module]
#generate an color palette for the functional annotation of each node
mod[[1]]$annotation <- as.factor(metadata[match(mod[[1]]$name, metadata$NAME),]$Annotation)
V(g)$color_border <- functions_colours[as.numeric(mod[[1]][match(V(g)$name, mod[[1]]$name),]$annotation)]

## Generate a weight column that includes the detected modules as an edge attribute
cl_for_layout <- com_mod
E(g)$weight <- edge.weights(cl_for_layout, g)

##### ATTENTION : THE CODE BELOW HAS TO BE RUN AT LEAST ONCE - UNCOMMENT ON FIRST RUN

# # create layout - run this line several times until you find a layout that you like
coord_graphopt <- layout_with_graphopt(g, mass = 20, charge = E(g)$weight)
igraph::plot.igraph(g,layout=coord_graphopt)
#After generating a good layout this can be saved for ploting all the other Mre11ID networks 
## and we store the layout coordinates to plot the other networks in the same way
coords_nodes <- coord_graphopt
row.names(coords_nodes) <- V(g)$name
coords_nodes
idxs <- match(V(g)$name, row.names(coords_nodes))
coord_graphopt <- coords_nodes[idxs,]
vsize <- igraph::degree(g)/vcount(g)

#Plot the networks with the module as the fill color, and the annotation as the color of the circle line
par(mar=c(0,0,0,0))
pdf(paste0('Mre11_protein-protein-network-whole.pdf'), width=30, height=30)
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('M = ', format(mod[[2]], digits=2)))
#
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.shape='circle', vertex.frame.color=V(g)$color_border, vertex.frame.width=4, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('whole network M = ', format(mod[[2]], digits = 2)), edge.arrow.size=.2)

legend('bottomleft', legend=levels(mod[[1]]$annotation), pt.lwd=4, cex=2, col=functions_colours, pch=21, pt.bg='white', title='function')

#recover the module results into a new table
cur_mod <- mod[[1]]
cur_mod$name <- as.character(cur_mod$name)
cur_mod$module <- cur_mod$module + 1
cur_mod$color <- custom_palette[cur_mod$module]
cur_mod$treatment <- 'IR8h'
cur_mod$degree <- degree(g, mode='all')[match(cur_mod$name, names(degree(g, mode='all')))]

#Plot the  Participation coefficient versus the degree 
par(mar=c(4,4,4,4))
ggplot(cur_mod, aes(degree, participation, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('degree-part-IR8h.pdf'), width = 12, height = 8)

# Plot the connectivity (within-module degree) versus the degree
ggplot(cur_mod, aes(degree, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('degree-conn-IR8h.pdf'), width = 12, height = 8)

#Plot the connectivity versus the participation
ggplot(cur_mod, aes(participation, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('part-conn-IR8h.pdf'), width = 12, height = 8)

#Print the results in a new table to merge them with the Node data frames
mod_overall <- cur_mod
write.csv(mod_overall,'./Mre11_IR8h_mod.csv')

###IR24h
metadata <- read.csv('./SensorID_Node_DataFrame.csv')#Load the Node metadata table
metadata<-subset(metadata,metadata$Sensor=="Mre11")#Isolate the node information corresponding to Mre11ID networks
df <- read.table('./Mre_IR24h_Edges.tsv', header=T) #Load the edge list table
g <- graph_from_data_frame(df, directed = FALSE) #generate an igraph objetct
g <- simplify(g, remove.multiple=TRUE, remove.loops=TRUE, edge.attr.comb = "first") #remove duplicated edges, and self-loops
edge.attributes(g) #Check no NAs in the edge attributes
adj_mat_g <- as_adjacency_matrix(g, attr = 'weight', sparse = FALSE) #convert the igraph object into an adjacency matrix
mod <- netcarto(adj_mat_g) #generate the modules using the netcarto algorithm
#Now generate a function to run the netcarto algorithm 20 times, to improve the modularity results from the first run
for(i in 1:20){
  temp_mod <- netcarto(adj_mat_g)
  if(temp_mod[[2]] > mod[[2]])
    mod <- temp_mod
}

# integrate the module results from netcarto to the igraph network 
com_mod <- make_clusters(
  g,
  membership = mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1,
  algorithm = 'netcarto',
  merges = NULL,
  modularity = mod[[2]]
)
#Generate a color palette for the modules
V(g)$module <- mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1
V(g)$color <- custom_palette[V(g)$module]
#generate an color palette for the functional annotation of each node
mod[[1]]$annotation <- as.factor(metadata[match(mod[[1]]$name, metadata$NAME),]$Annotation)
V(g)$color_border <- functions_colours[as.numeric(mod[[1]][match(V(g)$name, mod[[1]]$name),]$annotation)]

## Generate a weight column that includes the detected modules as an edge attribute
cl_for_layout <- com_mod
E(g)$weight <- edge.weights(cl_for_layout, g)

##### ATTENTION : THE CODE BELOW HAS TO BE RUN AT LEAST ONCE - UNCOMMENT ON FIRST RUN

# # create layout - run this line several times until you find a layout that you like
coord_graphopt <- layout_with_graphopt(g, mass = 20, charge = E(g)$weight)
igraph::plot.igraph(g,layout=coord_graphopt)
#After generating a good layout this can be saved for ploting all the other Mre11ID networks 
## and we store the layout coordinates to plot the other networks in the same way
coords_nodes <- coord_graphopt
row.names(coords_nodes) <- V(g)$name
coords_nodes
idxs <- match(V(g)$name, row.names(coords_nodes))
coord_graphopt <- coords_nodes[idxs,]
vsize <- igraph::degree(g)/vcount(g)

#Plot the networks with the module as the fill color, and the annotation as the color of the circle line
par(mar=c(0,0,0,0))
pdf(paste0('Mre11_protein-protein-network-whole.pdf'), width=30, height=30)
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('M = ', format(mod[[2]], digits=2)))
#
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.shape='circle', vertex.frame.color=V(g)$color_border, vertex.frame.width=4, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('whole network M = ', format(mod[[2]], digits = 2)), edge.arrow.size=.2)

legend('bottomleft', legend=levels(mod[[1]]$annotation), pt.lwd=4, cex=2, col=functions_colours, pch=21, pt.bg='white', title='function')

#recover the module results into a new table
cur_mod <- mod[[1]]
cur_mod$name <- as.character(cur_mod$name)
cur_mod$module <- cur_mod$module + 1
cur_mod$color <- custom_palette[cur_mod$module]
cur_mod$treatment <- 'IR24h'
cur_mod$degree <- degree(g, mode='all')[match(cur_mod$name, names(degree(g, mode='all')))]

#Plot the  Participation coefficient versus the degree 
par(mar=c(4,4,4,4))
ggplot(cur_mod, aes(degree, participation, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('degree-part-IR24h.pdf'), width = 12, height = 8)

# Plot the connectivity (within-module degree) versus the degree
ggplot(cur_mod, aes(degree, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('degree-conn-IR24h.pdf'), width = 12, height = 8)

#Plot the connectivity versus the participation
ggplot(cur_mod, aes(participation, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('part-conn-IR24h.pdf'), width = 12, height = 8)

#Print the results in a new table to merge them with the Node data frames
mod_overall <- cur_mod
write.csv(mod_overall,'./Mre11_IR24h_mod.csv')

###nonIR24
metadata <- read.csv('./SensorID_Node_DataFrame.csv')#Load the Node metadata table
metadata<-subset(metadata,metadata$Sensor=="Mre11")#Isolate the node information corresponding to Mre11ID networks
df <- read.table('./Mre_nonIR24_Edges.tsv', header=T) #Load the edge list table
g <- graph_from_data_frame(df, directed = FALSE) #generate an igraph objetct
g <- simplify(g, remove.multiple=TRUE, remove.loops=TRUE, edge.attr.comb = "first") #remove duplicated edges, and self-loops
edge.attributes(g) #Check no NAs in the edge attributes
adj_mat_g <- as_adjacency_matrix(g, attr = 'weight', sparse = FALSE) #convert the igraph object into an adjacency matrix
mod <- netcarto(adj_mat_g) #generate the modules using the netcarto algorithm
#Now generate a function to run the netcarto algorithm 20 times, to improve the modularity results from the first run
for(i in 1:20){
  temp_mod <- netcarto(adj_mat_g)
  if(temp_mod[[2]] > mod[[2]])
    mod <- temp_mod
}

# integrate the module results from netcarto to the igraph network 
com_mod <- make_clusters(
  g,
  membership = mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1,
  algorithm = 'netcarto',
  merges = NULL,
  modularity = mod[[2]]
)
#Generate a color palette for the modules
V(g)$module <- mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1
V(g)$color <- custom_palette[V(g)$module]
#generate an color palette for the functional annotation of each node
mod[[1]]$annotation <- as.factor(metadata[match(mod[[1]]$name, metadata$NAME),]$Annotation)
V(g)$color_border <- functions_colours[as.numeric(mod[[1]][match(V(g)$name, mod[[1]]$name),]$annotation)]

## Generate a weight column that includes the detected modules as an edge attribute
cl_for_layout <- com_mod
E(g)$weight <- edge.weights(cl_for_layout, g)

##### ATTENTION : THE CODE BELOW HAS TO BE RUN AT LEAST ONCE - UNCOMMENT ON FIRST RUN

# # create layout - run this line several times until you find a layout that you like
coord_graphopt <- layout_with_graphopt(g, mass = 20, charge = E(g)$weight)
igraph::plot.igraph(g,layout=coord_graphopt)
#After generating a good layout this can be saved for ploting all the other Mre11ID networks 
## and we store the layout coordinates to plot the other networks in the same way
coords_nodes <- coord_graphopt
row.names(coords_nodes) <- V(g)$name
coords_nodes
idxs <- match(V(g)$name, row.names(coords_nodes))
coord_graphopt <- coords_nodes[idxs,]
vsize <- igraph::degree(g)/vcount(g)

#Plot the networks with the module as the fill color, and the annotation as the color of the circle line
par(mar=c(0,0,0,0))
pdf(paste0('Mre11_protein-protein-network-whole.pdf'), width=30, height=30)
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('M = ', format(mod[[2]], digits=2)))
#
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.shape='circle', vertex.frame.color=V(g)$color_border, vertex.frame.width=4, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('whole network M = ', format(mod[[2]], digits = 2)), edge.arrow.size=.2)

legend('bottomleft', legend=levels(mod[[1]]$annotation), pt.lwd=4, cex=2, col=functions_colours, pch=21, pt.bg='white', title='function')

#recover the module results into a new table
cur_mod <- mod[[1]]
cur_mod$name <- as.character(cur_mod$name)
cur_mod$module <- cur_mod$module + 1
cur_mod$color <- custom_palette[cur_mod$module]
cur_mod$treatment <- 'nonIR24'
cur_mod$degree <- degree(g, mode='all')[match(cur_mod$name, names(degree(g, mode='all')))]

#Plot the  Participation coefficient versus the degree 
par(mar=c(4,4,4,4))
ggplot(cur_mod, aes(degree, participation, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('degree-part-nonIR24.pdf'), width = 12, height = 8)

# Plot the connectivity (within-module degree) versus the degree
ggplot(cur_mod, aes(degree, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('degree-conn-nonIR24.pdf'), width = 12, height = 8)

#Plot the connectivity versus the participation
ggplot(cur_mod, aes(participation, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')

ggsave(paste0('part-conn-nonIR24.pdf'), width = 12, height = 8)

#Print the results in a new table to merge them with the Node data frames
mod_overall <- cur_mod
write.csv(mod_overall,'./Mre11_nonIR24_mod.csv')




