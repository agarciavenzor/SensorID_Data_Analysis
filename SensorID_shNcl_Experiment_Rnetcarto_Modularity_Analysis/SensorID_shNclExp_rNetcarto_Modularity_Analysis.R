
### This R script contains the code required to perform the rnetcarto modularity analysis using the SensorID shNcl and shScr networks edge lists
## the results from the analysis were later manually verified, and plotted, and merged with previous results.

#Install the required packages
install.packages(c("rnetcarto","igraph","ggplot2"))
#Load the required packages
require(igraph)
require(rnetcarto)
require(ggplot2)

### rnetcarto modularity Analysis for the SensorID networks after silecing of Nucleolin (shNcl)

### set working directory
setwd("../SensorID_shNcl_Experiment_Rnetcarto_Modularity_Analysis")

### Generate a color palette with contrast for the resulting modules
custom_palette <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')
functions_colours <- rainbow(10)
?rainbow

###### ONLY FOR PLOTTING THE NETWORKS: For layout the networks using the detected modules to highlight the edges ######
edge.weights <- function(community, network, weight.within = .001, weight.between = .01) {
  bridges <- crossing(communities = community, graph = network)
  weights <- ifelse(test = bridges, yes = weight.between, no = weight.within)
  return(weights) 
}

### Sirt6ID shNcl experiment network modularity analysis ###
### Sirt6ID full nuclear network (it contain all the interactors found in the shScr and shNcl conditions, it serves to use a common layout for all the Sirt6ID networks)
## Load the required files for the analysis
metadata <- read.csv('./Sirt6ID_shNclExp_DataFrame.csv')#Open the Sirt6ID node metadata tables, containing the annotations for each protein
df <- read.delim('./Sirt6ID_shNclExp_FullNetwork_EdgeList.tsv', sep="\t",header=T) #Load the edge list of the Sirt6ID full Network
g <- graph_from_data_frame(df, directed = FALSE) # convert the edge list into an igraph network object
g <- simplify(g, remove.multiple=TRUE, remove.loops=TRUE, edge.attr.comb = "first") # Remove the duplicated edges and self-loops
edge.attributes(g) # Confirm that there is no NAs in the edge weights
adj_mat_g <- as_adjacency_matrix(g, attr = 'weight1', sparse = FALSE) # Transform the networks into an adjacency matrix
mod <- netcarto(adj_mat_g) # Run the netcarto modularity analysis algorithm 

# Now we run the netcarto modularity analysis 20 times to optimize the separation of modules
for(i in 1:20){
  temp_mod <- netcarto(adj_mat_g)
  if(temp_mod[[2]] > mod[[2]])
    mod <- temp_mod
}

# Translate the netcarto generated modules in a format compatible with igraph network objects
com_mod <- make_clusters(
  g,
  membership = mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1,
  algorithm = 'netcarto',
  merges = NULL,
  modularity = mod[[2]]
)

# Add the module and module colors as vertex attributes to the igraph networks
V(g)$module <- mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1
V(g)$color <- custom_palette[V(g)$module]
# add the functional annotations to the netcarto generated modules
mod[[1]]$annotation <- as.factor(metadata[match(mod[[1]]$name, metadata$NAME),]$Annotation)
V(g)$color_border <- functions_colours[as.numeric(mod[[1]][match(V(g)$name, mod[[1]]$name),]$annotation)]

## Generate the layout to map the modules into the network
cl_for_layout <- com_mod
E(g)$weight <- edge.weights(cl_for_layout, g)

##### ATTENTION : THE CODE BELOW HAS TO BE RUN AT LEAST ONCE - UNCOMMENT ON FIRST RUN
# # create layout - run this line several times until you find a layout that you like
coord_graphopt <- layout_with_graphopt(g, mass = 20, charge = E(g)$weight)
igraph::plot.igraph(g,layout=coord_graphopt)

# # and we store the layout coordinates to plot the other networks in the same way
coords_nodes <- coord_graphopt
row.names(coords_nodes) <- V(g)$name
coords_nodes
# match the nodes of the network with the coordinates in the layout
idxs <- match(V(g)$name, row.names(coords_nodes))
coord_graphopt <- coords_nodes[idxs,]
vsize <- igraph::degree(g)/vcount(g)

#Plot the network sith the module as fill color and the annotation as the perimeter color
par(mar=c(0,0,0,0))
pdf(paste0('Sirt6_protein-protein-network-whole.pdf'), width=30, height=30)
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('M = ', format(mod[[2]], digits=2)))
#
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.shape='circle', vertex.frame.color=V(g)$color_border, vertex.frame.width=4, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('whole network M = ', format(mod[[2]], digits = 2)), edge.arrow.size=.2)
legend('bottomleft', legend=levels(mod[[1]]$annotation), pt.lwd=4, cex=2, col=functions_colours, pch=21, pt.bg='white', title='function')

#Create a data table containing the module number, the treatment and the degree per node
cur_mod <- mod[[1]]
cur_mod$name <- as.character(cur_mod$name)
cur_mod$module <- cur_mod$module + 1

cur_mod$color <- custom_palette[cur_mod$module]

cur_mod$treatment <- 'whole'
cur_mod$degree <- degree(g, mode='all')[match(cur_mod$name, names(degree(g, mode='all')))]

#Plot the Participation (inter-module degree coefficient) versus the node degree (number of connections)
par(mar=c(4,4,4,4))
ggplot(cur_mod, aes(degree, participation, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')
#Save the plot
ggsave(paste0('degree-part-whole.pdf'), width = 12, height = 8)

#Plot the connectivity (within-module degree coefficient) versus the degree (number of connectios)
ggplot(cur_mod, aes(degree, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')
#Save the plot
ggsave(paste0('degree-conn-whole.pdf'), width = 12, height = 8)
# Plot the connectivity versus the participation (This plot allows an easy visualization of the connector and central nodes)
ggplot(cur_mod, aes(participation, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')
# Save the Plot
ggsave(paste0('part-conn-whole.pdf'), width = 12, height = 8)

mod_overall <- cur_mod

#write the results tables
write.csv(mod_overall,'./Sirt6ID_shNclExp_FullNetwork_overall_mod.csv')

### Sirt6ID shNcl experiment network modularity analysis ###
### Sirt6ID shScramble network 
## Load the required files for the analysis
metadata <- read.csv('./Sirt6ID_shNclExp_DataFrame.csv')#Open the Sirt6ID node metadata tables, containing the annotations for each protein
metadata<- subset(metadata, !metadata$shScr.Presence==0) #Subset the Sirt6ID datafram to keep only the results of the Sirt6ID shScrable condition
df <- read.delim('./Sirt6ID_shScrlvsCneg_Network_EdgeList.tsv', sep="\t",header=T) #Load the edge list of the Sirt6ID shScr Network
g <- graph_from_data_frame(df, directed = FALSE) # convert the edge list into an igraph network object
g <- simplify(g, remove.multiple=TRUE, remove.loops=TRUE, edge.attr.comb = "first") # Remove the duplicated edges and self-loops
edge.attributes(g) # Confirm that there is no NAs in the edge weights
adj_mat_g <- as_adjacency_matrix(g, attr = 'weight1', sparse = FALSE) # Transform the networks into an adjacency matrix
mod <- netcarto(adj_mat_g) # Run the netcarto modularity analysis algorithm 

# Now we run the netcarto modularity analysis 20 times to optimize the separation of modules
for(i in 1:20){
  temp_mod <- netcarto(adj_mat_g)
  if(temp_mod[[2]] > mod[[2]])
    mod <- temp_mod
}

# Translate the netcarto generated modules in a format compatible with igraph network objects
com_mod <- make_clusters(
  g,
  membership = mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1,
  algorithm = 'netcarto',
  merges = NULL,
  modularity = mod[[2]]
)

# Add the module and module colors as vertex attributes to the igraph networks
V(g)$module <- mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1
V(g)$color <- custom_palette[V(g)$module]
# add the functional annotations to the netcarto generated modules
mod[[1]]$annotation <- as.factor(metadata[match(mod[[1]]$name, metadata$NAME),]$Annotation)
V(g)$color_border <- functions_colours[as.numeric(mod[[1]][match(V(g)$name, mod[[1]]$name),]$annotation)]

## Generate the layout to map the modules into the network
cl_for_layout <- com_mod
E(g)$weight <- edge.weights(cl_for_layout, g)

##### ATTENTION : THE CODE BELOW HAS TO BE RUN AT LEAST ONCE - UNCOMMENT ON FIRST RUN
# # create layout - run this line several times until you find a layout that you like
coord_graphopt <- layout_with_graphopt(g, mass = 20, charge = E(g)$weight)
igraph::plot.igraph(g,layout=coord_graphopt)

# # and we store the layout coordinates to plot the other networks in the same way
coords_nodes <- coord_graphopt
row.names(coords_nodes) <- V(g)$name
coords_nodes
# match the nodes of the network with the coordinates in the layout
idxs <- match(V(g)$name, row.names(coords_nodes))
coord_graphopt <- coords_nodes[idxs,]
vsize <- igraph::degree(g)/vcount(g)

#Plot the network sith the module as fill color and the annotation as the perimeter color
par(mar=c(0,0,0,0))
pdf(paste0('Sirt6_protein-protein-network-shScr.pdf'), width=30, height=30)
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('M = ', format(mod[[2]], digits=2)))
#
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.shape='circle', vertex.frame.color=V(g)$color_border, vertex.frame.width=4, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('whole network M = ', format(mod[[2]], digits = 2)), edge.arrow.size=.2)
legend('bottomleft', legend=levels(mod[[1]]$annotation), pt.lwd=4, cex=2, col=functions_colours, pch=21, pt.bg='white', title='function')

#Create a data table containing the module number, the treatment and the degree per node
cur_mod <- mod[[1]]
cur_mod$name <- as.character(cur_mod$name)
cur_mod$module <- cur_mod$module + 1

cur_mod$color <- custom_palette[cur_mod$module]

cur_mod$treatment <- 'shScr'
cur_mod$degree <- degree(g, mode='all')[match(cur_mod$name, names(degree(g, mode='all')))]

#Plot the Participation (inter-module degree coefficient) versus the node degree (number of connections)
par(mar=c(4,4,4,4))
ggplot(cur_mod, aes(degree, participation, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')
#Save the plot
ggsave(paste0('degree-part-shScr.pdf'), width = 12, height = 8)

#Plot the connectivity (within-module degree coefficient) versus the degree (number of connectios)
ggplot(cur_mod, aes(degree, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')
#Save the plot
ggsave(paste0('degree-conn-shScr.pdf'), width = 12, height = 8)
# Plot the connectivity versus the participation (This plot allows an easy visualization of the connector and central nodes)
ggplot(cur_mod, aes(participation, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')
# Save the Plot
ggsave(paste0('part-conn-shScr.pdf'), width = 12, height = 8)

mod_overall <- cur_mod

#write the results tables
write.csv(mod_overall,'./Sirt6ID_shNclExp_shScr_Net_mod.csv')



### Sirt6ID shNcl experiment network modularity analysis ###
### Sirt6ID shNcl network 
## Load the required files for the analysis
metadata <- read.csv('./Sirt6ID_shNclExp_DataFrame.csv')#Open the Sirt6ID node metadata tables, containing the annotations for each protein
metadata<- subset(metadata, !metadata$shNcl.Presence==0) #Subset the Sirt6ID datafram to keep only the results of the Sirt6ID shNcl condition
df <- read.delim('./Sirt6ID_shNclvsCneg_Network_EdgeList.tsv', sep="\t",header=T) #Load the edge list of the Sirt6ID shNcl Network
g <- graph_from_data_frame(df, directed = FALSE) # convert the edge list into an igraph network object
g <- simplify(g, remove.multiple=TRUE, remove.loops=TRUE, edge.attr.comb = "first") # Remove the duplicated edges and self-loops
edge.attributes(g) # Confirm that there is no NAs in the edge weights
adj_mat_g <- as_adjacency_matrix(g, attr = 'weight1', sparse = FALSE) # Transform the networks into an adjacency matrix
mod <- netcarto(adj_mat_g) # Run the netcarto modularity analysis algorithm 

# Now we run the netcarto modularity analysis 20 times to optimize the separation of modules
for(i in 1:20){
  temp_mod <- netcarto(adj_mat_g)
  if(temp_mod[[2]] > mod[[2]])
    mod <- temp_mod
}

# Translate the netcarto generated modules in a format compatible with igraph network objects
com_mod <- make_clusters(
  g,
  membership = mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1,
  algorithm = 'netcarto',
  merges = NULL,
  modularity = mod[[2]]
)

# Add the module and module colors as vertex attributes to the igraph networks
V(g)$module <- mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1
V(g)$color <- custom_palette[V(g)$module]
# add the functional annotations to the netcarto generated modules
mod[[1]]$annotation <- as.factor(metadata[match(mod[[1]]$name, metadata$NAME),]$Annotation)
V(g)$color_border <- functions_colours[as.numeric(mod[[1]][match(V(g)$name, mod[[1]]$name),]$annotation)]

## Generate the layout to map the modules into the network
cl_for_layout <- com_mod
E(g)$weight <- edge.weights(cl_for_layout, g)

##### ATTENTION : THE CODE BELOW HAS TO BE RUN AT LEAST ONCE - UNCOMMENT ON FIRST RUN
# # create layout - run this line several times until you find a layout that you like
coord_graphopt <- layout_with_graphopt(g, mass = 20, charge = E(g)$weight)
igraph::plot.igraph(g,layout=coord_graphopt)

# # and we store the layout coordinates to plot the other networks in the same way
coords_nodes <- coord_graphopt
row.names(coords_nodes) <- V(g)$name
coords_nodes
# match the nodes of the network with the coordinates in the layout
idxs <- match(V(g)$name, row.names(coords_nodes))
coord_graphopt <- coords_nodes[idxs,]
vsize <- igraph::degree(g)/vcount(g)

#Plot the network sith the module as fill color and the annotation as the perimeter color
par(mar=c(0,0,0,0))
pdf(paste0('Sirt6_protein-protein-network-shNcl.pdf'), width=30, height=30)
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('M = ', format(mod[[2]], digits=2)))
#
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.shape='circle', vertex.frame.color=V(g)$color_border, vertex.frame.width=4, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('whole network M = ', format(mod[[2]], digits = 2)), edge.arrow.size=.2)
legend('bottomleft', legend=levels(mod[[1]]$annotation), pt.lwd=4, cex=2, col=functions_colours, pch=21, pt.bg='white', title='function')

#Create a data table containing the module number, the treatment and the degree per node
cur_mod <- mod[[1]]
cur_mod$name <- as.character(cur_mod$name)
cur_mod$module <- cur_mod$module + 1

cur_mod$color <- custom_palette[cur_mod$module]

cur_mod$treatment <- 'shNcl'
cur_mod$degree <- degree(g, mode='all')[match(cur_mod$name, names(degree(g, mode='all')))]

#Plot the Participation (inter-module degree coefficient) versus the node degree (number of connections)
par(mar=c(4,4,4,4))
ggplot(cur_mod, aes(degree, participation, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')
#Save the plot
ggsave(paste0('degree-part-shNcl.pdf'), width = 12, height = 8)

#Plot the connectivity (within-module degree coefficient) versus the degree (number of connectios)
ggplot(cur_mod, aes(degree, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')
#Save the plot
ggsave(paste0('degree-conn-shNcl.pdf'), width = 12, height = 8)
# Plot the connectivity versus the participation (This plot allows an easy visualization of the connector and central nodes)
ggplot(cur_mod, aes(participation, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')
# Save the Plot
ggsave(paste0('part-conn-shNcl.pdf'), width = 12, height = 8)

mod_overall <- cur_mod

#write the results tables
write.csv(mod_overall,'./Sirt6ID_shNclExp_shNcl_Net_mod.csv')

####################################
### Ku80ID shNcl experiment network modularity analysis ###
### Ku80ID full nuclear network (it contain all the interactors found in the shScr and shNcl conditions, it serves to use a common layout for all the Ku80ID networks)
## Load the required files for the analysis
metadata <- read.csv('./Ku80ID_shNclExp_DataFrame.csv')#Open the Ku80ID node metadata tables, containing the annotations for each protein
df <- read.delim('./Ku80ID_shNclExp_FullNetwork_EdgeList.tsv', sep="\t",header=T) #Load the edge list of the Ku80ID full Network
g <- graph_from_data_frame(df, directed = FALSE) # convert the edge list into an igraph network object
g <- simplify(g, remove.multiple=TRUE, remove.loops=TRUE, edge.attr.comb = "first") # Remove the duplicated edges and self-loops
edge.attributes(g) # Confirm that there is no NAs in the edge weights
adj_mat_g <- as_adjacency_matrix(g, attr = 'weight1', sparse = FALSE) # Transform the networks into an adjacency matrix
mod <- netcarto(adj_mat_g) # Run the netcarto modularity analysis algorithm 

# Now we run the netcarto modularity analysis 20 times to optimize the separation of modules
for(i in 1:20){
  temp_mod <- netcarto(adj_mat_g)
  if(temp_mod[[2]] > mod[[2]])
    mod <- temp_mod
}

# Translate the netcarto generated modules in a format compatible with igraph network objects
com_mod <- make_clusters(
  g,
  membership = mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1,
  algorithm = 'netcarto',
  merges = NULL,
  modularity = mod[[2]]
)

# Add the module and module colors as vertex attributes to the igraph networks
V(g)$module <- mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1
V(g)$color <- custom_palette[V(g)$module]
# add the functional annotations to the netcarto generated modules
mod[[1]]$annotation <- as.factor(metadata[match(mod[[1]]$name, metadata$NAME),]$Annotation)
V(g)$color_border <- functions_colours[as.numeric(mod[[1]][match(V(g)$name, mod[[1]]$name),]$annotation)]

## Generate the layout to map the modules into the network
cl_for_layout <- com_mod
E(g)$weight <- edge.weights(cl_for_layout, g)

##### ATTENTION : THE CODE BELOW HAS TO BE RUN AT LEAST ONCE - UNCOMMENT ON FIRST RUN
# # create layout - run this line several times until you find a layout that you like
coord_graphopt <- layout_with_graphopt(g, mass = 20, charge = E(g)$weight)
igraph::plot.igraph(g,layout=coord_graphopt)

# # and we store the layout coordinates to plot the other networks in the same way
coords_nodes <- coord_graphopt
row.names(coords_nodes) <- V(g)$name
coords_nodes
# match the nodes of the network with the coordinates in the layout
idxs <- match(V(g)$name, row.names(coords_nodes))
coord_graphopt <- coords_nodes[idxs,]
vsize <- igraph::degree(g)/vcount(g)

#Plot the network sith the module as fill color and the annotation as the perimeter color
par(mar=c(0,0,0,0))
pdf(paste0('Ku80_protein-protein-network-whole.pdf'), width=30, height=30)
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('M = ', format(mod[[2]], digits=2)))
#
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.shape='circle', vertex.frame.color=V(g)$color_border, vertex.frame.width=4, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('whole network M = ', format(mod[[2]], digits = 2)), edge.arrow.size=.2)
legend('bottomleft', legend=levels(mod[[1]]$annotation), pt.lwd=4, cex=2, col=functions_colours, pch=21, pt.bg='white', title='function')

#Create a data table containing the module number, the treatment and the degree per node
cur_mod <- mod[[1]]
cur_mod$name <- as.character(cur_mod$name)
cur_mod$module <- cur_mod$module + 1

cur_mod$color <- custom_palette[cur_mod$module]

cur_mod$treatment <- 'whole'
cur_mod$degree <- degree(g, mode='all')[match(cur_mod$name, names(degree(g, mode='all')))]

#Plot the Participation (inter-module degree coefficient) versus the node degree (number of connections)
par(mar=c(4,4,4,4))
ggplot(cur_mod, aes(degree, participation, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')
#Save the plot
ggsave(paste0('degree-part-whole.pdf'), width = 12, height = 8)

#Plot the connectivity (within-module degree coefficient) versus the degree (number of connectios)
ggplot(cur_mod, aes(degree, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')
#Save the plot
ggsave(paste0('degree-conn-whole.pdf'), width = 12, height = 8)
# Plot the connectivity versus the participation (This plot allows an easy visualization of the connector and central nodes)
ggplot(cur_mod, aes(participation, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')
# Save the Plot
ggsave(paste0('part-conn-whole.pdf'), width = 12, height = 8)

mod_overall <- cur_mod

#write the results tables
write.csv(mod_overall,'./Ku80ID_shNclExp_FullNetwork_overall_mod.csv')

### Ku80ID shNcl experiment network modularity analysis ###
### Ku80ID shScramble network 
## Load the required files for the analysis
metadata <- read.csv('./Ku80ID_shNclExp_DataFrame.csv')#Open the Ku80ID node metadata tables, containing the annotations for each protein
metadata<- subset(metadata, !metadata$shScr.Presence==0) #Subset the Ku80ID datafram to keep only the results of the Ku80ID shScrable condition
df <- read.delim('./Ku80ID_shScrlvsCneg_Network_EdgeList.tsv', sep="\t",header=T) #Load the edge list of the Ku80ID shScr Network
g <- graph_from_data_frame(df, directed = FALSE) # convert the edge list into an igraph network object
g <- simplify(g, remove.multiple=TRUE, remove.loops=TRUE, edge.attr.comb = "first") # Remove the duplicated edges and self-loops
edge.attributes(g) # Confirm that there is no NAs in the edge weights
adj_mat_g <- as_adjacency_matrix(g, attr = 'weight1', sparse = FALSE) # Transform the networks into an adjacency matrix
mod <- netcarto(adj_mat_g) # Run the netcarto modularity analysis algorithm 

# Now we run the netcarto modularity analysis 20 times to optimize the separation of modules
for(i in 1:20){
  temp_mod <- netcarto(adj_mat_g)
  if(temp_mod[[2]] > mod[[2]])
    mod <- temp_mod
}

# Translate the netcarto generated modules in a format compatible with igraph network objects
com_mod <- make_clusters(
  g,
  membership = mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1,
  algorithm = 'netcarto',
  merges = NULL,
  modularity = mod[[2]]
)

# Add the module and module colors as vertex attributes to the igraph networks
V(g)$module <- mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1
V(g)$color <- custom_palette[V(g)$module]
# add the functional annotations to the netcarto generated modules
mod[[1]]$annotation <- as.factor(metadata[match(mod[[1]]$name, metadata$NAME),]$Annotation)
V(g)$color_border <- functions_colours[as.numeric(mod[[1]][match(V(g)$name, mod[[1]]$name),]$annotation)]

## Generate the layout to map the modules into the network
cl_for_layout <- com_mod
E(g)$weight <- edge.weights(cl_for_layout, g)

##### ATTENTION : THE CODE BELOW HAS TO BE RUN AT LEAST ONCE - UNCOMMENT ON FIRST RUN
# # create layout - run this line several times until you find a layout that you like
coord_graphopt <- layout_with_graphopt(g, mass = 20, charge = E(g)$weight)
igraph::plot.igraph(g,layout=coord_graphopt)

# # and we store the layout coordinates to plot the other networks in the same way
coords_nodes <- coord_graphopt
row.names(coords_nodes) <- V(g)$name
coords_nodes
# match the nodes of the network with the coordinates in the layout
idxs <- match(V(g)$name, row.names(coords_nodes))
coord_graphopt <- coords_nodes[idxs,]
vsize <- igraph::degree(g)/vcount(g)

#Plot the network sith the module as fill color and the annotation as the perimeter color
par(mar=c(0,0,0,0))
pdf(paste0('Ku80_protein-protein-network-shScr.pdf'), width=30, height=30)
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('M = ', format(mod[[2]], digits=2)))
#
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.shape='circle', vertex.frame.color=V(g)$color_border, vertex.frame.width=4, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('whole network M = ', format(mod[[2]], digits = 2)), edge.arrow.size=.2)
legend('bottomleft', legend=levels(mod[[1]]$annotation), pt.lwd=4, cex=2, col=functions_colours, pch=21, pt.bg='white', title='function')

#Create a data table containing the module number, the treatment and the degree per node
cur_mod <- mod[[1]]
cur_mod$name <- as.character(cur_mod$name)
cur_mod$module <- cur_mod$module + 1

cur_mod$color <- custom_palette[cur_mod$module]

cur_mod$treatment <- 'shScr'
cur_mod$degree <- degree(g, mode='all')[match(cur_mod$name, names(degree(g, mode='all')))]

#Plot the Participation (inter-module degree coefficient) versus the node degree (number of connections)
par(mar=c(4,4,4,4))
ggplot(cur_mod, aes(degree, participation, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')
#Save the plot
ggsave(paste0('degree-part-shScr.pdf'), width = 12, height = 8)

#Plot the connectivity (within-module degree coefficient) versus the degree (number of connectios)
ggplot(cur_mod, aes(degree, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')
#Save the plot
ggsave(paste0('degree-conn-shScr.pdf'), width = 12, height = 8)
# Plot the connectivity versus the participation (This plot allows an easy visualization of the connector and central nodes)
ggplot(cur_mod, aes(participation, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')
# Save the Plot
ggsave(paste0('part-conn-shScr.pdf'), width = 12, height = 8)

mod_overall <- cur_mod

#write the results tables
write.csv(mod_overall,'./Ku80ID_shNclExp_shScr_Net_mod.csv')


### Ku80ID shNcl experiment network modularity analysis ###
### Ku80ID shNcl network 
## Load the required files for the analysis
metadata <- read.csv('./Ku80ID_shNclExp_DataFrame.csv')#Open the Ku80ID node metadata tables, containing the annotations for each protein
metadata<- subset(metadata, !metadata$shNcl.Presence==0) #Subset the Ku80ID datafram to keep only the results of the Ku80ID shNcl condition
df <- read.delim('./Ku80ID_shNclvsCneg_Network_EdgeList.tsv', sep="\t",header=T) #Load the edge list of the Ku80ID shNcl Network
g <- graph_from_data_frame(df, directed = FALSE) # convert the edge list into an igraph network object
g <- simplify(g, remove.multiple=TRUE, remove.loops=TRUE, edge.attr.comb = "first") # Remove the duplicated edges and self-loops
edge.attributes(g) # Confirm that there is no NAs in the edge weights
adj_mat_g <- as_adjacency_matrix(g, attr = 'weight1', sparse = FALSE) # Transform the networks into an adjacency matrix
mod <- netcarto(adj_mat_g) # Run the netcarto modularity analysis algorithm 

# Now we run the netcarto modularity analysis 20 times to optimize the separation of modules
for(i in 1:20){
  temp_mod <- netcarto(adj_mat_g)
  if(temp_mod[[2]] > mod[[2]])
    mod <- temp_mod
}

# Translate the netcarto generated modules in a format compatible with igraph network objects
com_mod <- make_clusters(
  g,
  membership = mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1,
  algorithm = 'netcarto',
  merges = NULL,
  modularity = mod[[2]]
)

# Add the module and module colors as vertex attributes to the igraph networks
V(g)$module <- mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1
V(g)$color <- custom_palette[V(g)$module]
# add the functional annotations to the netcarto generated modules
mod[[1]]$annotation <- as.factor(metadata[match(mod[[1]]$name, metadata$NAME),]$Annotation)
V(g)$color_border <- functions_colours[as.numeric(mod[[1]][match(V(g)$name, mod[[1]]$name),]$annotation)]

## Generate the layout to map the modules into the network
cl_for_layout <- com_mod
E(g)$weight <- edge.weights(cl_for_layout, g)

##### ATTENTION : THE CODE BELOW HAS TO BE RUN AT LEAST ONCE - UNCOMMENT ON FIRST RUN
# # create layout - run this line several times until you find a layout that you like
coord_graphopt <- layout_with_graphopt(g, mass = 20, charge = E(g)$weight)
igraph::plot.igraph(g,layout=coord_graphopt)

# # and we store the layout coordinates to plot the other networks in the same way
coords_nodes <- coord_graphopt
row.names(coords_nodes) <- V(g)$name
coords_nodes
# match the nodes of the network with the coordinates in the layout
idxs <- match(V(g)$name, row.names(coords_nodes))
coord_graphopt <- coords_nodes[idxs,]
vsize <- igraph::degree(g)/vcount(g)

#Plot the network sith the module as fill color and the annotation as the perimeter color
par(mar=c(0,0,0,0))
pdf(paste0('Ku80_protein-protein-network-shNcl.pdf'), width=30, height=30)
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('M = ', format(mod[[2]], digits=2)))
#
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.shape='circle', vertex.frame.color=V(g)$color_border, vertex.frame.width=4, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('whole network M = ', format(mod[[2]], digits = 2)), edge.arrow.size=.2)
legend('bottomleft', legend=levels(mod[[1]]$annotation), pt.lwd=4, cex=2, col=functions_colours, pch=21, pt.bg='white', title='function')

#Create a data table containing the module number, the treatment and the degree per node
cur_mod <- mod[[1]]
cur_mod$name <- as.character(cur_mod$name)
cur_mod$module <- cur_mod$module + 1

cur_mod$color <- custom_palette[cur_mod$module]

cur_mod$treatment <- 'shNcl'
cur_mod$degree <- degree(g, mode='all')[match(cur_mod$name, names(degree(g, mode='all')))]

#Plot the Participation (inter-module degree coefficient) versus the node degree (number of connections)
par(mar=c(4,4,4,4))
ggplot(cur_mod, aes(degree, participation, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')
#Save the plot
ggsave(paste0('degree-part-shNcl.pdf'), width = 12, height = 8)

#Plot the connectivity (within-module degree coefficient) versus the degree (number of connectios)
ggplot(cur_mod, aes(degree, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')
#Save the plot
ggsave(paste0('degree-conn-shNcl.pdf'), width = 12, height = 8)
# Plot the connectivity versus the participation (This plot allows an easy visualization of the connector and central nodes)
ggplot(cur_mod, aes(participation, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')
# Save the Plot
ggsave(paste0('part-conn-shNcl.pdf'), width = 12, height = 8)

mod_overall <- cur_mod

#write the results tables
write.csv(mod_overall,'./Ku80ID_shNclExp_shNcl_Net_mod.csv')

##############################################
### Mre11ID shNcl experiment network modularity analysis ###
### Mre11ID full nuclear network (it contain all the interactors found in the shScr and shNcl conditions, it serves to use a common layout for all the Mre11ID networks)
## Load the required files for the analysis
metadata <- read.csv('./Mre11ID_shNclExp_DataFrame.csv')#Open the Mre11ID node metadata tables, containing the annotations for each protein
df <- read.delim('./Mre11ID_shNclExp_FullNetwork_EdgeList.tsv', sep="\t",header=T) #Load the edge list of the Mre11ID full Network
g <- graph_from_data_frame(df, directed = FALSE) # convert the edge list into an igraph network object
g <- simplify(g, remove.multiple=TRUE, remove.loops=TRUE, edge.attr.comb = "first") # Remove the duplicated edges and self-loops
edge.attributes(g) # Confirm that there is no NAs in the edge weights
adj_mat_g <- as_adjacency_matrix(g, attr = 'weight1', sparse = FALSE) # Transform the networks into an adjacency matrix
mod <- netcarto(adj_mat_g) # Run the netcarto modularity analysis algorithm 

# Now we run the netcarto modularity analysis 20 times to optimize the separation of modules
for(i in 1:20){
  temp_mod <- netcarto(adj_mat_g)
  if(temp_mod[[2]] > mod[[2]])
    mod <- temp_mod
}

# Translate the netcarto generated modules in a format compatible with igraph network objects
com_mod <- make_clusters(
  g,
  membership = mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1,
  algorithm = 'netcarto',
  merges = NULL,
  modularity = mod[[2]]
)

# Add the module and module colors as vertex attributes to the igraph networks
V(g)$module <- mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1
V(g)$color <- custom_palette[V(g)$module]
# add the functional annotations to the netcarto generated modules
mod[[1]]$annotation <- as.factor(metadata[match(mod[[1]]$name, metadata$NAME),]$Annotation)
V(g)$color_border <- functions_colours[as.numeric(mod[[1]][match(V(g)$name, mod[[1]]$name),]$annotation)]

## Generate the layout to map the modules into the network
cl_for_layout <- com_mod
E(g)$weight <- edge.weights(cl_for_layout, g)

##### ATTENTION : THE CODE BELOW HAS TO BE RUN AT LEAST ONCE - UNCOMMENT ON FIRST RUN
# # create layout - run this line several times until you find a layout that you like
coord_graphopt <- layout_with_graphopt(g, mass = 20, charge = E(g)$weight)
igraph::plot.igraph(g,layout=coord_graphopt)

# # and we store the layout coordinates to plot the other networks in the same way
coords_nodes <- coord_graphopt
row.names(coords_nodes) <- V(g)$name
coords_nodes
# match the nodes of the network with the coordinates in the layout
idxs <- match(V(g)$name, row.names(coords_nodes))
coord_graphopt <- coords_nodes[idxs,]
vsize <- igraph::degree(g)/vcount(g)

#Plot the network sith the module as fill color and the annotation as the perimeter color
par(mar=c(0,0,0,0))
pdf(paste0('Mre11_protein-protein-network-whole.pdf'), width=30, height=30)
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('M = ', format(mod[[2]], digits=2)))
#
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.shape='circle', vertex.frame.color=V(g)$color_border, vertex.frame.width=4, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('whole network M = ', format(mod[[2]], digits = 2)), edge.arrow.size=.2)
legend('bottomleft', legend=levels(mod[[1]]$annotation), pt.lwd=4, cex=2, col=functions_colours, pch=21, pt.bg='white', title='function')

#Create a data table containing the module number, the treatment and the degree per node
cur_mod <- mod[[1]]
cur_mod$name <- as.character(cur_mod$name)
cur_mod$module <- cur_mod$module + 1

cur_mod$color <- custom_palette[cur_mod$module]

cur_mod$treatment <- 'whole'
cur_mod$degree <- degree(g, mode='all')[match(cur_mod$name, names(degree(g, mode='all')))]

#Plot the Participation (inter-module degree coefficient) versus the node degree (number of connections)
par(mar=c(4,4,4,4))
ggplot(cur_mod, aes(degree, participation, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')
#Save the plot
ggsave(paste0('degree-part-whole.pdf'), width = 12, height = 8)

#Plot the connectivity (within-module degree coefficient) versus the degree (number of connectios)
ggplot(cur_mod, aes(degree, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')
#Save the plot
ggsave(paste0('degree-conn-whole.pdf'), width = 12, height = 8)
# Plot the connectivity versus the participation (This plot allows an easy visualization of the connector and central nodes)
ggplot(cur_mod, aes(participation, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')
# Save the Plot
ggsave(paste0('part-conn-whole.pdf'), width = 12, height = 8)

mod_overall <- cur_mod

#write the results tables
write.csv(mod_overall,'./Mre11ID_shNclExp_FullNetwork_overall_mod.csv')

### Mre11ID shNcl experiment network modularity analysis ###
### Mre11ID shScramble network 
## Load the required files for the analysis
metadata <- read.csv('./Mre11ID_shNclExp_DataFrame.csv')#Open the Mre11ID node metadata tables, containing the annotations for each protein
metadata<- subset(metadata, !metadata$shScr.Presence==0) #Subset the Mre11ID datafram to keep only the results of the Mre11ID shScrable condition
df <- read.delim('./Mre11ID_shScrlvsCneg_Network_EsgeList.tsv', sep="\t",header=T) #Load the edge list of the Mre11ID shScr Network
g <- graph_from_data_frame(df, directed = FALSE) # convert the edge list into an igraph network object
g <- simplify(g, remove.multiple=TRUE, remove.loops=TRUE, edge.attr.comb = "first") # Remove the duplicated edges and self-loops
edge.attributes(g) # Confirm that there is no NAs in the edge weights
adj_mat_g <- as_adjacency_matrix(g, attr = 'weight1', sparse = FALSE) # Transform the networks into an adjacency matrix
mod <- netcarto(adj_mat_g) # Run the netcarto modularity analysis algorithm 

# Now we run the netcarto modularity analysis 20 times to optimize the separation of modules
for(i in 1:20){
  temp_mod <- netcarto(adj_mat_g)
  if(temp_mod[[2]] > mod[[2]])
    mod <- temp_mod
}

# Translate the netcarto generated modules in a format compatible with igraph network objects
com_mod <- make_clusters(
  g,
  membership = mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1,
  algorithm = 'netcarto',
  merges = NULL,
  modularity = mod[[2]]
)

# Add the module and module colors as vertex attributes to the igraph networks
V(g)$module <- mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1
V(g)$color <- custom_palette[V(g)$module]
# add the functional annotations to the netcarto generated modules
mod[[1]]$annotation <- as.factor(metadata[match(mod[[1]]$name, metadata$NAME),]$Annotation)
V(g)$color_border <- functions_colours[as.numeric(mod[[1]][match(V(g)$name, mod[[1]]$name),]$annotation)]

## Generate the layout to map the modules into the network
cl_for_layout <- com_mod
E(g)$weight <- edge.weights(cl_for_layout, g)

##### ATTENTION : THE CODE BELOW HAS TO BE RUN AT LEAST ONCE - UNCOMMENT ON FIRST RUN
# # create layout - run this line several times until you find a layout that you like
coord_graphopt <- layout_with_graphopt(g, mass = 20, charge = E(g)$weight)
igraph::plot.igraph(g,layout=coord_graphopt)

# # and we store the layout coordinates to plot the other networks in the same way
coords_nodes <- coord_graphopt
row.names(coords_nodes) <- V(g)$name
coords_nodes
# match the nodes of the network with the coordinates in the layout
idxs <- match(V(g)$name, row.names(coords_nodes))
coord_graphopt <- coords_nodes[idxs,]
vsize <- igraph::degree(g)/vcount(g)

#Plot the network sith the module as fill color and the annotation as the perimeter color
par(mar=c(0,0,0,0))
pdf(paste0('Mre11_protein-protein-network-shScr.pdf'), width=30, height=30)
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('M = ', format(mod[[2]], digits=2)))
#
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.shape='circle', vertex.frame.color=V(g)$color_border, vertex.frame.width=4, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('whole network M = ', format(mod[[2]], digits = 2)), edge.arrow.size=.2)
legend('bottomleft', legend=levels(mod[[1]]$annotation), pt.lwd=4, cex=2, col=functions_colours, pch=21, pt.bg='white', title='function')

#Create a data table containing the module number, the treatment and the degree per node
cur_mod <- mod[[1]]
cur_mod$name <- as.character(cur_mod$name)
cur_mod$module <- cur_mod$module + 1

cur_mod$color <- custom_palette[cur_mod$module]

cur_mod$treatment <- 'shScr'
cur_mod$degree <- degree(g, mode='all')[match(cur_mod$name, names(degree(g, mode='all')))]

#Plot the Participation (inter-module degree coefficient) versus the node degree (number of connections)
par(mar=c(4,4,4,4))
ggplot(cur_mod, aes(degree, participation, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')
#Save the plot
ggsave(paste0('degree-part-shScr.pdf'), width = 12, height = 8)

#Plot the connectivity (within-module degree coefficient) versus the degree (number of connectios)
ggplot(cur_mod, aes(degree, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')
#Save the plot
ggsave(paste0('degree-conn-shScr.pdf'), width = 12, height = 8)
# Plot the connectivity versus the participation (This plot allows an easy visualization of the connector and central nodes)
ggplot(cur_mod, aes(participation, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')
# Save the Plot
ggsave(paste0('part-conn-shScr.pdf'), width = 12, height = 8)

mod_overall <- cur_mod

#write the results tables
write.csv(mod_overall,'./Mre11ID_shNclExp_shScr_Net_mod.csv')


### Mre11ID shNcl experiment network modularity analysis ###
### Mre11ID shNcl network 
## Load the required files for the analysis
metadata <- read.csv('./Mre11ID_shNclExp_DataFrame.csv')#Open the Mre11ID node metadata tables, containing the annotations for each protein
metadata<- subset(metadata, !metadata$shNcl.Presence==0) #Subset the Mre11ID datafram to keep only the results of the Mre11ID shNcl condition
df <- read.delim('./Mre11ID_shNclvsCneg_Network_EdgeList.tsv', sep="\t",header=T) #Load the edge list of the Mre11ID shNcl Network
g <- graph_from_data_frame(df, directed = FALSE) # convert the edge list into an igraph network object
g <- simplify(g, remove.multiple=TRUE, remove.loops=TRUE, edge.attr.comb = "first") # Remove the duplicated edges and self-loops
edge.attributes(g) # Confirm that there is no NAs in the edge weights
adj_mat_g <- as_adjacency_matrix(g, attr = 'weight1', sparse = FALSE) # Transform the networks into an adjacency matrix
mod <- netcarto(adj_mat_g) # Run the netcarto modularity analysis algorithm 

# Now we run the netcarto modularity analysis 20 times to optimize the separation of modules
for(i in 1:20){
  temp_mod <- netcarto(adj_mat_g)
  if(temp_mod[[2]] > mod[[2]])
    mod <- temp_mod
}

# Translate the netcarto generated modules in a format compatible with igraph network objects
com_mod <- make_clusters(
  g,
  membership = mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1,
  algorithm = 'netcarto',
  merges = NULL,
  modularity = mod[[2]]
)

# Add the module and module colors as vertex attributes to the igraph networks
V(g)$module <- mod[[1]]$module[match(V(g)$name, mod[[1]]$name)] + 1
V(g)$color <- custom_palette[V(g)$module]
# add the functional annotations to the netcarto generated modules
mod[[1]]$annotation <- as.factor(metadata[match(mod[[1]]$name, metadata$NAME),]$Annotation)
V(g)$color_border <- functions_colours[as.numeric(mod[[1]][match(V(g)$name, mod[[1]]$name),]$annotation)]

## Generate the layout to map the modules into the network
cl_for_layout <- com_mod
E(g)$weight <- edge.weights(cl_for_layout, g)

##### ATTENTION : THE CODE BELOW HAS TO BE RUN AT LEAST ONCE - UNCOMMENT ON FIRST RUN
# # create layout - run this line several times until you find a layout that you like
coord_graphopt <- layout_with_graphopt(g, mass = 20, charge = E(g)$weight)
igraph::plot.igraph(g,layout=coord_graphopt)

# # and we store the layout coordinates to plot the other networks in the same way
coords_nodes <- coord_graphopt
row.names(coords_nodes) <- V(g)$name
coords_nodes
# match the nodes of the network with the coordinates in the layout
idxs <- match(V(g)$name, row.names(coords_nodes))
coord_graphopt <- coords_nodes[idxs,]
vsize <- igraph::degree(g)/vcount(g)

#Plot the network sith the module as fill color and the annotation as the perimeter color
par(mar=c(0,0,0,0))
pdf(paste0('Mre11_protein-protein-network-shNcl.pdf'), width=30, height=30)
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('M = ', format(mod[[2]], digits=2)))
#
plot(g, layout=coord_graphopt, vertex.size=vsize, vertex.shape='circle', vertex.frame.color=V(g)$color_border, vertex.frame.width=4, vertex.label.cex=1, vertex.label.font=2, vertex.label.color='black', main=paste0('whole network M = ', format(mod[[2]], digits = 2)), edge.arrow.size=.2)
legend('bottomleft', legend=levels(mod[[1]]$annotation), pt.lwd=4, cex=2, col=functions_colours, pch=21, pt.bg='white', title='function')

#Create a data table containing the module number, the treatment and the degree per node
cur_mod <- mod[[1]]
cur_mod$name <- as.character(cur_mod$name)
cur_mod$module <- cur_mod$module + 1

cur_mod$color <- custom_palette[cur_mod$module]

cur_mod$treatment <- 'shNcl'
cur_mod$degree <- degree(g, mode='all')[match(cur_mod$name, names(degree(g, mode='all')))]

#Plot the Participation (inter-module degree coefficient) versus the node degree (number of connections)
par(mar=c(4,4,4,4))
ggplot(cur_mod, aes(degree, participation, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')
#Save the plot
ggsave(paste0('degree-part-shNcl.pdf'), width = 12, height = 8)

#Plot the connectivity (within-module degree coefficient) versus the degree (number of connectios)
ggplot(cur_mod, aes(degree, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')
#Save the plot
ggsave(paste0('degree-conn-shNcl.pdf'), width = 12, height = 8)
# Plot the connectivity versus the participation (This plot allows an easy visualization of the connector and central nodes)
ggplot(cur_mod, aes(participation, connectivity, color=annotation, fill=as.factor(module))) + 
  geom_point(size=4, pch=21, stroke=2) + theme_bw() + ggtitle('whole') +
  geom_text(aes(label=ifelse( (!role %in% c('Peripheral', 'Ultra peripheral')), name,'')), size=3, color='black') +
  scale_fill_manual(values=custom_palette[1:length(unique(cur_mod$module))], name='Module') +
  scale_color_manual(values=functions_colours[1:length(as.numeric(cur_mod$annotation))], name='Function')
# Save the Plot
ggsave(paste0('part-conn-shNcl.pdf'), width = 12, height = 8)

mod_overall <- cur_mod

#write the results tables
write.csv(mod_overall,'./Mre11ID_shNclExp_shNcl_Net_mod.csv')




