
### This R Script uses the package Mfuzz to generate dynamic clusters 
#Install packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Mfuzz")
BiocManager::install("marray")

#Open the required librarys
require(Mfuzz)
require(marray)


#Set working directory
setwd("C:/Users/alfredog/Desktop/SensorID_Manuscript/GitHub_Scripts/Fuzzy_C_Means_Analysis")

### Sirt6ID Fuzzy C Means analysis using the FoldChange as a normalized measurement of protein intensity changes
#Load the Fold Change values and the metadata tables
DS<-read.csv("./Sirt6_Data.csv",row.names = "NAME")#Table with all the FoldChange values
DS<-as.matrix(DS)
DS<-ifelse(DS==0,1,DS)# The 0 fold changes were transformed to 1, since the standardization function of mfuzz will convert 0s to NAs
Anno<-read.csv("./Protein_Annotations.csv",row.names = "NAME")#Open the protein annotations for plotting functions
Anno<-Anno[rownames(Anno)%in%row.names(DS),]#subset the protein annotation to contain only those proteins included in the Sirt6ID dataset
Anno<-Anno[rownames(DS),]#reorder the rownames of the annotation data to match the order in the Fold Change matrix
is.na(match(rownames(DS),rownames(Anno)))#Make sure the name of the proteins is the same in the fold change matrix and in the Annotation data frame
colAnno<-read.csv("./ColAnnoDF.csv",row.names = "Sample") #Open the sample metadata
colAnno<-colAnno[rownames(colAnno)%in%colnames(DS),]
#Generate an expression dataset compatible with the mfuzz package
phenoData <- new("AnnotatedDataFrame", data=colAnno)
featData<-new("AnnotatedDataFrame", data=Anno)
Expr<- ExpressionSet(assayData = DS,
                     phenoData = phenoData,
                     featureData = featData,
                     annotation = "Sirt6ID LIMMA Analysis Fold Change")
#mfuzz Analysis

#Filter missing values and proteins with no variance in the experiment
Expr.r<-filter.NA(Expr,thres = 0.25)#Filter out proteins with less than 75% of its values missing
Expr.f<-fill.NA(Expr.r,mode = "knn")#Fill-in the NA values using the Knn method
tmp <- filter.std(Expr.f,min.std=0)#Remove those proteins with no variance in the experiment
#Standardisation: Scaling by mean centering at 0, and standard deviation scaling as 1
Expr.s <- standardise(Expr.f)#Standarize the NAs filled-in matrix
tmp <- filter.std(Expr.s,min.std=0)# Make sure the scaling of the matrix do not remove all the variance
#Recover the Standardized data into a new matrix
a<-Expr.s@assayData[["exprs"]]
write.table(a,file = "./Sirt6_Std.csv",sep = ",")
#Setting FCM clustering parameters
m1<-mestimate(Expr.s)# Calculate the fuzzifier m number to try to avoid the clustering of random data
m1 #1.638145
tmp  <- cselection(Expr.s,m=1.638145,crange=seq(2,10,1),repeats=5,visu=TRUE)#at cluster number above 8 the clustering generates empty clusters, hence we decided to use 5 clusters after using 7 and 6 clusters
#Generating clusters with m=1.63 and c=6
cl1 <- mfuzz(Expr.s,c=5,m=1.638145)
#Ploting
mfuzz.plot2(Expr.s,cl=cl1,mfrow=c(4,2),colo = "fancy",centre = TRUE,centre.col = "black",centre.lwd = 2,x11=FALSE,time.labels = c("nonIR5","IR5","IR30","IR2h","IR8h","IR24h","nonIR24h"))+
  mfuzzColorBar(col="fancy",main="Membership",cex.main=1)
#Overlap between clusters by PCA analysis
O1<-overlap(cl1)
overlap.plot(cl1,over=O1,thres=0.05)

#Retrieving fuzzy C means clusters and their components
clusters<-as.data.frame(cl1$cluster)
membership<-as.data.frame(cl1$membership)
Sirt6_Clusters<-cbind(clusters,membership)
Sirt6_Clusters$NAME<-rownames(Sirt6_Clusters)
Sirt6_Clusters
write.csv(Sirt6_Clusters,"./Sirt6_fuzzyCmeans_Clusters.csv",sep = ",",row.names = TRUE,col.names = colnames(Sirt6_Clusters))
#Retrieve the centers of each cluster
centers<-as.data.frame(cl1$centers)
centers
write.csv(centers,"./Sirt6_fuzzyCmeans_Centers.csv")
#Retrieve the size of each cluster
size<-as.data.frame(cl1$size)
write.csv(size,"./Sirt6_fuzzyCmeans_size.csv")

### After retrieving the Fuzzy C Means clusters, and the proteins grouped with each one, the cluster number was added to each protein manually in Excel to expand the results table.

#################################################
### Ku80ID Fuzzy C Means analysis using the FoldChange as a normalized measurement of protein intensity changes
#Load the Fold Change values and the metadata tables
DS<-read.csv("./Ku80_Data.csv",row.names = "NAME")#Table with all the FoldChange values
DS<-as.matrix(DS)
DS<-ifelse(DS==0,1,DS)# The 0 fold changes were transformed to 1, since the standardization function of mfuzz will convert 0s to NAs
Anno<-read.csv("./Protein_Annotations.csv",row.names = "NAME")#Open the protein annotations for plotting functions
Anno<-Anno[rownames(Anno)%in%row.names(DS),]#subset the protein annotation to contain only those proteins included in the Ku80ID dataset
Anno<-Anno[rownames(DS),]#reorder the rownames of the annotation data to match the order in the Fold Change matrix
is.na(match(rownames(DS),rownames(Anno)))#Make sure the name of the proteins is the same in the fold change matrix and in the Annotation data frame
colAnno<-read.csv("./ColAnnoDF.csv",row.names = "Sample") #Open the sample metadata
colAnno<-colAnno[rownames(colAnno)%in%colnames(DS),]
#Generate an expression dataset compatible with the mfuzz package
phenoData <- new("AnnotatedDataFrame", data=colAnno)
featData<-new("AnnotatedDataFrame", data=Anno)
Expr<- ExpressionSet(assayData = DS,
                     phenoData = phenoData,
                     featureData = featData,
                     annotation = "Ku80ID LIMMA Analysis Fold Change")

#mfuzz Analysis

#Filter missing values and proteins with no variance in the experiment
Expr.r<-filter.NA(Expr,thres = 0.25)#Filter out proteins with less than 75% of its values missing
Expr.f<-fill.NA(Expr.r,mode = "knn")#Fill-in the NA values using the Knn method
tmp <- filter.std(Expr.f,min.std=0)#Remove those proteins with no variance in the experiment
#Standardisation: Scaling by mean centering at 0, and standard deviation scaling as 1
Expr.s <- standardise(Expr.f)#Standarize the NAs filled-in matrix
tmp <- filter.std(Expr.s,min.std=0)# Make sure the scaling of the matrix do not remove all the variance
#Recover the Standardized data into a new matrix
a<-Expr.s@assayData[["exprs"]]
write.table(a,file = "./Ku80_Std.csv",sep = ",")
#Setting FCM clustering parameters
m1<-mestimate(Expr.s)# Calculate the fuzzifier m number to try to avoid the clustering of random data
m1 #1.648535
tmp  <- cselection(Expr.s,m=1.648535,crange=seq(2,10,1),repeats=5,visu=TRUE)#Interestigly with a fuzzifier m=1.62, we start getting empty clusters agter 7 cluster, hence after testing 6 and 5 number of clusters, 4 was the best option
#Generating clusters with m=1.62 and c=6
cl1 <- mfuzz(Expr.s,c=4,m=1.648535)
#Ploting
mfuzz.plot2(Expr.s,cl=cl1,mfrow=c(4,2),colo = "fancy",centre = TRUE,centre.col = "black",centre.lwd = 2,x11=FALSE,time.labels = c("nonIR5","IR5","IR30","IR2h","IR8h","IR24h","nonIR24h"))+
  mfuzzColorBar(col="fancy",main="Membership",cex.main=1)
#Overlap between clusters by PCA analysis
O1<-overlap(cl1)
overlap.plot(cl1,over=O1,thres=0.05)

#Retrieving fuzzy C means clusters and their components
clusters<-as.data.frame(cl1$cluster)
membership<-as.data.frame(cl1$membership)
Ku80_Clusters<-cbind(clusters,membership)
Ku80_Clusters$NAME<-rownames(Ku80_Clusters)
Ku80_Clusters
write.csv(Ku80_Clusters,"./Ku80_fuzzyCmeans_Clusters.csv",sep = ",",row.names = TRUE,col.names = colnames(Ku80_Clusters))
#Retrieve the centers of each cluster
centers<-as.data.frame(cl1$centers)
centers
write.csv(centers,"./Ku80_fuzzyCmeans_Centers.csv")
#Retrieve the size of each cluster
size<-as.data.frame(cl1$size)
write.csv(size,"./Ku80_fuzzyCmeans_size.csv")

### After retrieving the Fuzzy C Means clusters, and the proteins grouped with each one, the cluster number was added to each protein manually in Excel to expand the results table.


################
### Mre11ID Fuzzy C Means analysis using the FoldChange as a normalized measurement of protein intensity changes
#Load the Fold Change values and the metadata tables
DS<-read.csv("./Mre11_Data.csv",row.names = "NAME")#Table with all the FoldChange values
DS<-as.matrix(DS)
DS<-ifelse(DS==0,1,DS)# The 0 fold changes were transformed to 1, since the standardization function of mfuzz will convert 0s to NAs
Anno<-read.csv("./Protein_Annotations.csv",row.names = "NAME")#Open the protein annotations for plotting functions
Anno<-Anno[rownames(Anno)%in%row.names(DS),]#subset the protein annotation to contain only those proteins included in the Mre11ID dataset
Anno<-Anno[rownames(DS),]#reorder the rownames of the annotation data to match the order in the Fold Change matrix
is.na(match(rownames(DS),rownames(Anno)))#Make sure the name of the proteins is the same in the fold change matrix and in the Annotation data frame
colAnno<-read.csv("./ColAnnoDF.csv",row.names = "Sample") #Open the sample metadata
colAnno<-colAnno[rownames(colAnno)%in%colnames(DS),]
#Generate an expression dataset compatible with the mfuzz package
phenoData <- new("AnnotatedDataFrame", data=colAnno)
featData<-new("AnnotatedDataFrame", data=Anno)
Expr<- ExpressionSet(assayData = DS,
                     phenoData = phenoData,
                     featureData = featData,
                     annotation = "Mre11ID LIMMA Analysis Fold Change")

#mfuzz Analysis

#Filter missing values and proteins with no variance in the experiment
Expr.r<-filter.NA(Expr,thres = 0.25)#Filter out proteins with less than 75% of its values missing
Expr.f<-fill.NA(Expr.r,mode = "knn")#Fill-in the NA values using the Knn method
tmp <- filter.std(Expr.f,min.std=0)#Remove those proteins with no variance in the experiment
#Standardisation: Scaling by mean centering at 0, and standard deviation scaling as 1
Expr.s <- standardise(Expr.f)#Standarize the NAs filled-in matrix
tmp <- filter.std(Expr.s,min.std=0)# Make sure the scaling of the matrix do not remove all the variance
#Recover the Standardized data into a new matrix
a<-Expr.s@assayData[["exprs"]]
write.table(a,file = "./Mre11_Std.csv",sep = ",")
#Setting FCM clustering parameters
m1<-mestimate(Expr.s)# Calculate the fuzzifier m number to try to avoid the clustering of random data
m1 #1.627608
tmp  <- cselection(Expr.s,m=1.627608,crange=seq(2,100,1),repeats=5,visu=TRUE)#Interestigly with a fuzzifier m=1.62, we can obtain up to 100 non-empty clusters, hence we are going to attach to what we already saw in Sirt6 and Ku80 for 4 or 5 clusters max
#Generating clusters with m=1.62 and c=6
cl1 <- mfuzz(Expr.s,c=6,m=1.627608)
#Ploting
mfuzz.plot2(Expr.s,cl=cl1,mfrow=c(4,2),colo = "fancy",centre = TRUE,centre.col = "black",centre.lwd = 2,x11=FALSE,time.labels = c("nonIR5","IR5","IR30","IR2h","IR8h","IR24h","nonIR24h"))+
mfuzzColorBar(col="fancy",main="Membership",cex.main=1)
#Overlap between clusters by PCA analysis
O1<-overlap(cl1)
overlap.plot(cl1,over=O1,thres=0.05)

#Retrieving fuzzy C means clusters and their components
clusters<-as.data.frame(cl1$cluster)
membership<-as.data.frame(cl1$membership)
Mre11_Clusters<-cbind(clusters,membership)
Mre11_Clusters$NAME<-rownames(Mre11_Clusters)
Mre11_Clusters
write.csv(Mre11_Clusters,"./Mre11_fuzzyCmeans_Clusters.csv",sep = ",",row.names = TRUE,col.names = colnames(Mre11_Clusters))
#Retrieve the centers of each cluster
centers<-as.data.frame(cl1$centers)
centers
write.csv(centers,"./Mre11_fuzzyCmeans_Centers.csv")
#Retrieve the size of each cluster
size<-as.data.frame(cl1$size)
write.csv(size,"./Mre11_fuzzyCmeans_size.csv")

### After retrieving the Fuzzy C Means clusters, and the proteins grouped with each one, the cluster number was added to each protein manually in Excel to expand the results table.


