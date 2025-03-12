#Load the required packages
require(limma)
require(qvalue)
require(ggplot2)
require(EnhancedVolcano)

#Set up the working directory, please use the path where the raw data files were downloaded
setwd("C:/../Limma_t_Test_Analysis")

### Sirt6ID analysis
#Open the raw data files
dat<-read.csv("./Sirt6ID_RawData.csv",header = TRUE)#Open the raw data table, and coherce the column name as it is

#Review the structure of the dataset
dim(dat)
#[1] 1262   33
str(dat)

#define vector headlines
cha<- c("CN_1" , "CN_2" , "CN_3" , "CN_4" , "CN_5" , "CN_6" , "S6nonIR5_1" , "S6nonIR5_2" , "S6nonIR5_3" , "S6nonIR24_1" , "S6nonIR24_2" , "S6nonIR24_3" , "S6IR5_1" , "S6IR5_2" , "S6IR5_3" , "S6IR30_1" , "S6IR30_2" , "S6IR30_3" , "S6IR2_1" , "S6IR2_2" , "S6IR2_3" , "S6IR8_1" , "S6IR8_2" , "S6IR8_3" , "S6IR24_1" , "S6IR24_2" , "S6IR24_3")

#data procession download functions
source("http://www.biostat.jhsph.edu/~kkammers/software/eupa/source.functions.r")

#Readpeptides
#dat<- read.peptides(dat , cha)#Since we do not imput the peptide reads, but already the protein groups this one is not used
#scale the protein intensity measurements from the mass spectrometry results
dat<-quantify.proteins(dat , cha)# median polished log2 transformed values, as a sacaled intensity for each protein 

#boxplot: intensities of all channels after data processing and normalization
par(mfrow = c(1,1) , font.lab=2 , cex.lab=1.2 , font.axis=2 , cex.axis=1.2)
boxplot(dat[,1:length(cha)], ylim=c(-10,10) , main="boxplot normalized intensities")

##heatmap of the scaled intensities
require(ComplexHeatmap)
#generate the matrix for the heat map
mat<-as.matrix(mat)#Transform the dataset as matrix
mat<-mat[,7:33]
any(is.na(mat))#Check for NA values
mean(mat)# Mean
median(mat)# median
quantile(mat,c(0,0.25,0.5,0.75,1)) #Range of the dataset

#Color function
col_fun<-colorRamp2(c(min(mat),mean(mat),max(mat)),c("royalblue","aliceblue","tomato"))
col_fun(seq( 0, 15))

#Heatmap plot
Heatmap(mat , name = "Fold Change vs Background" ,  col = col_fun , cluster_rows = T , 
        cluster_columns = FALSE , clustering_distance_rows = "spearman" , clustering_distance_columns = "spearman" ,
        clustering_method_rows = "single" , clustering_method_columns = "single" , column_names_centered = TRUE , 
        column_names_gp = gpar(fontsize=15) , row_names_gp = gpar(fontsize=3) ,
        row_split = GO_Anno ,  row_gap = unit(1, "mm"))

### T-test nonIR5 vs Cneg
#define treatment and control groups for comparison
ct<-c("CN_1" , "CN_2" , "CN_3" , "CN_4" , "CN_5" , "CN_6")
trnonIR5<-c("S6nonIR5_1" , "S6nonIR5_2" , "S6nonIR5_3")

#define design according to syntax of limma package
designnonIR5<-model.matrix(~factor(c(1,1,1,1,1,1,2,2,2)))

#Limma statistical analysis
colnames(designnonIR5)<-c("intercept","Diff")
res_nonIR5.eb<-eb.fit(dat[,c(ct,trnonIR5)] , designnonIR5)
head(res_nonIR5.eb)

####### Export tables
write.table(res_nonIR5.eb, file = "DiffBind_S6nonIR5_Limma.csv" , sep = ",")

#Volcanoplot
#trying to correct xylim problems
rx <- c(-1,1)*max(abs(res_nonIR5.eb$logFC))*1.1
ry <- c(0, 11)

par(mfrow = c(1,2) , font.lab=2 , cex.lab=1.2 , font.axis=2 , cex.axis=1.2)
par(las=1 , xaxs="i" , yaxs="i")

#Ordinary pValues
plot(res_nonIR5.eb$logFC, -log10(res_nonIR5.eb$p.ord), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="gray", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="gray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of ordinary p-Values")

#Moderate pValues
plot(res_nonIR5.eb$logFC, -log10(res_nonIR5.eb$p.mod), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="lightgray", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="lightgray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of moderated p-Values")

#Enhanced volcano
library(EnhancedVolcano)

EnhancedVolcano(toptable = res_nonIR5.eb, lab = rownames(res_nonIR5.eb) , x = "logFC" , y = "p.ord" , xlim = c(min(res_nonIR5.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res_nonIR5.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res_nonIR5.eb[["p.ord"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot S6nonIR5" , subtitle = bquote(italic("Limma T-test")) , caption = paste0("total=" , nrow(res_nonIR5.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))
EnhancedVolcano(toptable = res_nonIR5.eb, lab = rownames(res_nonIR5.eb) , x = "logFC" , y = "p.mod" , xlim = c(min(res_nonIR5.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res_nonIR5.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res_nonIR5.eb[["p.mod"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot S6nonIR5" , subtitle = bquote(italic("Limma moderated T-test")) , caption = paste0("total=" , nrow(res_nonIR5.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))


### T-test nonIR24 vs Cneg
#define treatment and control groups for comparison
ct<-c("CN_1" , "CN_2" , "CN_3" , "CN_4" , "CN_5" , "CN_6")
trnonIR24<-c("S6nonIR24_1" , "S6nonIR24_2" , "S6nonIR24_3")

#define design according to syntax of limma package
designnonIR24<-model.matrix(~factor(c(1,1,1,1,1,1,2,2,2)))

#Limma statistical analysis
colnames(designnonIR24)<-c("intercept","Diff")
res_nonIR24.eb<-eb.fit(dat[,c(ct,trnonIR24)] , designnonIR24)
head(res_nonIR24.eb)

####### Export tables
write.table(res_nonIR24.eb, file = "DiffBind_S6nonIR24_Limma.csv" , sep = ",")

#Volcanoplot
#trying to correct xylim problems
rx <- c(-1,1)*max(abs(res_nonIR24.eb$logFC))*1.1
ry <- c(0, 16)

par(mfrow = c(1,2) , font.lab=2 , cex.lab=1.2 , font.axis=2 , cex.axis=1.2)
par(las=1 , xaxs="i" , yaxs="i")

# Ordinary pValues
plot(res_nonIR24.eb$logFC, -log10(res_nonIR24.eb$p.ord), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="gray", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="gray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of ordinary p-Values")

#Moderate pValues
plot(res_nonIR24.eb$logFC, -log10(res_nonIR24.eb$p.mod), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="lightgray", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="lightgray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of moderated p-Values")

#Enhanced volcano
library(EnhancedVolcano)

EnhancedVolcano(toptable = res_nonIR24.eb, lab = rownames(res_nonIR24.eb) , x = "logFC" , y = "p.ord" , xlim = c(min(res_nonIR24.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res_nonIR24.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res_nonIR24.eb[["p.ord"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot S6nonIR24" , subtitle = bquote(italic("Limma T-test")) , caption = paste0("total=" , nrow(res_nonIR24.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))
EnhancedVolcano(toptable = res_nonIR24.eb, lab = rownames(res_nonIR24.eb) , x = "logFC" , y = "p.mod" , xlim = c(min(res_nonIR24.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res_nonIR24.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res_nonIR24.eb[["p.mod"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot S6nonIR24" , subtitle = bquote(italic("Limma moderated T-test")) , caption = paste0("total=" , nrow(res_nonIR24.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))

### T-Test S6IR6 vs Negative controls
#S6IR5
#define treatment and control groups for comparison
ct<-c("CN_1" , "CN_2" , "CN_3" , "CN_4" , "CN_5" , "CN_6")
trIR5<-c("S6IR5_1" , "S6IR5_2" , "S6IR5_3")

#define design according to syntax of limma package
designIR5<-model.matrix(~factor(c(1,1,1,1,1,1,2,2,2)))

#Limma statistical analysis
colnames(designIR5)<-c("intercept","Diff")
res_IR5.eb<-eb.fit(dat[,c(ct,trIR5)] , designIR5)
head(res_IR5.eb)

####### Export tables
write.table(res_IR5.eb, file = "DiffBind_S6IR5_Limma.csv" , sep = ",")

#Volcanoplot
#trying to correct xylim problems
rx <- c(-1,1)*max(abs(res_IR5.eb$logFC))*1.1
ry <- c(0, 11)


par(mfrow = c(1,2) , font.lab=2 , cex.lab=1.2 , font.axis=2 , cex.axis=1.2)
par(las=1 , xaxs="i" , yaxs="i")

#Ordinary pValues
plot(res_IR5.eb$logFC, -log10(res_IR5.eb$p.ord), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="gray", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="gray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of ordinary p-Values")

#Moderate pValues
plot(res_IR5.eb$logFC, -log10(res_IR5.eb$p.mod), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="lightgray", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="lightgray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of moderated p-Values")

#Enhanced volcano
library(EnhancedVolcano)

EnhancedVolcano(toptable = res_IR5.eb, lab = rownames(res_IR5.eb) , x = "logFC" , y = "p.ord" , xlim = c(min(res_IR5.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res_IR5.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res_IR5.eb[["p.ord"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot S6IR5" , subtitle = bquote(italic("Limma T-test")) , caption = paste0("total=" , nrow(res_IR5.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))
EnhancedVolcano(toptable = res_IR5.eb, lab = rownames(res_IR5.eb) , x = "logFC" , y = "p.mod" , xlim = c(min(res_IR5.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res_IR5.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res_IR5.eb[["p.mod"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot S6IR5" , subtitle = bquote(italic("Limma moderated T-test")) , caption = paste0("total=" , nrow(res_IR5.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))

### T-Test S6IR30 vs Negative controls
#S6IR30
#define treatment and control groups for comparison
ct<-c("CN_1" , "CN_2" , "CN_3" , "CN_4" , "CN_5" , "CN_6")
trIR30<-c("S6IR30_1" , "S6IR30_2" , "S6IR30_3")

#define design according to syntax of limma package
designIR30<-model.matrix(~factor(c(1,1,1,1,1,1,2,2,2)))

#Limma statistical analysis
colnames(designIR30)<-c("intercept","Diff")
res_IR30.eb<-eb.fit(dat[,c(ct,trIR30)] , designIR30)
head(res_IR30.eb)


####### Export tables
write.table(res_IR30.eb, file = "DiffBind_S6IR30_Limma.csv" , sep = ",")

#Volcanoplot
#trying to correct xylim problems
rx <- c(-1,1)*max(abs(res_IR30.eb$logFC))*1.1
ry <- c(0, 14)


par(mfrow = c(1,2) , font.lab=2 , cex.lab=1.2 , font.axis=2 , cex.axis=1.2)
par(las=1 , xaxs="i" , yaxs="i")

#Ordinary pValue
plot(res_IR30.eb$logFC, -log10(res_IR30.eb$p.ord), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="gray", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="gray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of ordinary p-Values")

#Moderate pValues
plot(res_IR30.eb$logFC, -log10(res_IR30.eb$p.mod), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="lightgray", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="lightgray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of moderated p-Values")

#Enhanced volcano
library(EnhancedVolcano)

EnhancedVolcano(toptable = res_IR30.eb, lab = rownames(res_IR30.eb) , x = "logFC" , y = "p.ord" , xlim = c(min(res_IR30.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res_IR30.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res_IR30.eb[["p.ord"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot S6IR30" , subtitle = bquote(italic("Limma T-test")) , caption = paste0("total=" , nrow(res_IR30.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))
EnhancedVolcano(toptable = res_IR30.eb, lab = rownames(res_IR30.eb) , x = "logFC" , y = "p.mod" , xlim = c(min(res_IR30.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res_IR30.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res_IR30.eb[["p.mod"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot S6IR30" , subtitle = bquote(italic("Limma moderated T-test")) , caption = paste0("total=" , nrow(res_IR30.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))


### T-Test S6IR2h vs Negative Control
#S6IR2
#define treatment and control groups for comparison
ct<-c("CN_1" , "CN_2" , "CN_3" , "CN_4" , "CN_5" , "CN_6")
trIR2<-c("S6IR2_1" , "S6IR2_2" , "S6IR2_3")

#define design according to syntax of limma package
designIR2<-model.matrix(~factor(c(1,1,1,1,1,1,2,2,2)))

#Limma statistical analysis
colnames(designIR2)<-c("intercept","Diff")
res_IR2.eb<-eb.fit(dat[,c(ct,trIR2)] , designIR2)
head(res_IR2.eb)


####### Export tables
write.table(res_IR2.eb, file = "DiffBind_S6IR2_Limma.csv" , sep = ",")

#Volcanoplot
#trying to correct xylim problems
rx <- c(-1,1)*max(abs(res_IR2.eb$logFC))*1.1
ry <- c(0, 11)


par(mfrow = c(1,2) , font.lab=2 , cex.lab=1.2 , font.axis=2 , cex.axis=1.2)
par(las=1 , xaxs="i" , yaxs="i")

#Ordinary pValues
plot(res_IR2.eb$logFC, -log10(res_IR2.eb$p.ord), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="gray", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="gray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of ordinary p-Values")

#Moderate pValues
plot(res_IR2.eb$logFC, -log10(res_IR2.eb$p.mod), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="lightgray", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="lightgray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of moderated p-Values")

#Enhanced volcano
library(EnhancedVolcano)

EnhancedVolcano(toptable = res_IR2.eb, lab = rownames(res_IR2.eb) , x = "logFC" , y = "p.ord" , xlim = c(min(res_IR2.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res_IR2.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res_IR2.eb[["p.ord"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot S6IR2" , subtitle = bquote(italic("Limma T-test")) , caption = paste("total=" , nrow(res_IR2.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))
EnhancedVolcano(toptable = res_IR2.eb, lab = rownames(res_IR2.eb) , x = "logFC" , y = "p.mod" , xlim = c(min(res_IR2.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res_IR2.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res_IR2.eb[["p.mod"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot S6IR2" , subtitle = bquote(italic("Limma moderated T-test")) , caption = paste0("total=" , nrow(res_IR2.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))


###T-Test S6IR8h vs Negative controls
#S6IR8
#define treatment and control groups for comparison
ct<-c("CN_1" , "CN_2" , "CN_3" , "CN_4" , "CN_5" , "CN_6")
trIR8<-c("S6IR8_1" , "S6IR8_2" , "S6IR8_3")

#define design according to syntax of limma package
designIR8<-model.matrix(~factor(c(1,1,1,1,1,1,2,2,2)))

#Limma statistical analysis
colnames(designIR8)<-c("intercept","Diff")
res_IR8.eb<-eb.fit(dat[,c(ct,trIR8)] , designIR8)
head(res_IR8.eb)

####### Export tables
write.table(res_IR8.eb, file = "DiffBind_S6IR8_Limma.csv" , sep = ",")

#Volcanoplot
#trying to correct xylim problems
rx <- c(-1,1)*max(abs(res_IR2.eb$logFC))*1.1
ry <- c(0, 11)


par(mfrow = c(1,2) , font.lab=2 , cex.lab=1.2 , font.axis=2 , cex.axis=1.2)
par(las=1 , xaxs="i" , yaxs="i")

#Ordinary pValues
plot(res_IR8.eb$logFC, -log10(res_IR8.eb$p.ord), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="gray", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="gray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of ordinary p-Values")

#Moderate pValues
plot(res_IR8.eb$logFC, -log10(res_IR8.eb$p.mod), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="lightgray", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="lightgray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of moderated p-Values")

#Enhanced volcano
library(EnhancedVolcano)

EnhancedVolcano(toptable = res_IR8.eb, lab = rownames(res_IR8.eb) , x = "logFC" , y = "p.ord" , xlim = c(min(res_IR8.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res_IR8.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res_IR8.eb[["p.ord"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot S6IR8" , subtitle = bquote(italic("Limma T-test")) , caption = paste("total=" , nrow(res_IR8.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))
EnhancedVolcano(toptable = res_IR8.eb, lab = rownames(res_IR8.eb) , x = "logFC" , y = "p.mod" , xlim = c(min(res_IR8.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res_IR8.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res_IR8.eb[["p.mod"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot S6IR8" , subtitle = bquote(italic("Limma moderated T-test")) , caption = paste0("total=" , nrow(res_IR8.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))


### T-Test S6IR24h vs Negative controls
#S6IR24
#define treatment and control groups for comparison
ct<-c("CN_1" , "CN_2" , "CN_3" , "CN_4" , "CN_5" , "CN_6")
trIR24<-c("S6IR24_1" , "S6IR24_2" , "S6IR24_3")

#define design according to syntax of limma package
designIR24<-model.matrix(~factor(c(1,1,1,1,1,1,2,2,2)))

#Limma statistical analysis
colnames(designIR24)<-c("intercept","Diff")
res_IR24.eb<-eb.fit(dat[,c(ct,trIR24)] , designIR24)
head(res_IR24.eb)

####### Export tables
write.table(res_IR24.eb, file = "DiffBind_S6IR24_Limma.csv" , sep = ",")

#Volcanoplot
#trying to correct xylim problems
rx <- c(-1,1)*max(abs(res_IR24.eb$logFC))*1.1
ry <- c(0, 17)

par(mfrow = c(1,2) , font.lab=2 , cex.lab=1.2 , font.axis=2 , cex.axis=1.2)
par(las=1 , xaxs="i" , yaxs="i")

plot(res_IR24.eb$logFC, -log10(res_IR24.eb$p.ord), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="gray", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="gray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of ordinary p-Values")

#Moderate pValues
plot(res_IR24.eb$logFC, -log10(res_IR24.eb$p.mod), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="lightgray", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="lightgray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of moderated p-Values")

#Enhanced volcano
library(EnhancedVolcano)

EnhancedVolcano(toptable = res_IR24.eb, lab = rownames(res_IR24.eb) , x = "logFC" , y = "p.ord" , xlim = c(min(res_IR24.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res_IR24.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res_IR24.eb[["p.ord"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot S6IR24" , subtitle = bquote(italic("Limma T-test")) , caption = paste("total=" , nrow(res_IR24.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))
EnhancedVolcano(toptable = res_IR24.eb, lab = rownames(res_IR24.eb) , x = "logFC" , y = "p.mod" , xlim = c(min(res_IR24.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res_IR24.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res_IR24.eb[["p.mod"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot S6IR24" , subtitle = bquote(italic("Limma moderated T-test")) , caption = paste0("total=" , nrow(res_IR24.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))


##################################
###Ku80ID analysis
dat = read.csv("./Ku80ID_RawData.csv",header=TRUE)#Open the Ku80ID raw data table, and coherce the column names as it is
#Review the structure of the data
dim(dat)
str(dat)

#define vector headlines
cha<- c("CN_1" , "CN_2" , "CN_3" , "CN_4" , "CN_5" , "CN_6" , "KunonIR5_1" , "KunonIR5_2" , "KunonIR5_3" , "KunonIR24_1" , "KunonIR24_2" , "KunonIR24_3" , "KuIR5_1" , "KuIR5_2" , "KuIR5_3" , "KuIR30_1" , "KuIR30_2" , "KuIR30_3" , "KuIR2_1" , "KuIR2_2" , "KuIR2_3" , "KuIR8_1" , "KuIR8_2" , "KuIR8_3" , "KuIR24_1" , "KuIR24_2" , "KuIR24_3")

#data procession download functions
source("http://www.biostat.jhsph.edu/~kkammers/software/eupa/source.functions.r")

#Scale the protein intensity measurements from the mass spectrometry results
dat<-quantify.proteins(dat , cha)# median polished log2 transformed values, as a sacaled intensity for each protein 

#boxplot: intensities of all channels after data processing and normalization
par(mfrow = c(1,1) , font.lab=2 , cex.lab=1.2 , font.axis=2 , cex.axis=1.2)
boxplot(dat[,1:length(cha)], ylim=c(-10,10) , main="boxplot normalized intensities")

##heatmap of the scaled intensities
require(ComplexHeatmap)
#generate the matrix for the heat map
mat<-as.matrix(mat)#Transform the dat aset as matrix
mat<-mat[,7:33]
any(is.na(mat))#Check for NA values
mean(mat)# Mean
median(mat)# median
quantile(mat,c(0,0.25,0.5,0.75,1)) #Range of the data set

#Color function
col_fun<-colorRamp2(c(min(mat),mean(mat),max(mat)),c("royalblue","aliceblue","tomato"))
col_fun(seq( 0, 15))

#Heatmap plot
Heatmap(mat , name = "Fold Change vs Background" ,  col = col_fun , cluster_rows = T , 
        cluster_columns = FALSE , clustering_distance_rows = "spearman" , clustering_distance_columns = "spearman" ,
        clustering_method_rows = "single" , clustering_method_columns = "single" , column_names_centered = TRUE , 
        column_names_gp = gpar(fontsize=15) , row_names_gp = gpar(fontsize=3) ,
        row_split = GO_Anno ,  row_gap = unit(1, "mm"))

### T-test Ku80ID nonIR5 vs Cneg
#define design according to syntax of limma package
ct<-c("CN_1" , "CN_2" , "CN_3" , "CN_4" , "CN_5" , "CN_6")
trnonIR5<-c("KunonIR5_1" , "KunonIR5_2" , "KunonIR5_3")

#define design according to syntax of limma package
designnonIR5<-model.matrix(~factor(c(1,1,1,1,1,1,2,2,2)))

#Limma statistical analysis
colnames(designnonIR5)<-c("intercept","Diff")
res_nonIR5.eb<-eb.fit(dat[,c(ct,trnonIR5)] , designnonIR5)
head(res_nonIR5.eb)

####### Export tables
write.table(res_nonIR5.eb, file = "2DiffBind_KunonIR5_Limma.csv" , sep = ",")

#Volcanoplot
#trying to correct xylim problems
rx <- c(-1,1)*max(abs(res_nonIR5.eb$logFC))*1.1
ry <- c(0, 11)


par(mfrow = c(1,2) , font.lab=2 , cex.lab=1.2 , font.axis=2 , cex.axis=1.2)
par(las=1 , xaxs="i" , yaxs="i")

#Ordinary pValue
plot(res_nonIR5.eb$logFC, -log10(res_nonIR5.eb$p.ord), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="gray", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="gray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of ordinary p-Values")

#Moderate pValues
plot(res_nonIR5.eb$logFC, -log10(res_nonIR5.eb$p.mod), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="lightgray", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="lightgray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of moderated p-Values")

#Enhanced volcano
library(EnhancedVolcano)

EnhancedVolcano(toptable = res_nonIR5.eb, lab = rownames(res_nonIR5.eb) , x = "logFC" , y = "p.ord" , xlim = c(min(res_nonIR5.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res_nonIR5.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res_nonIR5.eb[["p.ord"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot KunonIR5" , subtitle = bquote(italic("Limma T-test")) , caption = paste0("total=" , nrow(res_nonIR5.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))
EnhancedVolcano(toptable = res_nonIR5.eb, lab = rownames(res_nonIR5.eb) , x = "logFC" , y = "p.mod" , xlim = c(min(res_nonIR5.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res_nonIR5.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res_nonIR5.eb[["p.mod"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot KunonIR5" , subtitle = bquote(italic("Limma moderated T-test")) , caption = paste0("total=" , nrow(res_nonIR5.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))

### T-Test Ku80UD nonIR24 vs Negative Controls
#define treatment and control groups for comparison
ct<-c("CN_1" , "CN_2" , "CN_3" , "CN_4" , "CN_5" , "CN_6")
trnonIR24<-c("KunonIR24_1" , "KunonIR24_2" , "KunonIR24_3")

#define design according to syntax of limma package
designnonIR24<-model.matrix(~factor(c(1,1,1,1,1,1,2,2,2)))

#Limma statistical analysis
colnames(designnonIR24)<-c("intercept","Diff")
res_nonIR24.eb<-eb.fit(dat[,c(ct,trnonIR24)] , designnonIR24)
head(res_nonIR24.eb)

####### Export tables
write.table(res_nonIR24.eb, file = "DiffBind_KunonIR24_Limma.csv" , sep = ",")

#Volcanoplot
#trying to correct xylim problems
rx <- c(-1,1)*max(abs(res_nonIR24.eb$logFC))*1.1
ry <- c(0, 11)


par(mfrow = c(1,2) , font.lab=2 , cex.lab=1.2 , font.axis=2 , cex.axis=1.2)
par(las=1 , xaxs="i" , yaxs="i")

#Ordinary pValue
plot(res_nonIR24.eb$logFC, -log10(res_nonIR24.eb$p.ord), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="gray", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="gray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of ordinary p-Values")

#Moderate pValues
plot(res_nonIR24.eb$logFC, -log10(res_nonIR24.eb$p.mod), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="lightgray", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="lightgray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of moderated p-Values")

#Enhanced volcano
library(EnhancedVolcano)

EnhancedVolcano(toptable = res_nonIR24.eb, lab = rownames(res_nonIR24.eb) , x = "logFC" , y = "p.ord" , xlim = c(min(res_nonIR24.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res_nonIR24.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res_nonIR24.eb[["p.ord"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot KunonIR24" , subtitle = bquote(italic("Limma T-test")) , caption = paste0("total=" , nrow(res_nonIR24.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))
EnhancedVolcano(toptable = res_nonIR24.eb, lab = rownames(res_nonIR24.eb) , x = "logFC" , y = "p.mod" , xlim = c(min(res_nonIR24.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res_nonIR24.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res_nonIR24.eb[["p.mod"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot KunonIR24" , subtitle = bquote(italic("Limma moderated T-test")) , caption = paste0("total=" , nrow(res_nonIR24.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))


### T-Test Ku80ID IR5 vs Negative controls
#define treatment and control groups for comparison
ct<-c("CN_1" , "CN_2" , "CN_3" , "CN_4" , "CN_5" , "CN_6")
trIR5<-c("KuIR5_1" , "KuIR5_2" , "KuIR5_3")

#define design according to syntax of limma package
designIR5<-model.matrix(~factor(c(1,1,1,1,1,1,2,2,2)))

#Limma statistical analysis
colnames(designIR5)<-c("intercept","Diff")
res_IR5.eb<-eb.fit(dat[,c(ct,trIR5)] , designIR5)
head(res_IR5.eb)


####### Export tables
write.table(res_IR5.eb, file = "DiffBind_KuIR5_Limma.csv" , sep = ",")

#Volcanoplot
#trying to correct xylim problems
rx <- c(-1,1)*max(abs(res_IR5.eb$logFC))*1.1
ry <- c(0, 11)


par(mfrow = c(1,2) , font.lab=2 , cex.lab=1.2 , font.axis=2 , cex.axis=1.2)
par(las=1 , xaxs="i" , yaxs="i")

#Ordinary pValues
plot(res_IR5.eb$logFC, -log10(res_IR5.eb$p.ord), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="gray", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="gray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of ordinary p-Values")

#Moderate pValues
plot(res_IR5.eb$logFC, -log10(res_IR5.eb$p.mod), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="lightgray", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="lightgray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of moderated p-Values")

#Enhanced volcano
library(EnhancedVolcano)

EnhancedVolcano(toptable = res_IR5.eb, lab = rownames(res_IR5.eb) , x = "logFC" , y = "p.ord" , xlim = c(min(res_IR5.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res_IR5.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res_IR5.eb[["p.ord"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot KuIR5" , subtitle = bquote(italic("Limma T-test")) , caption = paste0("total=" , nrow(res_IR5.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))
EnhancedVolcano(toptable = res_IR5.eb, lab = rownames(res_IR5.eb) , x = "logFC" , y = "p.mod" , xlim = c(min(res_IR5.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res_IR5.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res_IR5.eb[["p.mod"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot KuIR5" , subtitle = bquote(italic("Limma moderated T-test")) , caption = paste0("total=" , nrow(res_IR5.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))

### T-Test Ku80ID IR30 vs Negative controls
#define treatment and control groups for comparison
ct<-c("CN_1" , "CN_2" , "CN_3" , "CN_4" , "CN_5" , "CN_6")
trIR30<-c("KuIR30_1" , "KuIR30_2" , "KuIR30_3")

#define design according to syntax of limma package
designIR30<-model.matrix(~factor(c(1,1,1,1,1,1,2,2,2)))

#Limma statistical analysis
colnames(designIR30)<-c("intercept","Diff")
res_IR30.eb<-eb.fit(dat[,c(ct,trIR30)] , designIR30)
head(res_IR30.eb)

####### Export tables
write.table(res_IR30.eb, file = "DiffBind_KuIR30_Limma.csv" , sep = ",")

#Volcanoplot
#trying to correct xylim problems
rx <- c(-1,1)*max(abs(res_IR30.eb$logFC))*1.1
ry <- c(0, 11)


par(mfrow = c(1,2) , font.lab=2 , cex.lab=1.2 , font.axis=2 , cex.axis=1.2)
par(las=1 , xaxs="i" , yaxs="i")

#Ordinary pValues
plot(res_IR30.eb$logFC, -log10(res_IR30.eb$p.ord), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="gray", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="gray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of ordinary p-Values")

#Moderate pValues
plot(res_IR30.eb$logFC, -log10(res_IR30.eb$p.mod), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="lightgray", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="lightgray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of moderated p-Values")

#Enhanced volcano
library(EnhancedVolcano)

EnhancedVolcano(toptable = res_IR30.eb, lab = rownames(res_IR30.eb) , x = "logFC" , y = "p.ord" , xlim = c(min(res_IR30.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res_IR30.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res_IR30.eb[["p.ord"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot KuIR30" , subtitle = bquote(italic("Limma T-test")) , caption = paste0("total=" , nrow(res_IR30.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))
EnhancedVolcano(toptable = res_IR30.eb, lab = rownames(res_IR30.eb) , x = "logFC" , y = "p.mod" , xlim = c(min(res_IR30.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res_IR30.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res_IR30.eb[["p.mod"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot KuIR30" , subtitle = bquote(italic("Limma moderated T-test")) , caption = paste0("total=" , nrow(res_IR30.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))


### T-Test Ku80ID IR2h vs Negative Controls
#define treatment and control groups for comparison
ct<-c("CN_1" , "CN_2" , "CN_3" , "CN_4" , "CN_5" , "CN_6")
trIR2<-c("KuIR2_1" , "KuIR2_2" , "KuIR2_3")

#define design according to syntax of limma package
designIR2<-model.matrix(~factor(c(1,1,1,1,1,1,2,2,2)))

#Limma statistical analysis
colnames(designIR2)<-c("intercept","Diff")
res_IR2.eb<-eb.fit(dat[,c(ct,trIR2)] , designIR2)
head(res_IR2.eb)


####### Export tables
write.table(res_IR2.eb, file = "DiffBind_KuIR2_Limma.csv" , sep = ",")

#Volcanoplot
#trying to correct xylim problems
rx <- c(-1,1)*max(abs(res_IR2.eb$logFC))*1.1
ry <- c(0, 11)


par(mfrow = c(1,2) , font.lab=2 , cex.lab=1.2 , font.axis=2 , cex.axis=1.2)
par(las=1 , xaxs="i" , yaxs="i")

#Ordinary pValues
plot(res_IR2.eb$logFC, -log10(res_IR2.eb$p.ord), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="gray", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="gray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of ordinary p-Values")

#Moderate pValues
plot(res_IR2.eb$logFC, -log10(res_IR2.eb$p.mod), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="lightgray", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="lightgray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of moderated p-Values")

#Enhanced volcano
library(EnhancedVolcano)

EnhancedVolcano(toptable = res_IR2.eb, lab = rownames(res_IR2.eb) , x = "logFC" , y = "p.ord" , xlim = c(min(res_IR2.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res_IR2.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res_IR2.eb[["p.ord"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot KuIR2" , subtitle = bquote(italic("Limma T-test")) , caption = paste("total=" , nrow(res_IR2.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))
EnhancedVolcano(toptable = res_IR2.eb, lab = rownames(res_IR2.eb) , x = "logFC" , y = "p.mod" , xlim = c(min(res_IR2.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res_IR2.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res_IR2.eb[["p.mod"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot KuIR2" , subtitle = bquote(italic("Limma moderated T-test")) , caption = paste0("total=" , nrow(res_IR2.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))


### T-Test Ku80ID IR8h vs Negative Controls
#define treatment and control groups for comparison
ct<-c("CN_1" , "CN_2" , "CN_3" , "CN_4" , "CN_5" , "CN_6")
trIR8<-c("KuIR8_1" , "KuIR8_2" , "KuIR8_3")

#define design according to syntax of limma package
designIR8<-model.matrix(~factor(c(1,1,1,1,1,1,2,2,2)))

#Limma statistical analysis
colnames(designIR8)<-c("intercept","Diff")
res_IR8.eb<-eb.fit(dat[,c(ct,trIR8)] , designIR8)
head(res_IR8.eb)

####### Export tables
write.table(res_IR8.eb, file = "DiffBind_KuIR8_Limma.csv" , sep = ",")

#Volcanoplot
#trying to correct xylim problems
rx <- c(-1,1)*max(abs(res_IR2.eb$logFC))*1.1
ry <- c(0, 11)

par(mfrow = c(1,2) , font.lab=2 , cex.lab=1.2 , font.axis=2 , cex.axis=1.2)
par(las=1 , xaxs="i" , yaxs="i")

#Ordinary pValue
plot(res_IR8.eb$logFC, -log10(res_IR8.eb$p.ord), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="gray", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="gray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of ordinary p-Values")

#Moderate pValues
plot(res_IR8.eb$logFC, -log10(res_IR8.eb$p.mod), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="lightgray", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="lightgray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of moderated p-Values")

#Enhanced volcano
library(EnhancedVolcano)

EnhancedVolcano(toptable = res_IR8.eb, lab = rownames(res_IR8.eb) , x = "logFC" , y = "p.ord" , xlim = c(min(res_IR8.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res_IR8.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res_IR8.eb[["p.ord"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot KuIR8" , subtitle = bquote(italic("Limma T-test")) , caption = paste("total=" , nrow(res_IR8.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))
EnhancedVolcano(toptable = res_IR8.eb, lab = rownames(res_IR8.eb) , x = "logFC" , y = "p.mod" , xlim = c(min(res_IR8.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res_IR8.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res_IR8.eb[["p.mod"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot KuIR8" , subtitle = bquote(italic("Limma moderated T-test")) , caption = paste0("total=" , nrow(res_IR8.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))


### T-Test Ku80ID IR24h vs Negative controls
#define treatment and control groups for comparison
ct<-c("CN_1" , "CN_2" , "CN_3" , "CN_4" , "CN_5" , "CN_6")
trIR24<-c("KuIR24_1" , "KuIR24_2" , "KuIR24_3")

#define design according to syntax of limma package
designIR24<-model.matrix(~factor(c(1,1,1,1,1,1,2,2,2)))

#Limma statistical analysis
colnames(designIR24)<-c("intercept","Diff")
res_IR24.eb<-eb.fit(dat[,c(ct,trIR24)] , designIR24)
head(res_IR24.eb)

####### Export tables
write.table(res_IR24.eb, file = "DiffBind_KuIR24_Limma.csv" , sep = ",")

#Volcanoplot
#trying to correct xylim problems
rx <- c(-1,1)*max(abs(res_IR24.eb$logFC))*1.1
ry <- c(0, 11)

par(mfrow = c(1,2) , font.lab=2 , cex.lab=1.2 , font.axis=2 , cex.axis=1.2)
par(las=1 , xaxs="i" , yaxs="i")

#Ordinary pValue
plot(res_IR24.eb$logFC, -log10(res_IR24.eb$p.ord), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="gray", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="gray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of ordinary p-Values")

#Moderate pValues
plot(res_IR24.eb$logFC, -log10(res_IR24.eb$p.mod), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="lightgray", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="lightgray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of moderated p-Values")

#Enhanced volcano
library(EnhancedVolcano)

EnhancedVolcano(toptable = res_IR24.eb, lab = rownames(res_IR24.eb) , x = "logFC" , y = "p.ord" , xlim = c(min(res_IR24.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res_IR24.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res_IR24.eb[["p.ord"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot KuIR24" , subtitle = bquote(italic("Limma T-test")) , caption = paste("total=" , nrow(res_IR24.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))
EnhancedVolcano(toptable = res_IR24.eb, lab = rownames(res_IR24.eb) , x = "logFC" , y = "p.mod" , xlim = c(min(res_IR24.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res_IR24.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res_IR24.eb[["p.mod"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot KuIR24" , subtitle = bquote(italic("Limma moderated T-test")) , caption = paste0("total=" , nrow(res_IR24.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))

################
###### Mre11ID analysis
dat = read.csv("./Mre11ID_RawData.csv",header=TRUE)#Open the Ku80ID raw data table, and coherce the column names as it is
#Review the structure of the data
dim(dat)
str(dat)

#define vector headlines
cha<- c("CN_1" , "CN_2" , "CN_3" , "CN_4" , "CN_5" , "CN_6" , "MrenonIR5_1" , "MrenonIR5_2" , "MrenonIR5_3" , "MrenonIR24_1" , "MrenonIR24_2" , "MrenonIR24_3" , "MreIR5_1" , "MreIR5_2" , "MreIR5_3" , "MreIR30_1" , "MreIR30_2" , "MreIR30_3" , "MreIR2_1" , "MreIR2_2" , "MreIR2_3" , "MreIR8_1" , "MreIR8_2" , "MreIR8_3" , "MreIR24_1" , "MreIR24_2" , "MreIR24_3")

#data procession download functions
source("http://www.biostat.jhsph.edu/~kkammers/software/eupa/source.functions.r")

#Readpeptides
#dat<- read.peptides(dat , cha)#Since we do not imput the peptide reads, but already the protein groups this one is not used
#scale the protein intensity measurements from the mass spectrometry results
dat<-quantify.proteins(dat , cha)# median polished log2 transformed values, as a sacaled intensity for each protein 

#boxplot: intensities of all channels after data processing and normalization
par(mfrow = c(1,1) , font.lab=2 , cex.lab=1.2 , font.axis=2 , cex.axis=1.2)
boxplot(dat[,1:length(cha)], ylim=c(-10,10) , main="boxplot normalized intensities")

##heatmap of the scaled intensities
require(ComplexHeatmap)
#generate the matrix for the heat map
mat<-as.matrix(mat)#Transform the dataset as matrix
mat<-mat[,7:33]
any(is.na(mat))#Check for NA values
mean(mat)# Mean
median(mat)# median
quantile(mat,c(0,0.25,0.5,0.75,1)) #Range of the dataset

#Color function
col_fun<-colorRamp2(c(min(mat),mean(mat),max(mat)),c("royalblue","aliceblue","tomato"))
col_fun(seq( 0, 15))

#Heatmap plot
Heatmap(mat , name = "Fold Change vs Background" ,  col = col_fun , cluster_rows = T , 
        cluster_columns = FALSE , clustering_distance_rows = "spearman" , clustering_distance_columns = "spearman" ,
        clustering_method_rows = "single" , clustering_method_columns = "single" , column_names_centered = TRUE , 
        column_names_gp = gpar(fontsize=15) , row_names_gp = gpar(fontsize=3) ,
        row_split = GO_Anno ,  row_gap = unit(1, "mm"))

### T-Test Mre11ID nonIR5 vs Negative controls
#define treatment and control groups for comparison
ct<-c("CN_1" , "CN_2" , "CN_3" , "CN_4" , "CN_5" , "CN_6")
trnonIR5<-c("MrenonIR5_1" , "MrenonIR5_2" , "MrenonIR5_3")

#define design according to syntax of limma package
designnonIR5<-model.matrix(~factor(c(1,1,1,1,1,1,2,2,2)))

#Limma statistical analysis
colnames(designnonIR5)<-c("intercept","Diff")
res_nonIR5.eb<-eb.fit(dat[,c(ct,trnonIR5)] , designnonIR5)
head(res_nonIR5.eb)

####### Export tables
write.table(res_nonIR5.eb, file = "DiffBind_MrenonIR5_Limma.csv" , sep = ",")

#Volcanoplot
#trying to correct xylim problems
rx <- c(-1,1)*max(abs(res_nonIR5.eb$logFC))*1.1
ry <- c(0, 11)

par(mfrow = c(1,2) , font.lab=2 , cex.lab=1.2 , font.axis=2 , cex.axis=1.2)
par(las=1 , xaxs="i" , yaxs="i")

#Ordinary pValue
plot(res_nonIR5.eb$logFC, -log10(res_nonIR5.eb$p.ord), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="gray", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="gray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of ordinary p-Values")

#Moderate pValues
plot(res_nonIR5.eb$logFC, -log10(res_nonIR5.eb$p.mod), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="lightgray", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="lightgray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of moderated p-Values")

#Enhanced volcano
library(EnhancedVolcano)

EnhancedVolcano(toptable = res_nonIR5.eb, lab = rownames(res_nonIR5.eb) , x = "logFC" , y = "p.ord" , xlim = c(min(res_nonIR5.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res_nonIR5.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res_nonIR5.eb[["p.ord"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot MrenonIR5" , subtitle = bquote(italic("Limma T-test")) , caption = paste0("total=" , nrow(res_nonIR5.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))
EnhancedVolcano(toptable = res_nonIR5.eb, lab = rownames(res_nonIR5.eb) , x = "logFC" , y = "p.mod" , xlim = c(min(res_nonIR5.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res_nonIR5.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res_nonIR5.eb[["p.mod"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot MrenonIR5" , subtitle = bquote(italic("Limma moderated T-test")) , caption = paste0("total=" , nrow(res_nonIR5.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))


### T-test Mre11ID nonIR24 vs negative controls
#define treatment and control groups for comparison
ct<-c("CN_1" , "CN_2" , "CN_3" , "CN_4" , "CN_5" , "CN_6")
trnonIR24<-c("MrenonIR24_1" , "MrenonIR24_2" , "MrenonIR24_3")

#define design according to syntax of limma package
designnonIR24<-model.matrix(~factor(c(1,1,1,1,1,1,2,2,2)))

#Limma statistical analysis
colnames(designnonIR24)<-c("intercept","Diff")
res_nonIR24.eb<-eb.fit(dat[,c(ct,trnonIR24)] , designnonIR24)
head(res_nonIR24.eb)

####### Export tables
write.table(res_nonIR24.eb, file = "DiffBind_MrenonIR24_Limma.csv" , sep = ",")

#Volcanoplot
#trying to correct xylim problems
rx <- c(-1,1)*max(abs(res_nonIR24.eb$logFC))*1.1
ry <- c(0, 16)

par(mfrow = c(1,2) , font.lab=2 , cex.lab=1.2 , font.axis=2 , cex.axis=1.2)
par(las=1 , xaxs="i" , yaxs="i")

#Ordinary pValue
plot(res_nonIR24.eb$logFC, -log10(res_nonIR24.eb$p.ord), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="gray", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="gray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of ordinary p-Values")

#Moderate pValues
plot(res_nonIR24.eb$logFC, -log10(res_nonIR24.eb$p.mod), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="lightgray", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="lightgray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of moderated p-Values")

#Enhanced volcano
library(EnhancedVolcano)

EnhancedVolcano(toptable = res_nonIR24.eb, lab = rownames(res_nonIR24.eb) , x = "logFC" , y = "p.ord" , xlim = c(min(res_nonIR24.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res_nonIR24.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res_nonIR24.eb[["p.ord"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot MrenonIR24" , subtitle = bquote(italic("Limma T-test")) , caption = paste0("total=" , nrow(res_nonIR24.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))
EnhancedVolcano(toptable = res_nonIR24.eb, lab = rownames(res_nonIR24.eb) , x = "logFC" , y = "p.mod" , xlim = c(min(res_nonIR24.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res_nonIR24.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res_nonIR24.eb[["p.mod"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot MrenonIR24" , subtitle = bquote(italic("Limma moderated T-test")) , caption = paste0("total=" , nrow(res_nonIR24.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))


### T-Test Mre11ID IR5 vs Negative controls 
#define treatment and control groups for comparison
ct<-c("CN_1" , "CN_2" , "CN_3" , "CN_4" , "CN_5" , "CN_6")
trIR5<-c("MreIR5_1" , "MreIR5_2" , "MreIR5_3")

#define design according to syntax of limma package
designIR5<-model.matrix(~factor(c(1,1,1,1,1,1,2,2,2)))

#Limma statistical analysis
colnames(designIR5)<-c("intercept","Diff")
res_IR5.eb<-eb.fit(dat[,c(ct,trIR5)] , designIR5)
head(res_IR5.eb)

####### Export tables
write.table(res_IR5.eb, file = "DiffBind_MreIR5_Limma.csv" , sep = ",")

#Volcanoplot
#trying to correct xylim problems
rx <- c(-1,1)*max(abs(res_IR5.eb$logFC))*1.1
ry <- c(0, 11)

par(mfrow = c(1,2) , font.lab=2 , cex.lab=1.2 , font.axis=2 , cex.axis=1.2)
par(las=1 , xaxs="i" , yaxs="i")

#Ordinary pValue
plot(res_IR5.eb$logFC, -log10(res_IR5.eb$p.ord), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="gray", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="gray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of ordinary p-Values")

#Moderate pValues
plot(res_IR5.eb$logFC, -log10(res_IR5.eb$p.mod), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="lightgray", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="lightgray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of moderated p-Values")

#Enhanced volcano
library(EnhancedVolcano)

EnhancedVolcano(toptable = res_IR5.eb, lab = rownames(res_IR5.eb) , x = "logFC" , y = "p.ord" , xlim = c(min(res_IR5.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res_IR5.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res_IR5.eb[["p.ord"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot MreIR5" , subtitle = bquote(italic("Limma T-test")) , caption = paste0("total=" , nrow(res_IR5.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))
EnhancedVolcano(toptable = res_IR5.eb, lab = rownames(res_IR5.eb) , x = "logFC" , y = "p.mod" , xlim = c(min(res_IR5.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res_IR5.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res_IR5.eb[["p.mod"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot MreIR5" , subtitle = bquote(italic("Limma moderated T-test")) , caption = paste0("total=" , nrow(res_IR5.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))


### T-test Mre11ID IR30 vs Negative control
#define treatment and control groups for comparison
ct<-c("CN_1" , "CN_2" , "CN_3" , "CN_4" , "CN_5" , "CN_6")
trIR30<-c("MreIR30_1" , "MreIR30_2" , "MreIR30_3")

#define design according to syntax of limma package
designIR30<-model.matrix(~factor(c(1,1,1,1,1,1,2,2,2)))

#Limma statistical analysis
colnames(designIR30)<-c("intercept","Diff")
res_IR30.eb<-eb.fit(dat[,c(ct,trIR30)] , designIR30)
head(res_IR30.eb)

####### Export tables
write.table(res_IR30.eb, file = "DiffBind_MreIR30_Limma.csv" , sep = ",")

#Volcanoplot
#trying to correct xylim problems
rx <- c(-1,1)*max(abs(res_IR30.eb$logFC))*1.1
ry <- c(0, 14)

par(mfrow = c(1,2) , font.lab=2 , cex.lab=1.2 , font.axis=2 , cex.axis=1.2)
par(las=1 , xaxs="i" , yaxs="i")

#Ordinary pValue
plot(res_IR30.eb$logFC, -log10(res_IR30.eb$p.ord), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="gray", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="gray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of ordinary p-Values")

#Moderate pValues
plot(res_IR30.eb$logFC, -log10(res_IR30.eb$p.mod), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="lightgray", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="lightgray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of moderated p-Values")

#Enhanced volcano
library(EnhancedVolcano)

EnhancedVolcano(toptable = res_IR30.eb, lab = rownames(res_IR30.eb) , x = "logFC" , y = "p.ord" , xlim = c(min(res_IR30.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res_IR30.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res_IR30.eb[["p.ord"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot MreIR30" , subtitle = bquote(italic("Limma T-test")) , caption = paste0("total=" , nrow(res_IR30.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))
EnhancedVolcano(toptable = res_IR30.eb, lab = rownames(res_IR30.eb) , x = "logFC" , y = "p.mod" , xlim = c(min(res_IR30.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res_IR30.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res_IR30.eb[["p.mod"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot MreIR30" , subtitle = bquote(italic("Limma moderated T-test")) , caption = paste0("total=" , nrow(res_IR30.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))

### T-test Mre11ID IR2h vs Negative control
#define treatment and control groups for comparison
ct<-c("CN_1" , "CN_2" , "CN_3" , "CN_4" , "CN_5" , "CN_6")
trIR2<-c("MreIR2_1" , "MreIR2_2" , "MreIR2_3")

#define design according to syntax of limma package
designIR2<-model.matrix(~factor(c(1,1,1,1,1,1,2,2,2)))
designIR2

#Limma statistical analysis
colnames(designIR2)<-c("intercept","Diff")
res_IR2.eb<-eb.fit(dat[,c(ct,trIR2)] , designIR2)
head(res_IR2.eb)

row.names(res_IR2.eb)<-row.names(dat)

####### Export tables
write.table(res_IR2.eb, file = "DiffBind_MreIR2_Limma.csv" , sep = ",")

#Volcanoplot
#trying to correct xylim problems
rx <- c(-1,1)*max(abs(res_IR2.eb$logFC))*1.1
ry <- c(0, 11)


par(mfrow = c(1,2) , font.lab=2 , cex.lab=1.2 , font.axis=2 , cex.axis=1.2)
par(las=1 , xaxs="i" , yaxs="i")

#ordinary pvalues
plot(res_IR2.eb$logFC, -log10(res_IR2.eb$p.ord), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="gray", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="gray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of ordinary p-Values")

#Moderate pValues
plot(res_IR2.eb$logFC, -log10(res_IR2.eb$p.mod), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="lightgray", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="lightgray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of moderated p-Values")

#Enhanced volcano
library(EnhancedVolcano)

EnhancedVolcano(toptable = res_IR2.eb, lab = rownames(res_IR2.eb) , x = "logFC" , y = "p.ord" , xlim = c(min(res_IR2.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res_IR2.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res_IR2.eb[["p.ord"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot MreIR2" , subtitle = bquote(italic("Limma T-test")) , caption = paste("total=" , nrow(res_IR2.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))
EnhancedVolcano(toptable = res_IR2.eb, lab = rownames(res_IR2.eb) , x = "logFC" , y = "p.mod" , xlim = c(min(res_IR2.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res_IR2.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res_IR2.eb[["p.mod"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot MreIR2" , subtitle = bquote(italic("Limma moderated T-test")) , caption = paste0("total=" , nrow(res_IR2.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))


### T-Test Mre11ID IR8h vs Negative control
#define treatment and control groups for comparison
ct<-c("CN_1" , "CN_2" , "CN_3" , "CN_4" , "CN_5" , "CN_6")
trIR8<-c("MreIR8_1" , "MreIR8_2" , "MreIR8_3")

#define design according to syntax of limma package
designIR8<-model.matrix(~factor(c(1,1,1,1,1,1,2,2,2)))

#Limma statistical analysis
colnames(designIR8)<-c("intercept","Diff")
res_IR8.eb<-eb.fit(dat[,c(ct,trIR8)] , designIR8)
head(res_IR8.eb)


####### Export tables
write.table(res_IR8.eb, file = "DiffBind_MreIR8_Limma.csv" , sep = ",")

#Volcanoplot
#trying to correct xylim problems
rx <- c(-1,1)*max(abs(res_IR2.eb$logFC))*1.1
ry <- c(0, 11)


par(mfrow = c(1,2) , font.lab=2 , cex.lab=1.2 , font.axis=2 , cex.axis=1.2)
par(las=1 , xaxs="i" , yaxs="i")

#Ordinary pvalues
plot(res_IR8.eb$logFC, -log10(res_IR8.eb$p.ord), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="gray", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="gray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of ordinary p-Values")

#Moderate pValues
plot(res_IR8.eb$logFC, -log10(res_IR8.eb$p.mod), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="lightgray", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="lightgray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of moderated p-Values")

#Enhanced volcano
library(EnhancedVolcano)

EnhancedVolcano(toptable = res_IR8.eb, lab = rownames(res_IR8.eb) , x = "logFC" , y = "p.ord" , xlim = c(min(res_IR8.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res_IR8.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res_IR8.eb[["p.ord"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot MreIR8" , subtitle = bquote(italic("Limma T-test")) , caption = paste("total=" , nrow(res_IR8.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))
EnhancedVolcano(toptable = res_IR8.eb, lab = rownames(res_IR8.eb) , x = "logFC" , y = "p.mod" , xlim = c(min(res_IR8.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res_IR8.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res_IR8.eb[["p.mod"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot MreIR8" , subtitle = bquote(italic("Limma moderated T-test")) , caption = paste0("total=" , nrow(res_IR8.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))



### T-test Mre11ID IR24h vs Negative controls
#define treatment and control groups for comparison
ct<-c("CN_1" , "CN_2" , "CN_3" , "CN_4" , "CN_5" , "CN_6")
trIR24<-c("MreIR24_1" , "MreIR24_2" , "MreIR24_3")

#define design according to syntax of limma package
designIR24<-model.matrix(~factor(c(1,1,1,1,1,1,2,2,2)))

#Limma statistical analysis
colnames(designIR24)<-c("intercept","Diff")
res_IR24.eb<-eb.fit(dat1[,c(ct,trIR24)] , designIR24)
head(res_IR24.eb)

####### Export tables
write.table(res_IR24.eb, file = "DiffBind_MreIR24_Limma.csv" , sep = ",")

#Volcanoplot
#trying to correct xylim problems
rx <- c(-1,1)*max(abs(res_IR24.eb$logFC))*1.1
ry <- c(0, 17)


par(mfrow = c(1,2) , font.lab=2 , cex.lab=1.2 , font.axis=2 , cex.axis=1.2)
par(las=1 , xaxs="i" , yaxs="i")

#Ordinary pValues
plot(res_IR24.eb$logFC, -log10(res_IR24.eb$p.ord), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="gray", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="gray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of ordinary p-Values")

#Moderate pValues
plot(res_IR24.eb$logFC, -log10(res_IR24.eb$p.mod), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="lightgray", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="lightgray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of moderated p-Values")

#Enhanced volcano
library(EnhancedVolcano)

EnhancedVolcano(toptable = res_IR24.eb, lab = rownames(res_IR24.eb) , x = "logFC" , y = "p.ord" , xlim = c(min(res_IR24.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res_IR24.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res_IR24.eb[["p.ord"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot MreIR24" , subtitle = bquote(italic("Limma T-test")) , caption = paste("total=" , nrow(res_IR24.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))
EnhancedVolcano(toptable = res_IR24.eb, lab = rownames(res_IR24.eb) , x = "logFC" , y = "p.mod" , xlim = c(min(res_IR24.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res_IR24.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res_IR24.eb[["p.mod"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot MreIR24" , subtitle = bquote(italic("Limma moderated T-test")) , caption = paste0("total=" , nrow(res_IR24.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))


##### Heat Map with all the Fold Change Values
# The LIMMA T-Test result tables were manually filtered in excel to select only those interactors with a Fold Change above 1.5 and a moderate pValue below 0.05.
# The resulting protein lists with their Fold Change values were merged in excel to generate a single data frame containing the fold change values for the three DSB-sensors, and the 7 analyzed time points.
require(ComplexHeatmap)
#Big HeatMaps
setwd("C:/../Limma_t_Test_Analysis")#select working Directory
DF<- read.csv("./SensorID_FoldChange_Matrix.csv",row.names = "NAME")#Read Protein Annotations and foldChange values
colAnno<-read.csv("./ColAnnoDF.csv",row.names = 1)#Read sample annotations metadata
mat<-subset(DF,select = c(12:32),)#Create a matrix only with Fold Change Values
#review the structure of the resulting
str(mat)
#Z-Score scaling of the SensorID matrixes
matS6<-mat[,1:7]#Subset the Sirt6ID columns 
matS6<- t(scale(t(matS6),center = TRUE, scale=TRUE))# Z-score scaling of the fold change Sirt6 matrix in a row-wise method
matS6<-ifelse(is.na(matS6),0,matS6)#convert the NA values generating from scaling the 0s values to 0s
matKu<-mat[,8:14]#Subset the Sirt6ID columns 
matKu<- t(scale(t(matKu),center = TRUE, scale=TRUE))# Z-score scaling of the fold change Sirt6 matrix in a row-wise method
matKu<-ifelse(is.na(matKu),0,matKu)#convert the NA values generating from scaling the 0s values to 0s
matMre<-mat[,15:21]#Subset the Sirt6ID columns 
matMre<- t(scale(t(matMre),center = TRUE, scale=TRUE))# Z-score scaling of the fold change Sirt6 matrix in a row-wise method
matMre<-ifelse(is.na(matMre),0,matMre)#convert the NA values generating from scaling the 0s values to 0s
#Construct the entire matrix again
mat<-cbind(matS6,matKu,matMre)#Merging the three matrix to have the three SensorID matrixes in one singe matrix
mat<-as.matrix(mat)#Transform the dataset as matrix
#Review the structure of the matrix
any(is.na(mat))#Check for NA values
dim(mat) 
range(mat) 
mean(mat)
median(mat)
quantile(mat,c(0,0.25,0.5,0.75,1))

#Generate the column annotations for the heat map
colAnno<-read.csv("./ColAnnoDF.csv",row.names = 1)#Read sample annotations metadata
str(colAnno)
Col_Anno<-HeatmapAnnotation(df=colAnno, col = list(Treatment=c("IR"="red","nonIR"="blue"),
                                                   Sensor=c("Ku80"="seagreen1","Mre11"="steelblue1","Sirt6"="tomato")))

#Generate Annotation vectors for proteins to cluster the rows of the heat map based on their core functions or roles in DNA repair
DF$Annotation<-as.matrix(subset(DF,select = DNA_Repair))
DNArepair<-as.matrix(subset(DF,select = DNA_Repair))
ha_anno<-rowAnnotation(Annotation=as.factor(GO_Anno))#,col=Col_Anno)#annotation_legend_param=)#,labels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22"))
ha_DNArep<-rowAnnotation(DNA_Repair=as.factor(DNArepair))

#Color function for HeatMaps
require(circlize)
col_fun<-colorRamp2(c(min(mat),mean(mat),max(mat)),c("royalblue","aliceblue","tomato"))
col_fun(seq( 0, 15))  
#Heatmap script
Heatmap(mat , name = "Fold Change vs Background" ,  col = col_fun , cluster_rows = T , 
        cluster_columns = FALSE , clustering_distance_rows = "spearman" , clustering_distance_columns = "spearman" ,
        clustering_method_rows = "ward.D2" , clustering_method_columns = "single" , column_names_centered = TRUE , 
        column_names_gp = gpar(fontsize=15) , row_names_gp = gpar(fontsize=3) ,top_annotation = Col_Anno, show_column_names = F,
        row_split = GO_Anno ,  row_gap = unit(1, "mm") , 
        left_annotation = ha_anno)

