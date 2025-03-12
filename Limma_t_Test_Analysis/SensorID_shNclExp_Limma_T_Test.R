#Load the require packages
require(limma)
require(qvalue)
require(rlang)
require(devtools)
require(ggplot2)
require(ggrepel)
require(tidyverse)
require(dplyr)
require(tibble)

#Upload data
setwd("C:/../Limma_t_Test_Analysis")

#Analysis of Sirt6ID shScramble and shNcl StrepIP mass spectrometry results

#Load formated dataset as required for Limma analysis
dat<-read.csv("./SensorID_shNclExp_RawData_All.csv")
#Check datasets
dim(dat)
str(dat)
#Remove N/A from dataset
a<-is.na(dat)
unique(a)
dat<-na.omit(dat)
dim(dat)
#check for duplicate gene names
length(unique(dat$Protein.Group.Accessions))

#define vector headlines
cha<- c("CN_1" , "CN_2" , "CN_3" , "CN_4" , "CN_5" , "CN_6" ,"CN_7" , "CN_8" , "CN_9" ,
        "S6_shScr_1" , "S6_shScr_2" , "S6_shScr_3" , "S6_shNcl_1" , "S6_shNcl_2" , "S6_shNcl_3", 
        "Ku_shScr_1" , "Ku_shScr_2" , "Ku_shScr_3" , "Ku_shNcl_1" , "Ku_shNcl_2" , "Ku_shNcl_3",
        "Mre_shScr_1" , "Mre_shScr_2" , "Mre_shScr_3" , "Mre_shNcl_1" , "Mre_shNcl_2" , "Mre_shNcl_3")

#Load the packages for the LIMMA multiple T-Test
#data procession download functions
source("http://www.biostat.jhsph.edu/~kkammers/software/eupa/source.functions.r")

#Readpeptides (This is not necessary because we are already using the averaged protein intensities)
#dat<- read.peptides(dat , cha)

#Channel values for each protein are median polished log2 transformed values from the median value of all peptides belonging to each protein.
dat<-quantify.proteins(dat , cha)
dim(dat)

#boxplot: intensities of all channels after data processing and normalization
par(mfrow = c(1,1) , font.lab=1 , cex.lab=0.5 , font.axis=1 , cex.axis=0.8)
boxplot(dat[,1:length(cha)], ylim=c(-20,20) , main="boxplot normalized intensities")

### Sirt6ID shNcl Experiment LIMMA T-Test analysis

## Sirt6ID T-test of shScramble vs Negative Controls
#First lets subset the samples:
S6_Scr<-dat[,1:12]
S6_Scr
#define treatment and control groups for comparison
ct<-c("CN_1" , "CN_2" , "CN_3" , "CN_4" , "CN_5" , "CN_6" ,"CN_7" , "CN_8" , "CN_9" )
tr<-c("S6_shScr_1" , "S6_shScr_2" , "S6_shScr_3")

#define design according to syntax of limma package
design<-model.matrix(~factor(c(1,1,1,1,1,1,1,1,1,2,2,2)))
design

#Limma statistical analysis
colnames(design)<-c("intercept","Diff")
res.eb<-eb.fit(S6_Scr , design)
head(res.eb)
res.eb$NAME<-row.names(res.eb)
str(res.eb)

#Export dataset res.eb to an excel file

require(writexl)
write_xlsx(res.eb , "./Diff_S6_shScrVScneg.xlsx")

#Volcano plots for ordinary and moderate pValues
#Ordinary pValues
rx <- c(-1,1)*max(abs(res.eb$logFC))*1.1
ry <- c(0, 8)

par(mfrow = c(1,2) , font.lab=2 , cex.lab=1.2 , font.axis=2 , cex.axis=1.2)
par(las=1 , xaxs="i" , yaxs="i")

#Ordinary pValues
plot(res.eb$logFC, -log10(res.eb$p.ord), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="blue", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="blue", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of ordinary p-Values")

#Moderate pValues
plot(res.eb$logFC, -log10(res.eb$p.mod), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="lightgray", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="lightgray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of moderated p-Values")

#enhanced volcano plots
require("EnhancedVolcano")

EnhancedVolcano(toptable = res.eb, lab = rownames(res.eb) , x = "logFC" , y = "p.ord" , xlim = c(min(res.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res.eb[["p.ord"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot Sirt6ID shScramVsCneg" , subtitle = bquote(italic("Limma T-test")) , caption = paste0("total=" , nrow(res.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))
EnhancedVolcano(toptable = res.eb, lab = rownames(res.eb) , x = "logFC" , y = "p.mod" , xlim = c(min(res.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res.eb[["p.mod"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot Sirt6ID shScramVsCneg" , subtitle = bquote(italic("Limma moderated T-test")) , caption = paste0("total=" , nrow(res.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))

### Sirt6ID T-test of shNcl vs Negative Controls
#Subset the dataset
S6_Ncl<-dat[,c(1:9,13:15)]
S6_Ncl
#define treatment and control groups for comparison
ct<-c("CN_1" , "CN_2" , "CN_3" , "CN_4" , "CN_5" , "CN_6" ,"CN_7" , "CN_8" , "CN_9" )
tr<-c("S6_shNcl_1" , "S6_shNcl_2" , "S6_shNcl_3")

#define design according to syntax of limma package
design<-model.matrix(~factor(c(1,1,1,1,1,1,1,1,1,2,2,2)))
design

#Limma statistical analysis
colnames(design)<-c("intercept","Diff")
res.eb<-eb.fit(S6_Ncl , design)
head(res.eb)
res.eb$NAME<-row.names(res.eb)
str(res.eb)

#Export dataset res.eb to an excel file
write_xlsx(res.eb , "./Diff_S6_shNclVScneg.xlsx")

#Volcano plots for ordinary and moderate pValues
#Ordinary pValues
rx <- c(-1,1)*max(abs(res.eb$logFC))*1.1
ry <- c(0, 8)

par(mfrow = c(1,2) , font.lab=2 , cex.lab=1.2 , font.axis=2 , cex.axis=1.2)
par(las=1 , xaxs="i" , yaxs="i")

#Ordinary pValue
plot(res.eb$logFC, -log10(res.eb$p.ord), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="blue", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="blue", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of ordinary p-Values")

#Moderate pValues
plot(res.eb$logFC, -log10(res.eb$p.mod), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="lightgray", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="lightgray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of moderated p-Values")

#enhanced volcano plots

require(EnhancedVolcano)
EnhancedVolcano(toptable = res.eb, lab = rownames(res.eb) , x = "logFC" , y = "p.ord" , xlim = c(min(res.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res.eb[["p.ord"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot Sirt6ID shNclVsCneg" , subtitle = bquote(italic("Limma T-test")) , caption = paste0("total=" , nrow(res.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))
EnhancedVolcano(toptable = res.eb, lab = rownames(res.eb) , x = "logFC" , y = "p.mod" , xlim = c(min(res.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res.eb[["p.mod"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot Sirt6ID shNclVsCneg" , subtitle = bquote(italic("Limma moderated T-test")) , caption = paste0("total=" , nrow(res.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))

##############################
### Sirt6ID LIMMA T-test of shNcl vs shScramble
#Subset the data
S6_NclvsScr<-dat[,c(10:15)]
S6_NclvsScr
#define treatment and control groups for comparison
ct<-c("S6_shScr_1" , "S6_shScr_2" , "S6_shScr_3")
tr<-c("S6_shNcl_1" , "S6_shNcl_2" , "S6_shNcl_3")

#define design according to syntax of limma package
design<-model.matrix(~factor(c(1,1,1,2,2,2)))
design

#Limma statistical analysis
colnames(design)<-c("intercept","Diff")
res.eb<-eb.fit(S6_NclvsScr , design)
head(res.eb)
res.eb$NAME<-row.names(res.eb)
str(res.eb)

#Export dataset res.eb to an excel file
write_xlsx(res.eb , "./Diff_S6_shNclVSshScr.xlsx")

#Volcano plots for ordinary and moderate pValues
#Ordinary pValues
rx <- c(-1,1)*max(abs(res.eb$logFC))*1.1
ry <- c(0, 8)

par(mfrow = c(1,2) , font.lab=2 , cex.lab=1.2 , font.axis=2 , cex.axis=1.2)
par(las=1 , xaxs="i" , yaxs="i")

#Ordinary pValue
plot(res.eb$logFC, -log10(res.eb$p.ord), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="blue", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="blue", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of ordinary p-Values")

#Moderate pValues
plot(res.eb$logFC, -log10(res.eb$p.mod), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="lightgray", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="lightgray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of moderated p-Values")

#enhanced volcano plots
require(EnhancedVolcano)
EnhancedVolcano(toptable = res.eb, lab = rownames(res.eb) , x = "logFC" , y = "p.ord" , xlim = c(min(res.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res.eb[["p.ord"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot Sirt6ID shNcl_Vs_shScram" , subtitle = bquote(italic("Limma T-test")) , caption = paste0("total=" , nrow(res.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))
EnhancedVolcano(toptable = res.eb, lab = rownames(res.eb) , x = "logFC" , y = "p.mod" , xlim = c(min(res.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res.eb[["p.mod"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot Sirt6ID shNcl_Vs_shScram" , subtitle = bquote(italic("Limma moderated T-test")) , caption = paste0("total=" , nrow(res.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))


################################################
### Ku80ID shNcl Experiment LIMMA multiple T-Test analysis

### Ku80ID multiple T-Test analysis shScramble vs Negative Controls
#Subset the big dataset
Ku_Scr<-dat[,c(1:9,16:18)]
Ku_Scr
#T-test of shScramble vs Negative Controls
#define treatment and control groups for comparison
ct<-c("CN_1" , "CN_2" , "CN_3" , "CN_4" , "CN_5" , "CN_6" ,"CN_7" , "CN_8" , "CN_9"  )
tr<-c("Ku_shScr_1" , "Ku_shScr_2" , "Ku_shScr_3")

#define design according to syntax of limma package
design<-model.matrix(~factor(c(1,1,1,1,1,1,1,1,1,2,2,2)))
design

#Limma statistical analysis
colnames(design)<-c("intercept","Diff")
res.eb<-eb.fit(Ku_Scr, design)
head(res.eb)
res.eb$NAME<-row.names(res.eb)
str(res.eb)

#Export dataset res.eb to an excel file

write_xlsx(res.eb , "./Diff_Ku_shScrVScneg.xlsx")

#Volcano plots for ordinary and moderate pValues
#Ordinary pValues
rx <- c(-1,1)*max(abs(res.eb$logFC))*1.1
ry <- c(0, 8)

par(mfrow = c(1,2) , font.lab=2 , cex.lab=1.2 , font.axis=2 , cex.axis=1.2)
par(las=1 , xaxs="i" , yaxs="i")

#Ordinary pValues
plot(res.eb$logFC, -log10(res.eb$p.ord), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="blue", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="blue", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of ordinary p-Values")

#Moderate pValues
plot(res.eb$logFC, -log10(res.eb$p.mod), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="lightgray", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="lightgray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of moderated p-Values")

#enhanced volcano plots
require(EnhancedVolcano)
EnhancedVolcano(toptable = res.eb, lab = rownames(res.eb) , x = "logFC" , y = "q.ord" , xlim = c(min(res.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res.eb[["p.ord"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot Ku80ID shScramVsCneg" , subtitle = bquote(italic("Limma T-test")) , caption = paste0("total=" , nrow(res.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))
EnhancedVolcano(toptable = res.eb, lab = rownames(res.eb) , x = "logFC" , y = "p.mod" , xlim = c(min(res.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res.eb[["p.mod"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot Ku80ID shScramVsCneg" , subtitle = bquote(italic("Limma moderated T-test")) , caption = paste0("total=" , nrow(res.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))


### Ku80ID multiple T-test of shNcl vs Negative Controls
#Subset the big dataset
Ku_Ncl<-dat[,c(1:6,19:21)]
Ku_Ncl

#define treatment and control groups for comparison
ct<-c("CN_1" , "CN_2" , "CN_3" , "CN_4" , "CN_5" , "CN_6" ,"CN_7" , "CN_8" , "CN_9" )
tr<-c("Ku_shNcl_1" , "Ku_shNcl_2" , "Ku_shNcl_3")

#define design according to syntax of limma package
design<-model.matrix(~factor(c(1,1,1,1,1,1,2,2,2)))
design

#Limma statistical analysis
colnames(design)<-c("intercept","Diff")
res.eb<-eb.fit(Ku_Ncl , design)
head(res.eb)
res.eb$NAME<-row.names(res.eb)
str(res.eb)

#Export dataset res.eb to an excel file
write_xlsx(res.eb , "./Diff_Ku_shNclVScneg.xlsx")

#Volcano plots for ordinary and moderate pValues
#Ordinary pValues
rx <- c(-1,1)*max(abs(res.eb$logFC))*1.1
ry <- c(0, 8)

par(mfrow = c(1,2) , font.lab=2 , cex.lab=1.2 , font.axis=2 , cex.axis=1.2)
par(las=1 , xaxs="i" , yaxs="i")

#Ordinary pValue
plot(res.eb$logFC, -log10(res.eb$p.ord), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="blue", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="blue", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of ordinary p-Values")

#Moderate pValues
plot(res.eb$logFC, -log10(res.eb$p.mod), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="lightgray", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="lightgray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of moderated p-Values")

#enhanced volcano plots
require(EnhancedVolcano)
EnhancedVolcano(toptable = res.eb, lab = rownames(res.eb) , x = "logFC" , y = "q.ord" , xlim = c(min(res.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res.eb[["p.ord"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot Ku80ID shNclVsCneg" , subtitle = bquote(italic("Limma T-test")) , caption = paste0("total=" , nrow(res.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))
EnhancedVolcano(toptable = res.eb, lab = rownames(res.eb) , x = "logFC" , y = "p.mod" , xlim = c(min(res.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res.eb[["p.mod"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot Ku80ID shNclVsCneg" , subtitle = bquote(italic("Limma moderated T-test")) , caption = paste0("total=" , nrow(res.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))

############## 
### Ku80ID shNcl Exp multiple T-test of shNcl vs shScramble
#Subset the big dataset
Ku_ScrNcl<-dat[,16:21]
Ku_ScrNcl
#define treatment and control groups for comparison
ct<-c("Ku_shScr_1" , "Ku_shScr_2" , "Ku_shScr_3")
tr<-c("Ku_shNcl_1" , "Ku_shNcl_2" , "Ku_shNcl_3")

#define design according to syntax of limma package
design<-model.matrix(~factor(c(1,1,1,2,2,2)))
design

#Limma statistical analysis
colnames(design)<-c("intercept","Diff")
res.eb<-eb.fit(Ku_ScrNcl , design)
head(res.eb)
res.eb$NAME<-row.names(res.eb)
str(res.eb)

#Export dataset res.eb to an excel file
write_xlsx(res.eb , "./Diff_Ku_shNclVSshScr.xlsx")

#Volcano plots for ordinary and moderate pValues
#Ordinary pValues
rx <- c(-1,1)*max(abs(res.eb$logFC))*1.1
ry <- c(0, 8)

par(mfrow = c(1,2) , font.lab=2 , cex.lab=1.2 , font.axis=2 , cex.axis=1.2)
par(las=1 , xaxs="i" , yaxs="i")

#Ordinary pValues
plot(res.eb$logFC, -log10(res.eb$p.ord), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="blue", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="blue", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of ordinary p-Values")

#Moderate pValues
plot(res.eb$logFC, -log10(res.eb$p.mod), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="lightgray", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="lightgray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of moderated p-Values")

#enhanced volcano plots
require(EnhancedVolcano)
EnhancedVolcano(toptable = res.eb, lab = rownames(res.eb) , x = "logFC" , y = "p.ord" , xlim = c(min(res.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res.eb[["p.ord"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot Ku80ID shNcl_Vs_shScram" , subtitle = bquote(italic("Limma T-test")) , caption = paste0("total=" , nrow(res.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))
EnhancedVolcano(toptable = res.eb, lab = rownames(res.eb) , x = "logFC" , y = "p.mod" , xlim = c(min(res.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res.eb[["p.mod"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot Ku80ID shNcl_Vs_shScram" , subtitle = bquote(italic("Limma moderated T-test")) , caption = paste0("total=" , nrow(res.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))

############################
### Mre11ID shNcl Experiment LIMMA multiple T-Test analysis

### Mre11ID multiple T-Test analysis shScramble vs Negative Controls
Mre_Scr<-dat[,c(1:9,22:24)]
Mre_Scr
#define treatment and control groups for comparison
ct<-c("CN_1" , "CN_2" , "CN_3" , "CN_4" , "CN_5" , "CN_6" ,"CN_7" , "CN_8" , "CN_9"  )
tr<-c("Mre_shScr_1" , "Mre_shScr_2" , "Mre_shScr_3")

#define design according to syntax of limma package
design<-model.matrix(~factor(c(1,1,1,1,1,1,1,1,1,2,2,2)))
design

#Limma statistical analysis
colnames(design)<-c("intercept","Diff")
res.eb<-eb.fit(Mre_Scr , design)
head(res.eb)
res.eb$NAME<-row.names(res.eb)
str(res.eb)

#Export dataset res.eb to an excel file
write_xlsx(res.eb , "./Diff_Mre_shScrVScneg_9c.xlsx")

#Volcano plots for ordinary and moderate pValues
#Ordinary pValues
rx <- c(-1,1)*max(abs(res.eb$logFC))*1.1
ry <- c(0, 8)

par(mfrow = c(1,2) , font.lab=2 , cex.lab=1.2 , font.axis=2 , cex.axis=1.2)
par(las=1 , xaxs="i" , yaxs="i")

#Ordinary pValues
plot(res.eb$logFC, -log10(res.eb$p.ord), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="blue", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="blue", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of ordinary p-Values")

#Moderate pValues
plot(res.eb$logFC, -log10(res.eb$p.mod), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="lightgray", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="lightgray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of moderated p-Values")

#enhanced volcano plots
require(EnhancedVolcano)
EnhancedVolcano(toptable = res.eb, lab = rownames(res.eb) , x = "logFC" , y = "p.ord" , xlim = c(min(res.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res.eb[["p.ord"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot Mre11ID shScramVsCneg" , subtitle = bquote(italic("Limma T-test")) , caption = paste0("total=" , nrow(res.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))
EnhancedVolcano(toptable = res.eb, lab = rownames(res.eb) , x = "logFC" , y = "p.mod" , xlim = c(min(res.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res.eb[["p.mod"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot Mre11ID shScramVsCneg" , subtitle = bquote(italic("Limma moderated T-test")) , caption = paste0("total=" , nrow(res.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))

### Mre11ID multiple T-Test of shNcl vs Negative Controls 
#subset the original dataset
Mre_Ncl<-dat[,c(1:9,25:27)]
Mre_Ncl
#define treatment and control groups for comparison
ct<-c("CN_1" , "CN_2" , "CN_3" , "CN_4" , "CN_5" , "CN_6" ,"CN_7" , "CN_8" , "CN_9"  )
tr<-c("Mre_shNcl_1" , "Mre_shNcl_2" , "Mre_shNcl_3")

#define design according to syntax of limma package
design<-model.matrix(~factor(c(1,1,1,1,1,1,1,1,1,2,2,2)))
design

#Limma statistical analysis
colnames(design)<-c("intercept","Diff")
res.eb<-eb.fit(Mre_Ncl , design)
head(res.eb)
res.eb$NAME<-row.names(res.eb)
str(res.eb)

#Export dataset res.eb to an excel file
write_xlsx(res.eb , "./Diff_Mre_shNclVScneg_9c.xlsx")

#Volcano plots for ordinary and moderate pValues
#Ordinary pValues
rx <- c(-1,1)*max(abs(res.eb$logFC))*1.1
ry <- c(0, 8)

par(mfrow = c(1,2) , font.lab=2 , cex.lab=1.2 , font.axis=2 , cex.axis=1.2)
par(las=1 , xaxs="i" , yaxs="i")

#Ordinary pValues
plot(res.eb$logFC, -log10(res.eb$p.ord), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="blue", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="blue", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of ordinary p-Values")

#Moderate pValues
plot(res.eb$logFC, -log10(res.eb$p.mod), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="lightgray", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="lightgray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of moderated p-Values")

#enhanced volcano plots
require(EnhancedVolcano)
EnhancedVolcano(toptable = res.eb, lab = rownames(res.eb) , x = "logFC" , y = "p.ord" , xlim = c(min(res.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res.eb[["p.ord"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot Mre11ID shNclVsCneg" , subtitle = bquote(italic("Limma T-test")) , caption = paste0("total=" , nrow(res.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))
EnhancedVolcano(toptable = res.eb, lab = rownames(res.eb) , x = "logFC" , y = "p.mod" , xlim = c(min(res.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res.eb[["p.mod"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot Mre11ID shNclVsCneg" , subtitle = bquote(italic("Limma moderated T-test")) , caption = paste0("total=" , nrow(res.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))

### Mre11ID shNcl Experiment multiple T-Test shNcl vs shScramble condition
Mre_NclScr<-dat[,22:27]
Mre_NclScr
#define treatment and control groups for comparison
ct<-c("Mre_shScr_1" , "Mre_shScr_2" , "Mre_shScr_3")
tr<-c("Mre_shNcl_1" , "Mre_shNcl_2" , "Mre_shNcl_3")

#define design according to syntax of limma package
design<-model.matrix(~factor(c(1,1,1,2,2,2)))
design

#Limma statistical analysis
colnames(design)<-c("intercept","Diff")
res.eb<-eb.fit(Mre_NclScr , design)
head(res.eb)
res.eb$NAME<-row.names(res.eb)
str(res.eb)

#Export dataset res.eb to an excel file
write_xlsx(res.eb , "./Diff_Mre_shNclVSshScr_9c.xlsx")

#Volcano plots for ordinary and moderate pValues
#Ordinary pValues
rx <- c(-1,1)*max(abs(res.eb$logFC))*1.1
ry <- c(0, 8)

par(mfrow = c(1,2) , font.lab=2 , cex.lab=1.2 , font.axis=2 , cex.axis=1.2)
par(las=1 , xaxs="i" , yaxs="i")

#Ordinary pValue
plot(res.eb$logFC, -log10(res.eb$p.ord), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="blue", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="blue", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of ordinary p-Values")

#Moderate pValues
plot(res.eb$logFC, -log10(res.eb$p.mod), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="lightgray", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="lightgray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of moderated p-Values")

#enhanced volcano plots
require(EnhancedVolcano)
EnhancedVolcano(toptable = res.eb, lab = rownames(res.eb) , x = "logFC" , y = "p.ord" , xlim = c(min(res.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res.eb[["p.ord"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot Mre11ID shNcl_Vs_shScram" , subtitle = bquote(italic("Limma T-test")) , caption = paste0("total=" , nrow(res.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))
EnhancedVolcano(toptable = res.eb, lab = rownames(res.eb) , x = "logFC" , y = "p.mod" , xlim = c(min(res.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res.eb[["p.mod"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot Mre11ID shNcl_Vs_shScram" , subtitle = bquote(italic("Limma moderated T-test")) , caption = paste0("total=" , nrow(res.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))

######################HeatMap of the Mass Spec Data
dat<-read.csv("./SensorID_RawData_All.csv")#open the original dataset with the average intensities
rownames(dat)<-dat$Protein.Group.Accessions #Use the protein names as rownames
dat<-dat[,7:33]#Subset only numerical data
str(dat)
dat<-log(dat,2)#Log2 transform the data for normal distribution
mat<-as.matrix(dat)#convert the dataset into a matrix
mat<-t(scale(t(mat),center=TRUE,scale = TRUE))#Scale the matrix by Z-scoring in a row-wise manner (transpose the original matrix): center= x-mean, scale= x/sd
mat
mat<-na.omit(mat)#Remove missing values (there was only one, for PTRH1 protein)
a<-any(is.na(mat))#Check for NA values
a<-is.na(mat)
a
#Summary statistics
dim(mat) 
range(mat) 
mean(mat)
median(mat)
quantile(mat,c(0,0.25,0.5,0.75,1))

#Generate Column Annotation
colAnno<-read.csv("./ColAnnoDF_shNclExp.csv",row.names = 1)#Read sample annotations
Col_Anno<-HeatmapAnnotation(df=colAnno, col = list(Treatment=c("Cneg"="gray80","shScramble"="blue", "shNcl"="red3"),
                                                   Sensor=c("Ku80"="seagreen1","Mre11"="steelblue1","Sirt6"="tomato")))
str(colAnno)


#Opening required packages
require(ComplexHeatmap)
require(circlize)

#Color function for HeatMaps
col_fun<-colorRamp2(c(min(mat),median(mat),max(mat)),c("royalblue","aliceblue","tomato"))
col_fun(seq( 0, 15))  

#Heatmap plot
Heatmap(mat , name = "Protein Intenzities Z-Scores" ,  col = col_fun , cluster_rows = TRUE , 
        cluster_columns = FALSE , clustering_distance_rows = "spearman" , clustering_distance_columns = "pearson" ,
        clustering_method_rows = "ward.D" , clustering_method_columns = "average" , column_names_centered = TRUE , 
        column_names_gp = gpar(fontsize=15) , row_names_gp = gpar(fontsize=3) , show_column_names = T,top_annotation = Col_Anno)

###### Working with only the differentially bound proteins.
## The LIMMA multiple T-Test result tables were manually filtered in excel to keep only those proteins with a FoldChange above 1.5, and a moderate pValue below 0.05
## The list of differentially bound proteins against the negative controls is uploaded as a separate file
#Heat map from differentially bound proteins against the negative control
df<-read.csv("./SensorID_shNclList_Diff_Prots.csv")
str(df)
DiffProts<-df$NAME #Differentially bound proteins from the limma analysis
DiffProts
#Subbseting the matrix for only differentially bound proteins
mat_DiffProts<-mat[rownames(mat)%in%DiffProts,]
#
mat_DiffProts<-na.omit(mat_DiffProts)#Remove missing values (there was only one, for PTRH1 protein)
a<-any(is.na(mat_DiffProts))#Check for NA values
a<-is.na(mat_DiffProts)
a
#Summary statistics
dim(mat_DiffProts)
range(mat_DiffProts) 
mean(mat_DiffProts)
median(mat_DiffProts)
quantile(mat_DiffProts,c(0,0.25,0.5,0.75,1))

#Generate Column Annotation
colAnno<-read.csv("./ColAnnoDF_shNclExp.csv",row.names = 1)#Read sample annotations
Col_Anno<-HeatmapAnnotation(df=colAnno, col = list(Treatment=c("Cneg"="gray80","shScramble"="blue", "shNcl"="red3"),
                                                   Sensor=c("Ku80"="seagreen1","Mre11"="steelblue1","Sirt6"="tomato")))
str(colAnno)

#Opening required packages
require(ComplexHeatmap)
require(circlize)

#Color function for HeatMaps
col_fun<-colorRamp2(c(min(mat_DiffProts),median(mat_DiffProts),max(mat_DiffProts)),c("royalblue","aliceblue","tomato"))
col_fun(seq( 0, 15))  

#Heat map plot
Heatmap(mat_DiffProts , name = "Protein Intenzities Z-Scores" ,  col = col_fun , cluster_rows = TRUE , 
        cluster_columns = FALSE , clustering_distance_rows = "spearman" , clustering_distance_columns = "pearson" ,
        clustering_method_rows = "ward.D" , clustering_method_columns = "average" , column_names_centered = TRUE , 
        column_names_gp = gpar(fontsize=15) , row_names_gp = gpar(fontsize=3) , show_column_names = T,top_annotation = Col_Anno)



##################################################
####Limma T-Test Pairwise analysis of shNcl vs shScr with only proteins above the background level (Differentially bound proteins against the Negative controls).

#Upload data
#Analysis of Sirt6ID shScramble and shNcl StrepIP mass spec results
#Load formated dataset as required for Limma analysis
dat<-read.csv("./SensorID_shNclExp_RawData_All.csv")
#Check datasets
dim(dat)
str(dat)
#Remove N/A from dataset
a<-is.na(dat)
any(is.na(dat))
unique(a)
dat<-na.omit(dat)
dim(dat)
str(dat)

#Load the list of proteins above background level from previous Limma T-test against the negative control and with nuclear localization
NucProts<-read.csv("./Nuclear_Proteins.csv")
str(NucProts)
NucProts<-NucProts$UniprotID
#Subset the list of Nuclear Proteins
dat_nuc<-dat[dat$Protein.Descriptions%in%NucProts,]
#
a<-any(is.na(dat_nuc))#Check for NA values
dat_nuc<-na.omit(dat_nuc)#Remove missing values (there was only one, for PTRH1 protein)


#check for duplicate gene names
length(unique(dat_nuc$Protein.Group.Accessions))
#2348
#define vector headlines
#Channel values for each sample in the dataset.
cha<- c("CN_1" , "CN_2" , "CN_3" , "CN_4" , "CN_5" , "CN_6" ,"CN_7" , "CN_8" , "CN_9" ,
        "S6_shScr_1" , "S6_shScr_2" , "S6_shScr_3" , "S6_shNcl_1" , "S6_shNcl_2" , "S6_shNcl_3", 
        "Ku_shScr_1" , "Ku_shScr_2" , "Ku_shScr_3" , "Ku_shNcl_1" , "Ku_shNcl_2" , "Ku_shNcl_3",
        "Mre_shScr_1" , "Mre_shScr_2" , "Mre_shScr_3" , "Mre_shNcl_1" , "Mre_shNcl_2" , "Mre_shNcl_3")

#data procession download functions
source("http://www.biostat.jhsph.edu/~kkammers/software/eupa/source.functions.r")

#Channel values for each protein are median polished log2 transformed values from the median value of all peptides belonging to each protein.
dat_nuc<-quantify.proteins(dat_nuc , cha)
dim(dat_nuc) #[1] 1392   29
dat_nuc

#boxplot: intensities of all channels after data processing and normalization
par(mfrow = c(1,1) , font.lab=1 , cex.lab=0.5 , font.axis=1 , cex.axis=0.8)
boxplot(dat[,1:length(cha)], ylim=c(-20,20) , main="boxplot normalized intensities")

### Sirt6ID multiple T-test of shScramble vs Negative Controls
#Subset the data
S6_NclvsScr<-dat_nuc[,c(10:15)]
S6_NclvsScr
#define treatment and control groups for comparison
ct<-c("S6_shScr_1" , "S6_shScr_2" , "S6_shScr_3")
tr<-c("S6_shNcl_1" , "S6_shNcl_2" , "S6_shNcl_3")

#define design according to syntax of limma package
design<-model.matrix(~factor(c(1,1,1,2,2,2)))
design

#Limma statistical analysis
colnames(design)<-c("intercept","Diff")
res.eb<-eb.fit(S6_NclvsScr , design)
head(res.eb)
res.eb$NAME<-row.names(res.eb)
str(res.eb)

#Export dataset res.eb to an excel file
write_xlsx(res.eb , "./Diff_S6_shNclVSshScr_NuclearProtins.xlsx")

#Volcano plots for ordinary and moderate pValues
#Ordinary pValues
rx <- c(-1,1)*max(abs(res.eb$logFC))*1.1
ry <- c(0, 8)

par(mfrow = c(1,2) , font.lab=2 , cex.lab=1.2 , font.axis=2 , cex.axis=1.2)
par(las=1 , xaxs="i" , yaxs="i")

#Ordinary pValues
plot(res.eb$logFC, -log10(res.eb$p.ord), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="blue", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="blue", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of ordinary p-Values")

#Moderate pValues
plot(res.eb$logFC, -log10(res.eb$p.mod), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="lightgray", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="lightgray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of moderated p-Values")

#enhanced volcano plots
EnhancedVolcano(toptable = res.eb, lab = rownames(res.eb) , x = "logFC" , y = "p.ord" , xlim = c(min(res.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res.eb[["p.ord"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot Sirt6ID shNcl_Vs_shScram" , subtitle = bquote(italic("Limma T-test")) , caption = paste0("total=" , nrow(res.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))
EnhancedVolcano(toptable = res.eb, lab = rownames(res.eb) , x = "logFC" , y = "p.mod" , xlim = c(min(res.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res.eb[["p.mod"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot Sirt6ID shNcl_Vs_shScram" , subtitle = bquote(italic("Limma moderated T-test")) , caption = paste0("total=" , nrow(res.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))


### Ku80ID multiple T-test of shNcl vs shScramble
#Subset the big dataset
Ku_ScrNcl<-dat_nuc[,16:21]
Ku_ScrNcl
#define treatment and control groups for comparison
ct<-c("Ku_shScr_1" , "Ku_shScr_2" , "Ku_shScr_3")
tr<-c("Ku_shNcl_1" , "Ku_shNcl_2" , "Ku_shNcl_3")

#define design according to syntax of limma package
design<-model.matrix(~factor(c(1,1,1,2,2,2)))
design
#Limma statistical analysis
colnames(design)<-c("intercept","Diff")
res.eb<-eb.fit(Ku_ScrNcl , design)
head(res.eb)
res.eb$NAME<-row.names(res.eb)
str(res.eb)

#Export dataset res.eb to an excel file
write_xlsx(res.eb , "./Diff_Ku_shNclVSshScr_NuclearProteins.xlsx")

#Volcano plots for ordinary and moderate pValues
#Ordinary pValues
rx <- c(-1,1)*max(abs(res.eb$logFC))*1.1
ry <- c(0, 8)

par(mfrow = c(1,2) , font.lab=2 , cex.lab=1.2 , font.axis=2 , cex.axis=1.2)
par(las=1 , xaxs="i" , yaxs="i")

#Ordinary pValue
plot(res.eb$logFC, -log10(res.eb$p.ord), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="blue", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="blue", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of ordinary p-Values")

#Moderate pValues
plot(res.eb$logFC, -log10(res.eb$p.mod), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="lightgray", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="lightgray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of moderated p-Values")

#enhanced volcano plots
EnhancedVolcano(toptable = res.eb, lab = rownames(res.eb) , x = "logFC" , y = "p.ord" , xlim = c(min(res.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res.eb[["p.ord"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot Ku80ID shNcl_Vs_shScram" , subtitle = bquote(italic("Limma T-test")) , caption = paste0("total=" , nrow(res.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))
EnhancedVolcano(toptable = res.eb, lab = rownames(res.eb) , x = "logFC" , y = "p.mod" , xlim = c(min(res.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res.eb[["p.mod"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot Ku80ID shNcl_Vs_shScram" , subtitle = bquote(italic("Limma moderated T-test")) , caption = paste0("total=" , nrow(res.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))

### Mre11ID multiple T-test of shNcl vs shScramble
Mre_NclScr<-dat_nuc[,22:27]
Mre_NclScr
#define treatment and control groups for comparison
ct<-c("Mre_shScr_1" , "Mre_shScr_2" , "Mre_shScr_3")
tr<-c("Mre_shNcl_1" , "Mre_shNcl_2" , "Mre_shNcl_3")

#define design according to syntax of limma package
design<-model.matrix(~factor(c(1,1,1,2,2,2)))
design

#Limma statistical analysis
colnames(design)<-c("intercept","Diff")
res.eb<-eb.fit(Mre_NclScr , design)
head(res.eb)
res.eb$NAME<-row.names(res.eb)
str(res.eb)

#Export dataset res.eb to an excel file
write_xlsx(res.eb , "./Diff_Mre_shNclVSshScr_NuclearProteins.xlsx")

#Volcano plots for ordinary and moderate pValues
#Ordinary pValues
rx <- c(-1,1)*max(abs(res.eb$logFC))*1.1
ry <- c(0, 8)

par(mfrow = c(1,2) , font.lab=2 , cex.lab=1.2 , font.axis=2 , cex.axis=1.2)
par(las=1 , xaxs="i" , yaxs="i")

#Ordinary pValues
plot(res.eb$logFC, -log10(res.eb$p.ord), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="blue", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="blue", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of ordinary p-Values")

#Moderate pValues
plot(res.eb$logFC, -log10(res.eb$p.mod), pch=21, bg="Blue", cex=0.9, xlim=rx, ylim=ry, xaxt="n", xlab="Fold Change", ylab="-Log10 pVal")
abline(v=seq(-2,2,1), col="lightgray", lty="dotted")
abline(h=seq(0, ry[2] , 1), col="lightgray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("VolcanoPlot of moderated p-Values")

#enhanced volcano plots
EnhancedVolcano(toptable = res.eb, lab = rownames(res.eb) , x = "logFC" , y = "p.ord" , xlim = c(min(res.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res.eb[["p.ord"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot Mre11ID shNcl_Vs_shScram" , subtitle = bquote(italic("Limma T-test")) , caption = paste0("total=" , nrow(res.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))
EnhancedVolcano(toptable = res.eb, lab = rownames(res.eb) , x = "logFC" , y = "p.mod" , xlim = c(min(res.eb[["logFC"]], na.rm = TRUE) - 1.5 , max(res.eb[["logFC"]], na.rm = TRUE) +1.5) , ylim = c(0 , max(-log10(res.eb[["p.mod"]]) , na.rm = TRUE)+2) , xlab = "Fold Change" , ylab = "-Log10 pVal" , axisLabSize = 18 , title = "Volcano Plot Mre11ID shNcl_Vs_shScram" , subtitle = bquote(italic("Limma moderated T-test")) , caption = paste0("total=" , nrow(res.eb) , "variables") , pCutoff = 0.05 , FCcutoff = 1.5 , col = c("grey30" , "forestgreen" , "royalblue" , "red2"))

##Heatmaps
######################HeatMap of the Mass Spec Data
dat<-read.csv("./SensorID_shNclExp_RawData_All.csv")#open the original dataset with the average intensities
NucProts<-read.csv("./SensorID_shNclList_Diff_Prots.csv")
dat_nuc<-dat[dat$Protein.Descriptions%in%NucProts$UniprotID,]
row.names(dat_nuc)<-dat_nuc$Protein.Group.Accessions #Use the protein names as rownames
mat<-dat_nuc[,7:33]#Subset only numerical data
mat<-as.matrix(mat)#convert the dataset into a matrix
mat<-t(scale(t(mat),center=TRUE,scale = TRUE))#Scale the matrix by Z-scoring in a row-wise manner (transpose the original matrix): center= x-mean, scale= x/sd
str(mat)
mat<-na.omit(mat)#Remove missing values (there was only one, for PTRH1 protein)
a<-any(is.na(mat))#Check for NA values
a
#Summary statistics
dim(mat) 
range(mat) 
mean(mat)
median(mat)
quantile(mat,c(0,0.25,0.5,0.75,1))

#Opening required packages
library(ComplexHeatmap)
library(circlize)
#Generate Column Annotation
colAnno<-read.csv("./ColAnnoDF_shNclExp.csv",row.names = 1)#Read sample annotations
Col_Anno<-HeatmapAnnotation(df=colAnno, col = list(Treatment=c("Cneg"="gray80","shScramble"="blue", "shNcl"="red3"),
                                                   Sensor=c("Ku80"="seagreen1","Mre11"="steelblue1","Sirt6"="tomato")))
str(colAnno)

#Generate Annotation vectors for proteins
NucLoc<-as.factor(NucProts$CellCompartment)
NucLoc<-as.matrix(NucLoc)
NucAnno<-rowAnnotation(Nuclear_Localization=as.factor(NucLoc))#annotation_legend_param=)#,labels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22"))

#Color function for HeatMaps
col_fun<-colorRamp2(c(min(mat),median(mat),max(mat)),c("royalblue","aliceblue","tomato"))
col_fun(seq( 0, 15))  
#Heat map plot
Heatmap(mat , name = "Protein Intenzities Z-Scores" ,  col = col_fun , cluster_rows = TRUE , 
        cluster_columns = FALSE , clustering_distance_rows = "spearman" , clustering_distance_columns = "pearson" ,
        clustering_method_rows = "ward.D" , clustering_method_columns = "average" , column_names_centered = TRUE , 
        column_names_gp = gpar(fontsize=15) , row_names_gp = gpar(fontsize=3) , show_column_names = T,top_annotation = Col_Anno,
        row_split = NucLoc ,  row_gap = unit(1, "mm") , 
        left_annotation = NucAnno)

