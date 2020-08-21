#Author: Farnaz Fouladi
#Date: 08-17-2020
#Description: Ordination plot.

rm(list=ls())

#Libraries
library(vegan)

output<-"./output/"
taxa<-"Genus"
colors=c("red","blue","darkgreen","darkorange2","purple","hotpink","black","firebrick")
source("./Rcode/RYGB_IntegratedAnalysis/functions.R")


myT<-read.table(paste0(output,taxa,"_countTable_merged_log10.txt"),sep="\t",header = TRUE)
meta<-read.table(paste0(output,"metaData_merged.txt"),sep="\t",header = TRUE)
meta$TimeEachStudy<-paste0(meta$Study,"_",meta$time)

pdf(paste0(output,taxa,"_PCO.pdf"),width = 5,height = 10)
par(mfrow=c(2,1))

getPCO(myT,meta,"Study",names=levels(factor(meta$Study)),colors)
getPCO(myT,meta,"TimeEachStudy",names=levels(factor(meta$TimeEachStudy)),colors)

dev.off()

#PERMANOVA 
adonis(myT~factor(meta$time)*factor(meta$Study))

"                                      Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
factor(meta$time)                      1    0.8352 0.83523 11.2924 0.04014  0.001 ***
factor(meta$Study)                     3    6.4087 2.13625 28.8823 0.30798  0.001 ***
factor(meta$time):factor(meta$Study)   3    0.2518 0.08393  1.1347 0.01210  0.250    
Residuals                            180   13.3135 0.07396         0.63979           
Total                                187   20.8092                 1.00000           
---"


