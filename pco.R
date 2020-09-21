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
"Call:
adonis(formula = myT ~ factor(meta$time) * factor(meta$Study)) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

                                      Df SumsOfSqs MeanSqs F.Model      R2
factor(meta$time)                      1    0.7310 0.73104  9.1596 0.03296
factor(meta$Study)                     3    6.8119 2.27064 28.4503 0.30712
factor(meta$time):factor(meta$Study)   3    0.2710 0.09035  1.1320 0.01222
Residuals                            180   14.3660 0.07981         0.64770
Total                                187   22.1799                 1.00000
                                     Pr(>F)    
factor(meta$time)                     0.001 ***
factor(meta$Study)                    0.001 ***
factor(meta$time):factor(meta$Study)  0.277    
Residuals                                      
Tota"










