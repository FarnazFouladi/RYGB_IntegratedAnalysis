#Author: Farnaz Fouladi
#Date: 08-17-2020
#Description: Ordination plot.

rm(list=ls())

#Libraries
library(vegan)

set.seed(123)

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

                                      Df SumsOfSqs MeanSqs F.Model
factor(meta$time)                      1    0.7379 0.73786  9.1916
factor(meta$Study)                     3    6.7607 2.25357 28.0729
factor(meta$time):factor(meta$Study)   3    0.2781 0.09269  1.1546
Residuals                            180   14.4496 0.08028        
Total                                187   22.2262                
                                          R2 Pr(>F)    
factor(meta$time)                    0.03320  0.001 ***
factor(meta$Study)                   0.30418  0.001 ***
factor(meta$time):factor(meta$Study) 0.01251  0.219    
Residuals                            0.65012           
Total                                1.00000           
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1"








