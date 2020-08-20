#Author: Farnaz Fouladi
#Date: 08-17-2020
#Description: Heatmap

rm(list=ls())

#Libraries
library(ComplexHeatmap)

output<-"./output/"
taxa<-"Genus"
source("./Rcode/functions.R")
studies<-c("BS-1M","BS-6M","Assal-3M","Assal-1Y","Assal-2Y","Ilhan-6M","Ilhan-1Y","Afshar-6M")

myT.BS<-read.table(paste0(output,taxa,"_BS_MixedLinearModelResults.txt"),sep="\t",header = TRUE,row.names = 1)
myT.Assal<-read.table(paste0(output,taxa,"_Assal_MixedLinearModelResults.txt"),sep="\t",header = TRUE,row.names = 1)
myT.Ilhan<-read.table(paste0(output,taxa,"_Ilhan_MixedLinearModelResults.txt"),sep="\t",header = TRUE,row.names = 1)
myT.Afshar<-read.table(paste0(output,taxa,"_Afshar_MixedLinearModelResults.txt"),sep="\t",header = TRUE,row.names = 1)

myT.BS$logp_1M<-getlog10p(myT.BS$p1M,myT.BS$s1M)
myT.BS$logp_6M<-getlog10p(myT.BS$p6M,myT.BS$s6M)
myT.Assal$logp_3M<-getlog10p(myT.Assal$p3M,myT.Assal$s3M)
myT.Assal$logp_1Y<-getlog10p(myT.Assal$p1Y,myT.Assal$s1Y)
myT.Assal$logp_2Y<-getlog10p(myT.Assal$p2Y,myT.Assal$s2Y)
myT.Ilhan$logp_6M<-getlog10p(myT.Ilhan$p6M,myT.Ilhan$s6M)
myT.Ilhan$logp_1Y<-getlog10p(myT.Ilhan$p1Y,myT.Ilhan$s1Y)
myT.Afshar$logp_6M<-getlog10p(myT.Afshar$p6M,myT.Afshar$s6M)


bugs<-intersect(rownames(myT.BS),rownames(myT.Assal))
bugs1<-intersect(bugs,rownames(myT.Ilhan))
bugs2<-intersect(bugs1,rownames(myT.Afshar))

myT.BS_c<-myT.BS[bugs2,]
myT.Assal_c<-myT.Assal[bugs2,]
myT.Ilhan_c<-myT.Ilhan[bugs2,]
myT.Afshar_c<-myT.Afshar[bugs2,]

df<-matrix(c(myT.BS_c$logp_1M,myT.BS_c$logp_6M,
             myT.Assal_c$logp_3M,myT.Assal_c$logp_1Y,myT.Assal_c$logp_2Y,
             myT.Ilhan_c$logp_6M,myT.Ilhan_c$logp_1Y,
             myT.Afshar_c$logp_6M),ncol = length(bugs2),nrow=length(studies),byrow = TRUE)
rownames(df)<-studies
colnames(df)<-bugs2


#Heatmap
df<-t(df)
studies<-sapply(colnames(df), function(x){strsplit(x,"-")[[1]][1]})

col = list(Study = c("BS" = "yellow", "Assal" = "pink","Ilhan" = "blue","Afshar" = "grey"))

ha <- HeatmapAnnotation(
  Study = studies, col = col
)

pdf(paste0(output,taxa,"_HeatmapAndCluster.pdf"),height = 8)
Heatmap(df,top_annotation = ha,row_names_gp = gpar(fontsize = 6),heatmap_height = unit(20, "cm"),
        name = "log10 p-value")

dev.off()


