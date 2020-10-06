#Author: Farnaz Fouladi
#Date: 08-17-2020
#Description: Compare the metabolic pathways at each time point versus baseline using
#             mixed linear models and metagenomic data.

rm(list=ls())

#Libraries
library(nlme)

input<-"./input/humann2/BS/"
output<-"./output/MixedLinearModels/humann2/"

myT<-read.table(paste0(input,"humanN2_pathabundance_cpm.tsv"),
                sep="\t",header = TRUE,check.names = FALSE,comment.char = "",quote ="",row.names = 1)

meta<-read.table(paste0(input,"metagenomics_Metadata.txt"),
                 sep="\t",header = TRUE,row.names = 1)


myT<-as.data.frame(t(myT))

if(sum(rownames(myT)==rownames(meta))!=nrow(meta)) stop("Error!")

#Unstratified table
myT_unstratified<-myT[,-grep("|",colnames(myT),fixed = TRUE)]
dim(myT_unstratified) #135 493
meta$Timepoint<-factor(meta$Timepoint)

pval<-vector()
p1M<-vector()
p6M<-vector()
p1Y<-vector()
s1M<-vector()
s6M<-vector()
s1Y<-vector()
bugName<-vector()
index<-1

for (i in 1:ncol(myT_unstratified)){
  
  bug<-myT_unstratified[,i]
  
  if (mean(bug>0)>0.1){
    
    df<-data.frame(bug,meta)
    
    fit<-anova(lme(bug~Timepoint,method="REML",random=~1|PatientID,data=df))
    sm<-summary(lme(bug~Timepoint,method="REML",random=~1|PatientID,data=df))
    
    pval[index]<-fit$`p-value`[2]
    p1M[index]<-sm$tTable[2,5]
    p6M[index]<-sm$tTable[3,5]
    p1Y[index]<-sm$tTable[4,5]
    
    s1M[index]<-sm$tTable[2,1]
    s6M[index]<-sm$tTable[3,1]
    s1Y[index]<-sm$tTable[4,1]
    
    bugName[index]<-colnames(myT_unstratified)[i]
    index<-index+1
    
  }
}

df<-data.frame(bugName,pval,p1M,p6M,p1Y,s1M,s6M,s1Y)
df<-df[order(df$pval),]
df$Adjustedpval<-p.adjust(df$pval,method = "BH")
df$Adjustedp1M<-p.adjust(df$p1M,method = "BH")
df$Adjustedp6M<-p.adjust(df$p6M,method = "BH")
df$Adjustedp1Y<-p.adjust(df$p1Y,method = "BH")

write.table(df,paste0(output,"Pathway_BS_MixedLinearModelResults.txt"),sep="\t",row.names = FALSE,quote = FALSE)



