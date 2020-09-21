#Author: Farnaz Fouladi
#Date: 08-17-2020
#Description: Compare the metabolic pathways at each time point versus baseline using
#             mixed linear models and metagenomic data.

rm(list=ls())

#Libraries
library(nlme)

input<-"./input/"
output<-"./output/MixedLinearModels/humann2/"

myT<-read.table(paste0(input,"humann2/Palleja/humanN2_pathabundance_cpm.tsv"),
                sep="\t",header = TRUE,check.names = FALSE,comment.char = "",quote ="",row.names = 1)

#Unstrafied table
path_unstratified<-myT[-grep("|",rownames(myT),fixed = TRUE),]
dim(path_unstratified) #418  33
myT1<-as.data.frame(t(path_unstratified))
rownames(myT1)<-sapply(rownames(myT1), function(x){
  strsplit(x,"_")[[1]][1]
})

#Metadata
meta<-read.table(paste0(input,"palleja2016/metaData.txt"),sep=",",header = TRUE)
meta<-meta[meta$LibraryLayout=="PAIRED",]
meta$ID<-sapply(meta$Title,function(x){strsplit(x,"_")[[1]][1]})
meta$time<-factor(sapply(meta$Title,function(x){strsplit(x,"_")[[1]][2]}),levels = c("Baseline","3MO","1Y"))

#Ordering count table and metadata

myT2<-myT1[match(meta$Run,rownames(myT1)),]

if(sum(rownames(myT2)==meta$Run)!=nrow(meta)) stop("Error!")

pval<-vector()
p3M<-vector()
p1Y<-vector()
s3M<-vector()
s1Y<-vector()
bugName<-vector()
index<-1

for (i in 1:ncol(myT2)){
  
  bug<-myT2[,i]
  
  if (mean(bug>0)>0.1){
    
    df<-data.frame(bug,meta)
    
    fit<-anova(lme(bug~time,method="REML",random=~1|ID,data=df))
    sm<-summary(lme(bug~time,method="REML",random=~1|ID,data=df))
    
    pval[index]<-fit$`p-value`[2]
    p3M[index]<-sm$tTable[2,5]
    p1Y[index]<-sm$tTable[3,5]
    
    s3M[index]<-sm$tTable[2,1]
    s1Y[index]<-sm$tTable[3,1]
    
    
    bugName[index]<-colnames(myT2)[i]
    index<-index+1
    
  }
}

df<-data.frame(bugName,pval,p3M,p1Y,s3M,s1Y)
df<-df[order(df$pval),]
df$Adjustedpval<-p.adjust(df$pval,method = "BH")
df$Adjustedp3M<-p.adjust(df$p3M,method = "BH")
df$Adjustedp1Y<-p.adjust(df$p1Y,method = "BH")

write.table(df,paste0(output,"Pathway_Palleja_MixedLinearModelResults.txt"),sep="\t",row.names = FALSE,quote = FALSE)




