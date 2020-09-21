#Author: Farnaz Fouladi
#Date: 08-17-2020
#Description: Compare the metabolic pathways at each time point versus baseline using
#             mixed linear models and metagenomic data.

rm(list=ls())

#Libraries
library(nlme)

input<-"./input/"
output<-"./output/MixedLinearModels/humann2/"

myT<-read.table(paste0(input,"humann2/BS/humanN2_pathabundance_cpm.tsv"),
                sep="\t",header = TRUE,check.names = FALSE,comment.char = "",quote ="",row.names = 1)

meta<-read.table(paste0(input,"RYGB_BS_Metagenomics/Metadata_BS_RYGB.txt"),
                 sep="\t",comment.char = "",check.names = FALSE,quote = "",header = TRUE)


myT<-as.data.frame(t(myT))
rownames(myT)<-sapply(rownames(myT),function(x){
  if(substr(x,1,2)=="BS") return(substr(x,1,6))
  else  return(strsplit(x,"_")[[1]][1])
})

myT1<-myT[meta$Sample.ID,]
myT1_RYGB<-myT1[meta$Type.of.Surgery==1 & !is.na(meta$Type.of.Surgery),]
meta1<-meta[meta$Type.of.Surgery==1 & !is.na(meta$Type.of.Surgery),]

#Removing samples at 18 months after surgery (n=2)
myT2_RYGB<-myT1_RYGB[meta1$Timepoint!=18,]
meta2<-meta1[meta1$Timepoint!=18,]

#Removing samples that do not have baseline samples
myT3_RYGB<-myT2_RYGB[meta2$Participant_ID %in% (meta2$Participant_ID[meta2$Timepoint==0]),]
meta3<-meta2[meta2$Participant_ID %in% (meta2$Participant_ID[meta2$Timepoint==0]),]

#Removing duplicated ID ("BIO-2-013-1")
myT3_RYGB<-myT3_RYGB[!duplicated(meta3$Patient.ID.x),]
meta3<-meta3[!duplicated(meta3$Patient.ID.x),]

#Unstratified table
myT3_unstratified<-myT3_RYGB[,-grep("|",colnames(myT3_RYGB),fixed = TRUE)]
dim(myT3_unstratified) #135 481

meta3$Timepoint<-factor(meta3$Timepoint)

pval<-vector()
p1M<-vector()
p6M<-vector()
p1Y<-vector()
s1M<-vector()
s6M<-vector()
s1Y<-vector()
bugName<-vector()
index<-1

for (i in 1:ncol(myT3_unstratified)){
  
  bug<-myT3_unstratified[,i]
  
  if (mean(bug>0)>0.1){
    
    df<-data.frame(bug,meta3)
    
    fit<-anova(lme(bug~Timepoint,method="REML",random=~1|Participant_ID,data=df))
    sm<-summary(lme(bug~Timepoint,method="REML",random=~1|Participant_ID,data=df))
    
    pval[index]<-fit$`p-value`[2]
    p1M[index]<-sm$tTable[2,5]
    p6M[index]<-sm$tTable[3,5]
    p1Y[index]<-sm$tTable[4,5]
    
    s1M[index]<-sm$tTable[2,1]
    s6M[index]<-sm$tTable[3,1]
    s1Y[index]<-sm$tTable[4,1]
    
    bugName[index]<-colnames(myT3_unstratified)[i]
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



