#Author: Farnaz Fouladi
#Date: 08-17-2020
#Description: BS dataset:Compare the gut microbiome at each time point versus
#              baseline using mixed linear models. Count tables are classified by 
#             DADA2.

rm(list=ls())

#Libraries
library(nlme)

input<-"./input/RYGB_BS/"
output<-"./output/MixedLinearModels/"
taxa<-c("Phylum","Class","Order","Family","Genus","SV","Seq")

for (t in taxa){
  
  
  dada2<-read.table(paste0(input,t,"_norm_table.txt"),sep="\t",header = TRUE,row.names = 1,check.names = FALSE)
   
  #Selecting RYGB surgeries
  dada2_rygb<-dada2[dada2$Type.of.Surgery ==1 & !is.na(dada2$Type.of.Surgery),]
  
  #Removing samples that do not have baselines
  dada2_rygb<-dada2_rygb[dada2_rygb$Participant_ID %in% (dada2_rygb$Participant_ID[dada2_rygb$Timepoint==0]),]
  
  dada2_rygb$prepost<-sapply(dada2_rygb$Timepoint,function(x){if (x==0) return(0) else return(1)})
  
  if(t=="Genus")
    colnames(dada2_rygb)[colnames(dada2_rygb)=="Esherichica/Shigella"]<-"Escherichia/Shigella"
  
  #write table
  write.table(dada2_rygb,paste0(input,t,"_norm_table_updated.txt"),sep="\t",quote = FALSE)
  
  finishAbundanceIndex<-which(colnames(dada2_rygb)=="Patient.ID" )-1
  myT<-dada2_rygb[,1:finishAbundanceIndex]
  meta<-dada2_rygb[,(finishAbundanceIndex+1):ncol(dada2_rygb)]

  pval<-vector()
  p1M<-vector()
  p6M<-vector()
  s1M<-vector()
  s6M<-vector()
  bugName<-vector()
  index<-1
  
  for (i in 1:ncol(myT)){
    
    bug<-myT[,i]
    
    if (mean(bug>0)>0.1){
      
      df<-data.frame(bug,meta)
      
      fit<-anova(lme(bug~factor(Timepoint),method="REML",random=~1|Participant_ID,data=df))
      sm<-summary(lme(bug~factor(Timepoint),method="REML",random=~1|Participant_ID,data=df))
      
      pval[index]<-fit$`p-value`[2]
      p1M[index]<-sm$tTable[2,5]
      p6M[index]<-sm$tTable[3,5]
      s1M[index]<-sm$tTable[2,1]
      s6M[index]<-sm$tTable[3,1]
      bugName[index]<-colnames(myT)[i]
      index<-index+1
    }
  }
  
  df<-data.frame(bugName,pval,p1M,p6M,s1M,s6M)
  df<-df[order(df$pval),]
  df$Adjustedpval<-p.adjust(df$pval,method = "BH")
  df$Adjustedp1M<-p.adjust(df$p1M,method = "BH")
  df$Adjustedp6M<-p.adjust(df$p6M,method = "BH")
  
  write.table(df,paste0(output,t,"_BS_MixedLinearModelResults.txt"),sep="\t",row.names = FALSE,quote = FALSE)
}  
  
  