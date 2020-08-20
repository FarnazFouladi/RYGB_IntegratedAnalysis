#Author: Farnaz Fouladi
#Date: 08-17-2020
#Description: Compare the gut microbiome at each time point versus baseline.

rm(list=ls())

#Libraries
library(nlme)

input<-"./input/RYGB_BS/"
output<-"./output/MixedLinearModels/"
taxa<-c("Phylum","Class","Order","Family","Genus","SV")

for (t in taxa){
  
  
  dada2<-read.table(paste0(input,t,"_norm_table.txt"),sep="\t",header = TRUE,row.names = 1)
  surgeryType<-read.table(paste0(input,"TypeDateofSurgery-BiobehavioralR016-12-2020_FF.txt"),sep="\t",header = TRUE)
  dada2$Patient.ID<-sapply(as.character(dada2$Patient.ID),function(x){substr(x,1,9)})
  
  #Merge count table with type of surgery
  dada2_merge<-merge(dada2,surgeryType,by.x = "Patient.ID",by.y = "Patient.ID",all.x = TRUE,all.y = FALSE,sort = FALSE)
  
  #Selecting RYGB surgeries
  dada2_merge_rygb<-dada2_merge[dada2_merge$Type.of.Surgery==1 & !is.na(dada2_merge$Type.of.Surgery),]
  
  #Removing samples that do not have baselines
  dada2_merge_rygb<-dada2_merge_rygb[dada2_merge_rygb$Patient.ID %in% (dada2_merge_rygb$Patient.ID[dada2_merge_rygb$Timepoint==0]),]
  
  dada2_merge_rygb$prepost<-sapply(dada2_merge_rygb$Timepoint,function(x){if (x==0) return(0) else return(1)})
  
  #write table
  write.table(dada2_merge_rygb,paste0(input,t,"_norm_table_updated.txt"),sep="\t",quote = FALSE)
  
  finishAbundanceIndex<-which(colnames(dada2)=="Sample.ID")-1
  myT<-dada2_merge_rygb[,2:finishAbundanceIndex]
  meta<-dada2_merge_rygb[,(finishAbundanceIndex+1):ncol(dada2_merge_rygb)]

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
  
  
  
  
  