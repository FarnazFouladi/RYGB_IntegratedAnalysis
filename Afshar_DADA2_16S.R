#Author: Farnaz Fouladi
#Date: 08-17-2020
#Description: Afshar dataset:Compare the gut microbiome at each time point versus
#             baseline using mixed linear models.Count tables are classified by 
#             DADA2.

rm(list=ls())

#Libraries
library(nlme)

input<-"./input/RYGB_Afshar2018/"
output<-"./output/MixedLinearModels/"
taxa<-c("Phylum","Class","Order","Family","Genus","SV")

for (t in taxa){
  
  
  dada2<-read.table(paste0(input,t,"_norm_table.txt"),sep="\t",header = TRUE,row.names = 1,check.names = FALSE)
  
  dada2<-dada2[!is.na(dada2$ID),]
  dada2$time<-relevel(factor(dada2$time),ref="Pre")
  dada2$prepost<-sapply(dada2$time,function(x){if (x=="Pre") return(0) else return(1)})
  if(t=="Genus")
   colnames(dada2)[colnames(dada2)=="Esherichica/Shigella"]<-"Escherichia/Shigella"
  
  #write table
  write.table(dada2,paste0(input,t,"_norm_table_updated.txt"),sep="\t",quote = FALSE)
  
  finishAbundanceIndex<-which(colnames(dada2)=="Alias")-1
  myT<-dada2[,1:finishAbundanceIndex]
  meta<-dada2[,(finishAbundanceIndex+1):ncol(dada2)]
  
  pval<-vector()
  p6M<-vector()
  s6M<-vector()
  bugName<-vector()
  index<-1
  
  for (i in 1:ncol(myT)){
    
    bug<-myT[,i]
    
    if (mean(bug>0)>0.1){
      
      df<-data.frame(bug,meta)
      
      fit<-anova(lme(bug~time,method="REML",random=~1|ID,data=df))
      sm<-summary(lme(bug~time,method="REML",random=~1|ID,data=df))
      
      pval[index]<-fit$`p-value`[2]
      p6M[index]<-sm$tTable[2,5]
      s6M[index]<-sm$tTable[2,1]
      bugName[index]<-colnames(myT)[i]
      index<-index+1
    }
  }
  
  df<-data.frame(bugName,pval,p6M,s6M)
  df<-df[order(df$pval),]
  df$Adjustedpval<-p.adjust(df$pval,method = "BH")
  df$Adjustedp6M<-p.adjust(df$p6M,method = "BH")
 
  write.table(df,paste0(output,t,"_Afshar_MixedLinearModelResults.txt"),sep="\t",row.names = FALSE,quote = FALSE)
}  
