#Author: Farnaz Fouladi
#Date: 08-17-2020
#Description: Assal dataset:Compare the gut microbiome at each time point versus
#              baseline using mixed linear models.Count tables are classified by 
#             DADA2.

rm(list=ls())

#Libraries
library(nlme)

input<-"./input/RYGB_Assal2020/"
output<-"./output/MixedLinearModels/"
taxa<-c("Phylum","Class","Order","Family","Genus","SV","Seq")

for (t in taxa){
  
  dada2<-read.table(paste0(input,t,"_norm_table.txt"),sep="\t",header = TRUE,row.names = 1,check.names = FALSE)
  
  #Remove time point 3 years due to small sample size
  dada2<-dada2[dada2$time!="3Y",]
  dada2$prepost<-sapply(dada2$time,function(x){if (x=="Pre") return(0) else return(1)})
  
  if(t=="Genus")
    colnames(dada2)[colnames(dada2)=="Esherichica/Shigella"]<-"Escherichia/Shigella"
  
  #write table
  write.table(dada2,paste0(input,t,"_norm_table_updated.txt"),sep="\t",quote = FALSE)
  
  finishAbundanceIndex<-which(colnames(dada2)=="Assay.Type")-1
  myT<-dada2[,1:finishAbundanceIndex]
  meta<-dada2[,(finishAbundanceIndex+1):ncol(dada2)]
  meta$time<-factor(meta$time,levels = c("Pre","3M","1Y","2Y"))
  
  pval<-vector()
  p3M<-vector()
  p1Y<-vector()
  p2Y<-vector()
  s3M<-vector()
  s1Y<-vector()
  s2Y<-vector()
  bugName<-vector()
  index<-1
  
  for (i in 1:ncol(myT)){
    
    bug<-myT[,i]
    
    if (mean(bug>0)>0.1){
      
      df<-data.frame(bug,meta)
      
      fit<-anova(lme(bug~time,method="REML",random=~1|ID,data=df))
      sm<-summary(lme(bug~time,method="REML",random=~1|ID,data=df))
      
      pval[index]<-fit$`p-value`[2]
      p3M[index]<-sm$tTable[2,5]
      p1Y[index]<-sm$tTable[3,5]
      p2Y[index]<-sm$tTable[4,5]
      s3M[index]<-sm$tTable[2,1]
      s1Y[index]<-sm$tTable[3,1]
      s2Y[index]<-sm$tTable[4,1]
      
      bugName[index]<-colnames(myT)[i]
      index<-index+1
    }
  }
  
  df<-data.frame(bugName,pval,p3M,p1Y,p2Y,s3M,s1Y,s2Y)
  df<-df[order(df$pval),]
  df$Adjustedpval<-p.adjust(df$pval,method = "BH")
  df$Adjustedp3M<-p.adjust(df$p3M,method = "BH")
  df$Adjustedp1Y<-p.adjust(df$p1Y,method = "BH")
  df$Adjustedp2Y<-p.adjust(df$p2Y,method = "BH")
  
  write.table(df,paste0(output,t,"_Assal_MixedLinearModelResults.txt"),sep="\t",row.names = FALSE,quote = FALSE)
}  
