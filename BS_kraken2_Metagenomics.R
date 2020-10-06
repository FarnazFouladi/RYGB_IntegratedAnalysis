#Author: Farnaz Fouladi
#Date: 08-17-2020
#Description: Compare the gut microbiome at each time point versus baseline using
#             mixed linear models and metagenomic data.

rm(list=ls())

#Libraries
library(nlme)

input<-"./input/RYGB_BS_Metagenomics/"
output<-"./output/MixedLinearModels/Kraken2/"
taxa<-c("Phylum","Class","Order","Family","Genus","Species")


for (t in c("Phylum","Class","Order","Family","Genus","Species")){
  
  kraken<-read.table(paste0(input,"RYGB_BS_metagenomics_2020Oct01_taxaCount_norm_Log10_",t,".tsv"),
                     sep="\t",comment.char = "",check.names = FALSE,quote = "",row.names = 1,header = TRUE)
  
  meta<-read.table(paste0(input,"metagenomics_Metadata.txt"),
                   sep="\t",row.names = 1,header = TRUE)
  
  
  kraken<-kraken[match(rownames(meta),rownames(kraken)),]
  meta$Timepoint<-factor(meta$Timepoint)
  
  if(t=="Species"){
    myT<-cbind(kraken,meta)
    write.table(myT,paste0(input,"BSMetagenomics_speciesMetadata.txt"),sep="\t",quote = FALSE)
  }
  
  
  pval<-vector()
  p1M<-vector()
  p6M<-vector()
  p1Y<-vector()
  s1M<-vector()
  s6M<-vector()
  s1Y<-vector()
  bugName<-vector()
  index<-1
  
  for (i in 1:ncol(kraken)){
    
    bug<-kraken[,i]
    
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
      
      bugName[index]<-colnames(kraken)[i]
      index<-index+1
    }
  }
  
  df<-data.frame(bugName,pval,p1M,p6M,p1Y,s1M,s6M,s1Y)
  df<-df[order(df$pval),]
  df$Adjustedpval<-p.adjust(df$pval,method = "BH")
  df$Adjustedp1M<-p.adjust(df$p1M,method = "BH")
  df$Adjustedp6M<-p.adjust(df$p6M,method = "BH")
  df$Adjustedp1Y<-p.adjust(df$p1Y,method = "BH")
  
  write.table(df,paste0(output,t,"_BS_Kraken2_Metagenomics_MixedLinearModelResults.txt"),sep="\t",row.names = FALSE,quote = FALSE)
}
