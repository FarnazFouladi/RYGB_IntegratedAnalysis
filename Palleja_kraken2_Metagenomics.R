#Author: Farnaz Fouladi
#Date: 08-17-2020
#Description: Compare the gut microbiome at each time point versus baseline using
#             mixed linear models and metagenomic data.

rm(list=ls())

#Libraries
library(nlme)

input<-"./input/palleja2016/"
output<-"./output/MixedLinearModels/Kraken2/"
taxa<-c("Phylum","Class","Order","Family","Genus","Species")


for (t in taxa){
  
  
  kraken<-read.table(paste0(input,"RYGB_PRJEB12947_2020Aug16_taxaCount_norm_Log10_",t,".tsv"),
                     sep="\t",comment.char = "",check.names = FALSE,quote = "",
                     header = TRUE,row.names = 1)
  
  meta<-read.table(paste0(input,"metaData.txt"),sep=",",header = TRUE)
  meta<-meta[meta$LibraryLayout=="PAIRED",]
  meta$ID<-sapply(meta$Title,function(x){strsplit(x,"_")[[1]][1]})
  meta$time<-factor(sapply(meta$Title,function(x){strsplit(x,"_")[[1]][2]}),levels = c("Baseline","3MO","1Y"))
  
  #Ordering count table and metadata
  
  myT2<-kraken[match(meta$Run,rownames(kraken)),]
  
  if(sum(rownames(myT2)==meta$Run)!=nrow(meta)) stop("Error!")
  
  pval<-vector()
  p3M<-vector()
  p1Y<-vector()
  s3M<-vector()
  s1Y<-vector()
  bugName<-vector()
  index<-1
  
  
  if(t=="Species"){
    myT3<-cbind(myT2,meta)
    write.table(myT3,paste0(input,"PallejaMetagenomics_speciesMetadata.txt"),sep="\t",quote = FALSE)
  }
  
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
  
  write.table(df,paste0(output,t,"_Palleja_Kraken2_Metagenomics_MixedLinearModelResults.txt"),sep="\t",row.names = FALSE,quote = FALSE)
  
}
