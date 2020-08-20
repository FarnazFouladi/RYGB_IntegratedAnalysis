#Author: Farnaz Fouladi
#Date: 08-17-2020
#Description: Compare the gut microbiome at each time point versus baseline.

rm(list=ls())

#Libraries
library(nlme)

input<-"./input/RYGB_BS/Kraken2/"
output<-"./output/MixedLinearModels/Kraken2/"
taxa<-c("Phylum","Class","Order","Family","Genus","Species")


for (t in taxa){
  
  myT<-read.table(paste0(input,"RYGB_BS_2020Apr23_taxaCount_",t,".tsv"),
                  sep="\t",header = TRUE,row.names = 1,check.names = FALSE)
  
  rownames(myT)<-sapply(rownames(myT),function(x){strsplit(x,"_")[[1]][1]})
  rownames(myT)[c(103,104)]<-c("MMCC1","MMCC2")
  
  meta<-read.table(paste0(input,"2019.04.29_BS_16S_mapping_file_FF.txt"),sep="\t",header = TRUE)
  
  if(sum(rownames(myT)==meta$Sample.ID)!=nrow(meta)) stop("Error")
  
  #Remove duplicate samples
  myT<-myT[!duplicated(meta$Patient.ID),]
  meta<-meta[!duplicated(meta$Patient.ID),]
  
  #Remove low sequence depth and Normalize data
  
  myT1<-myT[rowSums(myT)>1000,]
  meta1<-meta[rowSums(myT)>1000,]
  average<-mean(rowSums(myT1))
  myT_relab<-sweep(myT1,1,rowSums(myT1),"/")
  myT_norm<-log10(myT_relab*average+1)
  
  surgeryType<-read.table("./input/RYGB_BS/TypeDateofSurgery-BiobehavioralR016-12-2020_FF.txt",sep="\t",header = TRUE)
  meta1$Patient.ID<-sapply(as.character(meta1$Patient.ID),function(x){substr(x,1,9)})
  meta_merge<-merge(meta1,surgeryType,by.x = "Patient.ID",by.y = "Patient.ID",all.x = TRUE,all.y = FALSE,sort = FALSE)
  meta_merge<-meta_merge[match(rownames(myT1),meta_merge$Sample.ID),]
  
  #Selecting RYGB surgeries
  meta_merge_rygb<-meta_merge[meta_merge$Type.of.Surgery==1 & !is.na(meta_merge$Type.of.Surgery),]
  meta_merge_rygb$prepost<-sapply(meta_merge_rygb$Timepoint,function(x){if (x==0) return(0) else return(1)})
  myT_rygb<-myT_norm[meta_merge$Type.of.Surgery==1 & !is.na(meta_merge$Type.of.Surgery),]
  
  #Removing samples that do not have baselines
  meta_merge_rygb1<-meta_merge_rygb[meta_merge_rygb$Patient.ID %in% (meta_merge_rygb$Patient.ID[meta_merge_rygb$Timepoint==0]),]
  myT_rygb1<-myT_rygb[meta_merge_rygb$Patient.ID %in% (meta_merge_rygb$Patient.ID[meta_merge_rygb$Timepoint==0]),]
  
  meta1<-meta_merge_rygb1
  myT1<-myT_rygb1
  
  if(t=="Species"){
    
    myT1$time<-meta1$Timepoint
    myT1$ID<-meta1$Patient.ID
    myT1$Sample<-meta1$Sample.ID
    write.table(myT1,paste0(input,"BS_speciesMetadata.txt"),sep="\t",quote = FALSE)
  }
 
  pval<-vector()
  p1M<-vector()
  p6M<-vector()
  s1M<-vector()
  s6M<-vector()
  bugName<-vector()
  index<-1
  
  for (i in 1:ncol(myT1)){
    
    bug<-myT1[,i]
    
    if (mean(bug>0)>0.1){
      
      df<-data.frame(bug,meta1)
      
      fit<-tryCatch({
        lme(bug~factor(Timepoint),method="REML",random=~1|Participant_ID,data=df)
      },
      error=function(e){cat("ERROR :",conditionMessage(e), "\n")
        return(NA)})
      
      if (is.na(fit)){
        
        pval[index]<-NA
        p1M[index]<-NA
        p6M[index]<-NA
        
        s1M[index]<-NA
        s6M[index]<-NA
        
        
      }else{
        
        fit<-anova(lme(bug~factor(Timepoint),method="REML",random=~1|Participant_ID,data=df))
        sm<-summary(lme(bug~factor(Timepoint),method="REML",random=~1|Participant_ID,data=df))
        
        pval[index]<-fit$`p-value`[2]
        p1M[index]<-sm$tTable[2,5]
        p6M[index]<-sm$tTable[3,5]
        
        s1M[index]<-sm$tTable[2,1]
        s6M[index]<-sm$tTable[3,1]
      }
      bugName[index]<-colnames(myT1)[i]
      index<-index+1
      
    }
  }
  
  df<-data.frame(bugName,pval,p1M,p6M,s1M,s6M)
  df<-df[order(df$pval),]
  df$Adjustedpval<-p.adjust(df$pval,method = "BH")
  df$Adjustedp1M<-p.adjust(df$p1M,method = "BH")
  df$Adjustedp6M<-p.adjust(df$p6M,method = "BH")
  
  write.table(df,paste0(output,t,"_BS_Kraken2_MixedLinearModelResults.txt"),sep="\t",row.names = FALSE,quote = FALSE)
}

