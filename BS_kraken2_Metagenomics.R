#Author: Farnaz Fouladi
#Date: 08-17-2020
#Description: Compare the gut microbiome at each time point versus baseline.

rm(list=ls())

#Libraries
library(nlme)

input<-"./input/RYGB_BS_Metagenomics/"
output<-"./output/MixedLinearModels/Kraken2/"
taxa<-c("Phylum","Class","Order","Family","Genus","Species")


for (t in c("Phylum","Class","Order","Family","Genus","Species")){
  
  
  kraken<-read.table(paste0(input,t,"_revisednames_Jan232020_taxaCount_norm_Log10.txt"),sep="\t",comment.char = "",check.names = FALSE,quote = "")
  
  kraken$SampleType<-sapply(rownames(kraken),function(x){
    if (x=="BS_BLANK" | x=="BS_DNA_BLANK" | x=="PL1H2O" | x=="PL2H2O" ) return("Blank")
    else if (x=="MMCC1" | x=="PL1MMC" | x=="PL2MMC" ) return("MMC")
    else return("Sample")
  })
  
  kraken$Run<-sapply(rownames(kraken),function(x){
    if(substr(x,1,2)=="BS" | x=="MMCC1") return("Run1")
    else return("Run2")
  })
  
  
  kraken$Site<-factor(sapply(as.character(kraken$patientID),function(x){
    if (x!= "positive.control" & x!= "water.blank" & x!= "DNA.iso.blank"){
      
      if (strsplit(x,"-")[[1]][2]=="1") return("1") 
      else if (strsplit(x,"-")[[1]][2]=="2") return("2")
      else return(NA)
    } else {
      return(NA)
    }
    
  }))
  
  kraken$timePoint<-factor(kraken$timePoint)
  
  #Removing BIO-2-014 (duplicates)
  kraken<-kraken[kraken$patientID!="BIO-2-014-0",]
  #Remove controls 
  myT1<-kraken[!is.na(kraken$Site),]
  #Selecting RYGB surgery
  surgeryType<-read.table(paste0(input,"TypeDateofSurgery-BiobehavioralR016-12-2020_FF.txt"),sep="\t",header = TRUE)
  myT1$patientID<-sapply(as.character(myT1$patientID),function(x){substr(x,1,9)})
  myT1<-merge(myT1,surgeryType,by.x = "patientID",by.y = "Patient.ID",all.x = TRUE,all.y = FALSE)
  myT1_RYGB<-myT1[myT1$Type.of.Surgery==1 & !is.na(myT1$Type.of.Surgery),]
  #Removing samples at 18 months after surgery (n=2)
  myT1_RYGB<-myT1_RYGB[myT1_RYGB$timePoint!=18,]
  #Removing samples that do not have baseline samples
  myT1_RYGB<-myT1_RYGB[myT1_RYGB$patientID %in% (myT1_RYGB$patientID[myT1_RYGB$timePoint==0]),]
  
  pval<-vector()
  p1M<-vector()
  p6M<-vector()
  p1Y<-vector()
  s1M<-vector()
  s6M<-vector()
  s1Y<-vector()
  bugName<-vector()
  index<-1
  finishAbundanceIndex<-which(colnames(myT1_RYGB)=="timePoint")-1
  myT2<-myT1_RYGB[,2:finishAbundanceIndex]
  meta<-myT1_RYGB[,(finishAbundanceIndex+1):ncol(myT1_RYGB)]
  
  meta$ID<-myT1_RYGB$patientID
  meta$time<-meta$timePoint
  
  if(t=="Species"){
    myT3<-cbind(myT2,meta)
    write.table(myT3,paste0(input,"BSMetagenomics_speciesMetadata.txt"),sep="\t",quote = FALSE)
  }
  
  for (i in 1:ncol(myT2)){
    
    bug<-myT2[,i]
    
    if (mean(bug>0)>0.1){
      
      df<-data.frame(bug,meta)
      
      fit<-anova(lme(bug~time,method="REML",random=~1|ID,data=df))
      sm<-summary(lme(bug~time,method="REML",random=~1|ID,data=df))
      
      pval[index]<-fit$`p-value`[2]
      p1M[index]<-sm$tTable[2,5]
      p6M[index]<-sm$tTable[3,5]
      p1Y[index]<-sm$tTable[4,5]
      
      s1M[index]<-sm$tTable[2,1]
      s6M[index]<-sm$tTable[3,1]
      s1Y[index]<-sm$tTable[4,1]
      
      bugName[index]<-colnames(myT2)[i]
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
