#Author: Farnaz Fouladi
#Date: 08-17-2020
#Description: Ilhan dataset:Compare the gut microbiome at each time point versus
#              baseline using mixed linear models. Count tables are classified by 
#             DADA2.


rm(list=ls())

#Libraries
library(nlme)

input<-"./input/RYGB_Ilhan2020/"
output<-"./output/MixedLinearModels/"
taxa<-c("Phylum","Class","Order","Family","Genus","SV")

for (t in taxa){
  
  dada2<-read.table(paste0(input,t,"_norm_table.txt"),sep="\t",header = TRUE,row.names = 1,check.names = FALSE)
  
  dada2_fecalRYGB<-dada2[dada2$env_material=="fecal" & (dada2$Group=="Baseline" |dada2$Group=="6M" | dada2$Group=="12M"), ]
  dada2_fecalRYGB$Group<-factor(dada2_fecalRYGB$Group,levels = c("Baseline","6M","12M"))
  dada2_fecalRYGB$prepost<-sapply(dada2_fecalRYGB$Group,function(x){if (x=="Baseline") return(0) else return(1)})
  
  if(t=="Genus")
    colnames(dada2_fecalRYGB)[colnames(dada2_fecalRYGB)=="Esherichica/Shigella"]<-"Escherichia/Shigella"
  
  #write table
  write.table(dada2_fecalRYGB,paste0(input,t,"_norm_table_updated.txt"),sep="\t",quote = FALSE)

  finishAbundanceIndex<-which(colnames(dada2_fecalRYGB)=="Assay.Type")-1
  myT<-dada2_fecalRYGB[,1:finishAbundanceIndex]
  meta<-dada2_fecalRYGB[,(finishAbundanceIndex+1):ncol(dada2_fecalRYGB)]
  
  pval<-vector()
  p6M<-vector()
  p1Y<-vector()
  s6M<-vector()
  s1Y<-vector()
  bugName<-vector()
  index<-1
  
  for (i in 1:ncol(myT)){
    
    bug<-myT[,i]
    
    if (mean(bug>0)>0.1){
      
      df<-data.frame(bug,meta)
      
      fit<-anova(lme(bug~Group,method="REML",random=~1|ID,data=df))
      sm<-summary(lme(bug~Group,method="REML",random=~1|ID,data=df))
      
      pval[index]<-fit$`p-value`[2]
      p6M[index]<-sm$tTable[2,5]
      p1Y[index]<-sm$tTable[3,5]
      s6M[index]<-sm$tTable[2,1]
      s1Y[index]<-sm$tTable[3,1]
      bugName[index]<-colnames(myT)[i]
      index<-index+1
    }
  }
  
  
  df<-data.frame(bugName,pval,p6M,p1Y,s6M,s1Y)
  df<-df[order(df$pval),]
  df$Adjustedpval<-p.adjust(df$pval,method = "BH")
  df$Adjustedp6M<-p.adjust(df$p6M,method = "BH")
  df$Adjustedp1Y<-p.adjust(df$p1Y,method = "BH")
 
  write.table(df,paste0(output,t,"_Ilhan_MixedLinearModelResults.txt"),sep="\t",row.names = FALSE,quote = FALSE)
}  
