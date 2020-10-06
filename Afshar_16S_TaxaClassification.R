#Author: Farnaz Fouladi
#Date: 10-01-2020
#Description: This script generates taxonomic tables

rm(list=ls())

path<-"./input/RYGB_Afshar2018/DADA2/"
metaData<-"./input/RYGB_Afshar2018/Metadata/metaData.txt"
source("./Rcode/RYGB_IntegratedAnalysis/functions.R")

dada<-read.table(paste0(path,"ForwardReads.txt"),sep="\t",header=TRUE)
taxa<-read.table(paste0(path,"taxForwardReads.txt"),sep="\t",header=TRUE)
meta<-read.table(metaData,sep="\t",header=TRUE,row.names = 1)

#Modify Escherichia name
taxa$Genus<-sapply(taxa$Genus,function(x){
  if (x=="Esherichica/Shigella" & !is.na(x)) return("Escherichia/Shigella")
  else return(x)})

meta$time<-sapply(as.character(meta$title),function(x){
  if (substr(x,1,1)=="P") return(strsplit(x,"_")[[1]][2])
  else return(NA)
})

meta$ID<-sapply(as.character(meta$title),function(x){
  if (substr(x,1,1)=="P") return(strsplit(x,"_")[[1]][1])
  else return(NA)
})

dada<-dada[match(rownames(meta),rownames(dada)),]

#Remove controls
dada<-dada[!is.na(meta$ID),]
meta1<-meta[!is.na(meta$ID),]

meta1$time<-relevel(factor(meta1$time),ref="Pre")
meta1$prepost<-sapply(meta1$time,function(x){if (x=="Pre") return(0) else return(1)})

for(t in c("Phylum","Class","Order","Family","Genus")){
  
  t1<-getTaxaTable(dada,taxa,t)
  t1_norm<-norm(t1)
  t1_normMeta<-cbind(t1_norm,meta1)
  write.table(t1_normMeta,paste0(path,t,"_norm_table.txt"),sep = "\t",row.names = TRUE,quote = FALSE)
}

#SV table
dada1<-norm(dada)
num<-c(1:nrow(taxa))
taxanomy<-apply(taxa,1,function(x){paste0(x[1],"_",x[2],"_",x[3],"_",x[4],"_",x[5],"_",x[6])})
taxanomy<-paste0(taxanomy,"_",num)
colnames(dada1)<-taxanomy
dada1_meta<-cbind(dada1,meta1)
write.table(dada1_meta,paste0(path,"SV_norm_table.txt"),sep="\t",row.names = TRUE,quote = FALSE)
