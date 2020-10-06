#Author: Farnaz Fouladi
#Date: 10-01-2020
#Description: This script generates taxonomic tables

rm(list=ls())

path<-"./input/RYGB_BS/DADA2/"
metaData<-"./input/RYGB_BS/Metadata/16S_metadata.txt"
source("./Rcode/RYGB_IntegratedAnalysis/functions.R")

dada<-read.table(paste0(path,"ForwardReads.txt"),sep="\t",header=TRUE)
taxa<-read.table(paste0(path,"taxForwardReads.txt"),sep="\t",header=TRUE)
meta<-read.table(metaData,sep="\t",header=TRUE,row.names = 1)

#Modify Escherichia name
taxa$Genus<-sapply(taxa$Genus,function(x){
  if (x=="Esherichica/Shigella" & !is.na(x)) return("Escherichia/Shigella")
  else return(x)})


dada<-dada[match(rownames(meta),rownames(dada)),]

#Write normalized taxonomic tables 
for (t in c("Phylum","Class","Order","Family","Genus")){
  
  t1<-getTaxaTable(dada,taxa,t)
  t1_norm<-norm(t1)
  t1_normMeta<-cbind(t1_norm,meta)
  write.table(t1_normMeta,paste0(path,t,"_norm_table.txt"),sep = "\t",row.names = TRUE,quote = FALSE)
}

#Write normalized sequence variant tables 
dada1<-norm(dada)
dada1_meta<-cbind(dada1,meta)
write.table(dada1_meta,paste0(path,"Seq_norm_table.txt"),sep="\t",row.names = TRUE,quote = FALSE)


#Write normalized sequence variant tables with their taxonomic classifications
num<-c(1:nrow(taxa))
taxanomy<-apply(taxa,1,function(x){paste0(x[1],"_",x[2],"_",x[3],"_",x[4],"_",x[5],"_",x[6])})
taxanomy<-paste0(taxanomy,"_",num)
colnames(dada1)<-taxanomy
dada1_meta<-cbind(dada1,meta)
write.table(dada1_meta,paste0(path,"SV_norm_table.txt"),sep="\t",row.names = TRUE,quote = FALSE)

#None-Normalized taxonomic tables
for (t in c("Phylum","Class","Order","Family","Genus")){
  
  t1<-getTaxaTable(dada,taxa,t)
  t1_Meta<-cbind(t1,meta)
  write.table(t1_Meta,paste0(path,t,"_count_table.txt"),sep = "\t",row.names = TRUE,quote = FALSE)
}

