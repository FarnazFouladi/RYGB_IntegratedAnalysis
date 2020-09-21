#Author: Farnaz Fouladi
#Date: 08-17-2020
#Description: Compare opportunistic pathogens between timepoints.

rm(list=ls())

#Libraries
library(ggplot2)
library(ggsignif)
library(reshape2)

input<-"./input/"
output<-"./output/"
taxa<-"Species"
bugs<-c("Klebsiella pneumoniae","Klebsiella oxytoca","Enterococcus faecalis","Enterobacter cloacae complex sp.",
        "Enterobacter cloacae", "Enterococcus faecium","Clostridium perfringens")

name="_speciesMetadata.txt"

########16S datasets
######BS study

study="BS"
myT<-read.table(paste0(input,"RYGB_BS/Kraken2/",study,name),sep="\t",header = TRUE,check.names = FALSE,comment.char = "",quote = "")

myT1<-myT[,colnames(myT) %in% bugs]
meta1<-data.frame(time=myT$Timepoint,ID=myT$Participant_ID,Sample=rownames(myT))

p<-matrix(NA,nrow = ncol(myT1),ncol=3,
          dimnames = list(colnames(myT1),c("01","06","16")))


for (bug in colnames(myT1)){
  if(mean(myT1[,bug]>0)>0.1){
    df<-data.frame(bug=myT1[,bug],meta1)
    df$time<-factor(df$time,levels = c(0,1,6))
    sm<-summary(lme(bug~time,method="REML",random=~1|ID,data=df))
    p[bug,1]<-sm$tTable[2,5]
    p[bug,2]<-sm$tTable[3,5]
    
    df$time<-factor(df$time,levels = c(1,0,6))
    sm<-summary(lme(bug~time,method="REML",random=~1|ID,data=df))
    p[bug,3]<-sm$tTable[3,5]
  }
}

padj<-matrix(p.adjust(as.numeric(p),method = "BH"),nrow =ncol(myT1),ncol=3,
             dimnames = list(colnames(myT1),c("01","06","16")))

write.table(padj,paste0(output,"BS_opportunisticPathogens.txt"),sep="\t")

######Assal study

study="Assal"
myT<-read.table(paste0(input,"RYGB_Assal2020/Kraken2/",study,name),sep="\t",header = TRUE,check.names = FALSE,comment.char = "",quote = "")
myT1<-myT[,colnames(myT) %in% bugs]

myT1<-myT[,colnames(myT) %in% bugs]
meta1<-data.frame(time=myT$time,ID=myT$ID,Sample=rownames(myT))

p<-matrix(NA,nrow =ncol(myT1),ncol=6,
          dimnames = list(colnames(myT1),
                          c("0-3M","0-1Y","0-2Y","3M-1Y","3M-2Y","1Y-2Y")))

for (bug in colnames(myT1)){
  
  if(mean(myT1[,bug]>0)>0.1){
    df<-data.frame(bug=myT1[,bug],meta1)
    df$time<-factor(df$time,levels = c("Pre","3M","1Y","2Y"))
    sm<-summary(lme(bug~time,method="REML",random=~1|ID,data=df))
    
    p[bug,1]<-sm$tTable[2,5]
    p[bug,2]<-sm$tTable[3,5]
    p[bug,3]<-sm$tTable[4,5]
  
    
    df$time<-factor(df$time,levels = c("3M","Pre","1Y","2Y"))
    sm<-summary(lme(bug~time,method="REML",random=~1|ID,data=df))
    p[bug,4]<-sm$tTable[3,5]
    p[bug,5]<-sm$tTable[4,5]

    df$time<-factor(df$time,levels = c("1Y","3M","Pre","2Y"))
    sm<-summary(lme(bug~time,method="REML",random=~1|ID,data=df))
    p[bug,6]<-sm$tTable[4,5]
    
  }
}

padj<-matrix(p.adjust(as.numeric(p),method = "BH"),nrow = ncol(myT1),ncol=6,
             dimnames = list(colnames(myT1),
             c("0-3M","0-1Y","0-2Y","3M-1Y","3M-2Y","1Y-2Y")))
write.table(padj,paste0(output,"Assal_opportunisticPathogens.txt"),sep="\t")

######Ilhan study

study="Ilhan"
myT<-read.table(paste0(input,"RYGB_Ilhan2020/Kraken2/",study,name),sep="\t",header = TRUE,check.names = FALSE,comment.char = "",quote = "")
myT1<-myT[,colnames(myT) %in% bugs]

myT1<-myT[,colnames(myT) %in% bugs]
meta1<-data.frame(time=myT$time,ID=myT$ID,Sample=rownames(myT))

p<-matrix(NA,nrow =ncol(myT1),ncol=3,
          dimnames = list(colnames(myT1),c("0-6M","0-12M","6M-12M")))

for (bug in colnames(myT1)){
  
  df<-data.frame(bug=myT1[,bug],meta1)
  df$time<-factor(df$time,levels = c("Baseline","6M","12M"))
  
  if(mean(myT1[,bug]>0)>0.1){
    sm<-summary(lme(bug~time,method="REML",random=~1|ID,data=df))
    
    p[bug,1]<-sm$tTable[2,5]
    p[bug,2]<-sm$tTable[3,5]
    
    df$time<-factor(df$time,levels = c("6M","Baseline","12M"))
    sm<-summary(lme(bug~time,method="REML",random=~1|ID,data=df))
    p[bug,3]<-sm$tTable[3,5]
  }
}

padj<-matrix(p.adjust(as.numeric(p),method = "BH"),nrow=ncol(myT1),ncol=3,
             dimnames = list(colnames(myT1),c("0-6M","0-12M","6M-12M")))

write.table(padj,paste0(output,"Ilhan_opportunisticPathogens.txt"),sep="\t")


######Afshar study

study="Afshar"
myT<-read.table(paste0(input,"RYGB_Afshar2018/Kraken2/",study,name),sep="\t",header = TRUE,check.names = FALSE,comment.char = "",quote = "")
myT1<-myT[,colnames(myT) %in% bugs]

myT1<-myT[,colnames(myT) %in% bugs]
meta1<-data.frame(time=myT$time,ID=myT$ID,Sample=rownames(myT))


p<-matrix(NA,nrow =ncol(myT1),ncol=1,
          dimnames = list(colnames(myT1),c("06")))

for (bug in colnames(myT1)){
  
  if(mean(myT1[,bug]>0)>0.1){
    df<-data.frame(bug=myT1[,bug],meta1)
    sm<-summary(lme(bug~time,method="REML",random=~1|ID,data=df))
    p[bug,1]<-sm$tTable[2,5]
    
  }
}

padj<-matrix(p.adjust(as.numeric(p),method = "BH"),nrow = ncol(myT1),ncol=1,
             dimnames = list(colnames(myT1),c("06")))

write.table(padj,paste0(output,"Afshar_opportunisticPathogens.txt"),sep="\t")

#######BSMetagenomics

study="BSMetagenomics"
myT<-read.table(paste0(input,"RYGB_BS_Metagenomics/",study,name),sep="\t",header = TRUE,check.names = FALSE,comment.char = "",quote = "")
myT1<-myT[,colnames(myT) %in% bugs]

myT1<-myT[,colnames(myT) %in% bugs]
meta1<-data.frame(time=myT$Timepoint,ID=myT$Participant_ID,Sample=rownames(myT))

p<-matrix(NA,nrow = ncol(myT1),ncol=6,
          dimnames = list(colnames(myT1),c("0-1","0-6","0-12","1-6","1-12","6-12")))

for (bug in colnames(myT1)){
  
  if(mean(myT1[,bug]>0)>0.1){
    df<-data.frame(bug=myT1[,bug],meta1)
    df$time<-factor(df$time,levels = c(0,1,6,12))
    sm<-summary(lme(bug~time,method="REML",random=~1|ID,data=df))
    
    p[bug,1]<-sm$tTable[2,5]
    p[bug,2]<-sm$tTable[3,5]
    p[bug,3]<-sm$tTable[4,5]
    
    df$time<-factor(df$time,levels = c(1,0,6,12))
    sm<-summary(lme(bug~time,method="REML",random=~1|ID,data=df))
    p[bug,4]<-sm$tTable[3,5]
    p[bug,5]<-sm$tTable[4,5]
    
    df$time<-factor(df$time,levels = c(6,0,1,12))
    sm<-summary(lme(bug~time,method="REML",random=~1|ID,data=df))
    p[bug,6]<-sm$tTable[4,5]
  }
}

padj<-matrix(p.adjust(as.numeric(p),method = "BH"),nrow =ncol(myT1),ncol=6,
             dimnames = list(colnames(myT1),c("0-1","0-6","0-12","1-6","1-12","6-12")))

write.table(padj,paste0(output,"BSMetagenomics_opportunisticPathogens.txt"),sep="\t")


#######PallejaMetagenomics

study="PallejaMetagenomics"
myT<-read.table(paste0(input,"palleja2016/",study,name),sep="\t",header = TRUE,check.names = FALSE,comment.char = "",quote = "")
myT1<-myT[,colnames(myT) %in% bugs]

myT1<-myT[,colnames(myT) %in% bugs]
meta1<-data.frame(time=myT$time,ID=myT$ID,Sample=rownames(myT))

p<-matrix(NA,nrow = ncol(myT1),ncol=3,
          dimnames = list(colnames(myT1),c("0-3","0-12","3-12")))

for (bug in colnames(myT1)){
  
  if(mean(myT1[,bug]>0)>0.1){
    df<-data.frame(bug=myT1[,bug],meta1)
    df$time<-factor(df$time,levels = c("Baseline","3MO","1Y"))
    sm<-summary(lme(bug~time,method="REML",random=~1|ID,data=df))
    
    p[bug,1]<-sm$tTable[2,5]
    p[bug,2]<-sm$tTable[3,5]
    
    df$time<-factor(df$time,levels = c("3MO","Baseline","1Y"))
    sm<-summary(lme(bug~time,method="REML",random=~1|ID,data=df))
    p[bug,3]<-sm$tTable[3,5]
  }
}

padj<-matrix(p.adjust(as.numeric(p),method = "BH"),nrow =ncol(myT1),ncol=3,
             dimnames = list(colnames(myT1),c("0-3","0-12","3-12")))

write.table(padj,paste0(output,"PallejaMetagenomics_opportunisticPathogens.txt"),sep="\t")

#####Plots
datasets<-c("BS","Assal","Afshar","Ilhan","BSMetagenomics","PallejaMetagenomics")
folders<-c("RYGB_BS/Kraken2/","RYGB_Assal2020/Kraken2/","RYGB_Afshar2018/Kraken2/",
           "RYGB_Ilhan2020/Kraken2/","RYGB_BS_Metagenomics/","palleja2016/")
myList<-list()
outdex<-1

for (study in datasets){
  
  myT<-read.table(paste0(input,folders[which(datasets==study)],study,name),
                  sep="\t",header = TRUE,check.names = FALSE,comment.char = "",quote = "")
  
  myT1<-myT[,colnames(myT) %in% bugs]
  
  if(study=="BS" | study=="BSMetagenomics"){
    myT1$time<-as.factor(myT$Timepoint)
    myT1$ID<-myT$Participant_ID
    myT1$Sample<-rownames(myT)
  }else{
    myT1$time<-as.factor(myT$time)
    myT1$ID<-myT$ID
    myT1$Sample<-rownames(myT)
  }
  
  
  if(study=="Assal")
    myT1$time<-factor(myT1$time,levels = c("Pre","3M","1Y","2Y"))
  else if (study=="Afshar")
    myT1$time<-factor(myT1$time,levels = c("Pre","Post"))
  else if (study=="Ilhan")
    myT1$time<-factor(myT1$time,levels = c("Baseline","6M","12M"))
  else if (study=="PallejaMetagenomics")
    myT1$time<-factor(myT1$time,levels = c("Baseline","3MO","1Y"))
  
  
  long<-melt(myT1,id.vars = c("Sample","time","ID"))
  
  
  col<-c("blue","darkorange2","darkgreen","hotpink","red","purple")
  col1<-col[long$time]
  
  theme_set(theme_classic(base_size = 16))
  plot<-ggplot(data=long, aes(x=variable,y=value))+
    geom_boxplot(aes(color = time), width = 0.5, size = 0.4,
                 position = position_dodge(0.8),outlier.shape = NA)+
    geom_jitter(aes(colour=time),shape=16, size=0.5,position = position_dodge(0.8))+
    labs(y=expression(Log[10] ~"normalized count"),x="")+scale_color_manual(values = col)+
    theme(axis.text.x = element_text(angle = 90))
  
  if(study=="BS")
    
    plot1<-plot+geom_signif(y_position=c(3.9,4.1,3.3,0.8), xmin=c(0.8,0.8,2.8,4.8), xmax=c(1,1.2,3,5),annotations = "",tip_length=0.03,vjust = 0.5,textsize =4)+
    labs(title = "BS-16S")+scale_color_manual(labels=c("Baseline","1 month","6 months"),values = col)+coord_flip()
  else if(study=="Assal")
    plot1<-plot+geom_signif(y_position=c(3.1,3.3,3.5,1.5), xmin=c(0.7,0.7,0.9,6.7), xmax=c(0.9,1.3,1.1,6.9),annotations = "",tip_length=0.03,vjust = 0.5,textsize =4)+
    coord_flip()+labs(title = "Assal-16S")+scale_color_manual(values = col,labels=c("Baseline","3 months","12 months","24 months"))
  else if (study=="Afshar")
    plot1<-plot+coord_flip()+labs(title = "Afshar-16S")+scale_color_manual(values = col,labels=c("Baseline","Post-surgery"))
  else if (study=="Ilhan")
    plot1<-plot+coord_flip()+labs(title = "Ilhan-16S")+scale_color_manual(values = col,labels=c("Baseline","6 months","12 months"))
  else if (study=="BSMetagenomics")
    plot1<-plot+geom_signif(y_position=c(6.1,4.6,4.8,4,4.2,4.4,6,6.2,6.4,5,5.2,5.4,6.7,6.8,7,6.7,6.8,7), xmin=c(0.7,1.7,1.7,2.7,2.7,2.7,3.7,3.7,3.8,4.7,4.7,4.7,5.7,5.7,5.7,6.7,6.7,6.7), xmax=c(1.05,1.9,2.3,2.9,3.1,3.2,4.05,4.2,4.05,4.9,5.1,5.2,5.9,6.1,6.2,6.9,7.1,7.2),annotations = "",tip_length=0.03,vjust = 0.5,textsize =4)+
    labs(title = "BS-Metagenomics")+coord_flip()+scale_color_manual(labels=c("Baseline","1 month","6 months","12 months"),values = col)
  else if (study=="PallejaMetagenomics")
    plot1<-plot+coord_flip()+scale_color_manual(labels=c("Baseline","3 months","12 months"),values = col)+
    geom_signif(y_position=c(4.2,4.5,4.2,4.5,4.6,6,6.3,7,7.3),xmin=c(1.7,1.7,2.7,2.7,3.7,5.7,5.7,6.7,6.7),xmax =c(2,2.2,3,3.2,4,6,6.2,7,7.2) ,annotations = "",tip_length=0.03,vjust = 0.5,textsize =4)+
    labs(title = "Palleja-Metagenomics")
  
  myList[[outdex]]<-plot1
  outdex<-outdex+1
  
}

pdf(paste0(output,"OportunisticPathogens.pdf"),width = 15,height = 7)
theme_set(theme_classic(base_size = 10))
grid.arrange(myList[[6]],myList[[5]],myList[[1]],myList[[2]],myList[[3]],myList[[4]],ncol=3,nrow=2)
dev.off()










