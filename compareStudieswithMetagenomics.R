#Author: Farnaz Fouladi
#Date: 08-17-2020
#Description: p-value versus p-value plots for 16S datasets

rm(list=ls())

#Libraries
library(ggplot2)
library(gridExtra)
library(ggsignif)
library(ggrepel)

output<-"./output/"
taxa<-"Genus"
source("./Rcode/RYGB_IntegratedAnalysis/functions.R")


#Metagenomics BS 1M 6M 1Y

#BS 1M 6M
#Assal 3M 1Y 2Y
#Ilhan 6M 1Y
#Afshar 6M


studies<-c("BS_Kraken2-p1M","BS_Kraken2-p6M","Assal_Kraken2-p3M","Assal_Kraken2-p1Y","Assal_Kraken2-p2Y","Ilhan_Kraken2-p6M","Ilhan_Kraken2-p1Y","Afshar_Kraken2-p6M")
studyNames<-c("BS-1 month","BS-6 months","Assal-3 months","Assal-1 year","Assal-2 years","Ilhan-6 months","Ilhan-1 year","Afshar-Post surgery")
metagenomicStudies<-c("BS_Kraken2_Metagenomics-p1M","BS_Kraken2_Metagenomics-p6M","BS_Kraken2_Metagenomics-p1Y")
metagenomicStudyNames<-c("BS-1 month","BS-6 months","BS-1 year")

r<-vector()
pval<-vector()
plotList<-list()
studyPairs<-vector()
compariosn<-vector()
compariosnTime<-vector()
index<-1

for (s in 1:length(metagenomicStudies)){
  
  for (s1 in 1:length(studies)){
    
    path<-paste0(output,"MixedLinearModels/Kraken2/")
    df<-compareStudies(path,taxa,strsplit(studies[s1],"-")[[1]][1],strsplit(metagenomicStudies[s],"-")[[1]][1],strsplit(studies[s1],"-")[[1]][2],strsplit(metagenomicStudies[s],"-")[[1]][2])
    r[index]<-correlationBetweenStudies(df)[[1]]
    pval[index]<-correlationBetweenStudies(df)[[2]]
    index<-index+1
    
  }
}
#Adjust p-values
pval<-p.adjust(pval,method = "BH")
count<-0

for (s in 1:length(metagenomicStudies)){
  
  for (s1 in 1:length(studies)){
    
    
    count<-count+1
    df<-compareStudies(path,taxa,strsplit(studies[s1],"-")[[1]][1],strsplit(metagenomicStudies[s],"-")[[1]][1],strsplit(studies[s1],"-")[[1]][2],strsplit(metagenomicStudies[s],"-")[[1]][2])
    plot<-plotPairwiseStudiesMetagenomics(df,paste0("16S ",studyNames[s1]),paste0("Metagenomics ",metagenomicStudyNames[s]),r[count],pval[count])
    plotList[[count]]<-plot
    studyPairs[count]<-paste0(metagenomicStudies[s],"_",studies[s1])
    if(strsplit(metagenomicStudies[s],"_")[[1]][1] == strsplit(studies[s1],"_")[[1]][1])
      compariosn[[count]]<-"Same study"
    else
      compariosn[[count]]<-"Different study"
    
    if(strsplit(metagenomicStudies[s],"-")[[1]][2] == strsplit(studies[s1],"-")[[1]][2])
      compariosnTime[[count]]<-"Same timepoint"
    else
      compariosnTime[[count]]<-"Different timepoint"
    
  }
}


df<-data.frame(studyPairs,compariosn,compariosnTime,pval,r)
write.table(df, paste0(output,taxa,"_PairwiseComaprison_Metagenomics.txt"),sep="\t",row.names = FALSE)

p1<-wilcox.test(df$r[df$compariosn=="Different study" & df$compariosnTime=="Different timepoint"],
                df$r[df$compariosn=="Different study" & df$compariosnTime=="Same timepoint"])$p.value

p2<-wilcox.test(df$r[df$compariosn=="Different study" & df$compariosnTime=="Different timepoint"],
                df$r[df$compariosn=="Same study" & df$compariosnTime=="Different timepoint"])$p.value


p3<-wilcox.test(df$r[df$compariosn=="Different study" & df$compariosnTime=="Different timepoint"],
                df$r[df$compariosn=="Same study" & df$compariosnTime=="Same timepoint"])$p.value

p4<-wilcox.test(df$r[df$compariosn=="Different study" & df$compariosnTime=="Same timepoint"],
                df$r[df$compariosn=="Same study" & df$compariosnTime=="Different timepoint"])$p.value

p5<-wilcox.test(df$r[df$compariosn=="Different study" & df$compariosnTime=="Same timepoint"],
                df$r[df$compariosn=="Same study" & df$compariosnTime=="Same timepoint"])$p.value

p6<-wilcox.test(df$r[df$compariosn=="Same study" & df$compariosnTime=="Different timepoint"],
                df$r[df$compariosn=="Same study" & df$compariosnTime=="Same timepoint"])$p.value

p7<-wilcox.test(df$r[df$compariosn=="Different study"],
                df$r[df$compariosn=="Same study"])$p.value

p.adjust(c(p1,p2,p3,p4,p5,p6,p7),method = "BH")

plot1<-ggplot(data=df,aes(x=compariosn,y=r))+
  geom_boxplot(aes(color=compariosnTime),position = position_dodge(0.6), width = 0.5, size = 0.4,outlier.shape = NA)+
  geom_jitter(aes(color=compariosnTime),size=0.8,position = position_dodge(0.6))+labs(x="",y="Spearman Coefficient",col="")+
  geom_signif(y_position = c(0.6,0.65),xmin = c(0.8,0.8),xmax=c(1.8,2.2),annotations = c("*","*"),tip_length=0.03,vjust = 0.5,textsize =4)



pdf(paste0(output,taxa,"_scatterPlots_ComparisonWithMetagenomics.pdf"),width = 10,height = 10)
theme_set(theme_classic(base_size = 9))
grid.arrange(plotList[[1]],plotList[[2]],plotList[[3]],
             plotList[[4]],plotList[[5]],plotList[[6]],
             plotList[[7]],plotList[[8]],ncol=3,nrow=3)

grid.arrange(plotList[[9]],plotList[[10]],plotList[[11]],
             plotList[[12]],plotList[[13]],plotList[[14]],
             plotList[[15]],plotList[[16]],ncol=3,nrow=3)

grid.arrange(plotList[[17]],plotList[[18]],plotList[[19]],
             plotList[[20]],plotList[[21]],plotList[[22]],
             plotList[[23]],plotList[[24]],ncol=3,nrow=3)


dev.off()

pdf(paste0(output,taxa,"_coefficientsFromScatterPlots_Metagenomics.pdf"),width = 7,height = 7)
theme_set(theme_classic(base_size = 18))
print(plot1)
dev.off()



