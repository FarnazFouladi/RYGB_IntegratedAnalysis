#Author: Farnaz Fouladi
#Date: 08-17-2020
#Description: p-value versus p-value plots for 16S datasets (Sequence variants)

rm(list=ls())

#Libraries
library(ggplot2)
library(gridExtra)
library(ggsignif)
library(ggrepel)

output<-"./output/"
taxa<-"Seq"
source("./Rcode/RYGB_IntegratedAnalysis/functions.R")


#BS 1M 6M
#Assal 3M 1Y 2Y
#Ilhan 6M 1Y
#Afshar 6M

studies<-c("BS-p1M","BS-p6M","Assal-p3M","Assal-p1Y","Assal-p2Y")
studyNames<-c("BS-1 month","BS-6 months","Assal-3 months","Assal-1 year","Assal-2 years")

r<-vector()
pval<-vector()
plotList<-list()
studyPairs<-vector()
compariosn<-vector()
compariosnTime<-vector()
index<-1

for (s in 1:length(studies)){
  
  if (s!=length(studies)){
    
    otherStudies<-c((s+1):length(studies))
    
    for (s1 in otherStudies){
      
      path<-paste0(output,"MixedLinearModels/")
      df<-compareStudies(path,taxa,strsplit(studies[s],"-")[[1]][1],strsplit(studies[s1],"-")[[1]][1],strsplit(studies[s],"-")[[1]][2],strsplit(studies[s1],"-")[[1]][2])
      
      r[index]<-correlationBetweenStudies(df)[[1]]
      pval[index]<-correlationBetweenStudies(df)[[2]]
      index<-index+1
      
    }
  }
}

pval<-p.adjust(pval,method = "BH")
count=0

for (s in 1:length(studies)){
  
  if (s!=length(studies)){
    otherStudies<-c((s+1):length(studies))
    
    for (s1 in otherStudies){
      
      count<-count+1
      df<-compareStudies(path,taxa,strsplit(studies[s],"-")[[1]][1],strsplit(studies[s1],"-")[[1]][1],strsplit(studies[s],"-")[[1]][2],strsplit(studies[s1],"-")[[1]][2])
      xlab=paste0(studyNames[s]," vs. baseline")
      ylab=paste0(studyNames[s1]," vs. baseline")
      plot<-plotPairwiseStudiesMetagenomics(df,xlab,ylab,r[count],pval[count])
      plotList[[count]]<-plot
      studyPairs[count]<-paste0(studies[s],"_",studies[s1])
      if(strsplit(studies[s],"-")[[1]][1] == strsplit(studies[s1],"-")[[1]][1])
        compariosn[[count]]<-"Same study"
      else
        compariosn[[count]]<-"Different study"
      
      if(strsplit(studies[s],"-")[[1]][2] == strsplit(studies[s1],"-")[[1]][2])
        compariosnTime[[count]]<-"Same timepoint"
      else
        compariosnTime[[count]]<-"Different timepoint"
    }
  }
  
}


df<-data.frame(studyPairs,compariosn,compariosnTime,pval,r)
write.table(df, paste0(output,taxa,"_PairwiseComaprison.txt"),sep="\t",row.names = FALSE)

#Wilcoxon test for comparing coefficients
p2<-wilcox.test(df$r[df$compariosn=="Different study" & df$compariosnTime=="Different timepoint"],
                df$r[df$compariosn=="Same study" & df$compariosnTime=="Different timepoint"])

#No significant difference between different time points
plot1<-ggplot(data=df,aes(x=compariosn,y=r))+
  geom_boxplot(aes(color=compariosnTime),position = position_dodge(0.6), width = 0.5, size = 0.4,outlier.shape = NA)+
  geom_jitter(aes(color=compariosnTime),size=0.8,position = position_dodge(0.6))+labs(x="",y="Spearman Coefficient",color="")

pdf(paste0(output,taxa,"_scatterPlots.pdf"),width = 10,height = 10)
theme_set(theme_classic(base_size = 9))
i=1
grid.arrange(plotList[[i]],plotList[[i+1]],plotList[[i+2]],
               plotList[[i+3]],plotList[[i+4]],plotList[[i+5]],
               plotList[[i+6]],plotList[[i+7]],plotList[[i+8]],ncol=3,nrow=3)

grid.arrange(plotList[[10]],ncol=3,nrow=3)
dev.off()


pdf(paste0(output,taxa,"_coefficientsFromScatterPlots.pdf"),width = 5,height = 5)
theme_set(theme_classic(base_size = 14))
print(plot1)
dev.off()

df_g<-read.table(paste0(output,"Genus_PairwiseComaprison.txt"),sep="\t",header = TRUE)
df_g$taxanomy<-rep("Genus",nrow(df_g))
df$taxanomy<-rep("SV",nrow(df))
df_all<-rbind(df_g,df)
df_all$study1<-sapply(df_all$studyPairs,function(x){
  strsplit(x,"-")[[1]][1]})


df_all$study2<-sapply(sapply(df_all$studyPairs,function(x){
  strsplit(x,"_")[[1]][2]}),function(i){
    strsplit(i,"-")[[1]][1]
  })

df_sub<-df_all[(df_all$study1=="BS" | df_all$study1=="Assal") & (df_all$study2=="BS" | df_all$study2=="Assal"),  ]
plot2<-ggplot(data=df_sub,aes(x=compariosn,y=r))+
  geom_boxplot(aes(color=taxanomy),position = position_dodge(0.6), width = 0.5, size = 0.4,outlier.shape = NA)+
  geom_jitter(aes(color=taxanomy),size=0.8,position = position_dodge(0.6))+labs(x="",y="Spearman Coefficient",color="")

pdf(paste0(output,"coefficientsFromScatterPlots_compareSVandGenus.pdf"),width = 5,height = 5)
theme_set(theme_classic(base_size = 14))
print(plot2)
dev.off()


p1<-wilcox.test(df_sub$r[df_sub$compariosn=="Different study" & df_sub$taxanomy=="Genus"],
                df_sub$r[df_sub$compariosn=="Different study" & df_sub$taxanomy=="SV"])


p2<-wilcox.test(df_sub$r[df_sub$compariosn=="Different study" & df_sub$taxanomy=="Genus"],
                df_sub$r[df_sub$compariosn=="Same study" & df_sub$taxanomy=="Genus"])

p3<-wilcox.test(df_sub$r[df_sub$compariosn=="Different study" & df_sub$taxanomy=="Genus"],
                df_sub$r[df_sub$compariosn=="Same study" & df_sub$taxanomy=="SV"])

p4<-wilcox.test(df_sub$r[df_sub$compariosn=="Different study" & df_sub$taxanomy=="SV"],
                df_sub$r[df_sub$compariosn=="Same study" & df_sub$taxanomy=="Genus"])

p5<-wilcox.test(df_sub$r[df_sub$compariosn=="Different study" & df_sub$taxanomy=="SV"],
                df_sub$r[df_sub$compariosn=="Same study" & df_sub$taxanomy=="SV"])

p6<-wilcox.test(df_sub$r[df_sub$compariosn=="Same study" & df_sub$taxanomy=="Genus"],
                df_sub$r[df_sub$compariosn=="Same study" & df_sub$taxanomy=="SV"])

adjustedp<-p.adjust(c(p1$p.value,p2$p.value,p3$p.value,p4$p.value,p5$p.value,p6$p.value),method = "BH")

