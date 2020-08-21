#Author: Farnaz Fouladi
#Date: 08-17-2020
#Description: Join table across studies.

rm(list=ls())

taxa<-"Genus"
pathToBS<-"./input/RYGB_BS/"
pathToAssal<-"./input/RYGB_Assal2020/"
pathToIlhan<-"./input/RYGB_Ilhan2020/"
pathToAfshar<-"./input/RYGB_Afshar2018/"
source("./Rcode/RYGB_IntegratedAnalysis/functions.R")
output<-"./output/"

normalized="log10"
normalized="relab"


###BS study
myT.BS<-read.table(paste0(pathToBS,taxa,"_norm_table_updated.txt"),sep="\t",header = TRUE)
myT.BS<-myT.BS[,-1] #Removing patient ID column

#### Assal study
myT.Assal<-read.table(paste0(pathToAssal,taxa,"_norm_table_updated.txt"),sep="\t",header = TRUE)

####Ilhan study
myT.Ilhan<-read.table(paste0(pathToIlhan,taxa,"_norm_table_updated.txt"),sep="\t",header = TRUE)

###Afshar study
myT.Afshar<-read.table(paste0(pathToAfshar,taxa,"_norm_table_updated.txt"),sep="\t",header = TRUE)


if(normalized=="relab"){
  
  #Convert log10 normalized count to relative abundance
  myT.BS_normalized<-getRelativeAbundance(myT.BS)
  myT.Assal_normalized<-getRelativeAbundance(myT.Assal)
  myT.Ilhan_normalized<-getRelativeAbundance(myT.Ilhan)
  myT.Afshar_normalized<-getRelativeAbundance(myT.Afshar)
  
} else {
  
  myT.BS_normalized<-getNormalizedCountTable(myT.BS)
  myT.Assal_normalized<-getNormalizedCountTable(myT.Assal)
  myT.Ilhan_normalized<-getNormalizedCountTable(myT.Ilhan)
  myT.Afshar_normalized<-getNormalizedCountTable(myT.Afshar)
  
}
meta.BS<-getMetaData(myT.BS)
meta.Assal<-getMetaData(myT.Assal)
meta.Ilhan<-getMetaData(myT.Ilhan)
meta.Afshar<-getMetaData(myT.Afshar)

#Get a list of all taxa in all count tables
bugs<-list()
bugs<-getListOfBug(myT.BS_normalized,"BS",bugs)
bugs<-getListOfBug(myT.Assal_normalized,"Assal",bugs)
bugs<-getListOfBug(myT.Ilhan_normalized,"Ilhan",bugs)
bugs<-getListOfBug(myT.Afshar_normalized,"Afshar",bugs)
  
namesOfTheList<-names(bugs)
namesOfTheBugs<-sapply(namesOfTheList,function(x){strsplit(x,"_")[[1]][1]})
uniqueBugs<-unique(namesOfTheBugs)

mat<-matrix(NA,nrow = sum(nrow(myT.BS),nrow(myT.Assal),nrow(myT.Ilhan),nrow(myT.Afshar)),
            ncol=length(uniqueBugs),
            dimnames =list(c(rownames(myT.BS),rownames(myT.Assal),rownames(myT.Ilhan),rownames(myT.Afshar)),
                           uniqueBugs ))

for (i in 1:ncol(mat)){
  
  bugName<-colnames(mat)[i]
  
  if (paste0(bugName,"_BS") %in% namesOfTheList){
    mat[1:nrow(myT.BS),i]<-bugs[[paste0(bugName,"_BS")]]
  }else{
    mat[1:nrow(myT.BS),i]<-0
  } 
    
  if (paste0(bugName,"_Assal") %in% namesOfTheList){
    mat[(nrow(myT.BS)+1):(nrow(myT.BS)+nrow(myT.Assal)),i]<-bugs[[paste0(bugName,"_Assal")]]
  }else{
    mat[(nrow(myT.BS)+1):(nrow(myT.BS)+nrow(myT.Assal)),i]<-0
  } 
  
  if (paste0(bugName,"_Ilhan") %in% namesOfTheList){
    mat[(nrow(myT.BS)+nrow(myT.Assal)+1):(nrow(myT.BS)+nrow(myT.Assal)+nrow(myT.Ilhan)),i]<-bugs[[paste0(bugName,"_Ilhan")]]
  }else{
    mat[(nrow(myT.BS)+nrow(myT.Assal)+1):(nrow(myT.BS)+nrow(myT.Assal)+nrow(myT.Ilhan)),i]<-0
  } 
  
  if (paste0(bugName,"_Afshar") %in% namesOfTheList){
    mat[(nrow(myT.BS)+nrow(myT.Assal)+nrow(myT.Ilhan)+1):nrow(mat),i]<-bugs[[paste0(bugName,"_Afshar")]]
  }else{
    mat[(nrow(myT.BS)+nrow(myT.Assal)+nrow(myT.Ilhan)+1):nrow(mat),i]<-0
  } 

}

df<-as.data.frame(mat)
meta_all<-data.frame(time=c(meta.BS$prepost,meta.Assal$prepost,meta.Ilhan$prepost,meta.Afshar$prepost),
                     timepoint=c(as.character(meta.BS$Timepoint),as.character(meta.Assal$time),as.character(meta.Ilhan$Group),as.character(meta.Afshar$time)),
                     ID=c(meta.BS$Participant_ID,meta.Assal$ID,meta.Ilhan$ID,meta.Afshar$ID),
                     Study=c(rep("BS",nrow(meta.BS)),rep("Assal",nrow(meta.Assal)),
                             rep("Ilhan",nrow(meta.Ilhan)),rep("Afshar",nrow(meta.Afshar))),
                     Sample_ID=c(rownames(myT.BS),rownames(myT.Assal),rownames(myT.Ilhan),rownames(myT.Afshar)))
rownames(meta_all)<-meta_all$Sample_ID


write.table(df,paste0(output,taxa,"_countTable_merged_",normalized,".txt"),sep="\t",quote = FALSE)
write.table(meta_all,paste0(output,"metaData_merged.txt"),sep="\t",quote = FALSE)









