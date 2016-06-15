#KR72
##for making heatmaps of PBS data
library(gplots)
?heatmap.2


#########################################
#### A) LOAD REQUIRED PACKAGES
#########################################

if (!require("gplots")){
  install.packages("gplots", dependencies=TRUE)
library(gplots)
  }
if (!require("RColorBrewer")){
  install.packages("RColorBrewer", dependencies=TRUE)
  library(RColorBrewer)
}

#########################################
#### B) READ IN THE DATA INTO R IN MATRIX FORMAT
#########################################
#get annotation table
setwd("/Volumes/deobhome/TarbellGroup/GroupMembers/Kameron/RNAseq Data analysis/")
annotation<-read.table("annotation key.txt", sep="\t",as.is=T,header=TRUE)
annotation<-annotation[,5:7]

#MAKE DESCRIPTIONS FILE
setwd("/Volumes/deobhome/TarbellGroup/GroupMembers/Kameron/RNAseq Data analysis/")
descriptions3<-read.table("descriptions key3.csv",sep=",",as.is=T)
descriptions3[1,]<-gsub("sample ID","sample.ID", descriptions3[1,])
descriptions3[1,]<-gsub("PCA title","PCA.title", descriptions3[1,])
descriptions3[1,]<-gsub("sort date","sort.date", descriptions3[1,])
descriptions3<-descriptions3[,1:10]
colnames(descriptions3)<-descriptions3[1,]
descriptions3<-descriptions3[-1,]
rownames(descriptions3)<-NULL
#descriptions3_nonKO_IDs<-as.numeric(descriptions3[48:141,"sample.ID"])
descriptions3$new_ID<-""
descriptions3$new_ID2<-""

#descriptions3$PCA.title<-factor(descriptions3$PCA.title)
descriptions3$PCA.title<-factor(descriptions3$PCA.title,levels=unique(descriptions3$PCA.title))
#descriptions3$sort.date<-as.factor(descriptions3$sort.date)
# descriptions3$PCA.title<-LETTERS[factor(descriptions3$PCA.title,labels=c(1:length(unique(descriptions3$PCA.title))))]

library(data.table)
for (x in 1:nrow(descriptions3)){
  #row.new_ID<-descriptions3[descriptions3$sample.ID==PCA_sample_IDs2[x],]
  row.new_ID<-descriptions3[x,]
  #n<-LETTERS[match(row.new_ID$PCA.title,unique(descriptions3$PCA.title))]
  n<-LETTERS[match(row.new_ID$PCA.title,levels(descriptions3$PCA.title))]
  #n<-LETTERS[labels(row.new_ID$PCA.title)]
  if (row.new_ID$condition=="PBS" & row.new_ID$strain %like% "B6") j=1
  if (row.new_ID$condition=="PBS" & row.new_ID$strain %like% "NOD") j=2
  if (row.new_ID$condition=="CpG" & row.new_ID$strain %like% "B6") j=3
  if (row.new_ID$condition=="CpG" & row.new_ID$strain %like% "NOD") j=4
  #for (q in 1:ncol(same_title_PCA_data)) {
  #xx_new_ID<-paste(LETTERS[n],".",j,sep="")
  xx_new_ID<-paste0(n,".",j)
  #descriptions3[descriptions3$sample.ID == PCA_sample_IDs2[x],"new_ID"]<-xx_new_ID
  descriptions3[x,"new_ID"]<-xx_new_ID
}
rm(xx_new_ID)

for (y in 1:length(levels(row.new_ID$PCA.title))){
  #slice by levels(descriptions3$PCA.titles)
  slice<-descriptions3[descriptions3$PCA.title==levels(descriptions3$PCA.title)[y],]
  slice$sort.date<-factor(slice$sort.date,levels=unique(slice$sort.date))
  #w=1
  for(w in 1:nrow(slice)){
    row.new_ID2<-slice[w,]
    m<-letters[match(row.new_ID2$sort.date,levels(slice$sort.date))]
    xx_new_ID2<-paste0(slice[w,"new_ID"],m)
    slice[w,"new_ID2"]<-xx_new_ID2
    descriptions3[,"new_ID2"][descriptions3[,"sample.ID"] == slice[w,"sample.ID"]]<-xx_new_ID2
    rm(xx_new_ID2)
  }
}

rm(list=c("row.new_ID","row.new_ID2","slice","j","m","n","w","x","y"))


#IMPORT NORMALIZED DATA 
setwd("/Volumes/deobhome/TarbellGroup/GroupMembers/Kameron/RNAseq Data analysis/vipul subset normalized")
master.file.list = list.files(pattern="*.txt", recursive = T)
master.file.list<-master.file.list[!grepl("old",master.file.list)]
file.names <- vector(mode="character",length=1)
#i=3
for (i in 1:length(master.file.list)){
  file.names[i]<-substr(master.file.list[i],1,nchar(master.file.list[i])-4)
  file.names[i]<-gsub(" ", "_", file.names[i])
  assign(file.names[i], read.table(master.file.list[i], sep="\t",as.is=T,header=TRUE))
} 
rm(list=c("i"))


library(dplyr)
library(data.table)
library(tidyr)

#CLEAN UP NORMALIZED DATA
#t=1
for (t in 1:length(file.names)){
  normalized<-get(file.names[t])
  colnames(normalized)<-gsub("X","",colnames(normalized))
  normalized<-normalized[4:nrow(normalized),]
  rownames(normalized)<-normalized[,"transcriptid"]
  normalized<-subset(normalized,select=-transcriptid)
  assign(paste0(file.names[t],"2"), normalized)
}

file.names2<-c(paste0(file.names,"2"))
rm(list=file.names)
rm(list=c("t","file.names","normalized"))



#########################################
#### C) CLEAN UP NORMALIZED DATA
#########################################

#Make list of sample descriptions of only PBS conditions that are NOT RAG-/-
sample.list <- as.data.frame(descriptions3 [ descriptions3$condition=="PBS" &  !grepl("Rag",descriptions3$strain) , ])

list.33D1<- sample.list["sample.ID"] [grepl("33D1", sample.list$subset) & sample.list$sort.date=="10.26.2015" ,]
list.MHCII<-sample.list["sample.ID"] [grepl("MHCII", sample.list$subset)& sample.list$sort.date=="01.05.2016" ,]
list.CD8a<-sample.list["sample.ID"] [grepl("CD8a", sample.list$subset)& sample.list$sort.date=="10.26.2015" ,]
list.pDC<- sample.list["sample.ID"] [grepl("pDC", sample.list$subset),]
list.Ly6c<-sample.list["sample.ID"] [grepl("Ly6C", sample.list$subset),]

?apply
map.33D1<- apply(`33D1_norm2` [,colnames(`33D1_norm2`) %in% list.33D1],1,as.numeric)
map.MHCII<- apply(MHCII_norm2[,colnames(MHCII_norm2) %in% list.MHCII],1,as.numeric)
map.CD8a<- apply(CD8a_norm2[,colnames(CD8a_norm2) %in% list.CD8a],1,as.numeric)
map.pDC<- apply(pDC_norm2[,colnames(pDC_norm2) %in% list.pDC],1,as.numeric)
map.Ly6c<- apply(Ly6c_norm2[,colnames(Ly6c_norm2) %in% list.Ly6c],1,as.numeric)


keep<-c("annotation","descriptions3","map.33D1","map.CD8a","map.MHCII","map.pDC","map.Ly6c")
rm(list=setdiff(ls(), keep))


#########################################
#### D) DEFINING COLOR PALETTE AND MAKING HEATMAP
#########################################
my_palette<- colorRampPalette(c("red","yellow","green"))(n=299)

col_breaks<-c(seq(-1,0,length=100),
              seq(0,0.8,length=100),
              seq(0.8,1,length=100))

heatmap.2(map.MHCII,
          cellnote=map.MHCII,
          main="Unstimulated Strain Differences",
          notecol="black",
          density.info="none",
          trace="none",
          margins=c(12,9),
          col=my_palette,
          breaks=col_breaks,
          dendrogram = "row",
          Colv="Rowv"
)




