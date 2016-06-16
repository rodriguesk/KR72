#KR72

setwd("/Volumes/deobhome/TarbellGroup/GroupMembers/Kameron/KR72")


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
map.33D1<- apply(`33D1_norm2` [,colnames(`33D1_norm2`) %in% list.33D1],2,as.numeric)
rownames(map.33D1)<-rownames(`33D1_norm2`)
map.MHCII<- apply(MHCII_norm2[,colnames(MHCII_norm2) %in% list.MHCII],2,as.numeric)
rownames(map.MHCII)<-rownames(MHCII_norm2)
map.CD8a<- apply(CD8a_norm2[,colnames(CD8a_norm2) %in% list.CD8a],2,as.numeric)
rownames(map.CD8a)<-rownames(CD8a_norm2)
map.pDC<- apply(pDC_norm2[,colnames(pDC_norm2) %in% list.pDC],2,as.numeric)
rownames(map.pDC)<-rownames(pDC_norm2)
map.Ly6c<- apply(Ly6c_norm2[,colnames(Ly6c_norm2) %in% list.Ly6c],2,as.numeric)
rownames(map.Ly6c)<-rownames(Ly6c_norm2)

keep<-c("annotation","descriptions3","map.33D1","map.CD8a","map.MHCII","map.pDC","map.Ly6c")
rm(list=setdiff(ls(), keep))





#########################################
#### D) SHRINK MATRICES BY KEEPING ONLY SIGNIFICANT DE GENES FOR EACH SUBSET
#########################################
#read in DE gene list for NOD PBS - B6g7 PBS comparison from Vipul's DE gene list found in KR68 (not KR65 and not KR71)
setwd("/Volumes/deobhome/TarbellGroup/GroupMembers/Kameron/KR68/DE gene lists/NOD - B6")
master.file.list = list.files(pattern="*.txt", recursive = T)
master.file.list<-master.file.list[!grepl("CpG",master.file.list)]
file.names <- vector(mode="character",length=1)
#i=3
for (i in 1:length(master.file.list)){
  file.names[i]<-substr(master.file.list[i],1,nchar(master.file.list[i])-4)
  file.names[i]<-gsub(" ", "_", file.names[i])
  assign(file.names[i], read.table(master.file.list[i], sep="\t",as.is=T,header=TRUE))
} 
rm(list=c("i"))

#ONLY KEEP DE GENES THAT ARE SIGNIFICANT (adj. P Value <0.05)
k=1
for (k in 1:length(file.names)){
  DEgenes<-get(file.names[k])
  assign(paste0("genes.",file.names[k]),DEgenes$SYMBOL[DEgenes$adj.P.Val<=0.05])
}
rm(list=c("k"))
rm(list=file.names,"master.file.list")

#SUBSET MATRICES BY THEIR RESPECTIVE FILTERED DE GENE LIST
r=1
for (r in 1:length(file.names)){
  map.final<-get(paste0("genes.",file.names[r]))
  
}
rm(list=c("r"))
map.final.33D1<- map.33D1[,][ rownames(map.33D1) %in% genes.33D1_PBS ,]
map.final.CD8a<- map.CD8a[,][ rownames(map.CD8a) %in% genes.CD8a_PBS ,]
map.final.Ly6c<- map.Ly6c[,][ rownames(map.Ly6c) %in% genes.Ly6C_PBS ,]
map.final.MHCII<- map.MHCII[,][ rownames(map.MHCII) %in% genes.MHCII_PBS ,]
map.final.pDC<- map.pDC[,][ rownames(map.pDC) %in% genes.pDC_PBS ,]



setwd("/Volumes/deobhome/TarbellGroup/GroupMembers/Kameron/KR72")
write.table(map.final.33D1, file="33D1.txt",sep="\t",col.names=NA,quote=FALSE)
?write.table
#########################################
#### E) PRODUCE LABELS FOR ROWS OF HEATMAP FROM GENE NAMES
#########################################

conversion.33D1<-data.frame(transcriptid_long=rownames(map.final.33D1))
names.33D1<-merge(annotation[,1:2],conversion.33D1,by.x="transcriptid_long",by.y="transcriptid_long",sort=FALSE)

conversion.CD8a<-data.frame(transcriptid_long=rownames(map.final.CD8a))
names.CD8a<-merge(annotation[,1:2],conversion.CD8a,by.x="transcriptid_long",by.y="transcriptid_long",sort=FALSE)

conversion.Ly6c<-data.frame(transcriptid_long=rownames(map.final.Ly6c))
names.Ly6c<-merge(annotation[,1:2],conversion.Ly6c,by.x="transcriptid_long",by.y="transcriptid_long",sort=FALSE)

conversion.pDC<-data.frame(transcriptid_long=rownames(map.final.pDC))
names.pDC<-merge(annotation[,1:2],conversion.pDC,by.x="transcriptid_long",by.y="transcriptid_long",sort=FALSE)

conversion.MHCII<-data.frame(transcriptid_long=rownames(map.final.MHCII))
names.MHCII<-merge(annotation[,1:2],conversion.MHCII,by.x="transcriptid_long",by.y="transcriptid_long",sort=FALSE)



#########################################
#### F) DEFINING COLOR PALETTE AND MAKING HEATMAP
#########################################
my_palette<- colorRampPalette(c("red","yellow","green"))(n=100)

# col_breaks<-c(seq(-1,0,length=100),
#               seq(0.008080808,0.8,length=100),
#               seq(0.802020202,1,length=100))



quartz()
heatmap.2(map.final.33D1,
          main="Unstimulated Strain Differences: 33D1",
          #cellnote=map.final.33D1,
          #notecex=0.2,
          #notecol="black",
          density.info="none",
          trace="none",
          # tracecol="black",
          #margins=c(30,70),
          # col=my_palette,
          #breaks=col_breaks,
          dendrogram = "both",
          labRow = names.33D1$gene.name,
          labCol=c(rep("B6g7 PBS",3),rep("NOD PBS",3)),#colnames(map.final.33D1),
          Colv=TRUE,
          Rowv=TRUE,
          key=TRUE,
          keysize=1,
          cexRow = 0.05,
          cexCol=0.7,
          offsetCol=0,
          offsetRow=0,
          srtCol=50
          #lhei=c(0.1,4),
          #lwid=c(1.5,2.0)
          #useRaster=TRUE
)

quartz()
heatmap.2(map.final.CD8a,
          main="Unstimulated Strain Differences: CD8a",
          #cellnote=map.final.33D1,
          #notecex=0.2,
          #notecol="black",
          density.info="none",
          trace="none",
          # tracecol="black",
          #margins=c(30,70),
          # col=my_palette,
          #breaks=col_breaks,
          dendrogram = "both",
          labRow = names.CD8a$gene.name,
          labCol=c(rep("B6g7 PBS",3),rep("NOD PBS",3)),#colnames(map.final.33D1),
          Colv=TRUE,
          Rowv=TRUE,
          key=TRUE,
          keysize=1,
          cexRow = 0.05,
          cexCol=0.7,
          offsetCol=0,
          offsetRow=0,
          srtCol=50
          #lhei=c(0.1,4),
          #lwid=c(1.5,2.0)
          #useRaster=TRUE
)

quartz()
heatmap.2(map.final.Ly6c,
          main="Unstimulated Strain Differences: Ly6c",
          #cellnote=map.final.33D1,
          #notecex=0.2,
          #notecol="black",
          density.info="none",
          trace="none",
          # tracecol="black",
          #margins=c(30,70),
          # col=my_palette,
          #breaks=col_breaks,
          dendrogram = "both",
          labRow = names.Ly6c$gene.name,
          labCol=c(rep("B6g7 PBS",3),rep("NOD PBS",3)),#colnames(map.final.33D1),
          Colv=TRUE,
          Rowv=TRUE,
          key=TRUE,
          keysize=1,
          cexRow = 0.05,
          cexCol=0.7,
          offsetCol=0,
          offsetRow=0,
          srtCol=50
          #lhei=c(0.1,4),
          #lwid=c(1.5,2.0)
          #useRaster=TRUE
)

quartz()
heatmap.2(map.final.pDC,
          main="Unstimulated Strain Differences: pDC",
          #cellnote=map.final.33D1,
          #notecex=0.2,
          #notecol="black",
          density.info="none",
          trace="none",
          # tracecol="black",
          #margins=c(30,70),
          # col=my_palette,
          #breaks=col_breaks,
          dendrogram = "both",
          labRow = names.pDC$gene.name,
          labCol=c(rep("B6g7 PBS",3),rep("NOD PBS",3)),#colnames(map.final.33D1),
          Colv=TRUE,
          Rowv=TRUE,
          key=TRUE,
          keysize=1,
          cexRow = 0.05,
          cexCol=0.7,
          offsetCol=0,
          offsetRow=0,
          srtCol=50
          #lhei=c(0.1,4),
          #lwid=c(1.5,2.0)
          #useRaster=TRUE
)

quartz()
heatmap.2(map.final.MHCII,
          main="Unstimulated Strain Differences: MHCII",
          #cellnote=map.final.33D1,
          #notecex=0.2,
          #notecol="black",
          density.info="none",
          trace="none",
          # tracecol="black",
          #margins=c(30,70),
          # col=my_palette,
          #breaks=col_breaks,
          dendrogram = "both",
          labRow = names.MHCII$gene.name,
          labCol=c(rep("B6g7 PBS",3),rep("NOD PBS",3)),#colnames(map.final.33D1),
          Colv=TRUE,
          Rowv=TRUE,
          key=TRUE,
          keysize=1,
          cexRow = 0.05,
          cexCol=0.7,
          offsetCol=0,
          offsetRow=0,
          srtCol=50
          #lhei=c(0.1,4),
          #lwid=c(1.5,2.0)
          #useRaster=TRUE
)


# heatmap.2(as.matrix(dat.bcd[,myorder]),Colv=FALSE,density.info="none",lhei=c(0.1,4),dendrogram="row",scale="row",RowSideColors=clustcol.height[clusters],col=hmcols,trace="none", margin=c(30,70), hclust=hclustfunc,distfun=distfunc,lwid=c(1.5,2.0),keysize=0.3)

dev.off()
?heatmap.2
