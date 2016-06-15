library(gplots)

install.packages('RColorBrewer')
library(RColorBrewer)

DC <- read.table("/Users/iwamurac2/heatmap/DC.txt", header=TRUE, row.names=1)

png('test-DC1.png', width = 2000, height = 2000, res = 300)

heatmap.2(as.matrix(DC), 
          col=brewer.pal(11, "RdBu"), 
          tracecol=NULL, 
          scale="none",
          dendrogram="both",
          Rowv = TRUE,
          Colv = TRUE, 
          key=TRUE, 
          density.info="none", 
          margin=c(5,4),
          main="heatmap2")

dev.off()
