rm(list = ls()) 
library(xlsx)

GO_overlap_testis <- read.xlsx("/Users/yaoyuelin/Desktop/Tissue GO/overlap_GO/Testis.xlsx",header = T,stringsAsFactors = F,sheetIndex = 1)
head(GO_overlap_testis)

GO_overlap_testis <- GO_overlap_testis[c(1,4,6),c(3,4,7)]
GO_overlap_testis
GO_overlap_testis$logFDR <- -log(GO_overlap_testis$p.adjust,10)




GO_cattle_testis <- read.xlsx("/Users/yaoyuelin/Desktop/Tissue GO/cattle-specific/cattle_Testis.xlsx",header = T,stringsAsFactors = F,sheetIndex = 1)

head(GO_cattle_testis)

GO_cattle_testis <- GO_cattle_testis[c(1,2,3),c(3,4,7)]
GO_cattle_testis$logFDR <- -log(GO_cattle_testis$p.adjust,10)


GO_cattle_testis


GO_human_testis <- read.xlsx("/Users/yaoyuelin/Desktop/Tissue GO/human-specific/human_Testis.xlsx",header = T,stringsAsFactors = F,sheetIndex = 1)

head(GO_human_testis,20)
GO_human_testis <- GO_human_testis[c(10,2,5),c(3,4,7)]
GO_human_testis$logFDR <- -log(GO_human_testis$p.adjust,10)
GO_human_testis

GO_cattle_testis
GO_human_testis
GO_overlap_testis

GO_rep <- rbind(GO_cattle_testis,GO_human_testis,GO_overlap_testis)
GO_rep$Category <- factor(rep(c("Cattle","Human","Overlap"),each=3),levels = c("Cattle","Human","Overlap"))
GO_rep$"Description" <- factor(GO_rep$"Description",levels = rev(GO_rep$"Description"))
GO_rep$Description
GO_rep
colnames(GO_rep)<-c("Description","GeneRatio","FDR","logFDR","Category")

GO_rep


library(xlsx)
library(ggplot2)
GO_rep$FDR

GO_rep
GO_rep




GO_rep$number <- factor(rev(1:nrow(GO_rep)))
## shorten the names of GO terms


tiff(file = "/Users/yaoyuelin/Desktop/ENRICHMENT.tiff",res = 300, width = 3300, height = 1200,compression = "lzw")

p<-ggplot(data=GO_rep, aes(x=factor(Description,level=Description), y=logFDR,fill=Category),split='Category') +scale_fill_manual(values=c("blue", "red", "orange"))+geom_bar(stat="identity")+xlab("")+ylab(expression(paste(-log[10],italic(FDR))))+theme(axis.text.x = element_text(size=12,color = "black"),axis.text.y = element_text(size=20,colour = "black"),legend.text =element_text(size=15),legend.title = element_text(size=15),strip.text =element_text(size=15),axis.title.x =element_text(size=15) )+facet_wrap( ~Category,nrow = 1)

p + coord_flip()+theme(panel.grid.major=element_line(colour=NA))
dev.off()
