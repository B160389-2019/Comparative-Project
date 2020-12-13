rm(list = ls())
setwd("~/Comparative/comparative_file/initial_dataset_cattle")
load("orthologous_cattle_human.Rdata") # ortholog information file
load("human10830_cattle4866_17315gene.Rdata") # human and cattle dataset
load("15696samples_17315genes.Rdata") # human and cattle combined dataset
library(ggplot2)

#################################################
tissue<-sort(unique(Meta_data_human$tissue_new))
Results <- array(NA,dim = c(length(tissue),1))
rownames(Results)<-tissue
colnames(Results)<-"expression_gene_number"
Results

## calculated the average gene expression in each tissue for human 
#####


for( i in 1: length(tissue)){
  tissue_sample <-Meta_data_human[Meta_data_human$tissue_new==tissue[i],]
  
  human_expression_tissue<-human_expression[match(tissue_sample$SAMPID,rownames(human_expression)),]
  human_expression_tissue[1:10,1:10]
  expression_gene_mean<-apply(human_expression_tissue,2,mean)
  Results[i,]<-length(expression_gene_mean[expression_gene_mean>1])
}





############################calculated the average gene expression in each tissue for cattle

tissue_cattle<-sort(unique(Meta_data_cattle$tissue_new))
Results_cattle <- array(NA,dim = c(length(tissue),1))
rownames(Results_cattle)<-tissue_cattle
colnames(Results_cattle)<-"expression_gene_number"
Results_cattle
#####


for( i in 1: length(tissue)){
  tissue_sample <-Meta_data_cattle[Meta_data_cattle$tissue_new==tissue_cattle[i],]
  
  cattle_expression_tissue<-Cattle_expression[match(tissue_sample$Sample,rownames(Cattle_expression)),]
  cattle_expression_tissue[1:10,1:10]
  expression_gene_mean<-apply(cattle_expression_tissue,2,mean)
  Results_cattle[i,]<-length(expression_gene_mean[expression_gene_mean>1])
}


colors<-c("#00008B","#008B8B","#8B008B","#8B0000","#828282","#00B2EE","#0000FF","#FF0000","#FF00FF","#FF8C00","#FFFF00","#006400","#00FF00",
          "#CD950C","#FF6A6A","#CCCC33","pink","#FF1493","#8B7355","#00F5FF")

Results
Results<-as.data.frame(Results)
Results_cattle<-as.data.frame(Results_cattle)
Results_cattle
Results<-Results/1000
Results_cattle<-Results_cattle/1000
combine<-data.frame(human=Results$expression_gene_number,cattle=Results_cattle$expression_gene_number)


##### calculated the correlation between the number of expressed gene between human and cattle tissue


library(ggpubr)
tiff("~/Desktop/human-cattle-tissue-gene-number.tiff",res = 300, width = 650, height = 600,compression = "lzw")
par(mar=c(5.1,6,6,10))
DTI1 <-ggscatter(combine, x = "human", y = "cattle",
                 
                 add = "reg.line", conf.int = TRUE, color =colors , size = 2, cor.coef.size = 4,
                 
                 cor.coef = TRUE,cor.method = "spearman", col = "deeppink", #label = rownames(comb),
                 
                 xlab = " ", ylab = " ",
                 
                 col.lab="red", cex.lab=0.8,cex.axis = 0.8)
DTI1+font("xlab", size =10 , color = "black")+
  font("ylab", size = 10, color = "black")+
  font("xy.text", size = 10, color = "black")

dev.off()



