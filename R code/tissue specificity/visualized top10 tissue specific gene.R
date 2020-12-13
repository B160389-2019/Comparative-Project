#BiocManager::install("ComplexHeatmap")
install.packages("viridis")
library(ComplexHeatmap)
library(viridis)
library(RColorBrewer)
library(circlize)
####################################

load("orthologous_cattle_human.Rdata")
setwd("/Users/yaoyuelin/Desktop/comparative_file/initial_dataset_cattle")
rm(list = ls())
load("human10830_cattle4866_17315gene.Rdata")
load("15696samples_17315genes.Rdata")
load("overlap_geneid.Rdata")

load("top10_gene.Rdata")
load("up_gene_cattle.Rdata")
load("overlap.Rdata")
load("up_gene_human.Rdata")
load("top10human.Rdata")

load("orthologous_cattle_human.Rdata")
######cattle_in cattle
length(unique(top10_gene_cattle$gene))
cattle_gene<-unique(top10_gene_cattle$gene)
TPM_gene<-Cattle_expression[,colnames(Cattle_expression)%in%unique(cattle_gene)]
dim(TPM_gene)
TPM_gene<-TPM_gene[,order(factor(colnames(TPM_gene),levels=unique(cattle_gene)))]
table(colnames(TPM_gene)==cattle_gene)
sub_info<-Meta_data_cattle[order(Meta_data_cattle$tissue_new),]
table(sub_info$Sample==rownames(TPM_gene))
TPM_gene<-TPM_gene[match(sub_info$Sample,rownames(TPM_gene)),]
TPM_gene<-t(log2(TPM_gene+0.25))
dim(TPM_gene)


library(pheatmap)
Tissue<-sub_info$tissue_new
annotation_mix = data.frame(Tissue)
col<-list(Tissue=c("#00008B","#008B8B","#8B008B","#8B0000","#828282","#00B2EE","#0000FF","#FF0000","#FF00FF","#FF8C00","#FFFF00","#006400","#00FF00",
                   "#CD950C","#FF6A6A","#CCCC33","pink","#FF1493","#8B7355","#00F5FF"))
names(col$Tissue) <- unique(sub_info$tissue_new)
rownames(annotation_mix)=colnames(TPM_gene)
length(rownames(annotation_mix))
length(colnames(TPM_gene))
breaksList = seq(0, 10, by = 1)
pheatmap(TPM_gene,scale = "row",cluster_rows  = F,
         annotation_col=annotation_mix,annotation_colors = col,
         cluster_cols  = F, show_rownames=F,show_colnames = F,color = colorRampPalette(c("white","firebrick3"))(10),
         main = '  ',legend = T,breaks = breaksList,
         file="/Users/yaoyuelin/Desktop/top10cattle_incattle.tiff")



##########cattle_in human

human_gene<-Orthologous_human_cattle$`Gene stable ID`[match(cattle_gene,Orthologous_human_cattle$`Cow gene stable ID`)]
TPM_gene<-human_expression[,colnames(human_expression)%in%unique(human_gene)]
dim(TPM_gene)
TPM_gene<-TPM_gene[,order(factor(colnames(TPM_gene),levels=unique(human_gene)))]
table(colnames(TPM_gene)==human_gene)
sub_info<-Meta_data_human[order(Meta_data_human$tissue_new),]
table(sub_info$SAMPID==rownames(TPM_gene))
TPM_gene<-TPM_gene[match(sub_info$SAMPID,rownames(TPM_gene)),]
table(sub_info$SAMPID==rownames(TPM_gene))
TPM_gene<-t(log2(TPM_gene+0.25))
dim(TPM_gene)

Tissue<-sub_info$tissue_new
annotation_mix = data.frame(Tissue)
col<-list(Tissue=c("#00008B","#008B8B","#8B008B","#8B0000","#828282","#00B2EE","#0000FF","#FF0000","#FF00FF","#FF8C00","#FFFF00","#006400","#00FF00",
                   "#CD950C","#FF6A6A","#CCCC33","pink","#FF1493","#8B7355","#00F5FF"))
names(col$Tissue) <- unique(sub_info$tissue_new)
rownames(annotation_mix)=colnames(TPM_gene)
length(rownames(annotation_mix))
length(colnames(TPM_gene))
breaksList = seq(0, 10, by = 1)
pheatmap(TPM_gene,scale = "row",cluster_rows  = F,
         annotation_col=annotation_mix,annotation_colors = col,
         cluster_cols  = F, show_rownames=F,show_colnames = F,color = colorRampPalette(c("white","firebrick3"))(10),
         main = '  ',legend = T,breaks = breaksList,
         file="/Users/yaoyuelin/Desktop/top10cattle_in_human.tiff")



######human_in cattle
length(unique(top10_gene_human$gene))
human_gene<-unique(top10_gene_human$gene)
TPM_gene<-human_expression[,colnames(human_expression)%in%unique(human_gene)]
dim(TPM_gene)
TPM_gene<-TPM_gene[,order(factor(colnames(TPM_gene),levels=unique(human_gene)))]
table(colnames(TPM_gene)==human_gene)
sub_info<-Meta_data_human[order(Meta_data_human$tissue_new),]
table(sub_info$SAMPID==rownames(TPM_gene))
TPM_gene<-TPM_gene[match(sub_info$SAMPID,rownames(TPM_gene)),]
table(sub_info$SAMPID==rownames(TPM_gene))
TPM_gene<-t(log2(TPM_gene+0.25))
dim(TPM_gene)


Tissue<-sub_info$tissue_new
annotation_mix = data.frame(Tissue)
col<-list(Tissue=c("#00008B","#008B8B","#8B008B","#8B0000","#828282","#00B2EE","#0000FF","#FF0000","#FF00FF","#FF8C00","#FFFF00","#006400","#00FF00",
                   "#CD950C","#FF6A6A","#CCCC33","pink","#FF1493","#8B7355","#00F5FF"))
names(col$Tissue) <- unique(sub_info$tissue_new)
rownames(annotation_mix)=colnames(TPM_gene)
length(rownames(annotation_mix))
length(colnames(TPM_gene))
breaksList = seq(0, 10, by = 1)
pheatmap(TPM_gene,scale = "row",cluster_rows  = F,
         annotation_col=annotation_mix,annotation_colors = col,
         cluster_cols  = F, show_rownames=F,show_colnames = F,color = colorRampPalette(c("white","firebrick3"))(10),
         main = '  ',legend = T,breaks = breaksList,
         file="/Users/yaoyuelin/Desktop/top10human_in_human.tiff")



##########human in cattle
human_gene<-unique(top10_gene_human$gene)
cattle_gene<-Orthologous_human_cattle$`Cow gene stable ID`[match(human_gene,Orthologous_human_cattle$`Gene stable ID`)]
length(cattle_gene)
TPM_gene<-Cattle_expression[,colnames(Cattle_expression)%in%unique(cattle_gene)]
dim(TPM_gene)
TPM_gene<-TPM_gene[,order(factor(colnames(TPM_gene),levels=unique(cattle_gene)))]
table(colnames(TPM_gene)==cattle_gene)
sub_info<-Meta_data_cattle[order(Meta_data_cattle$tissue_new),]
table(sub_info$Sample==rownames(TPM_gene))
TPM_gene<-TPM_gene[match(sub_info$Sample,rownames(TPM_gene)),]
table(sub_info$Sample==rownames(TPM_gene))
TPM_gene<-t(log2(TPM_gene+0.25))
dim(TPM_gene)

Tissue<-sub_info$tissue_new
annotation_mix = data.frame(Tissue)
col<-list(Tissue=c("#00008B","#008B8B","#8B008B","#8B0000","#828282","#00B2EE","#0000FF","#FF0000","#FF00FF","#FF8C00","#FFFF00","#006400","#00FF00",
                   "#CD950C","#FF6A6A","#CCCC33","pink","#FF1493","#8B7355","#00F5FF"))
names(col$Tissue) <- unique(sub_info$tissue_new)
rownames(annotation_mix)=colnames(TPM_gene)
length(rownames(annotation_mix))
length(colnames(TPM_gene))
breaksList = seq(0, 10, by = 1)
pheatmap(TPM_gene,scale = "row",cluster_rows  = F,
         annotation_col=annotation_mix,annotation_colors = col,
         cluster_cols  = F, show_rownames=F,show_colnames = F,color = colorRampPalette(c("white","firebrick3"))(10),
         main = '  ',legend = T,breaks = breaksList,
         file="/Users/yaoyuelin/Desktop/top10human_in_cattle.tiff")







