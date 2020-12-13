library(ggplot2)
file=dir("/Users/yaoyuelin/Desktop/results_new_deg/results_new_deg",full.names=T)
file
filenames<-substr(file,nchar("/Users/yaoyuelin/Desktop/results_new_deg/results_new_deg/")+1,nchar(file)-nchar(".sumstats.gz.results"))
filenames
names<-c("Adipose-C","Adipose-D","Adrenal-C","Adrenal-D",
         "Blood_Immune-C","Blood_Immune-D","Brain-C","Brain-D",
         "Heart-C","Heart-D","Kidney-C","Kidney-D","Large_intestine-C",
         "Large_intestine-D","Liver-C","Liver-D","Lung-C","Lung-D","Mammary-C",
         "Mammary-D","Muscle-C","Muscle-D","Ovary-C","Ovary-D","Pituitary-C","Pituitary-D",
         "Salivary_Gland-C","Salivary_Gland-D","Skin-C","Skin-D","Small_intestine-C","Small_intestine-D",
         "Spleen-C","Spleen-D","Stomach-C","Stomach-D","Testis-C","Testis-D","Uterus-C","Uterus-D")

Results <- array(NA,dim = c(length(filenames),length(names)))
rownames(Results)<-filenames
colnames(Results)<-names

Results
p_value<-array(NA,dim = c(length(filenames),length(names)))
rownames(p_value)<-filenames
colnames(p_value)<-names




for (i in 1:length(file)){
  filename<-filenames[i]
  
  results<-read.table(file[i],header = T)
  results$h_snp<-results$Prop._h2/results$Prop._SNPs
  Mydata_raw_FDR <- p.adjust(results$Enrichment_p,method = "BH")
  results$ad_p<-Mydata_raw_FDR
  Mydata_raw_FDR[Mydata_raw_FDR<=0.05] <- "*"
  Mydata_raw_FDR[Mydata_raw_FDR>0.05] <- " "
  Mydata_raw_FDR<-Mydata_raw_FDR[-1]
  results1<-results[-1,]
  index<-which(results1$h_snp<1)
  Mydata_raw_FDR[index]<- " "
  results1$star<-Mydata_raw_FDR
  Results[i,]<-scale(results1$h_snp-1)
  p_value[i,]<-results1$star
  
}





Classification<-rep(c("Conversed","Diverged"),20)

annotation_mix = data.frame(Classification)
annotation_mix
rownames(annotation_mix)=colnames(Results)


col<-list(Classification=c("#1b98e0","orange"))
names(col$Classification) <- unique(Classification)

summary(Results)
p_value

colnames(Results)
annotation_mix
library(pheatmap)
breaksList = seq(-1, 1, by = 0.1)
tiff(file = "/Users/yaoyuelin/Desktop/1GWAS_DEG.tiff",##reqiured to change
     res = 300, width = 3200, height = 3000,compression = "lzw")
pheatmap(Results,annotation_col = annotation_mix,annotation_colors = col, scale = "none",cluster_rows = T, display_numbers = p_value,fontsize_number = 13,color = colorRampPalette(c("white","firebrick3"))(20),breaks = breaksList,cellwidth = 10,cellheight = 10)
dev.off()


#######################3
library(ggplot2)
file=dir("/Users/yaoyuelin/Desktop/results_variance_new/results_variance_new",full.names=T)
file
filenames<-substr(file,nchar("/Users/yaoyuelin/Desktop/results_variance_new/results_variance_new/")+1,nchar(file)-nchar(".sumstats.gz.results"))
filenames
names<-c("Adipose-C","Adipose-D","Adrenal-C","Adrenal-D",
         "Blood_Immune-C","Blood_Immune-D","Brain-C","Brain-D",
         "Heart-C","Heart-D","Kidney-C","Kidney-D","Large_intestine-C",
         "Large_intestine-D","Liver-C","Liver-D","Lung-C","Lung-D","Mammary-C",
         "Mammary-D","Muscle-C","Muscle-D","Ovary-C","Ovary-D","Pituitary-C","Pituitary-D",
         "Salivary_Gland-C","Salivary_Gland-D","Skin-C","Skin-D","Small_intestine-C","Small_intestine-D",
         "Spleen-C","Spleen-D","Stomach-C","Stomach-D","Testis-C","Testis-D","Uterus-C","Uterus-D")

Results <- array(NA,dim = c(length(filenames),length(names)))
rownames(Results)<-filenames
colnames(Results)<-names

Results
p_value<-array(NA,dim = c(length(filenames),length(names)))
rownames(p_value)<-filenames
colnames(p_value)<-names




for (i in 1:length(file)){
  filename<-filenames[i]
  
  results<-read.table(file[i],header = T)
  results$h_snp<-results$Prop._h2/results$Prop._SNPs
  Mydata_raw_FDR <- p.adjust(results$Enrichment_p,method = "BH")
  results$ad_p<-Mydata_raw_FDR
  Mydata_raw_FDR[Mydata_raw_FDR<=0.05] <- "*"
  Mydata_raw_FDR[Mydata_raw_FDR>0.05] <- " "
  Mydata_raw_FDR<-Mydata_raw_FDR[-1]
  results1<-results[-1,]
  index<-which(results1$h_snp<1)
  Mydata_raw_FDR[index]<- " "
  results1$star<-Mydata_raw_FDR
  Results[i,]<-scale(results1$h_snp-1)
  p_value[i,]<-results1$star
  
}


Classification<-rep(c("Conversed","Diverged"),20)

annotation_mix = data.frame(Classification)
rownames(annotation_mix)=colnames(Results)


col<-list(Classification=c("#1b98e0","orange"))
names(col$Classification) <- unique(Classification)

Results
p_value

colnames(Results)
annotation_mix
library(pheatmap)
breaksList = seq(-1, 1, by =0.1)
tiff(file = "/Users/yaoyuelin/Desktop/1GWAS-variance.tiff",##reqiured to change
     res = 300, width = 3200, height = 3000,compression = "lzw")
pheatmap(Results,annotation_col = annotation_mix,annotation_colors = col, scale = "none",cluster_rows = T, display_numbers = p_value,fontsize_number = 13,color = colorRampPalette(c("white","firebrick3"))(20),breaks = breaksList,cellwidth = 10,cellheight = 10)
dev.off()


#######################
library(ggplot2)
file=dir("/Users/yaoyuelin/Desktop/results_cor_new/results_cor_new",full.names=T)
file
filenames<-substr(file,nchar("/Users/yaoyuelin/Desktop/results_cor_new/results_cor_new/")+1,nchar(file)-nchar(".sumstats.gz.results"))
filenames
names<-c("Adipose-D","Adipose-C","Adrenal-D","Adrenal-C",
         "Blood_Immune-D","Blood_Immune-C","Brain-D","Brain-C",
         "Heart-D","Heart-C","Kidney-D","Kidney-C","Large_intestine-D",
         "Large_intestine-C","Liver-D","Liver-C","Lung-D","Lung-C","Mammary-D",
         "Mammary-C","Muscle-D","Muscle-C","Ovary-D","Ovary-C","Pituitary-D","Pituitary-C",
         "Salivary_Gland-D","Salivary_Gland-C","Skin-D","Skin-C","Small_intestine-D","Small_intestine-C",
         "Spleen-D","Spleen-C","Stomach-D","Stomach-C","Testis-D","Testis-C","Uterus-D","Uterus-C")

Results <- array(NA,dim = c(length(filenames),length(names)))
rownames(Results)<-filenames
colnames(Results)<-names

Results
p_value<-array(NA,dim = c(length(filenames),length(names)))
rownames(p_value)<-filenames
colnames(p_value)<-names




for (i in 1:length(file)){
  filename<-filenames[i]
  
  results<-read.table(file[i],header = T)
  results$h_snp<-results$Prop._h2/results$Prop._SNPs
  Mydata_raw_FDR <- p.adjust(results$Enrichment_p,method = "BH")
  results$ad_p<-Mydata_raw_FDR
  Mydata_raw_FDR[Mydata_raw_FDR<=0.05] <- "*"
  Mydata_raw_FDR[Mydata_raw_FDR>0.05] <- " "
  Mydata_raw_FDR<-Mydata_raw_FDR[-1]
  results1<-results[-1,]
  index<-which(results1$h_snp<1)
  Mydata_raw_FDR[index]<- " "
  results1$star<-Mydata_raw_FDR
  Results[i,]<-scale(results1$h_snp-1)
  p_value[i,]<-results1$star
  
}


Classification<-rep(c("Diverged","Conserved"),20)

annotation_mix = data.frame(Classification)
annotation_mix
rownames(annotation_mix)=colnames(Results)


col<-list(Classification=c("orange","#1b98e0"))
names(col$Classification) <- unique(Classification)

Results
p_value

colnames(Results)
annotation_mix
library(pheatmap)
breaksList = seq(-1, 1, by = 0.1)
tiff(file = "/Users/yaoyuelin/Desktop/1gwas_cor.tiff",##reqiured to change
     res = 300, width = 3200, height = 3000,compression = "lzw")
pheatmap(Results,annotation_col = annotation_mix,annotation_colors = col, scale = "none",cluster_rows = T, display_numbers = p_value,fontsize_number = 13,color = colorRampPalette(c("white","firebrick3"))(20),breaks = breaksList,cellwidth = 10,cellheight = 10)
dev.off()


#######################3
library(ggplot2)
file=dir("/Users/yaoyuelin/Desktop/results_mad_new/results_mad_new",full.names=T)
file
filenames<-substr(file,nchar("/Users/yaoyuelin/Desktop/results_mad_new/results_mad_new/")+1,nchar(file)-nchar(".sumstats.gz.results"))
filenames
names<-c("Adipose-C","Adipose-V","Adrenal-C","Adrenal-V",
         "Blood_Immune-C","Blood_Immune-V","Brain-C","Brain-V",
         "Heart-C","Heart-V","Kidney-C","Kidney-V","Large_intestine-C",
         "Large_intestine-V","Liver-C","Liver-V","Lung-C","Lung-V","Mammary-C",
         "Mammary-V","Muscle-C","Muscle-V","Ovary-C","Ovary-V","Pituitary-C","Pituitary-V",
         "Salivary_Gland-C","Salivary_Gland-V","Skin-C","Skin-V","Small_intestine-C","Small_intestine-V",
         "Spleen-C","Spleen-V","Stomach-C","Stomach-V","Testis-C","Testis-V","Uterus-C","Uterus-V")

Results <- array(NA,dim = c(length(filenames),length(names)))
rownames(Results)<-filenames
colnames(Results)<-names

Results
p_value<-array(NA,dim = c(length(filenames),length(names)))
rownames(p_value)<-filenames
colnames(p_value)<-names




for (i in 1:length(file)){
  filename<-filenames[i]
  
  results<-read.table(file[i],header = T)
  results$h_snp<-results$Prop._h2/results$Prop._SNPs
  Mydata_raw_FDR <- p.adjust(results$Enrichment_p,method = "BH")
  results$ad_p<-Mydata_raw_FDR
  Mydata_raw_FDR[Mydata_raw_FDR<=0.05] <- "*"
  Mydata_raw_FDR[Mydata_raw_FDR>0.05] <- " "
  Mydata_raw_FDR<-Mydata_raw_FDR[-1]
  results1<-results[-1,]
  index<-which(results1$h_snp<1)
  Mydata_raw_FDR[index]<- " "
  results1$star<-Mydata_raw_FDR
  Results[i,]<-scale(results1$h_snp-1)
  p_value[i,]<-results1$star
  
}


Classification<-rep(c("Constant","Variable"),20)

annotation_mix = data.frame(Classification)
annotation_mix
rownames(annotation_mix)=colnames(Results)


col<-list(Classification=c("#1b98e0","orange"))
names(col$Classification) <- unique(Classification)

summary(Results)
p_value

colnames(Results)
annotation_mix
library(pheatmap)
breaksList = seq(-1, 1, by = 0.1)
tiff(file = "/Users/yaoyuelin/Desktop/1gwas_mad.tiff",##reqiured to change
     res = 300, width = 3200, height = 3000,compression = "lzw")
pheatmap(Results,annotation_col = annotation_mix,annotation_colors = col, scale = "none",cluster_rows = T, display_numbers = p_value,fontsize_number = 13,color = colorRampPalette(c("white","firebrick3"))(20),breaks = breaksList,cellwidth = 10,cellheight = 10)
dev.off()

##################################################
library(ggplot2)
file=dir("/Users/yaoyuelin/Desktop/results_tissue_specific_new/results_tissue_specific_new",full.names=T)
file
filenames<-substr(file,nchar("/Users/yaoyuelin/Desktop/results_tissue_specific_new/results_tissue_specific_new/")+1,nchar(file)-nchar(".sumstats.gz.results"))
filenames
names<-c("Adipose-N","Adipose-T","Adrenal-N","Adrenal-T",
         "Blood_Immune-N","Blood_Immune-T","Brain-N","Brain-T",
         "Heart-N","Heart-T","Kidney-N","Kidney-T","Large_intestine-N",
         "Large_intestine-T","Liver-N","Liver-T","Lung-N","Lung-T","Mammary-N",
         "Mammary-T","Muscle-N","Muscle-T","Ovary-N","Ovary-T","Pituitary-N","Pituitary-T",
         "Salivary_Gland-N","Salivary_Gland-T","Skin-N","Skin-T","Small_intestine-N","Small_intestine-T",
         "Spleen-N","Spleen-T","Stomach-N","Stomach-T","Testis-N","Testis-T","Uterus-N","Uterus-T")

Results <- array(NA,dim = c(length(filenames),length(names)))
rownames(Results)<-filenames
colnames(Results)<-names

Results
p_value<-array(NA,dim = c(length(filenames),length(names)))
rownames(p_value)<-filenames
colnames(p_value)<-names




for (i in 1:length(file)){
  filename<-filenames[i]
  
  results<-read.table(file[i],header = T)
  results$h_snp<-results$Prop._h2/results$Prop._SNPs
  Mydata_raw_FDR <- p.adjust(results$Enrichment_p,method = "BH")
  results$ad_p<-Mydata_raw_FDR
  Mydata_raw_FDR[Mydata_raw_FDR<=0.05] <- "*"
  Mydata_raw_FDR[Mydata_raw_FDR>0.05] <- " "
  Mydata_raw_FDR<-Mydata_raw_FDR[-1]
  results1<-results[-1,]
  index<-which(results1$h_snp<1)
  Mydata_raw_FDR[index]<- " "
  results1$star<-Mydata_raw_FDR
  Results[i,]<-scale(results1$h_snp-1)
  p_value[i,]<-results1$star
}


Classification<-rep(c("Non tissue specific","Tissue specific"),20)

annotation_mix = data.frame(Classification)
annotation_mix
rownames(annotation_mix)=colnames(Results)


col<-list(Classification=c("#1b98e0","orange"))
names(col$Classification) <- unique(Classification)

Results
p_value

colnames(Results)
annotation_mix
library(pheatmap)
breaksList = seq(-1, 1, by = 0.1)
tiff(file = "/Users/yaoyuelin/Desktop/1GWAS_tissue_specific.tiff",##reqiured to change
     res = 300, width = 3200, height = 3000,compression = "lzw")
pheatmap(Results,annotation_col = annotation_mix,annotation_colors = col, scale = "none",cluster_rows = T, display_numbers = p_value,fontsize_number = 13,color = colorRampPalette(c("white","firebrick3"))(20),breaks = breaksList,cellwidth = 10,cellheight = 10)
dev.off()

