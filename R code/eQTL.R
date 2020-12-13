rm(list = ls())
setwd("/Users/yaoyuelin/Desktop/comparative_file/initial_dataset_cattle")
load("orthologous_cattle_human.Rdata")
human_file<-dir("/Users/yaoyuelin/Desktop/GTEx_Analysis_v8_eQTL_independent",full.names = T)
cattle_file<-dir("/Users/yaoyuelin/Desktop/fine_map_eQTL",full.names = T)
human_file
cattle_file
humantissue<-substr(human_file,nchar("/Users/yaoyuelin/Desktop/GTEx_Analysis_v8_eQTL_independent/")+1,nchar(human_file)-nchar(".v8.independent_eqtls.txt"))
humantissue

cattletissue<-substr(cattle_file,nchar("/Users/yaoyuelin/Desktop/fine_map_eQTL/")+1,nchar(cattle_file)-nchar(".finemap.eQTL.txt"))
cattletissue

Results <- array(NA,dim = c(length(human_file),length(cattle_file)))
rownames(Results)<-humantissue
colnames(Results)<-cattletissue

Results1 <- array(NA,dim = c(length(human_file),length(cattle_file)))
rownames(Results1)<-humantissue
colnames(Results1)<-cattletissue
Results1



Results
Overlap_mat<-array(NA,dim = c(length(human_file),length(cattle_file)))
rownames(Overlap_mat)<-humantissue
colnames(Overlap_mat)<-cattletissue
human_egenes<-vector(length =length(human_file) )
names(human_egenes)<-humantissue
cattle_egenes<-vector(length =length(cattle_file) )
names(cattle_egenes)<-cattletissue
#######p-value
for(i in 1:length(human_file)){
  human1 <- read.table(human_file[i],header = T,stringsAsFactors = F)
  for(j in 1:length(cattle_file)){
    cattle1 <- read.table(cattle_file[j],header = T,stringsAsFactors = F)
    colnames(cattle1)<-c("Snp.id",	"Gene",	"findmapping_factor","Sample_number",	"shape1",	"shape2",	"dummy",	"ppval",	"bpval","dist",	"npval",	"slope",	"minor_allele_frequency",	"minor_allele_count",	"minor_allele_samples"   )
    
    tissue_cattle<-cattle1
    tissue_human<-human1
    tissue_human<-tissue_human[tissue_human$gene_id%in%Orthologous_human_cattle$`Gene stable ID version`,]
    tissue_cattle<-tissue_cattle[tissue_cattle$Gene%in%Orthologous_human_cattle$`Cow gene stable ID`,]
    
    #
    tissue_cattle<-tissue_cattle[tissue_cattle$npval<0.00001,]
    tissue_human<-tissue_human[tissue_human$pval_nominal<0.00001,]
    
    #
    human<-length(unique(tissue_human$gene_id))
    cattle<-length(unique(tissue_cattle$Gene))
    
    
    Orthologous<-Orthologous_human_cattle[Orthologous_human_cattle$`Cow gene stable ID`%in%tissue_cattle$Gene,]
    
    overlap<-intersect(tissue_human$gene_id,Orthologous$`Gene stable ID version`)
    Overlap<-length(overlap)
    
    phyper(Overlap-1 , human, 17315-human , cattle,lower.tail= FALSE)
    
    ########calculated 
    
    P_values <-     phyper(Overlap-1 , human, 17315-human , cattle,lower.tail= FALSE)
    Overlap_mat[i,j]<-Overlap
    Results[i,j] <- -log10(P_values) 
    human_egenes[i]<-human
    cattle_egenes[j]<-cattle
    Results1[i,j] <- P_values
    
    
  }
}













Mydata_raw_FDR <- p.adjust(Results1,method = "BH")
Mydata_raw_FDR[Mydata_raw_FDR<=0.05] <- "*"
#Mydata_raw_FDR[Mydata_raw_FDR>0.1&Mydata_raw_FDR<=0.2] <- "."
Mydata_raw_FDR[Mydata_raw_FDR>0.05] <- " "
Mydata_raw_m <- matrix(Mydata_raw_FDR,nrow = dim(Results1)[1],byrow = F)
Mydata_raw_m

Results



library(pheatmap)
breaksList = seq(0, 10, by = 1)
tiff(file = "/Users/yaoyuelin/Desktop/0.00001p-value.tiff",##reqiured to change
     res = 300, width = 2500, height = 2000,compression = "lzw")
pheatmap(Results,display_numbers = Mydata_raw_m,color = colorRampPalette(c("white","firebrick3"))(10),breaks = breaksList,cellwidth = 10,cellheight = 10)
dev.off()

