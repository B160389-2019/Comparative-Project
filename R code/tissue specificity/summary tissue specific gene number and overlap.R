
rm(list = ls())
setwd("~/comparative_file/initial_dataset_cattle")

load("human10830_cattle4866_17315gene.Rdata")
#load("15696samples_17315genes.Rdata")
#load("orthologous_cattle_human.Rdata")

##########3

Results<-array(data = NA,dim = c(20,2))
tissue_name<-unique(sort(Meta_data_cattle$tissue_new))
rownames(Results)<-tissue_name
colnames(Results)<-c("human_up","cattle_up","overlap")


logFC_cutoff<-1.5
for i in 1:20{
  tissue_cattle<-read.table(file = paste0("/Users/yaoyuelin/Desktop/tissue_specific_limma/cattle/",tissue_name[i],".txt")
  tissue_human<-read.table(file = paste0("/Users/yaoyuelin/Desktop/tissue_specific_limma/human/",tissue_name[i],".txt")
  tissue_cattle = as.factor(ifelse(tissue_cattle$adj.P.Val < 0.05 & abs(tissue_cattle$logFC) > logFC_cutoff,
                                                         ifelse(tissue_cattle$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
                           )
                           
  cattle_upgene<-tissue_cattle[tissue_cattle$change =='UP',]
  write.table(cattle_upgene,file=paste0("~/cattle-tissue-specific-gene/",tissue_name[i],".txt")
  Results[i,2]<- length(rownames(cattle_upgene))    
  
  tissue_human = as.factor(ifelse(tissue_human$adj.P.Val < 0.05 & abs(tissue_human$logFC) > logFC_cutoff,
                                   ifelse(tissue_human$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
  )
  
  human_upgene<-tissue_human[tissue_human$change =='UP',]
  write.table(human_upgene,file=paste0("~/human-tissue-specific-gene/",tissue_name[i],".txt")
              Results[i,1]<- length(rownames(human_upgene)) 
              
  cattle_trans<-Orthologous_human_cattle$`Gene stable ID`[Orthologous_human_cattle$`Cow gene stable ID`%in%rownames(cattle_upgene)]
  overlap<-intersect(cattle_trans,rownames(human_upgene))
  Results[i,3]<-length(overlap)
  write.table(overlap,file=paste0("~/overlap-tissue-specific-gene/",tissue_name[i],".txt")
              
}



