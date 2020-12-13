rm(list = ls()) 
setwd("/Users/yaoyuelin/Desktop/comparative_file/initial_dataset_cattle")
load("human10830_cattle4866_17315gene.Rdata")
load("15696samples_17315genes.Rdata")
load("orthologous_cattle_human.Rdata")
setwd("~/comparative_file/initial_dataset_cattle")

library(limma)


Meta_data_cattle$tissue_new[Meta_data_cattle$tissue_new=="Blood/Immune"]<-"Blood_Immune"
tissue_name<-unique(sort(Meta_data_cattle$tissue_new))

Meta_data_cattle_sort<-Meta_data_cattle[order(Meta_data_cattle$tissue_new),]
Cattle_expression<-t(Cattle_expression)
Cattle_expression_sort<-Cattle_expression[,match(Meta_data_cattle_sort$Sample,colnames(Cattle_expression))]
table(colnames(Cattle_expression_sort)==Meta_data_cattle_sort$Sample)
Cattle_expression_sort<-log2(Cattle_expression_sort+0.25)
table(colnames(Cattle_expression_sort)==Meta_data_cattle_sort$Sample)


dim(Meta_data_cattle_sort)

###########
for (i in 1:20){
  
  tissue_index<-which(Meta_data_cattle_sort$tissue_new==tissue_name[i])
  tissue_index
  df<-c(rep("others",tissue_index[1]-1),rep("tissue",length(tissue_index)),rep("others",dim(Meta_data_cattle_sort)[1]-length(tissue_index)-(tissue_index[1]-1)))
  table(df)
  
  df
  design <- model.matrix(~-1+factor(df))
  dim(design)
  colnames(design) <- c("Others","Tissue")
  table(design)
  rownames(design)=colnames(Cattle_expression_sort)
  design
  design[tissue_index[1],]
  contrastmatrix <- makeContrasts(Tissue-Others, levels=design)
  fit <- lmFit(Cattle_expression_sort, design) ##issue these commands to fit the model and make the contrasts
  fit2 <- contrasts.fit(fit, contrastmatrix)
  fit2 <- eBayes(fit2)
  fit2
  myresults <-topTable(fit2,coef=1, adjust="fdr", number=nrow(Cattle_expression_sort))
  write.table(myresults,file=paste0("~/tissue_specific/human/",tissue_name[i],".txt"),quote = F)
}


###################3
rm(list = ls()) 
setwd("/Users/yaoyuelin/Desktop/comparative_file/initial_dataset_cattle")
load("human10830_cattle4866_17315gene.Rdata")
load("15696samples_17315genes.Rdata")
load("orthologous_cattle_human.Rdata")
setwd("~/comparative_file/initial_dataset_cattle")






###########human
Meta_data_human$tissue_new[Meta_data_human$tissue_new=="Blood/Immune"]<-"Blood_Immune"
tissue_name<-unique(sort(Meta_data_human$tissue_new))
tissue_name
Meta_data_human_sort<-Meta_data_human[order(Meta_data_human$tissue_new),]
Meta_data_human_sort
human_expression<-t(Human_expression)
human_expression[1:10,1:10]
human_expression_sort<-human_expression[,match(Meta_data_human_sort$SAMPID,colnames(human_expression))]
human_expression_sort[1:10,1:10]
table(colnames(human_expression_sort)==Meta_data_human_sort$SAMPID)
human_expression_sort<-log2(human_expression_sort+0.25)
table(colnames(human_expression_sort)==Meta_data_human_sort$SAMPID)


dim(Meta_data_human_sort)
i=1
###########
for (i in 1:20){
  
  tissue_index<-which(Meta_data_human_sort$tissue_new==tissue_name[i])
  tissue_index
  df<-c(rep("others",tissue_index[1]-1),rep("tissue",length(tissue_index)),rep("others",dim(Meta_data_human_sort)[1]-length(tissue_index)-(tissue_index[1]-1)))
  table(df)
  
  df
  design <- model.matrix(~-1+factor(df))
  dim(design)
  colnames(design) <- c("Others","Tissue")
  table(design)
  rownames(design)=colnames(human_expression_sort)
  design
  design[tissue_index[1],]
  contrastmatrix <- makeContrasts(Tissue-Others, levels=design)
  fit <- lmFit(human_expression_sort, design) ##issue these commands to fit the model and make the contrasts
  fit2 <- contrasts.fit(fit, contrastmatrix)
  fit2 <- eBayes(fit2)
  fit2
  myresults <-topTable(fit2,coef=1, adjust="fdr", number=nrow(human_expression_sort))
  write.table(myresults,file=paste0("~/tissue_specific/human/",tissue_name[i],".txt"),quote = F)
}
