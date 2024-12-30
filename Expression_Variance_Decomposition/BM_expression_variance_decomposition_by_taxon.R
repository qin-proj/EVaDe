library(Seurat)
library(tidyverse)
pbmc=readRDS("bm_filtered.rds")
clu=unique(pbmc$cell_subset)
for(cl in clu){
vcell=pbmc@assays$RNA@data[,pbmc$cell_subset==cl & pbmc$species=="mouse"]
exp=vcell[,1:ncol(vcell)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data2=t(data)
sam=pbmc@meta.data[rownames(data2),]$sample
data3=data.frame(sam,data2)
gene_list0=colnames(data3)
gene_list=gene_list0[-which(gene_list0=="sam")]
aa <- c(0,0)
frame <- data.frame(aa)
rownames(frame)=c("sam","Residuals")
for(ge in gene_list){
        a1=aov(get(ge) ~ sam ,data3)
        a2=summary(a1)[[1]]["Mean Sq"]
        names(a2)[1] = ge
        frame[ge]=a2
}
write_csv(frame,file = paste0("mouse/", cl, ".csv"))
vcell=pbmc@assays$RNA@data[,pbmc$cell_subset==cl & pbmc$species=="nmr"]
exp=vcell[,1:ncol(vcell)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data2=t(data)
sam=pbmc@meta.data[rownames(data2),]$sample
data3=data.frame(sam,data2)
gene_list0=colnames(data3)
gene_list=gene_list0[-which(gene_list0=="sam")]
aa <- c(0,0)
frame <- data.frame(aa)
rownames(frame)=c("sam","Residuals")
for(ge in gene_list){
        a1=aov(get(ge) ~ sam ,data3)
        a2=summary(a1)[[1]]["Mean Sq"]
        names(a2)[1] = ge
        frame[ge]=a2
}
write_csv(frame,file = paste0("nmr/", cl, ".csv"))
}
