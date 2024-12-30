library(Seurat)
library(tidyverse)
pbmc=readRDS("pfc_rh.rds")
DefaultAssay(pbmc) <- "RNA"
clu=c("L2-3 IT", "L3-5 IT-1", "L3-5 IT-2", "L3-5 IT-3",  "L5-6 NP", "L6 CT", "L6 IT-1", "L6 IT-2", "L6B",
        "LAMP5 LHX6", "LAMP5 RELN", "VIP", "ADARB2 KCNG1", "SST", "PVALB", "PVALB ChC",
        "Astro","OPC","Oligo","Micro","Immune","Endo","PC","SMC","VLMC")
for(cl in clu){
vcell=pbmc@assays$RNA@data[,pbmc$subclass==cl]
exp=vcell[,1:ncol(vcell)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data2=t(data)
sam=pbmc@meta.data[rownames(data2),]$samplename
sp=pbmc@meta.data[rownames(data2),]$species
data3=data.frame(sam,data2)
gene_list0=colnames(data3)
gene_list=gene_list0[-which(gene_list0=="sam")]
data4=data.frame(sp,data3)
aa <- c(0,0,0)
frame <- data.frame(aa)
rownames(frame)=c("sp","sam","Residuals")
for(ge in gene_list){
	a1=aov(get(ge) ~ sp+sam ,data4)
	a2=summary(a1)[[1]]["Mean Sq"]
	names(a2)[1] = ge
	frame[ge]=a2
}
write_csv(frame,file = paste0("sum_frame_rh/", cl, ".csv"))
}
		
