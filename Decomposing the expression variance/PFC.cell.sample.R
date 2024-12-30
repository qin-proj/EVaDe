
library(Seurat)
library(tidyverse)
pbmc=readRDS("pfc.file")
class = unique(pbmc$subclass)
pbmc=pbmc[,(pbmc$species=="Rhesus" | pbmc$species=="Human")]
for(i in class) {
	cell=names(pbmc$cell_name[pbmc$subclass==i])
	if(length(cell)>1800){
	cha=length(cell)-1800
	w1=sample(cell,cha,replace=F)
	pbmc<- pbmc[,!colnames(pbmc) %in% w1]
	}
}
saveRDS(pbmc, file = "pfc_rh.rds")
