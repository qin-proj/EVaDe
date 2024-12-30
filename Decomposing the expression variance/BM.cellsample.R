library(Seurat)
pbmc=readRDS("../bm_data/all.rds")
class = unique(pbmc$cell_subset)
species_unique <- unique(pbmc$species)
for(i in class) {
	if(ncol(pbmc[,pbmc$cell_subset == i])>=3000){
	num=1
	cell1=colnames(pbmc[,pbmc$species == species_unique[1] & pbmc$cell_subset == i])
	cell2=colnames(pbmc[,pbmc$species == species_unique[2] & pbmc$cell_subset == i])
        num1=length(cell1)
	num2=length(cell2)
	if(num1>1500 & num2 > 1500){
        cha1=num1-1500
        w1=sample(cell1,cha1,replace=F)
        pbmc<- pbmc[,!colnames(pbmc) %in% w1]
	cha2=num2-1500
        w2=sample(cell2,cha2,replace=F)
        pbmc<- pbmc[,!colnames(pbmc) %in% w2]
        }
	if(num1<=1500 & num2 > 1500){
	cha2=num2-3000+num1
        w2=sample(cell2,cha2,replace=F)
        pbmc<- pbmc[,!colnames(pbmc) %in% w2]
	}
	if(num1>1500 & num2 <= 1500){
        cha1=num1-3000+num2
        w1=sample(cell1,cha1,replace=F)
        pbmc<- pbmc[,!colnames(pbmc) %in% w1]
        }
	}else{
	cell3=colnames(pbmc[,pbmc$cell_subset == i])
	pbmc<- pbmc[,!colnames(pbmc) %in% cell3]
	

	}
}
saveRDS(pbmc, file = "bm_filtered.rds")
