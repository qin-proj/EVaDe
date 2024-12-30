library(dplyr)
g01 <- read.csv("../nmr/gogene_all.csv",header=T)
pval1 <- g01[,'p1']
FDR1 <- p.adjust(pval1, method= "BH")
g02=cbind(g01, FDR1)
pval2 <- g01[,'p2']
FDR2 <- p.adjust(pval2, method= "BH")
g03=cbind(g02, FDR2)
pval3 <- g01[,'p3']
FDR3 <- p.adjust(pval3, method= "BH")
g04=cbind(g03, FDR3)
write.csv(g04,file="all-genes-fdr",row.names=F)
x01=g04[g04['FDR2']<0.05 & g04['Rps']<0 ,]
x02=g04[g04['FDR1']<0.05 & g04['Rpm']>0 ,]
x03=g04[g04['FDR3']<0.05 & g04['Rms']<0 ,]
filtered_rows <- anti_join(x01, rbind(x02, x03))
write.csv(filtered_rows,file="../nmr/key-genes-fdr",row.names=F)
