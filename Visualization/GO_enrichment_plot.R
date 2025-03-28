#Perform GO enrichment analysis and Visualization

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)
library(DOSE)
library(ggrepel)
g01 <- read.csv("key_exp02_genes",header=F)
gene1=unique(g01$V1)  # Key genes for enrichment analysis
g02 <- read.csv("all_exp02_genes",header=F)
gene2=unique(g02$V1)  # Background gene set

# Perform GO enrichment analysis
ego_ALL8 <- enrichGO(gene = gene1,keyType = "SYMBOL",universe=gene2,
OrgDb = org.Hs.eg.db,
ont = "ALL",
pAdjustMethod = "BH",
pvalueCutoff = 0.05,
qvalueCutoff = 0.05)
num_terms <- 18 # Number of terms to display
reversed_result <- ego_ALL8@result[rev(seq_len(nrow(ego_ALL8@result))),]
reversed_result <- reversed_result %>% mutate(Description = fct_inorder(Description))

# Create horizontal bar plot of enriched GO terms
p <- ggplot(reversed_result, aes(x = Description, y = Count, fill = p.adjust)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_gradient(low = "red", high = "blue", name = "Adjusted P-value") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 14, color = "black", hjust = 1), 
    axis.title.x = element_text(size = 14, color = "black"),
    legend.text = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 14, color = "black"),axis.title.y = element_blank()
  ) +
  labs(
    fill = "Adjusted P-value"
  )
ggsave("go_enrichment_plot.png", p, width = 12, height = 8, dpi=300)
