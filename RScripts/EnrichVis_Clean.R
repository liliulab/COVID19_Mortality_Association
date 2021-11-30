library("clusterProfiler")
library(enrichplot)
library(ggnewscale)
library(tidyverse)
library(plyr)
library(dplyr)
library(biomaRt)
library(simplifyEnrichment)

data(geneList, package="DOSE")

setwd("/Users/matthewlee/Projects/SARS-CoV-2")
dat = read.csv("Entrez.csv", header=TRUE)
de = dat$To

fc = dat$PCC
names(fc) = dat$To

ego <- data.frame(enrichGO(de, OrgDb = "org.Hs.eg.db", ont="ALL", readable=TRUE, pAdjustMethod = "BH",
                pvalueCutoff= 0.99, qvalueCutoff= 0.99))
ego.sub = ego[which(ego$p.adjust < .1), ]

write.csv(ego, file="GSEA.csv", row.names=FALSE)
write.csv(ego.sub, file="GSEA_sub.csv", row.names=FALSE)

pdf("CNET_Plot.pdf", height=12, width=15)
cnetplot(ego.sub, showCategory = 6, foldChange = fc)
dev.off()


# -------------------------------------------------- #
# -------------------------------------------------- #

# KEGG 

kegg = enrichKEGG(de, pAdjustMethod = "fdr", pvalueCutoff = 0.99, qvalueCutoff = 0.99)
kegg.sub = subset_enrichResult(kegg, c(1,2,7))
cnetplot(kegg.sub)



