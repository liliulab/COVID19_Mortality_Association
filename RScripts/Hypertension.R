# Hypertension

setwd("/Users/matthewlee/Projects/SARS-CoV-2")

library("affy")
library("illuminaHumanv1.db")
library("limma")
library("ggplot2")
library("ggrepel")
library("pheatmap")
library("RColorBrewer")
library("GEOquery")
library("Biobase")
library("biomaRt")
library("nortest")
library("hgu133plus2.db")
library("GSVA")
library("GSEABase")
library("ComplexHeatmap")
library("circlize")
("pd.hugene.1.0.st.v1")

dat = read.delim("Hypertension/Data/Hypertension_MicroArray.txt", header=TRUE, sep="\t")
conv = read.csv("Hypertension/Convert_HYP.csv", header=TRUE)
dat = cbind(conv, dat)
dat = dat[which(dat$TYPE == "mRNA"), ]
dat = subset(dat, select = c(GENE_SYMBOL, GSM3986167:GSM3986186))
meta = read.csv("Hypertension/Data/Hypertension_Meta.csv", header=TRUE)

control <- meta[meta$Hypertensive == "Normotensive", ]
hyp <- meta[meta$Hypertensive == "Hypertensive", ]

row.names(dat) <- dat$ID_REF

orderDat <- subset(dat, select = c(hyp$X.Sample_geo_accession, control$X.Sample_geo_accession))
dat = cbind(dat[,1], orderDat)
colnames(dat)[1] = "GeneID"

finalFrame = data.frame("GeneID" = dat$GeneID, "TValue" = NA, "PValue" = NA, "MeanX" = NA, "MeanY" = NA)

for(i in 1:nrow(dat)){
  temp = t.test(dat[i, 2:11], dat[i, 12:21])
  finalFrame[i, 2] = temp$statistic
  finalFrame[i, 3] = temp$p.value
  finalFrame[i, 4] = temp$estimate[1]
  finalFrame[i, 5] = temp$estimate[2]
}

finalFrame$PAdjust = p.adjust(finalFrame$PValue, method = "fdr")
length(which(finalFrame$PValue < 0.05))

finalFrame = finalFrame[complete.cases(finalFrame$GeneID), ]
finalFrame = finalFrame[order(finalFrame$PValue), ]
finalFrame = finalFrame[!duplicated(finalFrame$GeneID), ]
row.names(finalFrame) = finalFrame$GeneID

finalFrame$modlogFC = ifelse(finalFrame$PValue<.05, (finalFrame$MeanX/finalFrame$MeanY), 0)
finalFrame = subset(finalFrame, select = c(GeneID, TValue, PValue, PAdjust, MeanX, MeanY, modlogFC))
dim(finalFrame[which(finalFrame$PValue < .05), ])

write.csv(finalFrame, file="Hypertension/Output/Hypertension_Result.csv")





