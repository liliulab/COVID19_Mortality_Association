# Lupus

library("affy")
library("illuminaHumanv4.db")
library("illuminaHumanv3.db")
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

setwd("/Users/matthewlee/Projects/SARS-CoV-2/Lupus/Data")

dat <- read.delim("Lupus_MicroArray.txt", header=TRUE, stringsAsFactors=FALSE)
meta <- read.csv("Lupus_Meta.csv", header = TRUE, stringsAsFactors=FALSE)
meta <- meta[meta$X.Sample_characteristics_ch1.2 == "atherosclerosis status: no", ]

control <- meta[meta$X.Sample_characteristics_ch1.1 == "lupus status: control", ]
lup <- meta[meta$X.Sample_characteristics_ch1.1 == "lupus status: case", ]

row.names(dat) <- dat$ID_REF
dat <- dat[ , -1]
dat.orig = dat
dat[is.na(dat)] = 5.22

remove <- which(is.na(dat))
length(remove)
#dat <- dat[-remove, ]

dat[remove] = 5.22

dat <- subset(dat, select = c(control$X.Sample_geo_accession, lup$X.Sample_geo_accession))

designMatrix <- data.frame(Sample = colnames(dat), 
                           Status = c(replicate(14, "Control"), replicate(20, "Lupus")),
                           Experiment = c(1:14, 1:20))
with <- factor(designMatrix$Status, levels=c("Control", "Lupus"))
design <- model.matrix(~with)

fit <- lmFit(dat, design)
fit <- eBayes(fit, trend=TRUE, robust = TRUE)
results <- decideTests(fit)
summary(results)
limma.res <- topTable(fit, coef = "withLupus", sort = "none", n=Inf)
limma.res <- limma.res[row.names(limma.res) != "", ]

limma.res$logFC = 2^limma.res$logFC
limma.res$modlogFC = ifelse(limma.res$P.Value<.05, limma.res$logFC, 0)

limma.res = limma.res[-47232, ]

geneID <- data.frame(Gene=unlist(mget(x = row.names(limma.res), envir = illuminaHumanv4SYMBOL)))
limma.res$GeneID <- geneID$Gene

limma.res = limma.res[complete.cases(limma.res$GeneID), ]
limma.res = limma.res[order(limma.res$P.Value), ]
limma.res = limma.res[!duplicated(limma.res$GeneID), ]
limma.res <- subset(limma.res, select = c(GeneID, logFC:modlogFC))

setwd("/Users/matthewlee/Projects/SARS-CoV-2/Lupus/Output")
write.csv(limma.res, file="Lupus_Limma_Result.csv")














