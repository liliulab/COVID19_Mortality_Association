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
library("metaseqR")

setwd("/Users/matthewlee/Projects/SARS-CoV-2/Ethnicity/Data")

dat <- read.delim("Ethnic_MicroArray.txt", header=TRUE, stringsAsFactors=FALSE)
meta <- read.csv("Ethnic_Meta.csv", header = TRUE, stringsAsFactors=FALSE)
meta <- meta[meta$Vaccine == "vaccine: saline", ]

hispanic <- meta[meta$Ethnicity == "ethnicity: Hispanic", ]
white <- meta[meta$Ethnicity == "ethnicity: Not Hispanic or Latino" | meta$Ethnicity == "ethnicity: Non Hispanic or Latino", ]

row.names(dat) <- dat$ID_REF
dat <- dat[ , -1]
remove <- which(is.na(dat))
dat <- dat[-remove, ]

dat <- subset(dat, select = c(white$GEO, hispanic$GEO))


finalFrame = data.frame("ProbeID" = row.names(dat), "TValue" = NA, "PValue" = NA, "MeanX" = NA, "MeanY" = NA)

for(i in 1:nrow(dat)){
  tryCatch({
    temp = t.test(dat[i, 1:27], dat[i, 28:36])
    finalFrame[i, 2] = temp$statistic
    finalFrame[i, 3] = temp$p.value
    finalFrame[i, 4] = temp$estimate[1]
    finalFrame[i, 5] = temp$estimate[2]
  }, error=function(e){cat("ERROR :", "\n")})
}

finalFrame = finalFrame[order(finalFrame$PValue), ]
geneID <- data.frame(Gene=unlist(mget(x = finalFrame$ProbeID, envir = illuminaHumanv3SYMBOL)))
finalFrame$GeneID = geneID$Gene
finalFrame <- finalFrame[complete.cases(finalFrame$GeneID), ]
finalFrame <- finalFrame[!duplicated(finalFrame$GeneID), ]

finalFrame$modlogFC = ifelse(finalFrame$PValue<.05, (finalFrame$MeanX/finalFrame$MeanY), 0)
finalFrame$PAdjust <- p.adjust(finalFrame$PValue, method = 'fdr')

finalFrame = subset(finalFrame, select = c(GeneID, TValue, PValue, PAdjust, MeanX, MeanY, modlogFC))
row.names(finalFrame) <- finalFrame$GeneID
dim(finalFrame[which(finalFrame$PValue < .05), ])

#write.csv(finalFrame, file="Ethnicity_Result.csv")
