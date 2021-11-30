# Alzheimer's

setwd("/Users/matthewlee/Projects/SARS-CoV-2")

library("illuminaHumanv3.db")


dat <- read.delim("Alzheimer's/Data/Batch_1.txt", sep = "\t", header=TRUE, stringsAsFactors=FALSE)
meta <- read.csv("Alzheimer's/Data/AD_Meta.csv", header=TRUE, stringsAsFactors=FALSE)

names <- substring(colnames(dat), 2)
names[1] <- "ID_REF"
colnames(dat) <- names
cross <- intersect(meta$Sample_Title, colnames(dat))
meta <- meta[meta$Sample_Title == cross, ]

ctl <- meta[meta$Sample_Characteristics == "Control", ]
ad <- meta[meta$Sample_Characteristics == "AD", ]
meta <- rbind(ctl, ad)
dat <- subset(dat, select = c("ID_REF", meta$Sample_Title))
row.names(dat) <- dat$ID_REF
dat <- dat[ , -1]
remove <- which(is.na(dat))
dat <- dat[-remove, ]

finalFrame = data.frame("ProbeID" = row.names(dat), "TValue" = NA, "PValue" = NA, "MeanX" = NA, "MeanY" = NA)

for(i in 1:nrow(dat)){
  temp = t.test(dat[i, 105:249], dat[i, 1:104])
  finalFrame[i, 2] = temp$statistic
  finalFrame[i, 3] = temp$p.value
  finalFrame[i, 4] = temp$estimate[1]
  finalFrame[i, 5] = temp$estimate[2]
}

finalFrame$PAdjust = p.adjust(finalFrame$PValue, method = "fdr")

finalFrame = finalFrame[which(finalFrame$ProbeID != ""), ]
geneID <- data.frame(Gene=unlist(mget(x = finalFrame$ProbeID, envir = illuminaHumanv3SYMBOL)))
finalFrame$GeneID <- geneID$Gene
finalFrame <- finalFrame[complete.cases(finalFrame$GeneID), ]
finalFrame = finalFrame[order(finalFrame$PValue), ]
finalFrame = finalFrame[!duplicated(finalFrame$GeneID), ]
row.names(finalFrame) = finalFrame$GeneID


finalFrame$modlogFC = ifelse(finalFrame$PValue<.05, (finalFrame$MeanX/finalFrame$MeanY), 0)
finalFrame = subset(finalFrame, select = c(GeneID, TValue, PValue, PAdjust, MeanX, MeanY, modlogFC))

write.csv(finalFrame, file="Alzheimer's/Output/AD_Result.csv")


