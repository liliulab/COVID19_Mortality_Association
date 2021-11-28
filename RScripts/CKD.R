# CKD

setwd("/Users/matthewlee/Projects/SARS-CoV-2")

library("hgu133plus2.db")

dat = read.delim("CKD/Data/CKD_MicroArray.txt", header=TRUE, sep="\t")
meta = read.csv("CKD/Data/CKD_Meta.csv", header=TRUE)
meta = meta[grep('discovery', meta$X.Sample_title), ]

control <- meta[meta$X.Sample_characteristics_ch1 == "disease state: healthy control", ]
ckd <- meta[meta$X.Sample_characteristics_ch1 == "disease state: uremia", ]

row.names(dat) <- dat$ID_REF
dat <- dat[ , -1]
remove <- which(is.na(dat))
dat <- dat[-remove, ]

dat <- subset(dat, select = c(ckd$X.Sample_geo_accession, control$X.Sample_geo_accession))

finalFrame = data.frame("ProbeID" = row.names(dat), "TValue" = NA, "PValue" = NA, "MeanX" = NA, "MeanY" = NA)

for(i in 1:nrow(dat)){
  temp = t.test(dat[i, 1:63], dat[i, 64:83])
  finalFrame[i, 2] = temp$statistic
  finalFrame[i, 3] = temp$p.value
  finalFrame[i, 4] = temp$estimate[1]
  finalFrame[i, 5] = temp$estimate[2]
}

finalFrame$PAdjust = p.adjust(finalFrame$PValue, method = "fdr")

PROBES <- as.character(finalFrame$ProbeID)
OUT <- select(hgu133plus2.db, PROBES, c("SYMBOL"))
out <- OUT[-which(duplicated(OUT$PROBEID)), ]

finalFrame$GeneID <- out$SYMBOL
finalFrame = finalFrame[order(finalFrame$PValue), ]

finalFrame = finalFrame[complete.cases(finalFrame$GeneID), ]
finalFrame = finalFrame[!duplicated(finalFrame$GeneID), ]
row.names(finalFrame) = finalFrame$GeneID

finalFrame$modlogFC = ifelse(finalFrame$PValue<.05, (finalFrame$MeanX/finalFrame$MeanY), 0)
finalFrame = subset(finalFrame, select = c(GeneID, TValue, PValue, PAdjust, MeanX, MeanY, modlogFC))
head(finalFrame)
dim(finalFrame[which(finalFrame$PValue < .05), ])
#17688

write.csv(finalFrame, file="CKD/Output/CKD_Result.csv")





