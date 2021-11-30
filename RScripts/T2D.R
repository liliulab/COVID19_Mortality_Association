# T2D

setwd("/Users/matthewlee/Projects/SARS-CoV-2")

library("hgu133plus2.db")

dat <- read.csv("T2D/Data/T2D_MA.csv", header=TRUE, stringsAsFactors=FALSE)
meta <- read.csv("T2D/Data/T2D_Meta.csv", header = TRUE, stringsAsFactors=FALSE)
meta$Status = c(replicate(10,"Control"), replicate(12, "T2D"))

control <- meta[meta$Status == "Control", ]
t2d <- meta[meta$Status == "T2D", ]

row.names(dat) <- dat$ID_REF
dat <- dat[ , -1]
remove <- which(is.na(dat))
dat <- dat[-remove, ]

dat <- subset(dat, select = c(t2d$X.Sample_geo_accession, control$X.Sample_geo_accession))

finalFrame = data.frame("ProbeID" = row.names(dat), "TValue" = NA, "PValue" = NA, "MeanX" = NA, "MeanY" = NA)

for(i in 1:nrow(dat)){
  temp = t.test(dat[i, 1:12], dat[i, 13:22])
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
dim(finalFrame[which(finalFrame$PValue < .05), ])

write.csv(finalFrame, file="T2D/Output/T2D_Result.csv")





