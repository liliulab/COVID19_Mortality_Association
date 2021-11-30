# CLD

setwd("/Users/matthewlee/Projects/SARS-CoV-2")

library("hgu133plus2.db")

dat = read.delim("CLD/Data/CLD_MicroArray.txt", header=TRUE, sep="\t")
meta = read.csv("CLD/Data/CLD_Meta.csv", header=TRUE)

control <- meta[meta$X.Sample_characteristics_ch1 == "liver sample group: : Control", ]
cld <- meta[meta$X.Sample_characteristics_ch1 == "liver sample group: : Alcoholic hepatitis", ]

row.names(dat) <- dat$ID_REF
dat <- dat[ , -1]
remove <- which(is.na(dat))
dat <- dat[-remove, ]

dat <- subset(dat, select = c(cld$X.Sample_geo_accession, control$X.Sample_geo_accession))

finalFrame = data.frame("ProbeID" = row.names(dat), "TValue" = NA, "PValue" = NA, "MeanX" = NA, "MeanY" = NA)

for(i in 1:nrow(dat)){
  tryCatch({
    temp = t.test(dat[i, 1:15], dat[i, 16:22])
    finalFrame[i, 2] = temp$statistic
    finalFrame[i, 3] = temp$p.value
    finalFrame[i, 4] = temp$estimate[1]
    finalFrame[i, 5] = temp$estimate[2]
  }, error=function(cond){
    finalFrame[i, 2] = NA
    finalFrame[i, 3] = NA
    finalFrame[i, 4] = NA
    finalFrame[i, 5] = NA
  }
  )
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
#8323

write.csv(finalFrame, file="CLD/Output/CLD_Result.csv")





