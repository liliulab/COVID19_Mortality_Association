# Rheumatoid Arthritis

setwd("/Users/matthewlee/Projects/SARS-CoV-2")

library("hgu133plus2.db")

dat = read.delim("RA/Data/RA_MicroArray.txt", header=TRUE, sep="\t")
meta = read.csv("RA/Data/RA_Meta.csv", header=TRUE)

control <- meta[meta$X.Sample_characteristics_ch1 == "disease state: healthy control", ]

untreat = c("GSM2449643","GSM2449640","GSM2449647","GSM2449652","GSM2449667","GSM2449641",
            "GSM2449649","GSM2449656","GSM2449650","GSM2449655","GSM2449691","GSM2449648",
            "GSM2449658","GSM2449639","GSM2449653","GSM2449663","GSM2449642","GSM2449645",
            "GSM2449651","GSM2449646","GSM2449661","GSM2449654","GSM2449659","GSM2449660",
            "GSM2449644","GSM2449662","GSM2449664","GSM2449666","GSM2449665","GSM2449638",
            "GSM2449657","GSM2449792","GSM2449692","GSM2449706","GSM2449751","GSM2449754",
            "GSM2449720","GSM2449811","GSM2449793","GSM2449726","GSM2449766","GSM2449785",
            "GSM2449733","GSM2449739","GSM2449723")
ra = meta[which(meta$X.Sample_geo_accession %in% untreat), ]

row.names(dat) <- dat$ID_REF
dat <- dat[ , -1]
remove <- which(is.na(dat))
dat <- dat[-remove, ]

dat <- subset(dat, select = c(control$X.Sample_geo_accession, ra$X.Sample_geo_accession))

finalFrame = data.frame("ProbeID" = row.names(dat), "TValue" = NA, "PValue" = NA, "MeanX" = NA, "MeanY" = NA)

for(i in 1:nrow(dat)){
  temp = t.test(dat[i, 44:88], dat[i, 1:43])
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

write.csv(finalFrame, file="RA/Output/RA_Limma_Result.csv")


