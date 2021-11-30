# Obese

setwd("/Users/matthewlee/Projects/SARS-CoV-2")

library("hgu133plus2.db")

dat = read.delim("Asthma/Data/Asthma_MicroArray.txt", header=TRUE, sep="\t")
meta = read.csv("Asthma/Data/Asthma_Meta.csv", header=TRUE)

obese_control <- meta[meta$Disease_Status == "Obese_Control", ]
obese_asthma <- meta[meta$Disease_Status == "Obese_Asthma", ]
no_control <- meta[meta$Disease_Status == "NO_Control", ]
no_asthma <- meta[meta$Disease_Status == "NO_Asthma", ]

control <- rbind(no_asthma, no_control)
obese <- rbind(obese_asthma, obese_control)

row.names(dat) <- dat$ID_REF
dat <- dat[ , -1]
remove <- which(is.na(dat))
dat <- dat[-remove, ]

dat <- subset(dat, select = c(control$GEO_Accession, obese$GEO_Accession))

finalFrame = data.frame("ProbeID" = row.names(dat), "TValue" = NA, "PValue" = NA, "MeanX" = NA, "MeanY" = NA)

for(i in 1:nrow(dat)){
  temp = t.test(dat[i, 79:156], dat[i, 1:78])
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

write.csv(finalFrame, file="Obese/Output/Obese_Result.csv")





