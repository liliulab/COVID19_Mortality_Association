# HIV

setwd("/Users/matthewlee/Projects/SARS-CoV-2")

library("illuminaHumanv4.db")

hiv <- read.delim("HIV/Data/HIV_MicroArray.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE, row.names=1)
hiv.meta <- read.csv("HIV/Data/HIV_Meta.csv", header=TRUE, stringsAsFactors=FALSE)

hiv.meta.pos <- hiv.meta[hiv.meta$HIV_Status == 1, ]
hiv.meta.neg <- hiv.meta[hiv.meta$HIV_Status == 0, ]

hiv.pos <- subset(hiv, select = hiv.meta.pos$Sample)
hiv.neg <- subset(hiv, select = hiv.meta.neg$Sample)

hiv <- cbind(hiv.neg, hiv.pos)
remove <- which(is.na(hiv))
hiv <- hiv[-remove, ]

finalFrame = data.frame("ProbeID" = row.names(hiv), "TValue" = NA, "PValue" = NA, "MeanX" = NA, "MeanY" = NA)

for(i in 1:nrow(hiv)){
  temp = t.test(hiv[i, 61:248], hiv[i, 1:60])
  finalFrame[i, 2] = temp$statistic
  finalFrame[i, 3] = temp$p.value
  finalFrame[i, 4] = temp$estimate[1]
  finalFrame[i, 5] = temp$estimate[2]
}

finalFrame$PAdjust = p.adjust(finalFrame$PValue, method = "fdr")

finalFrame = finalFrame[which(finalFrame$ProbeID != ""), ]
geneID <- data.frame(Gene=unlist(mget(x = finalFrame$ProbeID, envir = illuminaHumanv4SYMBOL)))
finalFrame$GeneID <- geneID$Gene
finalFrame <- finalFrame[complete.cases(finalFrame$GeneID), ]
finalFrame = finalFrame[order(finalFrame$PValue), ]
finalFrame = finalFrame[!duplicated(finalFrame$GeneID), ]
row.names(finalFrame) = finalFrame$GeneID

finalFrame$modlogFC = ifelse(finalFrame$PValue<.05, (finalFrame$MeanX/finalFrame$MeanY) , 0)
finalFrame = subset(finalFrame, select = c(GeneID, TValue, PValue, PAdjust, MeanX, MeanY, modlogFC))
dim(finalFrame[which(finalFrame$PValue < .05), ])

write.csv(finalFrame, file="HIV/Output/HIV_Result.csv")





