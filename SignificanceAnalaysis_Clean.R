options(scipen=999)
setwd("/Users/matthewlee/Projects/SARS-CoV-2")

hR = read.csv("HazardRatios.csv", header=TRUE, stringsAsFactors=FALSE)
hR = hR[-c(9, 13), ]


alz = read.csv("Alzheimer's/Output/AD1_Limma_Result.csv", header=TRUE, stringsAsFactors=FALSE)
asthma = read.csv("Asthma/Output/Asthma_Result.csv", header=TRUE,stringsAsFactors=FALSE)
cad = read.csv("CAD/Output/CAD_Result.csv", header=TRUE, stringsAsFactors=FALSE)
ckd = read.csv("CKD/Output/CKD_Result.csv", header=TRUE, stringsAsFactors=FALSE)
cld = read.csv("CLD/Output/CLD_Result.csv", header=TRUE, stringsAsFactors=FALSE)
copd = read.csv("COPD/Output/COPD_Result.csv", header=TRUE, stringsAsFactors=FALSE)
ethnicity = read.csv("Ethnicity/Output/Ethnicity_Limma_Result.csv", header=TRUE, stringsAsFactors=FALSE)
colnames(ethnicity) = c("X", "GeneID", "logFC", "AveExpr", "t", "PValue", "AdjPValue", "B", "modlogFC")
hiv = read.csv("HIV/Output/HIV_Result.csv", header=TRUE, stringsAsFactors=FALSE)
lupus = read.csv("Lupus/Output/Lupus_Limma_Result.csv", header=TRUE, stringsAsFactors=FALSE)
colnames(lupus) = c("X", "GeneID", "logFC", "AveExpr", "t", "PValue", "AdjPValue", "B", "modlogFC")
obese = read.csv("Obese/Output/Obese_Result.csv", header=TRUE, stringsAsFactors=FALSE)
ond = read.csv("OND/Output/OND_Result.csv", header=TRUE, stringsAsFactors=FALSE)
ra = read.csv("RA/Output/RA_Limma_Result.csv", header=TRUE, stringsAsFactors=FALSE)
t1d = read.csv("T1D/Output/T1D_Result.csv", header=TRUE, stringsAsFactors=FALSE)

cc = c('GeneID', 'PValue', 'modlogFC')

com.list = list(alz[, cc], asthma[, cc], cad[, cc], ckd[,cc], cld[,cc], copd[,cc], ethnicity[,cc], 
                hiv[,cc], lupus[,cc], obese[,cc], ond[,cc], ra[,cc], t1d[,cc])

com.list = lapply(com.list, function(x) {
  x$PValue = ifelse(x$PValue > 0.05, 1, x$PValue)
  x$modlogFC = ifelse(x$modlogFC == 0, 1, x$modlogFC)
  return(unique(x))
})

com.merged = Reduce(function(x, y) merge(x, y, by='GeneID', all=T), com.list)
com.merged = com.merged[, c(1, grep('PValue', colnames(com.merged)), grep('modlogFC', colnames(com.merged)))]
colnames(com.merged) = c("GeneID", "AlzheimersPValue", "AsthmaPValue", "CADPValue", "CKDPValue", "CLDPValue", "COPDPValue", "EthnicityPValue", "HIVPValue",
                          "LupusPValue", "ObesePValue", "ONDPValue", "RAPValue", "T1DPValue", "AlzheimersFC", "AsthmaFC", "CADFC", "CKDFC", "CLDFC", "COPDFC", "EthnicityFC", "HIVFC",
                         "LupusFC", "ObeseFC", "ONDFC", "RAFC", "T1DFC")

s = apply(com.merged[, grep('FC', colnames(com.merged))], 1, function(x) length(which(x != 1)))
s4 = which(s >= 4)

finalCorrelation = data.frame(com.merged$GeneID, "TValue" = NA, "PValue" = NA, "Coefficient" = NA, sig.cnt=s)
colnames(finalCorrelation)[1] = "GeneID"


for(i in s3){
  temp = cor.test(hR$HazardRatio, as.numeric(com.merged[i, grep('FC', colnames(com.merged))]), method="pearson", use='pairwise.complete.obs')
  finalCorrelation[i,2] = temp$statistic
  finalCorrelation[i,3] = temp$p.value
  finalCorrelation[i,4] = temp$estimate
}

finalCorrelation = cbind(finalCorrelation, com.merged)
finalCorrelation = finalCorrelation[order(finalCorrelation$PValue), ]
finalCorrelation = finalCorrelation[order(finalCorrelation$PValue), ]
finalCorrelation = finalCorrelation[which(!is.na(finalCorrelation$PValue)), ]
finalCorrelation$AdjP = p.adjust(finalCorrelation$PValue, method = 'fdr')

sigCor = finalCorrelation[which(finalCorrelation$PValue < .01), ]

fileName = paste("Output/CorrelationResult.", Sys.Date(), ".csv", sep = "")
fileNameSig = paste("Output/CorrelationResultSig", Sys.Date(), ".csv", sep = "")
write.csv(finalCorrelation, file=fileName, row.names=FALSE)
write.csv(sigCor, file=fileNameSig, row.names=FALSE)


finalCorrelation.2 = data.frame(com.merged$GeneID, "TValue" = NA, "PValue" = NA, "Coefficient" = NA, sig.cnt=s)
colnames(finalCorrelation.2)[1] = "GeneID"

for(i in s3){
  a = data.frame(hr=hR$HazardRatio, fc=as.numeric(com.merged[i, grep('FC', colnames(com.merged))]))
  
  a$fc = log(a$fc)
  temp = cor.test(a$hr, a$fc, method="pearson", use='pairwise.complete.obs')
  finalCorrelation.2[i,2] = temp$statistic
  finalCorrelation.2[i,3] = temp$p.value
  finalCorrelation.2[i,4] = temp$estimate
}

finalCorrelation.2 = cbind(finalCorrelation.2, com.merged)
finalCorrelation.2 = finalCorrelation.2[order(finalCorrelation.2$PValue), ]
finalCorrelation.2 = finalCorrelation.2[order(finalCorrelation.2$PValue), ]
finalCorrelation.2 = finalCorrelation.2[which(!is.na(finalCorrelation.2$PValue)), ]
finalCorrelation.2$AdjP = p.adjust(finalCorrelation.2$PValue, method = 'fdr')
sigCor.2 = finalCorrelation.2[which(finalCorrelation.2$PValue < .01), ]

write.csv(sigCor.2, file="SigCor.2.csv", row.names=FALSE)

fileName = paste("Output/CorrelationResult.", Sys.Date(), ".csv", sep = "")
fileNameSig = paste("Output/CorrelationResultSig", Sys.Date(), ".csv", sep = "")
write.csv(finalCorrelation.2, file=fileName, row.names=FALSE)
write.csv(sigCor, file=fileNameSig, row.names=FALSE)



