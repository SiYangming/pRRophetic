### R code from vignette source 'vignetteOutline.Snw'

###################################################
### code chunk number 1: vignetteOutline.Snw:11-14
###################################################
library(pRRophetic)
library(ggplot2)
set.seed(12345)


###################################################
### code chunk number 2: vignetteOutline.Snw:23-24
###################################################
data("bortezomibData") #exprDataBortezomib, bortIndex, studyResponse and studyIndex


###################################################
### code chunk number 3: vignetteOutline.Snw:28-29
###################################################
pRRopheticQQplot("Bortezomib")


###################################################
### code chunk number 4: vignetteOutline.Snw:34-36
###################################################
cvOut <- pRRopheticCV("Bortezomib", cvFold=5, testExprData=exprDataBortezomib)
summary(cvOut)


###################################################
### code chunk number 5: vignetteOutline.Snw:40-41
###################################################
plot(cvOut)


###################################################
### code chunk number 6: vignetteOutline.Snw:45-51
###################################################
predictedPtype <- pRRopheticPredict(exprDataBortezomib, "Bortezomib", 
selection=1)
predictedPtype_blood <- pRRopheticPredict(exprDataBortezomib, "Bortezomib", 
"blood", selection=1)
predictedPtype_solid <- pRRopheticPredict(exprDataBortezomib, "Bortezomib",
"allSolidTumors", selection=1)


###################################################
### code chunk number 7: vignetteOutline.Snw:55-64
###################################################
t.test(predictedPtype[((studyResponse == "PGx_Responder = NR") & bortIndex)], 
predictedPtype[((studyResponse == "PGx_Responder = R") & bortIndex)], 
alternative="greater")
t.test(predictedPtype_blood[((studyResponse == "PGx_Responder = NR") & bortIndex)], 
predictedPtype_blood[((studyResponse == "PGx_Responder = R") & bortIndex)],
alternative="greater")
t.test(predictedPtype_solid[((studyResponse == "PGx_Responder = NR") & bortIndex)], 
predictedPtype_solid[((studyResponse == "PGx_Responder = R") & bortIndex)],
alternative="greater")


###################################################
### code chunk number 8: vignetteOutline.Snw:68-70
###################################################
allTissuesPpv <- getPPV(predResponders=predictedPtype[((studyResponse == "PGx_Responder = R") & bortIndex)], predNonResponders=predictedPtype[((studyResponse == "PGx_Responder = NR") & bortIndex)], drug="Bortezomib")
bloodPpv <- getPPV(predResponders=predictedPtype_blood[((studyResponse == "PGx_Responder = R") & bortIndex)], predNonResponders=predictedPtype_blood[((studyResponse == "PGx_Responder = NR") & bortIndex)], drug="Bortezomib", tissue="blood")


###################################################
### code chunk number 9: vignetteOutline.Snw:75-80
###################################################
df <- stack(list(NR=predictedPtype_blood[((studyResponse == "PGx_Responder = NR")
& bortIndex)], R=predictedPtype_blood[((studyResponse == "PGx_Responder = R") & 
bortIndex)]))
ggplot(data=df, aes(y=values, x=ind)) + geom_boxplot(alpha=.3, fill=c("#CC0033", "#006633")) + 
theme_bw() + ylab("Predicted Bortezomib Sensitivity") + xlab("Clinical Response")


###################################################
### code chunk number 10: vignetteOutline.Snw:86-87
###################################################
data(ccleData) #sensDataCcle, exprMatCcle


###################################################
### code chunk number 11: vignetteOutline.Snw:91-93
###################################################
cvOut_pd <- pRRopheticCV("PD.0325901", cvFold=5, testExprData=exprMatCcle)
summary(cvOut_pd)


###################################################
### code chunk number 12: vignetteOutline.Snw:97-98
###################################################
plot(cvOut_pd)


###################################################
### code chunk number 13: vignetteOutline.Snw:102-103
###################################################
predictedPtype_ccle <- pRRopheticPredict(exprMatCcle, "PD.0325901", selection=1)


###################################################
### code chunk number 14: vignetteOutline.Snw:108-115
###################################################
cellLinesWithCcleIc50s <- names(predictedPtype_ccle)[names(predictedPtype_ccle) %in%
sensDataCcle$CCLE.Cell.Line.Name]
predCcleOrd <- predictedPtype_ccle[names(predictedPtype_ccle)]
ccleActArea_pd <- -sensDataCcle$"ActArea"[sensDataCcle$Compound == "PD-0325901"]
names(ccleActArea_pd) <- sensDataCcle$"CCLE.Cell.Line.Name"[sensDataCcle$Compound ==
"PD-0325901"]
ccleActAreaord <- ccleActArea_pd[cellLinesWithCcleIc50s]


###################################################
### code chunk number 15: vignetteOutline.Snw:119-121
###################################################
cor.test(predictedPtype_ccle[cellLinesWithCcleIc50s], ccleActAreaord, 
method="spearman")


###################################################
### code chunk number 16: vignetteOutline.Snw:125-130
###################################################
df2 <- data.frame(predCcle=predictedPtype_ccle[cellLinesWithCcleIc50s], 
actAreaCcle=ccleActAreaord)
ggplot(data=df2, aes(y=predCcle, x=actAreaCcle)) + geom_point(alpha=0.5) + 
geom_smooth(method=lm) + theme_bw() + xlab("Measured Activity Area") +
ylab("Predicted Drug Sensitivity")


###################################################
### code chunk number 17: vignetteOutline.Snw:134-136
###################################################
predictedPtype_ccle_erlotinib <- pRRopheticLogisticPredict(exprMatCcle, "Erlotinib",
selection=1)


###################################################
### code chunk number 18: vignetteOutline.Snw:140-148
###################################################
cellLinesWithCcleIc50s <- 
names(predictedPtype_ccle_erlotinib)[names(predictedPtype_ccle_erlotinib) %in%
sensDataCcle$CCLE.Cell.Line.Name]
predCcleOrd <- predictedPtype_ccle_erlotinib[names(predictedPtype_ccle_erlotinib)]
ccleActArea_pd <- sensDataCcle$"ActArea"[sensDataCcle$Compound == "Erlotinib"]
names(ccleActArea_pd) <- sensDataCcle$"CCLE.Cell.Line.Name"[sensDataCcle$Compound ==
"Erlotinib"]
ccleActAreaord <- ccleActArea_pd[cellLinesWithCcleIc50s]


###################################################
### code chunk number 19: vignetteOutline.Snw:152-156
###################################################
resistant <- names(sort(ccleActAreaord))[1:55] #55 highly resistant cell lines.
sensitive <- names(sort(ccleActAreaord, decreasing=TRUE))[1:15] #15 sensitive
t.test(predictedPtype_ccle_erlotinib[resistant], 
predictedPtype_ccle_erlotinib[sensitive])


###################################################
### code chunk number 20: vignetteOutline.Snw:160-163
###################################################
boxplot(list(Resistant=predictedPtype_ccle_erlotinib[resistant], 
Sensitive=predictedPtype_ccle_erlotinib[sensitive]), pch=20, 
vertical=TRUE, method="jitter", ylab="Log-odds of sensitivity")


###################################################
### code chunk number 21: vignetteOutline.Snw:172-178
###################################################
trainExpr <- exprDataBortezomib[, (detailedResponse %in% c(1,2,3,4,5)) & 
studyIndex %in% c("studyCode = 25", "studyCode = 40")]
trainPtype <- detailedResponse[(detailedResponse %in% c(1,2,3,4,5)) & 
studyIndex %in% c("studyCode = 25", "studyCode = 40")]
testExpr <- exprDataBortezomib[, (detailedResponse %in% c(1,2,3,4,5)) & 
studyIndex %in% c("studyCode = 39")]


###################################################
### code chunk number 22: vignetteOutline.Snw:182-183
###################################################
ptypeOut <- calcPhenotype(trainExpr, trainPtype, testExpr, selection=1)


###################################################
### code chunk number 23: vignetteOutline.Snw:187-190
###################################################
testPtype <- detailedResponse[(detailedResponse %in% c(1,2,3,4,5)) & 
studyIndex %in% c("studyCode = 39")]
cor.test(testPtype, ptypeOut, alternative="greater")


###################################################
### code chunk number 24: vignetteOutline.Snw:195-197
###################################################
t.test(ptypeOut[testPtype %in% c(3,4,5)], ptypeOut[testPtype %in% c(1,2)], 
alternative="greater")


###################################################
### code chunk number 25: vignetteOutline.Snw:200-201
###################################################
sessionInfo()


