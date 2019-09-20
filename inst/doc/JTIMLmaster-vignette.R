## ---- include = TRUE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  fig.width=7, 
  fig.height=5, 
  fig.path='Figs/', 
  fig.align="left",
  warning=FALSE, 
  message=FALSE,
  comment = "#>"
)

## ----setup---------------------------------------------------------------
library(JTIMLmaster)
library(Biobase)
library(DESeq2)
library(edgeR)
library(ggplot2)
library(SmartSVA)
library(pheatmap)
library(data.table)
library(plyr)
library(dplyr)
library(DT)
library(car)
library(randomForest)
library(glmnet)
library(caret)
library(e1071)
library(VSURF)
library(neuralnet)
library(gbm)
library(verification)
library(pROC)

## ----Filtering Marine data set-------------------------------------------
## The original WTC raw data includes 25830 genes with 324 samples
## Filtering out genes having 0 sd and 0 interquantile 
sddat <- apply(WTC_raw, 1, sd)
sdidx <- which(sddat<=3 ) ## 7728 genes have <=3 SD. 
WTC_raw <- WTC_raw[-sdidx,]  ## 18102 genes left

## Filtering out genes having too extreme values and too small or too large variation
## Too small variance can inflate false positive and too large variance can inflate false negatives 
iqrdat <- apply( WTC_raw , 1 , function( x ) diff(quantile( x , c( 0.25 , 0.75 ) ) ) )
iqridx <- which(iqrdat==0)   ## n=2 genes have Inter Qualtile=0
WTC_raw <- WTC_raw[-iqridx,]   ## n=18100 genes left

## Too extreme cases 
iqrdat <- apply( WTC_raw , 1 , function( x ) diff(quantile( x , c( 0.25 , 0.75 ) ) ) )
iqridx <- c(which(iqrdat < quantile(iqrdat, 0.15)), which(iqrdat > quantile(iqrdat, 0.85)))  ## 6446 genes have <15% or >85% Inter Quantile. 
WTC_raw <- WTC_raw[-iqridx,]   ## n=12923 genes left

## ----setting DGE object--------------------------------------------------
feaD <- data.frame(gene=rownames(WTC_raw))
rownames(feaD) <- feaD$gene

phenoData <- new("AnnotatedDataFrame", data=WTC_pheno)
featureData <- new("AnnotatedDataFrame", data=feaD)
WTC_eset <- ExpressionSet(assayData=as.matrix(WTC_raw), phenoData=phenoData, featureData=featureData)

## ----Data QC-------------------------------------------------------------
## Data Normalization by TMM method ##
y <- DGEList(WTC_eset)
DGE <- calcNormFactors(y,method =c("TMM")) ## TMM = weighted trimmed mean of M-value
barplot(DGE$samples$lib.size, names=colnames(DGE),las=2)

# Add a title to the plot
title("Barplot of library sizes")   ## "id_187" "id_188" "id_190" "id_195" "id_204" "id_214" "id_217" having relatively low library size

pch <- c(0,1)
colors <- ifelse(pData(WTC_eset)$PTSD=="Never", "Blue", "Gray")
colors <- ifelse(pData(WTC_eset)$PTSD=="Past", "Orange", colors)

plotMDS(DGE, col=colors, main="MD plot(normalized)")

plotMD(WTC_raw, column = 1, main="MD plot(Before normalization, id_1)")
abline(h=0,col="red")
plotMD(WTC_raw, column = 7, main="MD plot(Before normalization, id_7) ")
abline(h=0,col="red")

plotMD(DGE, column = 1, main="MD plot(After normalization, id_1)")
abline(h=0,col="red")
plotMD(DGE, column = 7, main="MD plot(After normalization, id_7)")
abline(h=0,col="red")

## ----PCA and MDS---------------------------------------------------------
## PCA run
pcaOut <- prcomp(t(DGE$counts), center=TRUE, scale=TRUE)
dat <- data.frame(pcaOut$x, WTC_pheno)

pcaVariance <- round(unlist(lapply(1:length(pcaOut$sdev),function(i){pcaOut$sdev[i]^2/sum(pcaOut$sdev^2)})), digits=3)*100

ggplot(dat, aes(PC1, PC2, color=PTSD)) + geom_point(size=3) +
  xlab(paste0("PC1: ", pcaVariance[1], "% variance")) + ylab(paste0("PC2: ", pcaVariance[2], "% variance"))+ggtitle("PC1 vs PC2")

ggplot(dat, aes(PC1, PC3, color=PTSD)) + geom_point(size=3) +
  xlab(paste0("PC1: ", pcaVariance[1], "% variance")) + ylab(paste0("PC3: ", pcaVariance[3], "% variance")) +ggtitle("PC1 vs PC3")

WTC_pheno$PTSD <- factor(WTC_pheno$PTSD, levels = c("Never", "Current", "Past"))


## ----Feature Selection by DE analysis------------------------------------
design <- model.matrix(~-1 + PTSD, WTC_pheno) 
colnames(design) <- gsub(" ", "", colnames(design))
v <- voom(DGE, design, plot=FALSE)

groupCon <- makeContrasts(CurrentNever=PTSDCurrent-PTSDNever, 
                          CurrentPast=PTSDCurrent-PTSDPast, 
                          PastNever=PTSDNever-PTSDPast,
                          levels=design)
fit1 <- lmFit(v, design)

## Comparisons across three groups 
fit2 <- contrasts.fit(fit1, groupCon)
fit2 <- eBayes(fit2, trend=FALSE)

allOut_Current_to_Never <- topTable(fit2, number=nrow(v), coef="CurrentNever", sort="P")
allOut_Current_to_Past <- topTable(fit2, number=nrow(v), coef="CurrentPast", sort="P")
allOut_Past_to_Never <- topTable(fit2, number=nrow(v), coef="PastNever", sort="P")

## ----Vocano Plots--------------------------------------------------------
########################## Volcano Plot ###########################
with(allOut_Current_to_Never, plot(logFC, -log10(P.Value), pch=20, main="Volcano plot", cex=1.0, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~P~value)))
with(subset(allOut_Current_to_Never, adj.P.Val<0.05), points(logFC, -log10(P.Value), pch=20, col="red", cex=0.5))

with(allOut_Past_to_Never, plot(logFC, -log10(P.Value), pch=20, main="Volcano plot", cex=1.0, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~P~value)))
with(subset(allOut_Past_to_Never, adj.P.Val<0.05), points(logFC, -log10(P.Value), pch=20, col="red", cex=0.5))

with(allOut_Current_to_Past, plot(logFC, -log10(P.Value), pch=20, main="Volcano plot", cex=1.0, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~P~value)))
with(subset(allOut_Current_to_Past, adj.P.Val<0.05), points(logFC, -log10(P.Value), pch=20, col="red", cex=0.5))

sigGenesCoeff <- allOut_Current_to_Never[allOut_Current_to_Never$adj.P.Val < 0.05,]
varSelect <- rownames(sigGenesCoeff)

## ----heatmap-------------------------------------------------------------
sigGeneExp <- as.data.frame(v$E[(rownames(v$E) %in% varSelect),])
WTC_pheno_anno <- WTC_pheno[order(WTC_pheno$PTSD), c("PTSD"), drop=FALSE]
WTC_pheno_anno <- WTC_pheno_anno[which(WTC_pheno_anno$PTSD !="Past"),c("PTSD"), drop=FALSE]
sigGeneExp <- sigGeneExp[, match(rownames(WTC_pheno_anno), names(sigGeneExp))]

pheatmap(sigGeneExp[,names(sigGeneExp) %in% rownames(WTC_pheno_anno)], scale="row", cluster_rows=TRUE, cluster_cols=FALSE, annotation_col=WTC_pheno_anno)
pheatmap(sigGeneExp[,names(sigGeneExp) %in% rownames(WTC_pheno_anno)], scale="row", cluster_rows =TRUE, cluster_cols=TRUE, annotation_col=WTC_pheno_anno)

DT::datatable(sigGenesCoeff, caption = "Significant Genes between Never to Current PTSD")

## ----Model Training------------------------------------------------------
# https://cran.r-project.org/web/packages/VSURF/VSURF.pdf
dat <- data.frame(PTSD=WTC_pheno[, c("PTSD")], t(v$E[rownames(v$E)%in%varSelect,]))
dat <- dat[dat$PTSD !="Past", ]

dat$PTSD_binary <- as.factor(ifelse(dat$PTSD=="Current", 1, 0))
names(dat) <- gsub("-", "_", names(dat))

dat$PTSD <- NULL
set.seed(415)
K <- 10       ## Number of replicates 
n.case <- which(dat$PTSD_binary==1)
n.control <- which(dat$PTSD_binary==0)

train_rows <- lapply(1:K, function(x){c(sample(n.case, length(n.case), replace = TRUE), sample(n.control, length(n.control), replace = TRUE))})

depVar <- "PTSD_binary"

alphaLevel <- c(0.4, 0.7, 1)

ntree <- 100

### Random Forest ###
WTC_RFresult <- lapply(train_rows, function(x){RF.boot.err.Func(train=dat[x, ], test=dat[-x, ], depVar = depVar, ntree=ntree)})
names(WTC_RFresult) <- paste0("WTC_RF_iter", 1:K)

### Penalized Regression ###
WTC_PRresult <- lapply(train_rows, function(x){PR.boot.err.Func(x.train=dat[x, (names(dat)!=depVar)], y.train=dat[x, depVar], x.test=dat[-x, (names(dat)!=depVar)], y.test=dat[-x, depVar], alphaLevel=alphaLevel, family="binomial", type="class")})
names(WTC_PRresult) <- paste0("WTC_PR_iter", 1:K)

### SVM ###
WTC_SVMresult  <-  lapply(train_rows, function(x){SVM.boot.err.Func(train=dat[x, ], test=dat[-x, ], depVar= depVar, kernel="radial", cost=5)})
names(WTC_SVMresult) <- paste0("WTC_SVM_iter", 1:K)

### Neural Network ###
WTC_NNresult  <-  lapply(train_rows, function(x){NN.boot.err.Func(train=dat[x, ], test=dat[-x, ], depVar= depVar)})
names(WTC_NNresult) <- paste0("WTC_NN_iter", 1:K)

### Gredient Boosting (gbm R package) ###
WTC_GBMresult  <-  lapply(train_rows, function(x){GBM.boot.err.Func(train=dat[x, ], test=dat[-x, ], depVar= depVar, distribution="adaboost", n.trees = ntree)})
names(WTC_GBMresult) <- paste0("WTC_GBM_iter", 1:K)

## ----Building a model----------------------------------------------------
## Bootstrap sampling into relatively significant genes from WTC_raw's differential analysis
## Split them into Training and Testing
## Machine learners into Training set>> repeat it for 1000 times
## Identify candidate variables from Voting or Importance measures 

## Apply this model with the candidate predictors into Training set and Predict the outcome values in Testing set
## Repeat this for all replicates
## Calculate accuracy measures across Testing sets results. 
## Compare those accuracy measures for a few models
## Pick an optimal model 

## ----RF result summary---------------------------------------------------
WTC.RF.err.measure <- data.frame(AUC=unlist(lapply(WTC_RFresult, function(x){x$RF.err[1,"overall"]})),
                             Sensitivity= unlist(lapply(WTC_RFresult, function(x){x$RF.err[2,"overall"]})),
                             Specificity= unlist(lapply(WTC_RFresult, function(x){x$RF.err[3,"overall"]})),
                             Misclassification= unlist(lapply(WTC_RFresult, function(x){x$RF.err[4,"overall"]})))
WTC.RF.err.measure.melt <- melt(WTC.RF.err.measure)
names(WTC.RF.err.measure.melt) <- c("accuracy.measure", "value")
WTC.RF.err.measure.melt$accuracy.measure <- factor(WTC.RF.err.measure.melt$accuracy.measure, levels = c("AUC", "Sensitivity", "Specificity", "Misclassification"))

ggplot(WTC.RF.err.measure.melt) + geom_boxplot(aes(x=accuracy.measure, y=value*100, fill=accuracy.measure)) + ggtitle(paste0("Random Forest (K=", K, ", ntree=", ntree, ")")) + ylab("percentage(%)")

## ----Distribution of accuracy measures-----------------------------------
tempD <- ddply(WTC.RF.err.measure.melt, ~ accuracy.measure, summarise, mean= round(mean(value), digits = 4), sd=round(sd(value), digits = 4))
tempD$value <- paste0(tempD$mean, " (sd=", tempD$sd, ")")
tempD$mean <- tempD$sd <- NULL
DT::datatable(tempD, caption = paste0("Mean(SD) of Overall Accuracy Measures using RF (K=", K, ")"), rownames = FALSE)

## ----importance measures for RF------------------------------------------
all.important.measure <- lapply(WTC_RFresult, function(x){x$importance_measure[order(rownames(x$importance_measure)),]})
all.important.measure.DF <- do.call("cbind", all.important.measure)

all.important.measure.mean.decrease.acc.mean <- apply(all.important.measure.DF[, grep("Mean.Decrease.Accuracy", names(all.important.measure.DF))], 1, mean)
all.important.measure.mean.Decrease.Gini.mean <- apply(all.important.measure.DF[, grep("Mean.Decrease.Gini", names(all.important.measure.DF))], 1, mean)

important_measure_mean <- data.frame(Mean.Decrease.Accurary=round(all.important.measure.mean.decrease.acc.mean, 4), Mean.Decrease.Gini=round(all.important.measure.mean.Decrease.Gini.mean, 4))
important_measure_mean <- important_measure_mean[order(important_measure_mean$Mean.Decrease.Accurary, decreasing = TRUE),]

DT::datatable(important_measure_mean, caption = paste0("Mean of Importance Measures from Random Forest (K=", K, ")"))

## ----PR result-----------------------------------------------------------
## For each iteration, Penalized Regression run returns a matrix of ModelXerr.measures
## 1. Combine them across all K iterations. This returns a matrix of ModelxK dimension for each element
## 2. Plot the error.measure of AUC, Sensitivity, Sepecificity, and MC and faceted by Model
WTC.PR.err.measure <- list(AUC=do.call("rbind", lapply(WTC_PRresult, function(x){x$PR.err.vect["AUC.overall",]})),
                             Sensitivity= do.call("rbind", lapply(WTC_PRresult, function(x){x$PR.err.vect["sensitivity.overall",]})),
                             Specificity= do.call("rbind", lapply(WTC_PRresult, function(x){x$PR.err.vect["specificity.overall",]})),
                             Misclassification= do.call("rbind", lapply(WTC_PRresult, function(x){x$PR.err.vect["misclassification.overall",]})))

WTC.PR.err.measure <- lapply(WTC.PR.err.measure, function(x){as.data.frame(t(x), row.names = c("Elastic Net(alpha=0.4)", "Elastic Net(alpha=0.7)", "LASSO(alpha=1)"))})

WTC.PR.err.measure.DF <- do.call("cbind", WTC.PR.err.measure)
WTC.PR.err.measure.DF$Model <- rownames(WTC.PR.err.measure.DF)
WTC.PR.err.measure.DF.melt <- melt(WTC.PR.err.measure.DF, id="Model")
WTC.PR.err.measure.DF.melt$variable <- factor(unlist(lapply(strsplit(as.character(WTC.PR.err.measure.DF.melt$variable), "[.]"), function(x){x[[1]]})), levels=c("AUC", "Sensitivity", "Specificity", "Misclassification"))
WTC.PR.err.measure.DF.melt$Model <- factor(WTC.PR.err.measure.DF.melt$Model, levels=c("Elastic Net(alpha=0.4)", "Elastic Net(alpha=0.7)", "LASSO(alpha=1)"))

ggplot(WTC.PR.err.measure.DF.melt) + geom_boxplot(aes(x=variable, y=value*100,fill=variable)) + facet_grid(~ Model)+ ylab("percentage(%)") + ggtitle(paste0("Penalized Regressions (K=", K, ", Lasso/EN)")) + ylab("percentage(%)")

## 3. Mean and SD of error.measures for each Model
WTC.PR.err.measure.mean <- do.call("cbind", lapply(WTC.PR.err.measure, function(x){paste0(round(apply(x, 1, mean), 4), " (sd=", round(apply(x, 1, sd), 4), ")")}))
rownames(WTC.PR.err.measure.mean) <- c("Elastic Net(alpha=0.4)", "Elastic Net(alpha=0.7)", "LASSO(alpha=1)")
DT::datatable(WTC.PR.err.measure.mean, caption = paste0("Mean (SD) of Accuracy Measures using Penalized Regression, K=", K, ")"), rownames = TRUE)              

## ----Non zero coeffi for PR serum----------------------------------------
## 4. non zero coefficients and cound them how many times each variable selected
## For each model (alpha level), we order selected(non-zero coefficients) variables by their counts. 
WTC.PR.non0coeff <- list(alphaLevel_04=lapply(WTC_PRresult, function(x){x$non0coeff[[1]]}),
                     alphaLevel_07=lapply(WTC_PRresult, function(x){x$non0coeff[[2]]}),
                     alphaLevel_1=lapply(WTC_PRresult, function(x){x$non0coeff[[3]]}))
WTC.PR.non0coeff_variables <- list(alphaLevel_04=lapply(WTC_PRresult, function(x){names(x$non0coeff[[1]])}),
                     alphaLevel_07=lapply(WTC_PRresult, function(x){names(x$non0coeff[[2]])}),
                     alphaLevel_1=lapply(WTC_PRresult, function(x){names(x$non0coeff[[3]])}))
union_non0coeff <- data.frame(Variables=unique(unlist(lapply(WTC.PR.non0coeff_variables, function(x){unlist(x)}))))

union.all.non0coeff <- lapply(WTC.PR.non0coeff_variables, function(x){table(unlist(x))})
union.all.non0coeff <- lapply(union.all.non0coeff, function(x){merge(union_non0coeff, data.frame(Variables= names(x), Counts=as.vector(x)), by="Variables", all.x=TRUE)})
#union.all.non0coeff <- lapply(union.all.non0coeff, function(x){x[order(x$Counts, decreasing=TRUE),]})

union.all.non0coeff_all <- merge(merge(union.all.non0coeff[[1]], union.all.non0coeff[[2]], by="Variables"), union.all.non0coeff[[3]], by="Variables")
names(union.all.non0coeff_all) <- c("Variables", "alphaLevel=0.4", "alphaLevel=0.7", "alphaLevel=1")
union.all.non0coeff_all <- union.all.non0coeff_all[order(union.all.non0coeff_all$`alphaLevel=0.4`, decreasing = TRUE),]
union.all.non0coeff_all$Variables <- as.character(union.all.non0coeff_all$Variables)

DT::datatable(union.all.non0coeff_all[-1,], rownames = FALSE, caption = paste0("Selected Variables and its counts for each alpha (K=", K, ")"))

## Select top ranked genes having more than 70% out of K times selected in across all models ##
## Estimate these genes' coefficients in Logistic model and then use them as Polygenic score for the new testing data (Marine data)
polygenicGenes <- union.all.non0coeff_all$Variables[which(union.all.non0coeff_all$`alphaLevel=1`>=round(0.7*K))]
polygenicGenes <- polygenicGenes[polygenicGenes!="(Intercept)"]

## ----SVM result summary--------------------------------------------------
WTC.SVM.err.measure <- data.frame(AUC=unlist(lapply(WTC_SVMresult, function(x){x$SVM.err[1,"overall"]})),
                             Sensitivity= unlist(lapply(WTC_SVMresult, function(x){x$SVM.err[2,"overall"]})),
                             Specificity= unlist(lapply(WTC_SVMresult, function(x){x$SVM.err[3,"overall"]})),
                             Misclassification= unlist(lapply(WTC_SVMresult, function(x){x$SVM.err[4,"overall"]})))
WTC.SVM.err.measure.melt <- melt(WTC.SVM.err.measure)
names(WTC.SVM.err.measure.melt) <- c("accuracy.measure", "value")
WTC.SVM.err.measure.melt$accuracy.measure <- factor(WTC.SVM.err.measure.melt$accuracy.measure, levels = c("AUC", "Sensitivity", "Specificity", "Misclassification"))

ggplot(WTC.SVM.err.measure.melt) + geom_boxplot(aes(x=accuracy.measure, y=value*100, fill=accuracy.measure)) + ggtitle(paste0("Support Vector Machine (K=", K, ", Radial Kernal)")) + ylab("percentage(%)")

## ----SVM accuracy measures-----------------------------------------------
tempD <- ddply(WTC.SVM.err.measure.melt, ~ accuracy.measure, summarise, mean= round(mean(value), digits = 4), sd=round(sd(value), digits = 4))
tempD$value <- paste0(tempD$mean, " (sd=", tempD$sd, ")")
tempD$mean <- tempD$sd <- NULL
DT::datatable(tempD, caption = paste0("Mean(SD) of Overall Accuracy Measures using SVM (K=", K, ")"), rownames = FALSE)

## ----Neural Network result summary---------------------------------------
WTC.NN.err.measure <- data.frame(AUC=unlist(lapply(WTC_NNresult, function(x){x$NN.err[1,"overall"]})),
                             Sensitivity= unlist(lapply(WTC_NNresult, function(x){x$NN.err[2,"overall"]})),
                             Specificity= unlist(lapply(WTC_NNresult, function(x){x$NN.err[3,"overall"]})),
                             Misclassification= unlist(lapply(WTC_NNresult, function(x){x$NN.err[4,"overall"]})))
WTC.NN.err.measure.melt <- melt(WTC.NN.err.measure)
names(WTC.NN.err.measure.melt) <- c("accuracy.measure", "value")
WTC.NN.err.measure.melt$accuracy.measure <- factor(WTC.NN.err.measure.melt$accuracy.measure, levels = c("AUC", "Sensitivity", "Specificity", "Misclassification"))

ggplot(WTC.NN.err.measure.melt) + geom_boxplot(aes(x=accuracy.measure, y=value*100, fill=accuracy.measure)) + ggtitle(paste0("Neural Network (K=", K, ", 3 layers)")) + ylab("percentage(%)")

## ----NN  accuracy measures-----------------------------------------------
tempD <- ddply(WTC.NN.err.measure.melt, ~ accuracy.measure, summarise, mean= round(mean(value), digits = 4), sd=round(sd(value), digits = 4))
tempD$value <- paste0(tempD$mean, " (sd=", tempD$sd, ")")
tempD$mean <- tempD$sd <- NULL
DT::datatable(tempD, caption = paste0("Mean(SD) of Overall Accuracy Measures using SVM (K=", K, ")"), rownames = FALSE)

## ----Gredient Boosting result summary------------------------------------
GBM.err.measure <- data.frame(AUC=unlist(lapply(WTC_GBMresult, function(x){x$GBM.err[1,"overall"]})),
                             Sensitivity= unlist(lapply(WTC_GBMresult, function(x){x$GBM.err[2,"overall"]})),
                             Specificity= unlist(lapply(WTC_GBMresult, function(x){x$GBM.err[3,"overall"]})),
                             Misclassification= unlist(lapply(WTC_GBMresult, function(x){x$GBM.err[4,"overall"]})))
GBM.err.measure.melt <- melt(GBM.err.measure)
names(GBM.err.measure.melt) <- c("accuracy.measure", "value")
GBM.err.measure.melt$accuracy.measure <- factor(GBM.err.measure.melt$accuracy.measure, levels = c("AUC", "Sensitivity", "Specificity", "Misclassification"))

ggplot(GBM.err.measure.melt) + geom_boxplot(aes(x=accuracy.measure, y=value*100, fill=accuracy.measure)) + ggtitle(paste0("Gredient Boosting (K=", K, ")")) + ylab("percentage(%)")

## ----GBM  accuracy measures----------------------------------------------
tempD <- ddply(GBM.err.measure.melt, ~ accuracy.measure, summarise, mean= round(mean(value), digits = 4), sd=round(sd(value), digits = 4))
tempD$value <- paste0(tempD$mean, " (sd=", tempD$sd, ")")
tempD$mean <- tempD$sd <- NULL
DT::datatable(tempD, caption = paste0("Mean(SD) of Overall Accuracy Measures using SVM (K=", K, ")"), rownames = FALSE)

## ----Preprocessing the Marine data---------------------------------------
tdat <- Marine_raw            ## In total, n=27974 genes
sddat <- apply(tdat, 1, sd)
sdidx <- which(sddat>3)    ## 11585 genes removed
tdat <- tdat[sdidx,]        ## 16389 genes left 

## Filtering out genes having too extreme values and too small or too large variation
iqrdat <- apply( tdat , 1 , function( x ) diff(quantile( x , c( 0.25 , 0.75 ) ) ) )
iqridx <- which(iqrdat!=0)  ## 3723 genes removed 
tdat <- tdat[iqridx,]       ## 19736 gene left

iqrdat <- apply( tdat , 1 , function( x ) diff(quantile( x , c( 0.25 , 0.75 ) ) ) )
iqridx <- c(which(iqrdat < quantile(iqrdat, 0.15)), which(iqrdat > quantile(iqrdat, 0.85)))  ## 5431 genes have <15% or >85% Inter Quantile. 
tdat <- tdat[-iqridx,]   ## n=14305 genes left

feaD <- data.frame(gene=rownames(tdat))
rownames(feaD) <- feaD$gene

phenoData <- new("AnnotatedDataFrame", data=Marine_pheno)
featureData <- new("AnnotatedDataFrame", data=feaD)
esett <- ExpressionSet(assayData=as.matrix(tdat), phenoData=phenoData, featureData=featureData)

y=DGEList(esett)
DGE=calcNormFactors(y,method =c("TMM"))
barplot(DGE$samples$lib.size,names=colnames(DGE),las=2)

## ----PCA and MDS for Marine----------------------------------------------
## PCA run
pcaOut <- prcomp(t(DGE$counts), center=TRUE, scale=TRUE)
Marine_pheno_dat <- data.frame(pcaOut$x, Marine_pheno)
Marine_pheno_dat$case <- as.factor(Marine_pheno_dat$case)

pcaVariance <- round(unlist(lapply(1:length(pcaOut$sdev),function(i){pcaOut$sdev[i]^2/sum(pcaOut$sdev^2)})), digits=3)*100

ggplot(Marine_pheno_dat, aes(PC1, PC2, color=case)) + geom_point(size=3) +
  xlab(paste0("PC1: ", pcaVariance[1], "% variance")) + ylab(paste0("PC2: ", pcaVariance[2], "% variance"))+ggtitle("PC1 vs PC2 by Case")

ggplot(Marine_pheno_dat, aes(PC1, PC2, color=time)) + geom_point(size=3) +
  xlab(paste0("PC1: ", pcaVariance[1], "% variance")) + ylab(paste0("PC2: ", pcaVariance[2], "% variance")) +ggtitle("PC1 vs PC2 by Time")


Marine_pheno$PTSD <- ifelse(Marine_pheno$case=="1", "Yes", "No")
Marine_pheno$PTSD <- factor(Marine_pheno$PTSD, levels = c("No", "Yes"))


## ----Post Deployment data------------------------------------------------
## Split data into two Pre and Post data sets 
######################################### Post Deplyment data ############################
postMarine_pheno <- Marine_pheno[Marine_pheno$time=="Post",]
postMarine_raw <- Marine_raw[, rownames(postMarine_pheno)]

postTdat <- postMarine_raw     ## n=27974 genes 
postsddat <- apply(postTdat, 1, sd)
postsdidx <- which(postsddat!=0)        ## n=4916 genes having 0 sd
postTdat <- postTdat[postsdidx,]           ## 23058 genes left

## Filtering out genes having too extreme values and too small or too large variation
iqrdat <- apply( postTdat , 1 , function( x ) diff(quantile( x , c( 0.25 , 0.75 ) ) ) )
iqridx <- which(iqrdat >3)     ## n=3338 genes having 0 Interquantile
postTdat <- postTdat[iqridx,]          ## n=19720 gene left 

iqrdat <- apply( postTdat , 1 , function( x ) diff(quantile( x , c( 0.25 , 0.75 ) ) ) )
iqridx <- c(which(iqrdat < quantile(iqrdat, 0.15)), which(iqrdat > quantile(iqrdat, 0.85)))  ## 5480 genes have <15% or >85% Inter Quantile. 
postTdat <- postTdat[-iqridx,]   ## n=14240 genes left

feaD <- data.frame(gene=rownames(postTdat))
rownames(feaD) <- feaD$gene

phenoData <- new("AnnotatedDataFrame", data=postMarine_pheno)
featureData <- new("AnnotatedDataFrame", data=feaD)
postEsett <- ExpressionSet(assayData=as.matrix(postTdat), phenoData=phenoData, featureData=featureData)

posty=DGEList(postEsett)
postDGE=calcNormFactors(posty,method =c("TMM"))

design <- model.matrix(~ + case, postMarine_pheno)
postv <- voom(postDGE, design, plot=TRUE)

# postfit1 <- lmFit(postv, design)
# postfit1 <- eBayes(postfit1, trend=FALSE)
# 
# PostAllOut <- topTable(postfit1, number=nrow(postv), coef=1, sort="P",adjust.method = "BY")
# postsigallOut <- PostAllOut[PostAllOut$adj.P.Val < 0.05,]


## ----Only Post case------------------------------------------------------
## Select top ranked Genes from non-zero coefficients from Penalized regression model result 
polyGenesEN <- polygenicGenes
polyGenesEN <- polyGenesEN[polyGenesEN%in%rownames(postv$E)]   ## 21 out of 29 important genes from WTC data are included in Maine data
  
polyGenesDat <- cbind(WTC_pheno[, c("PTSD"), drop=FALSE], (t(v$E[rownames(v$E)%in%polyGenesEN,])))
polyGenesDat <- polyGenesDat[polyGenesDat$PTSD!="Past",]

polyGenesDat$PTSD <- as.factor(ifelse(polyGenesDat$PTSD=="Never", "No", "Yes"))

## ----logistic model fitting----------------------------------------------
polyGenesDatlogit <- glm(PTSD ~ ., data = polyGenesDat, family = "binomial")

MarinePolygenicDat<- cbind(Marine_pheno[rownames(postv$targets), c("time","PTSD")], t(postv$E[rownames(postv$E)%in%polyGenesEN,]))
MarinePolygenicDat$time <- NULL
MarinePolygenicDat$sampleID <- rownames(MarinePolygenicDat)

MarinePolygenic_Post <- predict(polyGenesDatlogit, newdata=MarinePolygenicDat, type="response")
MarinePolygenic_Post <- data.frame(polyScore=MarinePolygenic_Post, sampleID=names(MarinePolygenic_Post))
MarinePolygenic_Post <- merge(MarinePolygenic_Post, MarinePolygenicDat, by="sampleID")

MarinePolygenic_Post.melt <- melt(MarinePolygenic_Post, id=c("sampleID", "PTSD", "polyScore"))

ggplot(MarinePolygenic_Post.melt) + geom_boxplot(aes(x=PTSD, y=polyScore, fill=PTSD))  +
  ggtitle("Polygenic Expression Score for Marine by PTSD") + ylab("Polygenic Expression Score")

MarinePolyROC <- pROC::roc(MarinePolygenic_Post$PTSD, MarinePolygenic_Post$polyScore)
#roc1 <- pROC::ggroc(MarinePolyROC) 

plot(MarinePolyROC , col = 1, lty = 2, main = "Marine ROC curve", xlab="False Positive Rate (1-Specificity)", ylab="True Positive Rate (Sensitivity)")
text(1, 0.8, labels = paste0("AUC: ", round(MarinePolyROC$auc[1], digits = 3)))


## ----Random Forest model fitting-----------------------------------------
polyGenesDatRF <- randomForest(PTSD ~ ., data = polyGenesDat, importance = TRUE, ntree = ntree)
Marine_Post <- predict(polyGenesDatRF, MarinePolygenicDat)
Marine_Post <- data.frame(Marine_Pred=Marine_Post, sampleID=names(Marine_Post))
Marine_Post <- merge(Marine_Post, MarinePolygenicDat, by="sampleID")
Marine_Post$Marine_Pred <- factor(ifelse(Marine_Post$Marine_Pred=="No", 0, 1), ordered = TRUE)
Marine_Post$PTSD <- as.factor(ifelse(Marine_Post$PTSD=="No", 0, 1))

MarinePolyROC <- pROC::roc(response=Marine_Post$PTSD, predictor=Marine_Post$Marine_Pred)
#roc1 <- pROC::ggroc(MarinePolyROC) 

plot(MarinePolyROC , col = 1, lty = 2, main = "Marine ROC curve", xlab="False Positive Rate (1-Specificity)", ylab="True Positive Rate (Sensitivity)")
text(1, 0.8, labels = paste0("AUC: ", round(MarinePolyROC$auc[1], digits = 3)))


