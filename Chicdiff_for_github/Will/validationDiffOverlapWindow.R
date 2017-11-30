##Count overlap of enhancers in differential other ends
setwd("/bi/group/sysgen/chicdiff/")
source("Pipeline/00ConfigMonVsTCD4ActWindowCheck.R")

library(AnnotationHub)
library(rtracklayer)
library(caTools)

strict = FALSE

##1) get the regulatory build annotation

regBuildAnno <- fread("zcat Resources/regBuild/RegBuild_BLUEPRINT_annotations.txt.gz")
##discard unnecessary information
regBuildAnno <- regBuildAnno[CellType %in% c("Mon", "nCD4"),]
setkey(regBuildAnno, "BaitID", "OtherEndID")

##reduce only to regions observed in mon/cd4
sel <- regBuildAnno[,sum(SignifInteraction) > 0, by=c("BaitID", "OtherEndID")] ##TRUEs are required
sel <- sel[V1 == TRUE,]
sel[,V1:=NULL]
setkey(sel, "BaitID", "OtherEndID")
regBuildAnno <- regBuildAnno[sel]

##Function assessing if there is a differential enhancer
#TRUE if so
checkEnh <- function(IDs, acts, strict)
{
  if("." == IDs[1])
  {
    return(FALSE) ##no enhancers
  }else{
    #cat(IDs, "\n")
    #cat(acts, "\n")
    type <- strsplit(IDs[1], split = c("[_,0-9]+"))[[1]]
    annotAct <- strsplit(acts, split = ",")
    
    annotAct <- lapply(annotAct, as.integer)
    annotAct <- as.data.table(annotAct)
    annotAct <- cbind(annotAct, type)
    
    annotAct <- annotAct[type %in% c("distal","tss","proximal"),]
    
    ##narrow it down
    if(strict)
    {
      annotAct <- annotAct[(V1 == V2 & V3 == V4 & V2 != V3 & V2 < 2 & V3 < 2),] ##1,1,0,0 or 0,0,1,1
      # annotAct <- annotAct[(V1 == V2 & V3 == V4 & V2 != V3 & V1%%2 == 1 & V3%%2 == 1,] ##1,1,3,3 or 3,3,1,1
    } else {
      annotAct <- annotAct[(V1 == V2 & V3 == V4 & V2 != V3 & V2 < 2 & V3 < 2)|(((V1 + V2 + V3 + V4) %in% c(1,3)) & V1 < 2 & V2 < 2 & V3 < 2 & V4 < 2),] ##1,1,0,0 or 0,0,1,1 or just one (1,0,0,0; etc)
    }
    return(nrow(annotAct) > 0)}
}

toCheck <- regBuildAnno[,checkEnh(OtherEndRegFeatIDs, OtherEndRegFeatActivity, strict), by=c("BaitID", "OtherEndID")]
##subset reg annotation
toCheck <- toCheck[V1 == TRUE,]
toCheck[,V1:=NULL]
setkey(toCheck, BaitID, OtherEndID)

##0 - dead, 1 - active, 2 - poised, 3 - Polycomb-repressed
##Start with dead<->active (0<->1) and count how many regions we get.

##2) get the calls
weightedpvalues <- readRDS(file=file.path(filePath, "06Weighted.Rds"))
out <- readRDS(file=file.path(filePath, "05outPilot.Rds"))
out.null <- readRDS(file=file.path(filePath, "05outStandardDESeq.Rds"))
out$pvalue.null <- out.null$pvalue
out$padj.null <- out.null$padj
selwpv <- order(weightedpvalues$regionID)
weightedpvaluesRID <- weightedpvalues[selwpv,]
out$weighted_pvalue <- weightedpvaluesRID$weighted_pvalue

#   ##add pvals from t-test
# pvals <- readRDS(file=file.path(filePath, "05ExpectedTtestPvals.Rds"))
# out$expDiffPval <- pvals[out$regionID]

##3) check which regionIDs have a diff enhancer inside them
IRref <- IRanges(start = unique(toCheck$OtherEndID), width=1)
IRquery <- IRanges(start = out$minOE, end = out$maxOE)

temp <- findOverlaps(IRref, IRquery)
outTags <- unique(as.data.frame(temp)$subjectHits)
out[,isDiffEnh := regionID %in% outTags]
out.null[,isDiffEnh := regionID %in% outTags]

saveRDS(out, "preroccurveoutww_CheckLenient.Rds")
saveRDS(out.null, "preroccurveoutnullww_CheckLenient.Rds")
# 
# ##4) roc curve
out <- readRDS("preroccurveoutww_CheckLenient.Rds")
out.null <- readRDS("preroccurveoutnullww_CheckLenient.Rds")

sumTRUE <- sum(out$isDiffEnh)
sumFALSE <- sum(!out$isDiffEnh)
sel <- order(out$pvalue)
out.ordered <- out[sel,]
xvals <- cumsum(!out.ordered$isDiffEnh)/sumFALSE
yvals <- cumsum(out.ordered$isDiffEnh)/sumTRUE
plot(xvals, yvals, type="l", col="blue",
     main = "ROC curve for differential enhancers",
     xlab="False positive rate",
     ylab="True positive rate")
abline(0,1)
sel.null <- order(out$pvalue.null)
out.ordered.null <- out[sel.null,]
xvals.null <- cumsum(!out.ordered.null$isDiffEnh)/sumFALSE
yvals.null <- cumsum(out.ordered.null$isDiffEnh)/sumTRUE
lines(xvals.null, yvals.null, col="red")
orderedrows <- order(out$weighted_pvalue)
out.orderedW <- out[orderedrows,]
xvalsW <- cumsum(!out.orderedW$isDiffEnh)/sumFALSE
yvalsW <- cumsum(out.orderedW$isDiffEnh)/sumTRUE
lines(xvalsW, yvalsW, col="orange")
legend("bottomright", pch=19, col = c("red", "blue", "orange"), c("Null DESeq2 model", "Norm factor model", "Norm factor model + weighting"))

##get AUCs
trapz(xvals, yvals)
trapz(xvals.null, yvals.null)

##Will's test area

##No. TRUE/ No. FALSE plots 
#Linear vs linear
orderedrows <- order(out$pvalue)
out.ordered <- out[orderedrows,]
yvalues <- cumsum(head(out.ordered$isDiffEnh, 10000))/cumsum(head(!out.ordered$isDiffEnh, 10000))
xvalues <- 1:10000
plot(xvalues, yvalues, type="l", col="blue",
     main = "No. TRUE/No. FALSE vs No. called",
     xlab="No. called",
     ylab="No. TRUE/No. FALSE")
orderedrows.null <- order(out$pvalue.null)
out.ordered.null <- out[orderedrows.null,]
yvalues.null <- cumsum(head(out.ordered.null$isDiffEnh, 10000))/cumsum(head(!out.ordered.null$isDiffEnh, 10000))
xvalues.null <- 1:10000
lines(xvalues.null, yvalues.null, col="red")
legend("topright", pch=19, col = c("red", "blue"), c("Null DESeq2 model", "Norm factor model"))

#Log vs linear
orderedrows <- order(out$weighted_pvalue)
out.ordered <- out[orderedrows,]
yvalues <- cumsum(out.ordered$isDiffEnh)/cumsum(!out.ordered$isDiffEnh)
xvalues <- log(1:length(out.ordered$isDiffEnh))
plot(xvalues, yvalues, type="l", col="blue",
     main = "Odds vs log(rank) with weighting - Lenient",
     xlab="log(rank)",
     ylab="No. TRUE/No. FALSE")
orderedrows.null <- order(out$pvalue.null)
out.ordered.null <- out[orderedrows.null,]
yvalues.null <- cumsum(out.ordered.null$isDiffEnh)/cumsum(!out.ordered.null$isDiffEnh)
xvalues.null <- log(1:length(out.ordered.null$isDiffEnh))
lines(xvalues.null, yvalues.null, col="red")
legend("topright", pch=19, col = c("red", "blue"), c("Null DESeq2 model", "Norm factor model"))

#Log vs log
orderedrows <- order(out$pvalue)
out.ordered <- out[orderedrows,]
yvalues <- log(cumsum(head(out.ordered$isDiffEnh, 10000))/cumsum(head(!out.ordered$isDiffEnh, 10000)))
xvalues <- log(1:10000)
plot(xvalues, yvalues, type="l", col="blue",
     main = "log(No. TRUE/No. FALSE) vs No. called",
     xlab="log(No. called)",
     ylab="log(No. TRUE/No. FALSE)")
orderedrows.null <- order(out$pvalue.null)
out.ordered.null <- out[orderedrows.null,]
yvalues.null <- log(cumsum(head(out.ordered.null$isDiffEnh, 10000))/cumsum(head(!out.ordered.null$isDiffEnh, 10000)))
xvalues.null <- log(1:10000)
lines(xvalues.null, yvalues.null, col="red")
legend("topright", pch=19, col = c("red", "blue"), c("Null DESeq2 model", "Norm factor model"))

##Percentage plots
#Linear vs linear
orderedrows <- order(out$pvalue)
out.ordered <- out[orderedrows,]
xvalues <- 1:10000
yvalues <- cumsum(head(out.ordered$isDiffEnh, 10000))/xvalues
orderedrows.null <- order(out$pvalue.null)
out.ordered.null <- out[orderedrows.null,]
xvalues.null <- 1:10000
yvalues.null <- cumsum(head(out.ordered.null$isDiffEnh, 10000))/xvalues.null
plot(xvalues.null, yvalues.null, type="l", col="red",
     main = "Fraction TRUE in no. called so far vs No. called",
     xlab="No. called",
     ylab="Fraction TRUE")

lines(xvalues, yvalues, col="blue")
legend("topright", pch=19, col = c("red", "blue"), c("Null DESeq2 model", "Norm factor model"))

#Log vs linear
orderedrows <- order(out$pvalue)
out.ordered <- out[orderedrows,]
xvalues <- log(1:50000)
yvalues <- cumsum(head(out.ordered$isDiffEnh, 50000))/(1:50000)
plot(xvalues, yvalues, type="l", col="blue",
     main = "log(Fraction TRUE) vs No. called",
     xlab="log(No. called)",
     ylab="Fraction TRUE")
orderedrows.null <- order(out$pvalue.null)
out.ordered.null <- out[orderedrows.null,]
xvalues.null <- log(1:50000)
yvalues.null <- cumsum(head(out.ordered.null$isDiffEnh, 50000))/(1:50000)
lines(xvalues.null, yvalues.null, col="red")
legend("topright", pch=19, col = c("red", "blue"), c("Null DESeq2 model", "Norm factor model"))

#Log vs log
orderedrows <- order(out$pvalue)
out.ordered <- out[orderedrows,]
xvalues <- log(1:10000)
yvalues <- log(cumsum(head(out.ordered$isDiffEnh, 10000))/(1:10000))
plot(xvalues, yvalues, type="l", col="blue",
     main = "log(Fraction TRUE) vs No. called",
     xlab="log(No. called)",
     ylab="log(Fraction TRUE)")
orderedrows.null <- order(out$pvalue.null)
out.ordered.null <- out[orderedrows.null,]
xvalues.null <- log(1:10000)
yvalues.null <- log(cumsum(head(out.ordered.null$isDiffEnh, 10000))/(1:10000))
lines(xvalues.null, yvalues.null, col="red")
legend("topright", pch=19, col = c("red", "blue"), c("Null DESeq2 model", "Norm factor model"))

##5) get some regionIDs that are called with normFactors,
##   but not with DESeq2, and have diffEnh.
out[,interesting := (!is.na(padj) & padj < 0.05) & (is.na(padj.null) | padj.null > 0.05) & isDiffEnh]
table(out$interesting)
test <- out[(interesting),]
sel <- order(test$padj)
test[sel,]
##top - 39004, 75260, 66185, 39062, 66193

# ##6) can we predict diff enhancer status based on diffl expectation?
# plot(out.ordered$expDiffPval[1:10000],
#      col=ifelse(out.ordered$isDiffEnh, "blue", "black"),
#      pch=ifelse(out.ordered$isDiffEnh, ".", "X"))
