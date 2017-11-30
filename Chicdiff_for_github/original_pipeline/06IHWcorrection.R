library(IHW)
library(data.table)
library(ggplot2)

##Retrieve the appropriate output files first

getOutput <- function(myfilePath)
{
  test <- readRDS(file.path(myfilePath, "05outPilot.Rds"))
  RU.recast <- readRDS(file.path(myfilePath,"04recast.Rds"))
  RU.distances <- RU.recast[,list(avDist = mean(distSign)),by="regionID"]
  
  test$avDist <- RU.distances$avDist
  
  ##Comparison against genuinely uniform p-vals
  test$uniform <- runif(nrow(test))
  test$shuff <- sample(test$pvalue)
  
  ##Convert
  as.data.frame(test)
}

source("Pipeline/00ConfigMonVsTCD4ActWindow.R")
out <- getOutput(filePath)
setDT(out)
out.control <- getOutput(file.path(filePath, "control"))

##Train weights on the control sample
ihwRes <- ihw(pvalue ~ abs(avDist),  data = out.control, alpha = 0.05)
plot(ihwRes)
plot(ihwRes, what = "decisionboundary")

##---------------------------------------------------------------------
##Learn distance dependency

test <- ihwRes@df
setDT(test)
test[,log(covariate), by="group"]

distLookup <- test[,list(
  avgLogDist=mean(log(covariate)),
  minLogDist=min(log(covariate)),
  maxLogDist=max(log(covariate))
    ), by="group"]
distLookup <- distLookup[!is.na(group),]
setkey(distLookup, group)
if(distLookup[,!identical(as.integer(group), seq_along(group))])
{
  stop("Assumption violated")
}
w <- ihwRes@weights
distLookup[,group := as.integer(group)]
distLookup$avWeights <- rowSums(w)/ncol(w)
distLookup$minLogDist[1] = 0
distLookup$maxLogDist[nrow(distLookup)] = Inf

##--------------------------------------------------------------
##Apply to test data

out[,avgLogDist := log(abs(avDist))]
breaks <- (c(distLookup$minLogDist, Inf) + c(0, distLookup$maxLogDist))/2
out$group <- as.integer(cut(log(abs(out$avDist)), breaks))
##TODO: check that group & avDist correlate
setkey(distLookup, group)
setkey(out, group)

out <- merge(out, distLookup[,c("group", "avWeights"),with=FALSE], all.x=TRUE, by="group")
out$weight <- out$avWeights/mean(out$avWeights) ##renormalize
out[,weighted_pvalue := pvalue/weight]
saveRDS(out, file=file.path(filePath, "06Weighted.Rds"))
