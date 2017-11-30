##By default, DESeq2 uses sample-specific scaling factors.
##Idea: Allow bait-specific scaling factors. Show this improves model fit.
##      May lead to a drop in the dispersion.

library(data.table)
library(DESeq2)

logit <- function(p) {p/1-p}
expit <- function(x) {1/(1+exp(-x))}
geoMean <- function(x, na.rm=FALSE) {
  if(na.rm){
    exp(mean(log(x[!is.na(x)])))
  }else{
    exp(mean(log(x)))
  }
}

##Input data:

fragData <- readRDS(file.path(filePath, "04recast.Rds"))
setkey(fragData, otherEndID)


#impute=TRUE
#if(impute)
#{
##Option1) impute (i.e. s_j: NA -> 1)
fragData.impute <- copy(fragData)
fragData.impute[,s_j.impute := geoMean(s_j, na.rm = TRUE), by=c("baitID", "otherEndID")]
fragData.impute[is.na(s_j), s_j:=s_j.impute]
fragData.impute[,FullMean.impute := geoMean(FullMean, na.rm = TRUE), by=c("baitID", "otherEndID")]
fragData.impute[is.na(FullMean), FullMean:=FullMean.impute]
#} else {
##Option2) discard (i.e. is.na(s_j) -> N = NA)
fragData[is.na(s_j),N := as.integer(NA)]
##cheeky attempt to bypass this
fragData[is.na(s_j),N := 0]
fragData[is.na(s_j),s_j := 0]
#}
##NOTE: Standard normalization factors from DESeq2 used for certain sites later.

##construct an appropriate DESeq object

regionData <- fragData[,list(
  N=sum(N),
  avDist=(min(distSign)+max(distSign))/2,
  Bmean=sum(Bmean),
  Tmean=sum(Tmean),
  FullMean=sum(FullMean),
  s_j=s_j[1]
),by=c("baitID", "regionID", "sample")]
regionData.impute <- fragData.impute[,list(
  N=sum(N),
  avDist=(min(distSign)+max(distSign))/2,
  Bmean=sum(Bmean),
  Tmean=sum(Tmean),
  FullMean=sum(FullMean),
  s_j=s_j[1]
),by=c("baitID", "regionID", "sample")]

nSamples <- length(unique(regionData$sample))
       
regionDataMatrix <- matrix(regionData.impute$N,
                           ncol=nSamples,
                           byrow = TRUE)
colnames(regionDataMatrix) <- unique(regionData.impute$sample)
colData <- data.frame(condition = fragData$condition[1:ncol(regionDataMatrix)])
dds <- DESeqDataSetFromMatrix(countData = regionDataMatrix,
                              colData = colData,
                              design = ~ condition)

dds.nullModel <- estimateSizeFactors(dds)
dds.M3 <- copy(dds.nullModel)

nullSizeFactors <- sizeFactors(dds.nullModel)

##model 1) standard model
dds.nullModel <- estimateDispersions(dds.nullModel)
dds.nullModel <- nbinomWaldTest(dds.nullModel)

# ##model 2) multiplicative size factors: s_jk = s_j*s_k
# normFactorsSk <- matrix(rep(sizeFactors(dds), nrow(regionData)/nSamples), ncol=nSamples, byrow = TRUE)
# normFactorsSj <- matrix(regionData.impute$s_j, ncol=nSamples, byrow = TRUE)
# 
# normFactors <- normFactorsSj*normFactorsSk
# normFactors <- normFactors / exp(rowMeans(log(normFactors))) ##<- normalize to get geometric mean = 1
# normalizationFactors(dds) <- normFactors
# 
# dds <- estimateDispersions(dds)
# dds <- nbinomWaldTest(dds)

##model 3) use FullMean instead of s_j
##should make little difference at small distance, but should dampen out effect of s_j at large distance
normFactorsFull <- matrix(regionData.impute$FullMean, ncol=nSamples, byrow = TRUE)

normFactorsM3 <- normFactorsFull
normFactorsM3 <- normFactorsM3 / exp(rowMeans(log(normFactorsM3))) ##<- normalize to get geometric mean = 1
##remove NA rows
selNA <- apply(normFactorsM3, 1, function(x){any(is.na(x))})
normFactorsM3[selNA,] <- rep(nullSizeFactors, each=sum(selNA))
normalizationFactors(dds.M3) <- normFactorsM3

dds.M3 <- estimateDispersions(dds.M3)
dds.M3 <- nbinomWaldTest(dds.M3)

##get annotation information -----------


##other ends:
rmap <- fread(rmapfile)
colnames(rmap) <- c("OEchr", "OEstart", "OEend", "otherEndID")

RU <- readRDS(file.path(filePath, "01regionUniverse.Rds"))
RUsummary <- RU[,
  list(
    baitID=baitID[1],
    minOE=min(otherEndID),
    maxOE=max(otherEndID)
  ), by = "regionID"]
annoData <- merge(RUsummary, rmap[,c("otherEndID", "OEchr", "OEstart"), with=FALSE], by.x="minOE", by.y="otherEndID")
annoData <- merge(annoData, rmap[,c("otherEndID", "OEend"), with=FALSE], by.x="maxOE", by.y="otherEndID")

#baits
colnames(rmap) <- c("baitchr", "baitstart", "baitend", "baitID")
annoData <- merge(annoData, rmap, by.x="baitID", by.y="baitID")


setkey(annoData, regionID)
stopifnot(identical(1:nrow(annoData), annoData$regionID))

# out <- results(dds.M3)
# saveRDS(out, file=file.path(filePath, "05outPilot.Rds"))
# 
# out.nullModel <- results(dds.nullModel)
# saveRDS(out.nullModel, file=file.path(filePath, "05outStandardDESeq.Rds"))

out <- cbind(
  as.data.table(as.data.frame(results(dds.M3))),
  annoData
)
saveRDS(out, file=file.path(filePath, "05outPilot.Rds"))

out.nullModel <- cbind(
  as.data.table(as.data.frame(results(dds.nullModel))),
  annoData
)
saveRDS(out.nullModel, file=file.path(filePath, "05outStandardDESeq.Rds"))
