##Putting CHiCAGO results together to get regions large enough for diff binding

##Functions

getFile <- function(directory, suffix)
{
  candidates <- dir(directory, suffix)
  if(length(candidates) == 1)
  {
    out <- file.path(directory, candidates)
    return(out)
  }else{
    stop("Found ", length(candidates), " possible files. Require 1.")
  }
}

runDiff <- function(x)
{
  n <- length(x)
  if(n >= 2)
  {
    x[2:n] - x[1:(n-1)]
  } else {
    return(NULL)
  }
}

##e.g. xSeed <- c(12,16,25,40); s <- 5

expand <- function(xSeed, s)
{
  ##Return any integer that is within s of a "seed" in x
  temp <- IRanges(start=xSeed, width=1)
  temp <- temp + s
  reduce(temp)
}

expandAvoidBait <- function(bait, oe, s)
{
  ##Function usually returns (oe-s):(oe+s) - range around oe
  ##But if bait too close to oe, stop the range short.
  if(abs(bait - oe) > s+1)
  {
    return((oe - s):(oe + s))
  } else if (oe > bait) {
    return((bait + 2L):(oe + s))
  } else if (oe < bait){
    return((oe - s):(bait - 2L))
  } else {
    stop("Invalid parameters bait=",bait," oe=",oe," s=",s)
  }
}

##1) Read in CHiCAGO results for conditions of interest (subset for now)

x <- fread(peakMatrix)

sel <- which(colnames(x) %in% targetColumns)
x <- x[,c(1:11, sel), with=FALSE]

##Get any rows where at least one score > 5.
sel <- rep(FALSE, nrow(x))
for(cl in targetColumns)
{
  sel <- sel | x[,get(cl) > 5]
}
x <- x[sel,]

x <- x[!is.na(dist),] ## <- FILTER OUT TRANS INTERACTIONS - Assumed in (*)

##2) Expand the "point universe" to get "region universe" ------------------
if(RUmergeMode == "union")
{
  ##Stategy 2a): Fill in gaps with size < s.
  ##Reduce to small subset for now:
  baits <- unique(x$baitID)
  baits <- head(baits, 5000)
  x <- x[baitID %in% baits,]

  regionUniverse <- vector("list", length(baits))
  setkey(x, baitID)
  for(i in seq_along(baits))
  {
    bait <- baits[i]
    x.sel <- x[J(bait),]
    regionUniverse[[i]] <- expand(x.sel$oeID, RUexpand) ##RUexpand should be set in 00[...].R
  }
  
  #Construct RU.DT: a data.table containing the region universe
  temp <- sapply(regionUniverse, function(x){sum(width(x))})
  RU.DT <- data.table(
    baitID=rep(baits, temp)
  )
  RU.DT$otherEndID <- unlist(
    lapply(regionUniverse, unlist)
  )
  temp <- unlist(sapply(regionUniverse, width))
  RU.DT$regionID <- rep(1:length(temp), temp)
  
  ##delete self-interaction/re-ligation
  #sum(RU.DT$baitID == RU.DT$otherEndID)
  RU.DT <- RU.DT[abs(baitID - otherEndID) > 1,]
  
} else if (RUmergeMode == "window") {
  ##Stategy 2b): Just expand calls by s in each direction
  ##Reduce to small subset for now:
  baits <- unique(x$baitID)
  baits <- head(baits, 5000)
  x <- x[baitID %in% baits,]
  
  #RU.DT <- data.table(
  #  baitID=rep(baits, each=2*s + 1)
  #)
  RU.DT <- x[,c("baitID", "oeID"),with=FALSE]
  RU.DT$regionID <- 1:nrow(RU.DT)
  RU.DT <- RU.DT[,
                 list(otherEndID = expandAvoidBait(baitID, oeID, RUexpand)), ##RUexpand should be set in 00[...].R
                 by=c("baitID","regionID")
                ]
} else {
  stop("invalid RUmergeMode")
}

##3a) Final filters: Require that the average CHiCAGO score in a region is sufficiently high -------------

##FIXME

## b) Ensure regions remain within the genome
maxfrag <- max(fread(rmapfile)$V4)
RU.DT <- RU.DT[otherEndID <= maxfrag,]

## c) Ensure that regions stay on the correct chromosome (*)

## get chr info
rmap <- fread(rmapfile)
colnames(rmap) <- c("chr", "start", "end", "ID")
setkey(RU.DT, otherEndID)
setkey(rmap, ID)
RU.DT <- rmap[,c("chr", "ID"),with=FALSE][RU.DT]
setnames(RU.DT, "ID", "otherEndID")
setkey(RU.DT, baitID)
RU.DT <- rmap[,c("chr", "ID"),with=FALSE][RU.DT]
setnames(RU.DT, "ID", "baitID")

## keep only cis interactions
RU.DT <- RU.DT[chr == i.chr,]
RU.DT[,chr:=NULL]
RU.DT[,i.chr:=NULL]

##4) Output the region universe --------------------------------------

#saveRDS(RU.DT, file = "Merging/regionUniverse_Simple.Rds")
saveRDS(RU.DT, file = file.path(filePath, "01regionUniverse.Rds"))