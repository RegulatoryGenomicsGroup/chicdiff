##Given a region universe, collect the number of reads in each region from the chinputs.

##Inputs and parameters:

##requires targetChs

RUfile <- file.path(filePath, "01regionUniverse.Rds")

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

##1) Read in the region universe --------------------

RU <- readRDS(RUfile)
baits <- sort(unique(RU$baitID)) ##baits to get information from

##2) Collect all of the read count information ----------------

countData <- vector("list", length(targetChFiles))
gc()
for(i in 1:length(countData))
{
  #message(i, "\n")
  x <- fread(targetChFiles[i])
  setkey(x, baitID)
  #message(key(x), " ")
  x <- x[J(baits),]
  setkey(x, baitID) ##Required for countData[4:7], where key is dropped - not sure why
  #message(key(x), " \n")
  countData[[i]] <- x
  rm(x)
}

##2a) For each region, count the number of reads in each replicate -----------------
##(FIXME: could be optimized a lot)

x <- RU
setkey(x, baitID, otherEndID)
for(i in 1:length(countData))
{
  columnName <- paste0("N.", names(targetChFiles)[i])
  temp <- countData[[i]][,c("baitID","otherEndID","N"), with=FALSE]
  setkey(temp, baitID, otherEndID)
  x <- merge(x, temp, all.x=TRUE)
  
  x[is.na(N), N:=0]
  setnames(x, old="N", new=columnName)
}

##3) Collect distance -----------------

##3a) Distance ------------
rmap <- Chicago:::.readRmap(list(rmapfile=rmapfile))
colnames(rmap) <- c("chr","start","end","fragID")
setkey(rmap, fragID)
rmap[,midpoint:=round(0.5*(start+end))]
rmap[,start:=NULL]
rmap[,end:=NULL]
x <- merge(x, rmap, by.x="otherEndID", by.y="fragID")
x <- merge(x, rmap, by.x="baitID", by.y="fragID")
setkey(x, baitID, otherEndID)

x[, distSign :=
    ifelse(chr.x == chr.y,
           midpoint.x - midpoint.y,
           NA)]
x[, c("chr.x","midpoint.x","chr.y","midpoint.y"):=NULL]

saveRDS(x, file = file.path(filePath, "02CountOut.Rds"))