##Synthesise the count + interaction parameter data together

##inputs:

##requires targetRDAs
x <- readRDS(file.path(filePath, "02CountOut.Rds"))
params <- readRDS(file.path(filePath, "03interactionParameters.Rds"))

setkey(x, baitID, otherEndID)

for(i in 1:length(targetRDAFiles))
{
  s <- names(targetRDAFiles[i])
  temp <- params[[i]][,c("baitID", "otherEndID", "s_j", "Bmean", "Tmean", "score"), with=FALSE]
  temp <- unique(temp)
  temp[,FullMean:=Bmean+Tmean]
  #temp[,Bmean:=NULL]
  #temp[,Tmean:=NULL]
  for(nm in colnames(temp))
    if(!nm %in% c("baitID", "otherEndID"))
      {
        setnames(temp, nm, paste0(nm, ".", s))
      }
  
  setkey(temp, baitID, otherEndID)
  
  x <- merge(x, temp, all.x=TRUE)
  #x <- merge(x, temp, all.x=TRUE, allow.cartesian=TRUE)
}

##rework table to have columns sample, N, FullMean, etc
## (rather than N.[sample], FullMean.[sample])
myCols <- c("N", "s_j", "Bmean", "Tmean", "score", "FullMean")
getCols <- function(col) {grep(paste0("^", col), colnames(x), perl=TRUE)}
sels <- lapply(myCols, getCols)

recast = melt(x, measure = sels, value.name = myCols,
                    variable.name="sample")
recast$sample <- names(targetRDAFiles)[recast$sample]

recast$condition <- rep(
    rep(names(targetRDAs), sapply(targetRDAs, length))
  , each=nrow(x))


setkey(recast, regionID)


saveRDS(recast, file = file.path(filePath, "04recast.Rds"))
