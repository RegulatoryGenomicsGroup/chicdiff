##Given a region universe, collect Bmean + Tmean.
##Note: Some baits missing (filtered out from <250 reads)

##Inputs and parameters:

##requires targetRDAs
RUfile <- file.path(filePath, "01regionUniverse.Rds")

##Gutted version of Chicago function
estimateDistFun <- function (x, settings=defaultSettings()) {
  # Take the "refBinMean" column of the data x as f(d_b)
  # then interpolate & extrapolate to get f(d).
  
  # Get f(d_b)
  setkey(x, distbin, refBinMean)
  f.d <- unique(x, by=key(x))[is.na(refBinMean)==FALSE][, c("distbin", "refBinMean"), with=FALSE]
  f.d <- f.d[order(refBinMean, decreasing=TRUE)]
  
  setDF(f.d) # f.d is tiny, so no need to bother with it being a data.table
  f.d$midpoint <- seq(from=round(settings$binsize/2), by=settings$binsize, length.out=nrow(f.d))
  
  obs.min <- log(min(f.d$midpoint))
  obs.max <- log(max(f.d$midpoint))
  
    ##Spline - Cubic fit over observed interval, linear fit elsewhere, assume continuity of f(d) & f'(d).
    distFunParams <- list(method="cubic")
    
    ##cubic fit (quadratic not immensely different TBH)
    f.d.cubic <- lm(log(refBinMean) ~ log(midpoint) + I(log(midpoint)^2) + I(log(midpoint)^3), data = f.d)
    fit <- f.d.cubic$coefficients
    distFunParams[["cubicFit"]] <- fit
    
    ##Extrapolation: Fit the "head" and "tail" of f using continuity
    distFunParams[["obs.min"]] <- log(min(f.d$midpoint))
    distFunParams[["obs.max"]] <- log(max(f.d$midpoint))
    
    beta <- fit[2] + 2*fit[3]*c(obs.min, obs.max) + 3*fit[4]*(c(obs.min, obs.max)^2)
    alpha <- fit[1] + (fit[2] - beta)*c(obs.min, obs.max) + fit[3]*c(obs.min, obs.max)^2 + fit[4]*c(obs.min, obs.max)^3
    
    distFunParams[["head.coef"]] <- c(alpha[1], beta[1])
    distFunParams[["tail.coef"]] <- c(alpha[2], beta[2])
    

  distFunParams
}

RU <- readRDS(RUfile)
RU[,regionID:=NULL]
RU <- unique(RU)

##prep output
outData <- vector("list", length(targetRDAFiles))
dispersions <- numeric(length(targetRDAFiles))

##1) Read in each file...
for(i in 1:length(targetRDAFiles))
{
  message("File ",i, " of ", length(targetRDAFiles))
  file <- targetRDAFiles[i]
  load(file) ##this produces the object 'x'
  
  ##1a) collect the dispersion
  dispersions[i] <- attributes(x)$dispersion
  
  ##1b) collect the Bmean, Tmean information that is present
  setDT(x)
  setkey(x, baitID, otherEndID)
  setkey(RU, baitID, otherEndID)
  out <- x[RU, c("baitID", "otherEndID","distSign","Bmean", "Tmean", "score"), with=FALSE]
  ##(note that score=NA occurs when N=0 - hence replace score=0)
  
  ##1c) recalculate distSign (need to do this for control regions)
  if(any(is.na(out$distSign)))
  {
    rmap <- Chicago:::.readRmap(list(rmapfile=rmapfile))
    colnames(rmap) <- c("chr", "start", "end", "ID")
    setkey(rmap, "ID")
    temp <- merge(RU, rmap, all.x=TRUE, by.x="baitID", by.y="ID")
    temp <- merge(temp, rmap, all.x=TRUE, by.x="otherEndID", by.y="ID")
    setkey(temp, baitID, otherEndID)
    temp[,distSign:=round(((start.y + end.y)-(start.x + end.x))/2)]
    if(any(abs(temp$distSign - out$distSign) > 1, na.rm=TRUE))
    {
      stop("Error calculating distances.")
    }
    out$distSign <- temp$distSign
  }

  ##1d) bait annotation: collect s_js, tblb
  ##Note that s_j can be NA (for baits that were filtered out!)
  temp <- x[,list(s_j=s_j[1], tblb=tblb[1]),by="baitID"]
  setkey(temp, baitID)
  setkey(out, baitID)
  out <- merge(out, temp, all.x=TRUE)
    
  ##1e) OE annotation: collect s_is, tlb
  ##s_i not allowed to be NA - just assume 1.
  ##(Could be an issue if too many OEs were filtered out
  ## => s_i should be 0)
  temp <- x[,list(s_i=s_i[1], tlb=tlb[1]),by="otherEndID"]
  setkey(temp, otherEndID)
  setkey(out, otherEndID)
  out <- merge(out, temp, all.x=TRUE)
  out[is.na(s_i),s_i:=1] ## <-- assumption
  setkey(out, baitID, otherEndID)
  
  ##1f) reconstruct Tmean information
  setkey(x, tlb, tblb)
  setkey(out, tlb, tblb)
  temp <- x[,Tmean[1],by=c("tblb", "tlb")]
  setnames(temp, "V1", "Tmean")
  out[,Tmean:=NULL]
  out <- merge(out, temp, all.x=TRUE)
  
  ##1g) impute some of the missing Tmeans
  ##    missing tlb suggests that the other end was rarely observed, so assume lowest
  ##    Tmean corresponding to the appropriate tblb.
  templowest <- temp[,c(Tmean = min(Tmean)),by="tblb"]
  tempDictionary <- templowest$V1
  names(tempDictionary) <- templowest$tblb
  out[is.na(tlb) & !is.na(tblb), Tmean := tempDictionary[tblb]]
  
  ##2) Reconstruct the distance function from distbin info
  distFunParams <- estimateDistFun(x)
  
  ##3) Reconstruct Bmean information (But!! When s_j = NA,
  ##   ensure that Bmean comes out as NA too.)
  out <- Chicago:::.estimateBMean(out, distFunParams)
  out[is.na(s_j),Bmean:=NA]
  setkey(out, baitID, otherEndID)

  outData[[i]] <- out
}

saveRDS(outData, file = file.path(filePath, "03interactionParameters.Rds"))
saveRDS(dispersions, file = file.path(filePath, "03dispersionParameters.Rds"))