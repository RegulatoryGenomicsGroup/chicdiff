## Some useful functions as I come up with them (possibly all total shit):
## + now me trying to separate the pipeline into functions


#------------------------------------------------------------------------------#

readAndFilterPeakMatrix <- function(peakMatrix, targetColumns = c(), score = 5){
## Essentially reads in the peak matrix, taking the target columns if specified (else just takes everything) and filters for rows where at least
##one score is > 5 (or whatever is specified) and filters out trans interactions. 
  if(length(targetColumns) != 0){
    if(length(targetColumns) == length(targetRDSorRDAFiles)){
      suppressWarnings(x <- fread(peakMatrix)) #Suppressing an irrelevant warning - but probably a good idea to not do this just in case
      sel <- which(colnames(x) %in% targetColumns)
      x <- x[,c(1:11, sel), with=FALSE]
      sel <- rep(FALSE, nrow(x))
      for(cl in targetColumns)
      {
        sel <- sel | x[,get(cl) > score] ##Get any rows where at least one score > 5.
      }
      x <- x[sel,]
    } else{
      message("The number of target columns needs to be equal to the number of provided data files")
      stop
    }
  } else if(length(targetColumns) == 0){
    suppressWarnings(x <- fread(peakMatrix)) #Suppressing an annoying warning - but probably a good idea to not do this
    sel <- rep(FALSE, nrow(x))
    for(cl in c(12:length(x)))
    {
      sel <- sel | x[,get(colnames(x)[cl]) > score] #This might be the stupidest way I could have fixed a problem 
    }
    x <- x[sel,] ##Get any rows where at least one score > 5.
  }
  x <- x[!is.na(dist),] ## <- FILTER OUT TRANS INTERACTIONS - Assumed in (*)
}

#-------------------------------Basically getRegionUniverse----------------------------------#

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

getRegionUniverse <- function(RUmergeMode = "window", RUexpand = 5L, rmapfile = rmapfile, peakMatrix = peakMatrix, 
                              targetColumns = targetColumns, score = 5, saveRDS = FALSE, suffix = ""){
 
  
  
  x <- readAndFilterPeakMatrix(peakMatrix = peakMatrix, targetColumns = targetColumns, score = score)
  
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
  if(saveRDS == TRUE){
    saveRDS(RU.DT, paste0("./RegionUniverse", suffix, ".Rds"))
  }
  RU.DT
}

#-------------------------------Basically getControlRegionUniverse----------------------------------#

##function to sample appropriately from baitIDs on a given chromosome
getRandomBaitID <- function(ch, dists, bmap, rmap)
{
  vIDs <- bmap[chr == ch,ID]
  maxID <- rmap[chr == ch, max(ID)]
  minID <- rmap[chr == ch, min(ID)]
  unlist(lapply(dists,
                function(d)
                {
                  temp <- (vIDs <= maxID - d) | (vIDs >= minID + d)
                  sample(vIDs[temp], size=1)
                }
  ))
}

#------------------------------------------------------------------------------#

##Function to perform shuffle
shuffle <- function(x, bmap, rmap)
{
  ##chose new baits
  ##keep the chromosomes the same, randomly select new baits
  
  ##whether we will go right or left (can be overridden if it makes placement impossible)
  x$goRight <- rbinom(nrow(x),size=1,prob=0.5) == 1
  
  ##pick a bait that will give us a valid pair
  for(ch in unique(x[,CHROM]))
  {
    x[CHROM == ch, baitIDshuff := getRandomBaitID(ch, fragdist + WIDTH, bmap, rmap)]
  }
  
  ##pick other ends the correct distance away
  x[, STARTshuff := ifelse(goRight, baitIDshuff + fragdist, baitIDshuff - fragdist - WIDTH)]
  x[, STOPshuff := STARTshuff + WIDTH]
  
  ##verify other ends are on correct chromosome...
  temp2 <- rmap[,c("chr", "ID"), with=FALSE]
  x <- merge(x, temp2, all.x=TRUE, by.x="STARTshuff", by.y="ID")
  x <- merge(x, temp2, all.x=TRUE, by.x="STOPshuff", by.y="ID")
  x[, notOk := is.na(chr.x) | is.na(chr.y) | chr.x != CHROM | chr.y != CHROM]
  ##...and that OEIDs are positive!
  x[STOPshuff < 1 | STARTshuff < 1, notOk := TRUE]
  x$chr.x <- NULL
  x$chr.y <- NULL
  
  ##Fix the broken other ends
  x[(notOk), goRight:=!goRight]
  
  ##pick other ends the correct distance away
  x[, STARTshuff := ifelse(goRight, baitIDshuff + fragdist, baitIDshuff - fragdist - WIDTH)]
  x[, STOPshuff := STARTshuff + WIDTH]
  
  ##verify other ends are on correct chromosome...
  temp2 <- rmap[,c("chr", "ID"), with=FALSE]
  x <- merge(x, temp2, all.x=TRUE, by.x="STARTshuff", by.y="ID")
  x <- merge(x, temp2, all.x=TRUE, by.x="STOPshuff", by.y="ID")
  x[(notOk), notOk := is.na(chr.x) | is.na(chr.y) | chr.x != CHROM | chr.y != CHROM]
  ##...and that OEIDs are positive!
  x[STOPshuff < 1 | STARTshuff < 1, notOk := TRUE]
  x$chr.x <- NULL
  x$chr.y <- NULL
  
  if(any(x$notOk)) {stop("Random selection failed")}
  
  x
}

#------------------------------------------------------------------------------#

## RU is the output from getRegionUniverse()

getControlRegionUniverse <- function(RU, rmapfile = rmapfile, baitmapfile = baitmapfile, saveRDS = FALSE, 
                                     suffix = ""){
  ##Return to short format
  RUshort <- RU[,list(START=min(otherEndID), STOP=max(otherEndID)),by=c("baitID", "regionID")]
  
  ##Scramble enhancer locations
  ##get chromosome lengths
  rmap <- Chicago:::.readRmap(list(rmapfile=rmapfile))
  bmap <- Chicago:::.readBaitmap(list(baitmapfile=baitmapfile))
  colnames(rmap) <- c("chr","start","stop","ID")
  colnames(bmap) <- c("chr","start","stop","ID","gene")
  chrlengths <- rmap[,list(maxwidth=max(ID)),by="chr"]
  chrlengths <- chrlengths[!(chr %in% c("MT")),]
  #chrlengths <- chrlengths[!(chr %in% c("MT", "Y")),]
  RUshort <- merge(RUshort,
                   rmap[,c("chr","ID"),with=FALSE],
                   by.x="baitID",
                   by.y="ID")
  setnames(RUshort, "chr","CHROM")
  
  RUshort[,WIDTH:=STOP-START]
  
  ##reshuffle the enhancers
  #RUshort$CHROM <- substr(RUshort$CHROM, 4, nchar(RUshort$CHROM))
  RUshort <- merge(RUshort, chrlengths, by.x = "CHROM", by.y="chr")
  RUshort[,fragdist:=pmin(abs(START - baitID), abs(STOP - baitID))]
  
  
  #RUshortBackup <- copy(RUshort)
  RUshort <- shuffle(RUshort, bmap=bmap, rmap=rmap)
  
  ##Return to long format
  RUshuffled <- RUshort[,c("baitIDshuff", "STARTshuff", "STOPshuff"),with=FALSE]
  RUshuffled$regionID <- 1:nrow(RUshuffled)
  names(RUshuffled) <- c("baitID","start","end","regionID")
  RUcontrol <- RUshuffled[,list(otherEndID = start:end),by=c("baitID", "regionID")]
  if(saveRDS == TRUE){
    saveRDS(RUcontrol, paste0("./ControlRegionUniverse", suffix, ".Rds"))
  }
  RUcontrol
}

#-------------------------------Function combining 02, 03 and 04----------------------------------#
#Function should collect the count data from the chinputs, get the dispersions and 
#interaction parameters and then merge to give full data for each region

readRDSorRDA <- function(file){
  
  suppressWarnings(tmp <- try(load(file), silent = TRUE))
  if("try-error" %in% class(tmp)){
    #message("Not an RDA. Checking whether RDS.")
    suppressWarnings(tmp <- try(output <- readRDS(file), silent = TRUE))
    if("try-error" %in% class(tmp)){
      message("Cannot read file. Must be an RDA or RDA file.")
    } else{
      #message("Detected RDS")
      output
    }
  } else{
    #message("Detected RDA")
    eval(parse(text = get(ls()[ls() != "file"])))
  }
}

#------------------------------------------------------------------------------#

sreadRDSorRDA <- function(file){
  if(tolower(tools::file_ext(file)) == "rda"){
    load(file)
    x
    # eval(parse(text = get(ls()[ls() != "file"]))) if 'x' isn't the object produced by the RDA
  } else if(tolower(tools::file_ext(file)) == "rds"){
    output <- readRDS(file)
    output
  } else{
    message("Cannot determine whether input is RDA or RDA. Please provide file extension.")
  }
}

#------------------------------------------------------------------------------#

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
#I don't think this function is actually ever used but I have kept it for the moment 

#------------------------------------------------------------------------------#

#This will probably need another name as a very similar function of the same name appears in CHiCAGO
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

#---------------------------First the version which acts on (non)control individually-------------#

getFullRegionData <- function(RU, targetChFiles = targetChFiles, targetRDSorRDAFiles = targetRDSorRDAFiles, 
                              targetRDSorRDAs = targetRDSorRDAs, rmapfile = rmapfile, saveRDS = FALSE, suffix = ""){
  
  #Moved the 02getCounts.R script part AFTER the 03getInteractionParameters.R as I think it made more sense with the if statement
  #regarding whether or not the count data was coming from chinputs or the RDSs (and avoids having to read the RDSs in twice)
  
  
  
  ## 03getInteractionParameters
  
  RU_IntParams <- copy(RU)
  RU_IntParams[,regionID:=NULL]
  RU_IntParams <- unique(RU_IntParams)
  
  ##prep output
  outData <- vector("list", length(targetRDSorRDAFiles))
  dispersions <- numeric(length(targetRDSorRDAFiles))
  tempForCounts <- vector("list", length(targetRDSorRDAFiles))
  
  ##1) Read in each file...
  for(i in 1:length(targetRDSorRDAFiles))
  {
    message("File ",i, " of ", length(targetRDSorRDAFiles))
    file <- targetRDSorRDAFiles[i]
    #load(file) ##this produces the object 'x' --old 
    x <- readRDSorRDA(file)
    ##1a) collect the dispersion
    #dispersions[i] <- attributes(x)$dispersion --old
    
    # The following if statement might break if the class of the data being read being 'chicagoData'...
    #...is not a necessary and sufficient condition for it to be in the form assumed below. I also have...
    #...used setDT for non-chicagoData type because I can and the original did - not sure this is the correct decision to make though.
    
    if("chicagoData" %in% class(x)){
      dispersions[i] <- x@params$dispersion
      x <- as.data.table(x@x)
    } else{
      dispersions[i] <- attributes(x)$dispersion
      setDT(x)
    }
    
    setkey(x, baitID, otherEndID)
    setkey(RU_IntParams, baitID, otherEndID)
    out <- x[RU_IntParams, c("baitID", "otherEndID","distSign","Bmean", "Tmean", "score"), with=FALSE]
    ##(note that score=NA occurs when N=0 - hence replace score=0)
    
    ##1c) recalculate distSign (need to do this for control regions)
    if(any(is.na(out$distSign)))
    {
      rmap <- Chicago:::.readRmap(list(rmapfile=rmapfile))
      colnames(rmap) <- c("chr", "start", "end", "ID")
      setkey(rmap, "ID")
      temp <- merge(RU_IntParams, rmap, all.x=TRUE, by.x="baitID", by.y="ID")
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
    
    #This just keeps all of the loaded RDAs if we still need them for the count data because we don't have chinputs
    if(length(targetChFiles) == 0){  #my lame place holder conditional for checking whether chinputs have been provided. Easily changed. 
      x <- x[, c("baitID", "otherEndID", "N"), with = FALSE]
      columnNames <- paste0("N.", names(targetRDSorRDAFiles))
      setnames(x, old="N", new=columnNames[i])
      setkey(x, baitID, otherEndID)
      tempForCounts[[i]] <- x
      rm(x)
    }
  }
  saveRDS(outData, paste0("03interactionsParameters", suffix, ".Rds"))   ##REMEMBER TO REMOVE THESE TWO SAVES
  saveRDS(dispersions, paste0("03dispersionParamters", suffix, ".Rds"))
  baits <- sort(unique(RU$baitID)) ##baits to get information from
  
  if(length(targetChFiles) == 0){
    gc()
    mergedFiles <- Reduce(merge, tempForCounts)
    #Including the potentially pointless extra steps below so as to introduce as little deviation from
    #Jonathan's original code for the time being. Probably best way to understand why I have done what
    #I have is to look at the original 02getCounts.R script and compare
    countData <- vector("list", length(targetRDSorRDAFiles))
    for(i in seq_along(countData)){
      countData[[i]] <- mergedFiles[, c("baitID", "otherEndID", columnNames[i]), with = FALSE]
      setkey(countData[[i]], baitID)
      countData[[i]] <- countData[[i]][J(baits),] #issue that rds may not contain all baitid otherendid pairs which are in RU
      rows2remove <- which(is.na(countData[[i]]$otherEndID)) #I'm assuming that otherEndIDs should never be NA here
      countData[[i]] <- countData[[i]][-rows2remove, ] #Annoyingly, I can't find a way to remove rows by reference so I have this clumsy, memory hungry fix instead
      setkey(countData[[i]], baitID)
    }
    
    x <- RU
    setkey(x, baitID, otherEndID)
    
    for(i in seq_along(countData))
    {
      #columnName <- paste0("N.", names(targetChFiles)[i])
      temp <- countData[[i]] #[,c("baitID","otherEndID","N"), with=FALSE]
      setnames(temp, old = columnNames[i], new = "N")
      setkey(temp, baitID, otherEndID)
      x <- merge(x, temp, all.x=TRUE)
      x[is.na(N), N:=0]
      setnames(x, old="N", new=columnNames[i])
    }
  }
  
  #Recall that outData and dispersions are our useful outputs from this part above 
  
  ## 02getCounts.R
  
  ##Given a region universe, collect the number of reads in each region from the chinputs.
  
  ##Inputs and parameters:
  
  ##requires targetChs
  
  ##2) Collect all of the read count information ----------------
  if(length(targetChFiles) != 0){
    
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
  CountOut <- x
  saveRDS(CountOut, paste0("02CountOut", suffix, ".Rds")) ##REMEMBER TO REMOVE
  #04mergeData
  
  #x <- readRDS(file.path(filePath, "02CountOut.Rds")) #CountOut
  #params <- readRDS(file.path(filePath, "03interactionParameters.Rds")) #outData
  
  setkey(CountOut, baitID, otherEndID)
  
  for(i in 1:length(targetRDSorRDAFiles))
  {
    s <- names(targetRDSorRDAFiles[i])
    temp <- outData[[i]][,c("baitID", "otherEndID", "s_j", "Bmean", "Tmean", "score"), with=FALSE]
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
    
    CountOut <- merge(CountOut, temp, all.x=TRUE)
    #x <- merge(x, temp, all.x=TRUE, allow.cartesian=TRUE)
  }
  
  ##rework table to have columns sample, N, FullMean, etc
  ## (rather than N.[sample], FullMean.[sample])
  myCols <- c("N", "s_j", "Bmean", "Tmean", "score", "FullMean")
  getCols <- function(col) {grep(paste0("^", col), colnames(CountOut), perl=TRUE)}
  sels <- lapply(myCols, getCols)
  
  recast = melt(CountOut, measure = sels, value.name = myCols,
                variable.name="sample")
  recast$sample <- names(targetRDSorRDAFiles)[recast$sample]
  
  recast$condition <- rep(
    rep(names(targetRDSorRDAs), sapply(targetRDSorRDAs, length))
    , each=nrow(CountOut))
  
  
  setkey(recast, regionID)
  
  if(saveRDS == TRUE){
    saveRDS(recast, paste0("./FullRegionData", suffix, ".Rds"))
  }
  recast
}

#------------Below outdated as relies on chinputs still---------------Now the version which acts on both (in a rather inelegant way)-------------# 

getFullRegionData2 <- function(RU, RUcontrol, targetChFiles = targetChFiles, targetRDSorRDAFiles = targetRDSorRDAFiles, 
                              targetRDSorRDAs = targetRDSorRDAs, rmapfile = rmapfile, saveRDS = FALSE, suffix = ""){
  
  baits <- sort(unique(RU$baitID)) 
  baitsControl <- sort(unique(RUcontrol$baitID))
  
  ##2) Collect all of the read count information ----------------
  
  countData <- vector("list", length(targetChFiles))
  countDataControl <- vector("list", length(targetChFiles))
  gc()
  
  for(i in 1:length(countData))
  {
    #message(i, "\n")
    x <- fread(targetChFiles[i])
    setkey(x, baitID)
    #message(key(x), " ")
    y <- x[J(baitsControl),]
    x <- x[J(baits),]
    setkey(x, baitID) ##Required for countData[4:7], where key is dropped - not sure why
    setkey(y, baitID)
    #message(key(x), " \n")
    countData[[i]] <- x
    countDataControl[[i]] <- y
    rm(x)
    rm(y)
  }
  #Added some new lines above to create countDataControl with (hopefully) minimal redundancy.
  
  ##2a) For each region, count the number of reads in each replicate -----------------
  ##(FIXME: could be optimized a lot)
  
  CountOut <- RU 
  setkey(CountOut, baitID, otherEndID)
  for(i in 1:length(countData))
  {
    columnName <- paste0("N.", names(targetChFiles)[i])
    temp <- countData[[i]][,c("baitID","otherEndID","N"), with=FALSE]
    setkey(temp, baitID, otherEndID)
    CountOut <- merge(CountOut, temp, all.x=TRUE)
    
    CountOut[is.na(N), N:=0]
    setnames(CountOut, old="N", new=columnName)
  }
  
  #Repeating the above code for the control
  CountOutControl <- RUcontrol 
  setkey(CountOutControl, baitID, otherEndID)
  for(i in 1:length(countDataControl))
  {
    columnName <- paste0("N.", names(targetChFiles)[i])
    temp <- countData[[i]][,c("baitID","otherEndID","N"), with=FALSE]
    setkey(temp, baitID, otherEndID)
    CountOutControl <- merge(CountOutControl, temp, all.x=TRUE)
    
    CountOutControl[is.na(N), N:=0]
    setnames(CountOutControl, old="N", new=columnName)
  }
  
  ##3) Collect distance -----------------
  
  ##3a) Distance ------------
  rmap <- Chicago:::.readRmap(list(rmapfile=rmapfile))
  colnames(rmap) <- c("chr","start","end","fragID")
  setkey(rmap, fragID)
  rmap[,midpoint:=round(0.5*(start+end))]
  rmap[,start:=NULL]
  rmap[,end:=NULL]
  
  CountOut <- merge(CountOut, rmap, by.x="otherEndID", by.y="fragID")
  CountOut <- merge(CountOut, rmap, by.x="baitID", by.y="fragID")
  setkey(CountOut, baitID, otherEndID)
  
  CountOut[, distSign :=
             ifelse(chr.x == chr.y,
                    midpoint.x - midpoint.y,
                    NA)]
  CountOut[, c("chr.x","midpoint.x","chr.y","midpoint.y"):=NULL]
  
  #Just repeating the above code for the control
  CountOutControl <- merge(CountOutControl, rmap, by.x="otherEndID", by.y="fragID")  
  CountOutControl <- merge(CountOutControl, rmap, by.x="baitID", by.y="fragID")
  setkey(CountOutControl, baitID, otherEndID)
  
  CountOutControl[, distSign :=
                    ifelse(chr.x == chr.y,
                           midpoint.x - midpoint.y,
                           NA)]
  CountOutControl[, c("chr.x","midpoint.x","chr.y","midpoint.y"):=NULL]
  
  #CountOut used to be called 'x' and CountOutControl was called 'y'
  #saveRDS(CountOut, file = file.path(filePath, "02CountOut.Rds"))
  #saveRDS(CountOutControl, file = file.path(filePathControl, "02CountOut.Rds"))
  
  ## 03getInteractionParameters.R         13.15   29/08/17
  
  ##Given a region universe, collect Bmean + Tmean.
  ##Note: Some baits missing (filtered out from <250 reads)
  
  ##Inputs and parameters:
  
  ##requires targetRDAs
  
  #RUfile <- file.path(filePath, "01regionUniverse.Rds")   #Not necessary as we still have unaltered RU.DT 
  
  ##A couple of functions for reading in either RDS or RDA files. The first one basically just checks whether...
  #...trying to load() or readRDS() the file produces an error and then does whichever does not. The second...
  #...one just checks the file extension and then uses the correct function. 
  
  #RU <- readRDS(RUfile) RU.DT still exists - but would have to be altered in the next two lines, so:
  RU_IntParams <- copy(RU.DT) 
  RU_IntParams[,regionID:=NULL]
  RU_IntParams <- unique(RU_IntParams)
  
  RUcontrol_IntParams <- copy(RUcontrol)
  RUcontrol_IntParams[,regionID:=NULL]
  RUcontrol_IntParams <- unique(RUcontrol_IntParams)
  #Past this point, any mention of RU or RUcontrol in the 03 section of this function is a mistake. Should be RU_IntParams etc.
  
  ##prep output
  outData <- vector("list", length(targetRDSorRDAFiles))
  outDataControl <- vector("list", length(targetRDSorRDAFiles)) #New to collect control data
  dispersions <- numeric(length(targetRDSorRDAFiles)) #only need one of these
  
  ##1) Read in each file...
  for(i in 1:length(targetRDSorRDAFiles))
  {
    message("File ",i, " of ", length(targetRDSorRDAFiles))
    file <- targetRDSorRDAFiles[i]
    #load(file) ##this produces the object 'x' --old 
    x <- readRDSorRDA(file)
    ##1a) collect the dispersion
    #dispersions[i] <- attributes(x)$dispersion --old
    
    # The following if statement might break if the class of the data being read being 'chicagoData'...
    #...is not a necessary and sufficient condition for it to be in the form assumed below. I also have...
    #...used setDT for non-chicagoData type because I can. 
    
    if("chicagoData" %in% class(x)){
      dispersions[i] <- x@params$dispersion
      x <- as.data.table(x@x)
    } else{
      dispersions[i] <- attributes(x)$dispersion
      setDT(x)
    }
    
    ##1b) collect the Bmean, Tmean information that is present
    #setDT(x) --old (now included in the above if statement)
    setkey(x, baitID, otherEndID)
    setkey(RU_IntParams, baitID, otherEndID)
    setkey(RUcontrol_IntParams, baitID, otherEndID)
    out <- x[RU_IntParams, c("baitID", "otherEndID","distSign","Bmean", "Tmean", "score"), with=FALSE]
    outControl <- x[RUcontrol_IntParams, c("baitID", "otherEndID","distSign","Bmean", "Tmean", "score"), with=FALSE]
    #New lines above to produce outControl in the same way as out sis produced
    
    ##(note that score=NA occurs when N=0 - hence replace score=0)
    
    ##1c) recalculate distSign (need to do this for control regions)
    if(any(is.na(out$distSign)))
    {
      rmap <- Chicago:::.readRmap(list(rmapfile=rmapfile))
      colnames(rmap) <- c("chr", "start", "end", "ID")
      setkey(rmap, "ID")
      temp <- merge(RU_IntParams, rmap, all.x=TRUE, by.x="baitID", by.y="ID")
      temp <- merge(temp, rmap, all.x=TRUE, by.x="otherEndID", by.y="ID")
      setkey(temp, baitID, otherEndID)
      temp[,distSign:=round(((start.y + end.y)-(start.x + end.x))/2)]
      if(any(abs(temp$distSign - out$distSign) > 1, na.rm=TRUE))
      {
        stop("Error calculating distances.")
      }
      out$distSign <- temp$distSign
    }
    #Code below just a duplicate of that above but should be adapted for the control outputs. The >1 is just so...
    #...discrepencies due to rounding aren't a problem.
    
    if(any(is.na(outControl$distSign)))
    {
      rmap <- Chicago:::.readRmap(list(rmapfile=rmapfile))
      colnames(rmap) <- c("chr", "start", "end", "ID")
      setkey(rmap, "ID")
      temp <- merge(RUcontrol_IntParams, rmap, all.x=TRUE, by.x="baitID", by.y="ID")
      temp <- merge(temp, rmap, all.x=TRUE, by.x="otherEndID", by.y="ID")
      setkey(temp, baitID, otherEndID)
      temp[,distSign:=round(((start.y + end.y)-(start.x + end.x))/2)]
      if(any(abs(temp$distSign - outControl$distSign) > 1, na.rm=TRUE))
      {
        stop("Error calculating distances.")
      }
      outControl$distSign <- temp$distSign
    }
    
    ##1d) bait annotation: collect s_js, tblb
    ##Note that s_j can be NA (for baits that were filtered out!)
    temp <- x[,list(s_j=s_j[1], tblb=tblb[1]),by="baitID"]
    setkey(temp, baitID)
    setkey(out, baitID)
    setkey(outControl, baitID)
    out <- merge(out, temp, all.x=TRUE)
    outControl <- merge(outControl, temp, all.x = TRUE)
    #Added the new lines for outControl
    
    ##1e) OE annotation: collect s_is, tlb
    ##s_i not allowed to be NA - just assume 1.
    ##(Could be an issue if too many OEs were filtered out
    ## => s_i should be 0)
    temp <- x[,list(s_i=s_i[1], tlb=tlb[1]),by="otherEndID"]
    setkey(temp, otherEndID)
    setkey(out, otherEndID)
    setkey(outControl, otherEndID)
    
    out <- merge(out, temp, all.x=TRUE)
    outControl <- merge(outControl, temp, all.x = TRUE)
    
    out[is.na(s_i),s_i:=1] ## <-- assumption
    setkey(out, baitID, otherEndID)
    
    outControl[is.na(s_i),s_i:=1] ## <-- assumption
    setkey(outControl, baitID, otherEndID)
    
    ##1f) reconstruct Tmean information
    setkey(x, tlb, tblb)
    setkey(out, tlb, tblb)
    setkey(outControl, tlb, tblb)
    temp <- x[,Tmean[1],by=c("tblb", "tlb")]
    setnames(temp, "V1", "Tmean")
    out[,Tmean:=NULL]
    outControl[,Tmean:=NULL]
    out <- merge(out, temp, all.x=TRUE)
    outControl <- merge(outControl, temp, all.x=TRUE)
    
    ##1g) impute some of the missing Tmeans
    ##    missing tlb suggests that the other end was rarely observed, so assume lowest
    ##    Tmean corresponding to the appropriate tblb.
    templowest <- temp[,c(Tmean = min(Tmean)),by="tblb"]
    tempDictionary <- templowest$V1
    names(tempDictionary) <- templowest$tblb
    out[is.na(tlb) & !is.na(tblb), Tmean := tempDictionary[tblb]]
    outControl[is.na(tlb) & !is.na(tblb), Tmean := tempDictionary[tblb]]
    
    ##2) Reconstruct the distance function from distbin info
    distFunParams <- estimateDistFun(x)
    
    ##3) Reconstruct Bmean information (But!! When s_j = NA,
    ##   ensure that Bmean comes out as NA too.)
    out <- Chicago:::.estimateBMean(out, distFunParams)
    out[is.na(s_j),Bmean:=NA]
    setkey(out, baitID, otherEndID)
    
    outControl <- Chicago:::.estimateBMean(outControl, distFunParams)
    outControl[is.na(s_j),Bmean:=NA]
    setkey(outControl, baitID, otherEndID)
    
    outData[[i]] <- out
    outDataControl[[i]] <- outControl
  }
  
  #saveRDS(outData, file = file.path(filePath, "03interactionParameters.Rds"))
  #saveRDS(outDataControl, file = file.path(filePathControl, "03interactionParameters.Rds"))
  
  #saveRDS(dispersions, file = file.path(filePath, "03dispersionParameters.Rds")) #This only needs to be produced once and at the moment I have it going into the Resources directory (and not into the control directory within)
  
  # 04mergeData.R
  
  ##Synthesise the count + interaction parameter data together
  
  ##inputs:
  
  ##requires targetRDAs 
  #x <- readRDS(file.path(filePath, "02CountOut.Rds"))    ##CountOut and CountOutControl
  #params <- readRDS(file.path(filePath, "03interactionParameters.Rds"))  ##outData and outDataControl
  
  #Don't need these lines above. Have two options here. Continue with CountOut (instead of x) and outData 
  #(instead of params) or define x and params as copies of the aforementioned objects.  Will do the first for now
  #but keep in mind this should probably be changed.
  
  setkey(CountOut, baitID, otherEndID)
  setkey(CountOutControl, baitID, otherEndID)
  
  for(i in 1:length(targetRDSorRDAFiles))
  {
    s <- names(targetRDSorRDAFiles[i])
    temp <- outData[[i]][,c("baitID", "otherEndID", "s_j", "Bmean", "Tmean", "score"), with=FALSE]
    tempControl <- outDataControl[[i]][,c("baitID", "otherEndID", "s_j", "Bmean", "Tmean", "score"), with=FALSE]
    temp <- unique(temp)
    tempControl <- unique(tempControl)
    temp[,FullMean:=Bmean+Tmean]
    tempControl[,FullMean:=Bmean+Tmean]
    #temp[,Bmean:=NULL]
    #temp[,Tmean:=NULL]
    
    for(nm in colnames(temp))
      if(!nm %in% c("baitID", "otherEndID"))
      {
        setnames(temp, nm, paste0(nm, ".", s))
      }
    #Just a repeat of the above but for control unfortunately. Also hoping this simple copy doesn't break anything.
    for(nm in colnames(tempControl))
      if(!nm %in% c("baitID", "otherEndID"))
      {
        setnames(tempControl, nm, paste0(nm, ".", s))
      }
    
    setkey(temp, baitID, otherEndID)
    setkey(tempControl, baitID, otherEndID)
    
    CountOut <- merge(CountOut, temp, all.x=TRUE)
    CountOutControl <- merge(CountOutControl, tempControl, all.x=TRUE)
    #x <- merge(x, temp, all.x=TRUE, allow.cartesian=TRUE)
  }
  
  ##rework table to have columns sample, N, FullMean, etc
  ## (rather than N.[sample], FullMean.[sample])
  myCols <- c("N", "s_j", "Bmean", "Tmean", "score", "FullMean")
  getCols <- function(col) {grep(paste0("^", col), colnames(CountOut), perl=TRUE)} #colnames(CountOut) and colnames(CountOutControl) will always be identical - so I have just chosen one
  sels <- lapply(myCols, getCols)
  
  recast = melt(CountOut, measure = sels, value.name = myCols,
                variable.name="sample")
  recast$sample <- names(targetRDSorRDAFiles)[recast$sample]
  
  recast$condition <- rep(
    rep(names(targetRDSorRDAs), sapply(targetRDSorRDAs, length))
    , each=nrow(CountOut))
  
  # Again - just a copy and paste for the control here 
  recastControl = melt(CountOutControl, measure = sels, value.name = myCols,
                       variable.name="sample")
  recastControl$sample <- names(targetRDSorRDAFiles)[recastControl$sample]
  
  recastControl$condition <- rep(
    rep(names(targetRDSorRDAs), sapply(targetRDSorRDAs, length))
    , each=nrow(CountOutControl))
  
  setkey(recast, regionID)
  setkey(recastControl, regionID)
  
  if(saveRDS == TRUE){
    saveRDS(recast, paste0("./FullRegionData", suffix, ".Rds"))
    saveRDS(recast, paste0("./FullControlRegionData", suffix, ".Rds"))
  }
  #### Don't really know what would be the best output for this - so as a place holder I am doing the following:
  FullRegionData <- list(recast, recastControl)
}

#-------------------------------Basically DESeq2 wrapper function - for individual input----------------------------------#

logit <- function(p) {p/1-p}
expit <- function(x) {1/(1+exp(-x))}
geoMean <- function(x, na.rm=FALSE) {
  if(na.rm){
    exp(mean(log(x[!is.na(x)])))
  }else{
    exp(mean(log(x)))
  }
}

#------------------------------------------------------------------------------#

DESeq2Wrap <- function(RU, FullRegionData, rmapfile = rmapfile, saveRDS = FALSE, suffix = ""){
  ##By default, DESeq2 uses sample-specific scaling factors.
  ##Idea: Allow bait-specific scaling factors. Show this improves model fit.
  ##      May lead to a drop in the dispersion.
  
  library(data.table) ##Just keeping these in for the moment to remind me to look at this
  library(DESeq2)
  
  ##Input data:
  
  fragData <- copy(FullRegionData) ## As otherwise it is altered by reference below and the IHW part breaks
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
  
  #RU <- readRDS(file.path(filePath, "01regionUniverse.Rds"))
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
  
  out.nullModel <- cbind(
    as.data.table(as.data.frame(results(dds.nullModel))),
    annoData
  )
  
  if(saveRDS == TRUE){
    saveRDS(out, paste0("./05outPilot", suffix, ".Rds"))
    saveRDS(out.nullModel, paste0("./05outStandardDESeq", suffix, ".Rds"))
  }
  out ##Currently have the norm factor model as the output of this wrapper - but presumably we want the 'naive' DESeq2 output as well 
}

#------------------------------------------------------------------------------#

IHWcorrection <- function(DESeqOut, FullRegionData, DESeqOutControl, FullControlRegionData, saveRDS = TRUE,
                          DiagPlot = TRUE, suffix = ""){
  library(IHW)
  library(data.table)
  library(ggplot2)
  library(hexbin) #Became necessary (for some reason) when saving the graphs
  
  out <- DESeqOut ##DESeqOut and DESeqOutControl
  RU.recast <- FullRegionData ##FullRegionData and FullControlRegionData
  RU.distances <- RU.recast[,list(avDist = mean(distSign)),by="regionID"]
  
  out$avDist <- RU.distances$avDist
  
  ##Comparison against genuinely uniform p-vals
  out$uniform <- runif(nrow(out))
  out$shuff <- sample(out$pvalue)
  
  ##Convert
  as.data.frame(out) #is this line still necessary? I will check and report back.
  setDT(out)
  
  out.control <- DESeqOutControl ##DESeqOut and DESeqOutControl
  RU.recastControl <- FullControlRegionData ##FullRegionData and FullControlRegionData
  RU.distancesControl <- RU.recastControl[,list(avDist = mean(distSign)),by="regionID"]
  
  out.control$avDist <- RU.distancesControl$avDist
  
  ##Comparison against genuinely uniform p-vals
  out.control$uniform <- runif(nrow(out.control))
  out.control$shuff <- sample(out.control$pvalue)
  
  ##Convert
  as.data.frame(out.control) #This line is necessary even if the one for out may not be
  
  ##Train weights on the control sample
  ihwRes <- ihw(pvalue ~ abs(avDist),  data = out.control, alpha = 0.05)
  
  if(DiagPlot == TRUE){
    plot(ihwRes)
    ggsave("IHWweightPlot.png", device = "png", path = "./")
    plot(ihwRes, what = "decisionboundary")
    ggsave("IHWdecisionBoundaryPlot.png", device = "png", path = "./")
  }
 
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
  if(saveRDS == TRUE){
    saveRDS(out, paste0("06Weighted", suffix, ".Rds"))
  }
  out
}