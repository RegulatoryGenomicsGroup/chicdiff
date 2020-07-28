## Settings functions -----------------

defaultChicdiffSettings <- function()
{
  list(
    inputfiles= NA,
    peakfiles= NA,
    chicagoData= NA,
    countData= NA,
    rmapfile= NA,
    targetColumns= NA,
    baitmapfile= NA,
    RUexpand= 5L, 
    score= 5,
    norm="combined", 
    theta = NULL,
    theta_grid = seq(0,1,0.25),
    saveAuxData = FALSE,
    parallel = FALSE,
    device = "png",
    printMemory = FALSE,
    outprefix = ""    
  )
}

### This is where we now set the defaults
### Order of priority:
### settings override settings from settingsFile
### both override chicdiff.settings

setChicdiffExperiment = function(designDir="", chicagoData = NA, countData = NA, peakfiles = NA, 
                                 outprefix = "test", settings=list(), settingsFile=NULL, inputfiles = NA,
                                 chicdiff.settings=defaultChicdiffSettings())
{
  
  if(designDir == "" & identical(settings, list()) & is.null(settingsFile) & identical(chicdiff.settings, defaultChicdiffSettings()))
  {
    stop("Design not specified. Please specify a design directory or design files.")
  }
  chicdiff.settings = .updateChicdiffSettings(designDir, chicagoData, countData, peakfiles, outprefix,
                                              settings, settingsFile, inputfiles, chicdiff.settings)
  
  
  saveRDS(chicdiff.settings, paste0(chicdiff.settings[["outprefix"]], "_settings.Rds"))
  
  chicdiff.settings
  
}

.updateChicdiffSettings = function(designDir, chicagoData, countData, peakfiles, outprefix,
                                   settings, settingsFile, inputfiles, chicdiff.settings, 
                                   updateDesign=FALSE){
  modSettings = vector("list")
  
  if(!is.null(settingsFile)){
    
    message(paste0("Reading custom experimental settings from ", settingsFile, "..."))
    
    # http://stackoverflow.com/questions/6602881/text-file-to-list-in-r
    sf <- scan(settingsFile, what="", sep="\n")
    modSettings <- strsplit(sf, "[[:space:]]+")
    names(modSettings) <- sapply(modSettings, `[[`, 1)
    modSettings <- lapply(modSettings, `[`, -1)
    
    # convert numerically defined settings to numbers
    # do the same for logical settings
    # suppressing "NAs introduced by coercion" warnings
    suppressWarnings({
      modSettings <- lapply(modSettings, function(s){
        num_s = as.numeric(s);
        if(!is.na(num_s)) { return (num_s) };
        bool_s = as.logical(s);
        if(!is.na(bool_s)) { return(bool_s) };
        s 
      })
    })
  }
  
  for (s in names(settings)){
    modSettings[[s]] = settings[[s]]
  }
  
  for (s in names(modSettings)){
    chicdiff.settings[[s]] = modSettings[[s]]
  }
  
  chicdiff.settings[["outprefix"]] = outprefix
  chicdiff.settings[["chicagoData"]] = chicagoData
  chicdiff.settings[["countData"]] = countData
  
  if(any(is.na(chicagoData)) | any(is.na(countData))){
    if(is.na(chicdiff.settings[["inputfiles"]])){
      chicdiff.settings[["inputfiles"]] = inputfiles
    }else{
      if (!file.exists(chicdiff.settings[["inputfiles"]])){
        stop(paste("No config file found at the specified location", chicdiff.settings[["inputfiles"]]))
      }
    }
  } else{
    if(!(all.equal(names(chicagoData), names(countData)))){
      stop("Conditions for the RDS/RDA files and chinputs must be the same")
    }
    temp <- unlist(chicagoData)
    temp_2 <- unlist(countData)
    if(!(length(temp) == length(temp_2))){
      stop("Must provide the same number of RDS/RDA files as chinputs")
    }
  }
  
  if(!is.na(inputfiles)){
    input_lists <- .makeTargetFilesList(inputfiles)
    chicagoData <- input_lists[[1]]
    countData <- input_lists[[2]]
  }
  
  chicdiff.settings[["peakfiles"]] = peakfiles
  
  if(is.na(chicdiff.settings[["peakfiles"]])){
    stop("No peak files provided")
  } else if(!file.exists(chicdiff.settings[["peakfiles"]])){
    stop(paste("No peak files found at the specified location", chicdiff.settings[["peakfiles"]]))
  }
  
  targetColumns <- .getTargetColumns(chicagoData)
  chicdiff.settings[["targetColumns"]] = targetColumns
  
  # peakfile_columns <- colnames(.multimerge(peakfiles))
  
  #if(!all(targetColumns %in% peakfile_columns)){
  #  stop("All specified columns must be present in the peak files")  
  #}
  
  if(is.na(chicdiff.settings[["baitmapfile"]]) | (updateDesign & !is.null(designDir))){
    chicdiff.settings[["baitmapfile"]] = .locateFile("<baitmapfile>.baitmap", designDir, "\\.baitmap$")
  }else{
    if (!file.exists(chicdiff.settings[["baitmapfile"]])){
      stop(paste("No baitmap file found at the specified location", chicdiff.settings[["baitmapfile"]]))
    }
  }
  
  if(is.na(chicdiff.settings[["rmapfile"]]) | (updateDesign & !is.null(designDir))){
    chicdiff.settings[["rmapfile"]] = .locateFile("<rmapfile>.rmap", designDir, "\\.rmap$")
  }else{
    if (!file.exists(chicdiff.settings[["rmapfile"]])){
      stop(paste("No rmap file found at the specified location", chicdiff.settings[["rmapfile"]]))
    }
  }
  
  chicdiff.settings[["norm"]] = tolower(chicdiff.settings[["norm"]])
  if (length(chicdiff.settings[["norm"]])>1){ stop ("Parameter error: Only one normalisation method can be specified at a time") }
  if (!chicdiff.settings[["norm"]] %in% c("standard", "fullmean", "combined")){
    stop ("Parameter error: normalisation method should be one of 'standard', 'fullmean', 'combined'")
  }  
  
  message("Checking the design files...")
  
  baitmapfile = Chicago:::.readBaitmap(chicdiff.settings)
  
  rmapfile = Chicago:::.readRmap(chicdiff.settings)
  if (ncol(rmapfile)<4){
    stop("There are fewer columns in the rmap file than expected. This file should have 4 columns, listing the genomic coordinates and IDs for each restriction fragment.")
  }
  if (ncol(rmapfile)>4){
    stop("There are more columns in the rmap file than expected. This file should have 4 columns, listing the genomic coordinates and IDs for each restriction fragment. Check that rmapfile and baitmap files aren't swapped round\n")
  }
  
  if (any(duplicated(rmapfile[[4]]))){
    stop(paste("Duplicated fragment IDs found in rmapfile (listed below):\n", 
               paste(rmapfile[[4]][duplicated(rmapfile[[4]])], collapse=",")))
  }
  
  chicdiff.settings
}

# ----------

.getTargetColumns <- function(chicagoData){
  
  conditions <- names(chicagoData)
  targetColumns <- vector("list", length(conditions))
  for(i in seq_along(conditions)){
    temp <- conditions[[i]]
    targetColumns[[i]] <- names(chicagoData[[temp]])
  }
  targetColumns <- unlist(targetColumns)
  if(is.null(targetColumns)){
    targetColumns <- conditions
  } 
  okay = (length(targetColumns) == length(conditions) | length(targetColumns) == length(unlist(chicagoData)))
  if(!okay){ stop("Peak file columns incorrectly specified") }
  
  targetColumns
}

.strsplit_unlist <- function(x){
  temp <- strsplit(x, split= ",")
  unlist(temp)
}

.makeTargetFilesList <- function(config){
  temp <- fread(config, header = FALSE)
  temp_1 <- temp[, paste(V2, collapse = ","), by = V1]
  temp_2 <- temp[, paste(V3, collapse = ","), by = V1]
  
  names(temp_1) <- c("condition", "files")
  names(temp_2) <- c("condition", "files")
  
  temp_1_2 <- lapply(temp_1[,files], .strsplit_unlist)
  temp_2 <- lapply(temp_2[,files], .strsplit_unlist)
  
  names(temp_1_2) <- temp_1[, condition]
  names(temp_2) <- temp_1[, condition]
  return(list(temp_1_2, temp_2))
}

#-------------------------------Auxilliary functions-----------------------------------------------#

.multimerge <- function(peakFiles, targetColumns)
{
  peakDTs <- lapply(peakFiles, fread)
  peakMatrix <- Reduce(function(x,y)
    merge(x, y, 
          by=c("baitChr","baitStart","baitEnd","baitID",         
               "baitName","oeChr","oeStart","oeEnd",          
               "oeID","oeName","dist"), all=T), peakDTs)
  
  peakMatrix
}

readAndFilterPeakMatrix <- function(peakFiles, targetColumns, chicagoData, conditions, score, outprefix=""){
  ## Essentially reads in the peak matrix, taking the target columns if specified (else just takes everything) and filters for rows where at least
  ##one score is > 5 (or whatever is specified) and filters out trans interactions. 
  
  if(length(peakFiles) > 1){
    x <- .multimerge(peakFiles, targetColumns)
  } else {
    x <- fread(peakFiles)
  }
  
  peakfile_columns <- colnames(x)
  
  if(!all(targetColumns %in% peakfile_columns)){
    stop("All specified targetColumns must be present in the peak file(s)")  
  }
  
  all_baits <- unique(x$baitID)
  sel <- which(colnames(x) %in% targetColumns)
  x <- x[,c(1:11, sel), with=FALSE]
  sel <- rep(FALSE, nrow(x))
  for(cl in targetColumns)
  {
    sel <- sel | x[,get(cl) > score & !is.na(get(cl))] ##Get any rows where at least one score > 5.
  }
  x <- x[sel,]
  
  if(length(targetColumns) > length(conditions)){
    
    sel2 <- rep(TRUE, nrow(x))
    for(i in seq_along(conditions)){
      temp <- conditions[[i]]
      cols <- colnames(x)[colnames(x) %in% names(chicagoData[[temp]])]
      
      sel2 <- sel2 & rowSums(x[,!is.na(.SD[,cols, with = FALSE])]) >= 2 ##Remove any rows where there isn't at least 2 replicates for each condition with non-NA values
    } 
    
    x <- x[sel2]
  }
  
  x <- x[!is.na(dist),] ## <- FILTER OUT TRANS INTERACTIONS - Assumed in (*)
  x <- x[!(oeID == (baitID + 1) | oeID == (baitID - 1)),] ## FILTER OUT DIRECTLY ADJACENT INTERACTIONS
  
  filtered_baits <- all_baits[which(!(all_baits %in% unique(x$baitID)))]
  filtered_baits <- list(filtered_baits)
  fwrite(filtered_baits, paste0(outprefix, "_filteredBaits.txt"))
  
  x
}

.printMemory <- function(printMemory){
  if(printMemory == TRUE){
    print(gc(reset=TRUE))
  }
}

.locateFile = function(what, where, pattern){
  message("Locating ", what, " in ", where, "...")
  filename = list.files(where, pattern)
  
  if (length(filename)!=1){
    stop(paste0("Could not unambigously locate a ", what, " in ", where, ". Please specify explicitly in settings\n"))
  }
  
  message("Found ", filename)
  file.path(where, filename)
}

#--------------------------------------Pipeline ---------------------------------------------#

#Pipeline - many of these arguments have defaults, but are assigned here explicitly for clarity.

chicdiffPipeline <- function(chicdiff.settings)
{
  
  message("\n*** Running getRegionUniverse\n")
  RU <- getRegionUniverse(chicdiff.settings)  
  
  .printMemory(chicdiff.settings[["printMemory"]])
  
  message("\n*** Running getControlRegionUniverse\n")
  RUcontrol <- getControlRegionUniverse(chicdiff.settings, RU)
  
  .printMemory(chicdiff.settings[["printMemory"]])
  
  message("\n*** Running getFullRegionData\n")
  FullRegionData <- getFullRegionData(chicdiff.settings, RU, RUcontrol, suffix = "")
  
  .printMemory(chicdiff.settings[["printMemory"]])
  
  message("\n*** Running DESeq2Wrap for FullRegion\n")
  DESeqOut <- DESeq2Wrap(chicdiff.settings, RU, FullRegionData[[1]])
  
  .printMemory(chicdiff.settings[["printMemory"]])
  
  message("\n*** Running DESeq2Wrap for FullControlRegion\n")
  if(chicdiff.settings[["norm"]]=="combined" & 
       is.null(attributes(DESeqOut)$theta) & 
       is.null(chicdiff.settings[["theta"]])){
    warning("Normalisation weight theta is not defined and its inference may fail on control regions")
  }
  
  DESeqOutControl <- DESeq2Wrap(chicdiff.settings, RUcontrol, FullRegionData[[2]], 
                                suffix = "Control", theta = attributes(DESeqOut)$theta) 
  
  .printMemory(chicdiff.settings[["printMemory"]])
  
  message("\n*** Running IHWcorrection\n")
  output <- IHWcorrection(chicdiff.settings, DESeqOut, FullRegionData[[1]], 
                          DESeqOutControl, FullRegionData[[2]], countput = FullRegionData[[3]]) 
  
  .printMemory(chicdiff.settings[["printMemory"]])
  
  print(chicdiff.settings)
  print(sessionInfo())
  
  output
  
}

#-------------------------------Pipeline functions----------------------------------#

####------------------------------- getRegionUniverse-------------------------------#

.expandAvoidBait <- function(bait, oe, s)
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

getRegionUniverse <- function(chicdiff.settings, suffix = ""){
  
  RUexpand = chicdiff.settings[["RUexpand"]]
  rmapfile = chicdiff.settings[["rmapfile"]]
  peakFiles = chicdiff.settings[["peakfiles"]]
  targetColumns = chicdiff.settings[["targetColumns"]]
  chicagoData = chicdiff.settings[["chicagoData"]]
  score = chicdiff.settings[["score"]]
  saveRDS = chicdiff.settings[["saveAuxData"]]
  outprefix = chicdiff.settings[["outprefix"]]
  
  conditions <- names(chicagoData)
  
  x <- readAndFilterPeakMatrix(peakFiles = peakFiles, score = score, targetColumns = targetColumns, chicagoData=chicagoData, conditions = conditions, outprefix = outprefix)
  
  ## Expand the "point universe" to get "region universe" by "window" mode------------------
  
  #baits <- unique(x$baitID)
  #baits <- head(baits, 5000)
  #x <- x[baitID %in% baits,]
  
  ##Just expand calls by s in each direction
  
  RU.DT <- x[,c("baitID", "oeID"),with=FALSE]
  RU.DT$regionID <- 1:nrow(RU.DT)
  RU.DT <- RU.DT[,
                 list(otherEndID = .expandAvoidBait(baitID, oeID, RUexpand)), 
                 by=c("baitID","regionID")
                 ]  
  
  ## Ensure regions remain within the genome
  rmap <- fread(rmapfile)
  maxfrag <- max(rmap$V4)
  RU.DT <- RU.DT[otherEndID <= maxfrag,]
  
  ## Ensure that regions stay on the correct chromosome (*)
  
  ## get chr info
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
    saveRDS(RU.DT, paste0(outprefix, "_RegionUniverse", suffix, ".Rds"))
  }
  
  RU.DT
}

#-----------------------Functions needed for getControlRegionUniverse ----------------------------#

giveOneSeed <- function(bait, dist, min, max){
  ifelse(bait + dist < min, bait - dist, ifelse(bait + dist > max, bait - dist, bait + dist))   
}

giveDists <- function(bait, min, max, std){
  dist <- round(rnorm(n = 1, mean = 0, sd = std))
  
  onChr_notBait <- (((bait + abs(dist)) < max) | ((bait - abs(dist)) > min)) & (dist != 0)
  
  while(!onChr_notBait){
    dist <- round(rnorm(n = 1, mean = 0, sd = std))
    onChr_notBait <- (((bait + abs(dist)) < max) | ((bait - abs(dist)) > min)) & (dist != 0)
  }
  dist
}

giveManySeeds <- function(bait, min, max, std){
  dists <- unlist(lapply(bait, giveDists, min = min, max = max, std = std))
  mapply(giveOneSeed, bait, dists, MoreArgs = list(min = min, max = max))  
}

#-------------------------------------------getControlRegionUniverse----------------------------------------#

getControlRegionUniverse <- function(chicdiff.settings, RU){
  
  RUexpand <- chicdiff.settings[["RUexpand"]]
  saveRDS <- chicdiff.settings[["saveAuxData"]]
  outprefix <- chicdiff.settings[["outprefix"]]
  bmap <- fread(chicdiff.settings[["baitmapfile"]])
  rmap <- fread(chicdiff.settings[["rmapfile"]])
  colnames(rmap) <- c("chr", "start", "end", "ID")
  colnames(bmap) <- c("chr", "start", "end", "ID", "baitname")
  
  RUshorter <- merge(RU, rmap[,c("chr", "ID")], by.x = "baitID", by.y = "ID")
  max_contacts <- RUshorter[, list(max_contact = max(abs(baitID - otherEndID))),by = "chr"]
  
  RUcontrol <- as.data.table(sample(bmap$ID, size = length(unique(RU$regionID)), replace = TRUE))
  setnames(RUcontrol, "V1", "baitID")
  RUcontrol <- merge(bmap[chr %in% max_contacts$chr, c("chr", "ID")], RUcontrol, by.x = "ID", by.y = "baitID")
  
  for(ch in unique(RUcontrol[,chr])){
    
    std <- max_contacts[chr == ch, max_contact]/3
    max <- rmap[chr == ch, max(ID)]
    min <- rmap[chr == ch, min(ID)]
    
    RUcontrol[chr == ch, seed := giveManySeeds(ID, min = min, max = max, std = std)]
  }
  
  colnames(RUcontrol) <- c("baitID", "chr", "oeID")
  setkey(RUcontrol, baitID, oeID)
  RUcontrol[,regionID:=1:nrow(RUcontrol)]
  
  RUcontrol <- RUcontrol[,
                         list(otherEndID = .expandAvoidBait(baitID, oeID, RUexpand)), 
                         by=c("baitID","regionID")
                         ]  
  
  ## Ensure regions remain within the genome
  maxfrag <- max(rmap[,ID])
  RUcontrol <- RUcontrol[otherEndID <= maxfrag,]
  
  ## Ensure regions stay on the correct chromosome
  setkey(RUcontrol, otherEndID)
  setkey(rmap, ID)
  RUcontrol <- rmap[,c("chr", "ID"),with=FALSE][RUcontrol]
  setnames(RUcontrol, "ID", "otherEndID")
  setkey(RUcontrol, baitID)
  RUcontrol <- rmap[,c("chr", "ID"),with=FALSE][RUcontrol]
  setnames(RUcontrol, "ID", "baitID")
  
  ## keep only cis interactions
  RUcontrol <- RUcontrol[chr == i.chr,]
  RUcontrol[,chr:=NULL]
  RUcontrol[,i.chr:=NULL]
  
  if(saveRDS){
    saveRDS(RUcontrol, paste0(outprefix, "_ControlRegionUniverse.Rds"))
  }
  
  return(RUcontrol)
}

#-------------------------------Function combining 02, 03 and 04 R scripts ----------------------------------#
#Function should collect the count data from the chinputs, get the dispersions and 
#interaction parameters and then merge to give full data for each region

readRDSorRDA <- function(file){
  
  tmpEnv <- new.env()
  suppressWarnings(tmp <- try(load(file, envir = tmpEnv), silent = TRUE))
  if("try-error" %in% class(tmp)){
    #message("Not an RDA. Checking whether RDS.")
    suppressWarnings(tmp <- try(output <- readRDS(file), silent = TRUE))
    if("try-error" %in% class(tmp)){
      message("Cannot read file. Must be an RDS or RDA file.")
    } else{
      #message("Detected RDS")
      output
    }
  } else{
    #message("Detected RDA")
    get(ls(envir=tmpEnv)[1], envir=tmpEnv)
  }
}

#------------------------------------------------------------------------------#

.chicEstimateDistFun <- function (x, settings=Chicago::defaultSettings()) {
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

#--------------------------------------Non-parallel version--------------------------------------#

getFullRegionData1 <- function(chicdiff.settings, RU, is_control = FALSE, suffix = ""){
  
  countData = chicdiff.settings[["countData"]]
  chicagoData = chicdiff.settings[["chicagoData"]]
  rmapfile = chicdiff.settings[["rmapfile"]]
  saveRDS = chicdiff.settings[["saveAuxData"]]
  outprefix = chicdiff.settings[["outprefix"]]
  
  
  targetRDSorRDAFiles <- unlist(chicagoData)
  targetChFiles <- unlist(countData)
  
  ## 03getInteractionParameters
  
  RU_IntParams <- copy(RU)
  RU_IntParams[,regionID:=NULL]
  RU_IntParams <- unique(RU_IntParams)
  
  ##prep output
  outData <- vector("list", length(targetRDSorRDAFiles))
  dispersions <- numeric(length(targetRDSorRDAFiles))
  tempForCounts <- vector("list", length(targetRDSorRDAFiles))
  
  if(!is_control){
    
    conditions <- rep(names(chicagoData), sapply(chicagoData, length))
    countput <- lapply(unique(conditions), assign, 
                       value = vector("list", length(targetRDSorRDAFiles)))
    names(countput) <- unique(conditions)
  }
  
  ##1) Read in each file...
  for(i in 1:length(targetRDSorRDAFiles))
  {
    message("\nReading Chicago dataset ",i, " of ", length(targetRDSorRDAFiles), " : ", names(targetRDSorRDAFiles)[i])
    file <- targetRDSorRDAFiles[i]
    #load(file) ##this produces the object 'x' --old 
    x <- readRDSorRDA(file)
    ##1a) collect the dispersion
    
    if("chicagoData" %in% class(x)){
      dispersions[i] <- x@params$dispersion
      x <- as.data.table(x@x)
    } else{
      dispersions[i] <- attributes(x)$dispersion
      setDT(x)
    }
    
    # message("Finished reading")
    
    
    ##1b) collect the Bmean, Tmean information that is present
    
    message("Collecting Bmean and Tmean")
    
    setkey(x, baitID, otherEndID)
    setkey(RU_IntParams, baitID, otherEndID)
    out <- x[RU_IntParams, c("baitID", "otherEndID","distSign","Bmean", "Tmean", "score"), with=FALSE]
    ##(note that score=NA occurs when N=0 - hence replace score=0)
    
    
    ##1c) recalculate distSign (need to do this for control regions)
    message("Recalculating distSign")
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
    message("Collecting bait parameters (s_j and tblb)")
    temp <- x[,list(s_j=s_j[1], tblb=tblb[1]),by="baitID"]
    setkey(temp, baitID)
    setkey(out, baitID)
    out <- merge(out, temp, all.x=TRUE)
    
    ##1e) OE annotation: collect s_is, tlb
    ##s_i not allowed to be NA - just assume 1.
    ## => s_i should be 0)
    message("Collecting other end parameters (s_i and tlb)")
    temp <- x[,list(s_i=s_i[1], tlb=tlb[1]),by="otherEndID"]
    setkey(temp, otherEndID)
    setkey(out, otherEndID)
    out <- merge(out, temp, all.x=TRUE)
    out[is.na(s_i),s_i:=1] 
    setkey(out, baitID, otherEndID)
    
    
    ##1f) reconstruct Tmean information
    message("Collecting Tmean")
    setkey(x, tlb, tblb)
    setkey(out, tlb, tblb)
    temp <- x[,Tmean[1],by=c("tblb", "tlb")]
    setnames(temp, "V1", "Tmean")
    out[,Tmean:=NULL]
    out <- merge(out, temp, all.x=TRUE)
    
    ##1g) impute some of the missing Tmeans
    ##    missing tlb suggests that the other end was rarely observed, so assume lowest
    ##    Tmean corresponding to the appropriate tblb.
    message("Imputing missing Tmeans")
    templowest <- temp[,c(Tmean = min(Tmean)),by="tblb"]
    tempDictionary <- templowest$V1
    names(tempDictionary) <- templowest$tblb
    out[is.na(tlb) & !is.na(tblb), Tmean := tempDictionary[tblb]]
    
    ##2) Reconstruct the distance function from distbin info
    message("Reconstructing the distance function")
    distFunParams <- .chicEstimateDistFun(x)
    
    ##3) Reconstruct Bmean information (But!! When s_j = NA,
    ##   ensure that Bmean comes out as NA too.)
    message("Reconstructing Bmean")
    out <- Chicago:::.estimateBMean(out, distFunParams)
    out[is.na(s_j),Bmean:=NA]
    setkey(out, baitID, otherEndID)
    
    outData[[i]] <- out
    
    
    if(!is_control){
      rmap_copy <- copy(rmap)
      rmap_copy[, chr := NULL]
      rmap_copy[, midpoint := (start + end)/2]
      rmap_copy[, c("start", "end") := NULL]
      setnames(rmap_copy, "ID", "otherEndID")
      
      x <- x[!is.na(x[,distSign])]
      if("newScore" %in% names(x)){
        x[, names(x)[!(names(x) %in% c("baitID", "otherEndID", "N", "Bmean", "newScore"))] := NULL]
        setnames(x, "newScore", "score")
      } else{
        x[, names(x)[!(names(x) %in% c("baitID", "otherEndID", "N", "Bmean", "score"))] := NULL]
      }
      
      gc()
      x <- merge(x, rmap_copy, by= "otherEndID")
      gc()
      
      #message("Assigning conditions")
      for(n in 1:length(countput)){
        if(conditions[[i]] == names(countput)[[n]]){
          countput[[n]][[i]] <- x
        }
      }
      
      
    }
    
    #This just keeps all of the loaded RDAs if we still need them for the count data because we don't have chinputs
    
    if(is.null(countData)){ 
      
      x <- x[, c("baitID", "otherEndID", "N"), with = FALSE]
      columnNames <- paste0("N.", names(targetRDSorRDAFiles))
      setnames(x, old="N", new=columnNames[i])
      setkey(x, baitID, otherEndID)
      tempForCounts[[i]] <- x
      rm(x)
    }
    
    
    
  }
  
  message("")
  
  if(!is_control){
    message("Saving counts\n")
    
    for(i in 1:length(countput)){
      countput[[i]] <- rbindlist(countput[[i]])
      gc()
      countput[[i]] <- countput[[i]][,list(Nav = mean(N), Bav = mean(Bmean), score = max(score), midpoint = midpoint[1]), by = c("baitID", "otherEndID")]
      gc()
      countput[[i]][, condition := rep(names(countput)[[i]], nrow(countput[[i]]))]
    }
    
    countput <- rbindlist(countput)
    setnames(countput, "midpoint", "oeID_mid")
    
    saveRDS(countput, paste0(outprefix, "_countput.Rds"))
  }
  
  baits <- sort(unique(RU$baitID)) ##baits to get information from
  
  if(is.null(countData)){
    message("Reconstructing countData")
    
    gc()
    mergedFiles <- Reduce(merge, tempForCounts)
    
    countData <- vector("list", length(targetRDSorRDAFiles))
    for(i in seq_along(countData)){
      countData[[i]] <- mergedFiles[, c("baitID", "otherEndID", columnNames[i]), with = FALSE]
      setkey(countData[[i]], baitID)
      countData[[i]] <- countData[[i]][J(baits),]
      setkey(countData[[i]], baitID)
    }
    
    x <- RU
    setkey(x, baitID, otherEndID)
    
    message("Merging countData")
    
    for(i in seq_along(countData))
    {
      
      temp <- countData[[i]] 
      setnames(temp, old = columnNames[i], new = "N")
      setkey(temp, baitID, otherEndID)
      x <- merge(x, temp, all.x=TRUE)
      x[is.na(N), N:=0]
      setnames(x, old="N", new=columnNames[i])
      
    }
    
    # message("no chinput case finished")
    
  }
  
  #Recall that outData and dispersions are our useful outputs from this part above 
  
  ## 02getCounts.R
  
  ##Given a region universe, collect the number of reads in each region from the chinputs.
  
  ##Inputs and parameters:
  
  ##requires countData
  
  ##2) Collect all of the read count information ----------------
  if(!is.null(countData)){
    
    message("")
    countData <- vector("list", length(targetChFiles))
    gc()
    for(i in 1:length(countData))
    {
      message("Reading count data for ", names(targetChFiles)[i])
      x <- fread(targetChFiles[i])
      setkey(x, baitID)
      
      x <- x[J(baits),]
      setkey(x, baitID) ##Required for countData[4:7], where key is dropped - not sure why
      
      countData[[i]] <- x
      rm(x)
      
      
      gc()
    }
    
    ##2a) For each region, count the number of reads in each replicate -----------------
    
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
      
      #message("count the number of reads in each replicate for CountOut")
      
    }
    
  }
  
  ##3) Collect distance -----------------
  
  ##3a) Distance ------------
  
  message("\nCollecting distances for CountOut")
  
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
  
  message("Processing count data")
  
  setkey(CountOut, baitID, otherEndID)
  
  for(i in 1:length(targetRDSorRDAFiles))
  {
    s <- names(targetRDSorRDAFiles[i])
    temp <- outData[[i]][,c("baitID", "otherEndID", "s_j", "Bmean", "Tmean", "score"), with=FALSE]
    temp <- unique(temp)
    temp[,FullMean:=Bmean+Tmean]
    
    for(nm in colnames(temp))
      if(!nm %in% c("baitID", "otherEndID"))
      {
        setnames(temp, nm, paste0(nm, ".", s))
      }
    
    setkey(temp, baitID, otherEndID)
    
    CountOut <- merge(CountOut, temp, all.x=TRUE)
    
    #message("count + interaction parameter data together")
    
  }
  
  myCols <- c("N", "s_j", "Bmean", "Tmean", "score", "FullMean")
  getCols <- function(col) {grep(paste0("^", col), colnames(CountOut), perl=TRUE)}
  sels <- lapply(myCols, getCols)
  
  recast = melt(CountOut, measure = sels, value.name = myCols,
                variable.name="sample")
  recast$sample <- names(targetRDSorRDAFiles)[recast$sample]
  
  recast$condition <- rep(
    rep(names(chicagoData), sapply(chicagoData, length))
    , each=nrow(CountOut))
  
  
  setkey(recast, regionID)
  
  #message("recast")
  
  
  if(saveRDS == TRUE){
    message("Saving FullRegionData")
    if(!is_control){
      saveRDS(recast, paste0(outprefix, "_FullRegionData", suffix, ".Rds"))
    }else{
      saveRDS(recast, paste0(outprefix, "_FullControlRegionData", suffix, ".Rds"))
    }
  }
  
  if(!is_control){
    resultsList <- vector("list", length = 3)
    resultsList[[1]] <- recast
    resultsList[[3]] <- countput
    resultsList
  }else{
    recast
  }
  
}


#--------------------------- The version which reads significant and control interactions in parallel-------------#

getFullRegionData2 <- function(chicdiff.settings, RU, RUcontrol, suffix = ""){
  
  countData = chicdiff.settings[["countData"]]
  chicagoData = chicdiff.settings[["chicagoData"]]
  targetColumns = chicdiff.settings[["targetColumns"]]
  rmapfile = chicdiff.settings[["rmapfile"]]
  saveRDS = chicdiff.settings[["saveAuxData"]]
  outprefix = chicdiff.settings[["outprefix"]]
  
  targetRDSorRDAFiles <- unlist(chicagoData)
  targetChFiles <- unlist(countData)
  
  ## 03getInteractionParameters.R  
  
  ##Given a region universe, collect Bmean + Tmean  
  ##Inputs and parameters:
  ##requires targetRDAs
  
  RU_IntParams <- copy(RU) 
  RU_IntParams[,regionID:=NULL]
  RU_IntParams <- unique(RU_IntParams)
  
  RUcontrol_IntParams <- copy(RUcontrol)
  RUcontrol_IntParams[,regionID:=NULL]
  RUcontrol_IntParams <- unique(RUcontrol_IntParams)
  
  ##prep output
  outData <- vector("list", length(targetRDSorRDAFiles))
  outDataControl <- vector("list", length(targetRDSorRDAFiles)) #New to collect control data
  dispersions <- numeric(length(targetRDSorRDAFiles)) #only need one of these
  tempForCounts <- vector("list", length(targetRDSorRDAFiles))
  
  conditions <- rep(names(chicagoData), sapply(chicagoData, length))
  countput <- lapply(unique(conditions), assign, 
                     value = vector("list", length(targetRDSorRDAFiles)))
  names(countput) <- unique(conditions)
  
  ##1) Read in each file...
  for(i in 1:length(targetRDSorRDAFiles))
  { 
    message("\nReading Chicago dataset ",i, " of ", length(targetRDSorRDAFiles), " : ", names(targetRDSorRDAFiles)[i])
    file <- targetRDSorRDAFiles[i]
    x <- readRDSorRDA(file)
    ##1a) collect the dispersion
    
    if("chicagoData" %in% class(x)){
      dispersions[i] <- x@params$dispersion
      x <- as.data.table(x@x)
    } else{
      dispersions[i] <- attributes(x)$dispersion
      setDT(x)
    }
    
    #message("Finished reading")
    
    
    ##1b) collect the Bmean, Tmean information that is present
    message("Collecting Bmean and Tmean")
    
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
      message("Recalculating distSign for significant interactions")
      
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
    
    
    if(any(is.na(outControl$distSign)))
    {
      message("Recalculating distSign for control interactions")
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
    message("Collecting bait parameters (s_j and tblb)")
    temp <- x[,list(s_j=s_j[1], tblb=tblb[1]),by="baitID"]
    setkey(temp, baitID)
    setkey(out, baitID)
    setkey(outControl, baitID)
    out <- merge(out, temp, all.x=TRUE)
    outControl <- merge(outControl, temp, all.x = TRUE)
    
    ##1e) OE annotation: collect s_is, tlb
    ##s_i not allowed to be NA - just assume 1.
    ##(Could be an issue if too many OEs were filtered out => s_i should be 0)
    message("Collecting other end parameters (s_i and tlb)")
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
    
    #message("collected s_is and tlb")
    
    
    ##1f) reconstruct Tmean information
    message("Collecting Tmean")
    setkey(x, tlb, tblb)
    setkey(out, tlb, tblb)
    setkey(outControl, tlb, tblb)
    temp <- x[,Tmean[1],by=c("tblb", "tlb")]
    setnames(temp, "V1", "Tmean")
    out[,Tmean:=NULL]
    outControl[,Tmean:=NULL]
    out <- merge(out, temp, all.x=TRUE)
    outControl <- merge(outControl, temp, all.x=TRUE)
    
    #message("collected s_is and tlb")
    
    
    ##1g) impute some of the missing Tmeans
    ##    missing tlb suggests that the other end was rarely observed, so assume lowest
    ##    Tmean corresponding to the appropriate tblb.
    message("Imputing missing Tmeans")
    templowest <- temp[,c(Tmean = min(Tmean)),by="tblb"]
    tempDictionary <- templowest$V1
    names(tempDictionary) <- templowest$tblb
    out[is.na(tlb) & !is.na(tblb), Tmean := tempDictionary[tblb]]
    outControl[is.na(tlb) & !is.na(tblb), Tmean := tempDictionary[tblb]]
    
    ##2) Reconstruct the distance function from distbin info
    message("Reconstructing the distance function")
    distFunParams <- .chicEstimateDistFun(x)
    
    ##3) Reconstruct Bmean information (But!! When s_j = NA,
    ##   ensure that Bmean comes out as NA too.)
    message("Reconstructing Bmean")
    
    out <- Chicago:::.estimateBMean(out, distFunParams)
    out[is.na(s_j),Bmean:=NA]
    setkey(out, baitID, otherEndID)
    
    outControl <- Chicago:::.estimateBMean(outControl, distFunParams)
    outControl[is.na(s_j),Bmean:=NA]
    
    setkey(outControl, baitID, otherEndID)
    
    outData[[i]] <- out
    outDataControl[[i]] <- outControl
    
    
    #This keeps all of the loaded RDAs if we still need them for the count data because we don't have chinputs
    rmap_copy <- copy(rmap)
    rmap_copy[, chr := NULL]
    rmap_copy[, midpoint := (start + end)/2]
    rmap_copy[, c("start", "end") := NULL]
    setnames(rmap_copy, "ID", "otherEndID")
    
    x <- x[!is.na(x[,distSign])]
    if("newScore" %in% names(x)){
      x[, names(x)[!(names(x) %in% c("baitID", "otherEndID", "N", "Bmean", "newScore"))] := NULL]
      setnames(x, "newScore", "score")
    } else{
      x[, names(x)[!(names(x) %in% c("baitID", "otherEndID", "N", "Bmean", "score"))] := NULL]
    }
    
    x <- merge(x, rmap_copy, by= "otherEndID")
    
    #message("Assigning conditions")
    for(n in 1:length(countput)){
      if(conditions[[i]] == names(countput)[[n]]){
        countput[[n]][[i]] <- x
      }
    }
    
    message("")
    
    if(is.null(countData)){ 
      x <- x[, c("baitID", "otherEndID", "N"), with = FALSE]
      columnNames <- paste0("N.", names(targetRDSorRDAFiles))
      setnames(x, old="N", new=columnNames[i])
      setkey(x, baitID, otherEndID)
      tempForCounts[[i]] <- x
      rm(x)
      
    }
    gc()
  }
  
  for(i in 1:length(countput)){
    countput[[i]] <- rbindlist(countput[[i]])
    gc()
    countput[[i]] <- countput[[i]][,list(Nav = mean(N), Bav = mean(Bmean), score = max(score), midpoint = midpoint[1]), by = c("baitID", "otherEndID")]
    gc()
    countput[[i]][, condition := rep(names(countput)[[i]], nrow(countput[[i]]))]
  }
  
  countput <- rbindlist(countput)
  setnames(countput, "midpoint", "oeID_mid")
  
  message("Saving counts\n")
  saveRDS(countput, paste0(outprefix, "_countput.Rds"))
  
  
  # outData, outDataControl are our useful outputs from this part above  
  
  
  baits <- sort(unique(RU$baitID)) 
  baitsControl <- sort(unique(RUcontrol$baitID))
  
  ## 02getCounts.R  
  ## Given a region universe, collect the number of reads in each region
  
  ##2-Rds) If Chinputs are not given: Collect all of the read count information from Rds
  
  if(is.null(countData)){
    gc()
    message("Reconstructing countData")
    
    mergedFiles <- Reduce(merge, tempForCounts)
    countData <- vector("list", length(targetRDSorRDAFiles))
    countDataControl <- vector("list", length(targetRDSorRDAFiles))
    
    for(i in seq_along(countData)){
      countData[[i]] <- mergedFiles[, c("baitID", "otherEndID", columnNames[i]), with = FALSE]
      setkey(countData[[i]], baitID)
      countData[[i]] <- countData[[i]][J(baits),] #issue that rds may not contain all baitid otherendid pairs which are in RU
      setkey(countData[[i]], baitID)
      
    }
    
    for(i in seq_along(countDataControl)){
      countDataControl[[i]] <- mergedFiles[, c("baitID", "otherEndID", columnNames[i]), with = FALSE]
      setkey(countDataControl[[i]], baitID)
      countDataControl[[i]] <- countDataControl[[i]][J(baitsControl),] #issue that rds may not contain all baitid otherendid pairs which are in RU
      setkey(countDataControl[[i]], baitID)
      
    }
    
    CountOut <- RU
    CountOutControl <- RUcontrol
    setkey(CountOut, baitID, otherEndID)
    setkey(CountOutControl, baitID, otherEndID)
    
    
    message("Merging countData")
    for(i in seq_along(countData))
    {
      temp <- countData[[i]] #[,c("baitID","otherEndID","N"), with=FALSE]
      setnames(temp, old = columnNames[i], new = "N")
      setkey(temp, baitID, otherEndID)
      CountOut <- merge(CountOut, temp, all.x=TRUE)
      CountOut[is.na(N), N:=0]
      setnames(CountOut, old="N", new=columnNames[i])
      
      #message("countData info merged")
      
      
    }
    
    for(i in seq_along(countDataControl))
    {
      temp <- countDataControl[[i]] #[,c("baitID","otherEndID","N"), with=FALSE]
      setnames(temp, old = columnNames[i], new = "N")
      setkey(temp, baitID, otherEndID)
      CountOutControl <- merge(CountOutControl, temp, all.x=TRUE)
      CountOutControl[is.na(N), N:=0]
      setnames(CountOutControl, old="N", new=columnNames[i])
      
      #message("countDatacontrol info merged")     
    }
    #message("no chinput case finished")
    
  }
  
  ##2-Chi) If Chinputs are given: Collect all of the read count information ----------------
  if(!is.null(countData)){
    
    countData <- vector("list", length(targetChFiles))
    countDataControl <- vector("list", length(targetChFiles))
    gc()
    
    for(i in 1:length(countData)){
      message("Reading count data for ", names(targetChFiles)[i])
      
      x <- fread(targetChFiles[i])
      setkey(x, baitID)
      
      y <- x[J(baitsControl),]
      x <- x[J(baits),]
      
      setkey(x, baitID) 
      setkey(y, baitID)
      
      countData[[i]] <- x
      countDataControl[[i]] <- y
      rm(x)
      rm(y)
      
      #message("Finished reading of chinputs")
      
      gc()
    }
    
    ## For each region, count the number of reads in each replicate -----------------
    
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
      
      #message("count the number of reads in each replicate for CountOut")
      
    }
    
    #Repeating the above code for the control
    CountOutControl <- RUcontrol 
    setkey(CountOutControl, baitID, otherEndID)
    for(i in 1:length(countDataControl))
    {
      columnName <- paste0("N.", names(targetChFiles)[i])
      temp <- countDataControl[[i]][,c("baitID","otherEndID","N"), with=FALSE]
      setkey(temp, baitID, otherEndID)
      CountOutControl <- merge(CountOutControl, temp, all.x=TRUE)
      
      CountOutControl[is.na(N), N:=0]
      setnames(CountOutControl, old="N", new=columnName)
      
      #message("count the number of reads in each replicate for CountOutControl")
      
    }
  }
  
  ##3) Collect distance -----------------
  
  ##3a) Distance ------------
  
  message("\nCollecting distances for CountOut")
  
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
  
  #message("collected distances for CountOut")
  
  
  #Just repeating the above code for the control
  CountOutControl <- merge(CountOutControl, rmap, by.x="otherEndID", by.y="fragID")  
  CountOutControl <- merge(CountOutControl, rmap, by.x="baitID", by.y="fragID")
  setkey(CountOutControl, baitID, otherEndID)
  
  CountOutControl[, distSign :=
                    ifelse(chr.x == chr.y,
                           midpoint.x - midpoint.y,
                           NA)]
  CountOutControl[, c("chr.x","midpoint.x","chr.y","midpoint.y"):=NULL]
  
  #message("collected distances for CountOutcontrol") 
  
  # 04mergeData.R
  
  ##Synthesise the count + interaction parameter data together
  
  message("Processing count data")
  
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
    
    for (nm in colnames(temp))
    { 
      if(!nm %in% c("baitID", "otherEndID"))
      {
        setnames(temp, nm, paste0(nm, ".", s))
      }}
    
    for (nmc in colnames(tempControl))
    {  
      if(!nmc %in% c("baitID", "otherEndID"))
      {
        setnames(tempControl, nmc, paste0(nmc, ".", s))
      }}
    
    setkey(temp, baitID, otherEndID)
    setkey(tempControl, baitID, otherEndID)
    
    CountOut <- merge(CountOut, temp, all.x=TRUE)
    CountOutControl <- merge(CountOutControl, tempControl, all.x=TRUE)
    
    #message("count + interaction parameter data together")
    
    gc()
    
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
    rep(names(chicagoData), sapply(chicagoData, length))
    , each=nrow(CountOut))
  
  setkey(recast, regionID)
  
  #message("recast")
  
  
  # Again - just a copy and paste for the control here 
  recastControl = melt(CountOutControl, measure = sels, value.name = myCols,
                       variable.name="sample")
  recastControl$sample <- names(targetRDSorRDAFiles)[recastControl$sample]
  
  recastControl$condition <- rep(
    rep(names(chicagoData), sapply(chicagoData, length))
    , each=nrow(CountOutControl))
  
  setkey(recastControl, regionID)
  
  #message("recastControl")
  
  
  if(saveRDS == TRUE){
    message("Saving FullRegionData")
    saveRDS(recast, paste0(outprefix, "_FullRegionData", suffix, ".Rds"))
    message("Saving FullControlRegionData")
    saveRDS(recastControl, paste0(outprefix, "_FullControlRegionData", suffix, ".Rds"))
    
  }
  
  #message("reported FullRegionDatalist")
  
  
  return(list(recast, recastControl, countput))
}

#-------------------------------The full getFullRegionDataWrapper----------------------------------------#

getFullRegionData <- function(chicdiff.settings, RU, RUcontrol, suffix = ""){
  
  parallel = chicdiff.settings[["parallel"]]
  
  if(parallel){
    
    FullRegionData <- getFullRegionData2(chicdiff.settings = chicdiff.settings, RU = RU, RUcontrol = RUcontrol, suffix = "")
    
  } else{
    
    message("Reading data for significant interactions")
    FullRegionData <- getFullRegionData1(chicdiff.settings = chicdiff.settings, RU = RU, is_control = FALSE, suffix = "")
    message("\nReading data for control interactions")
    FullRegionData[[2]] <- getFullRegionData1(chicdiff.settings = chicdiff.settings, RU = RUcontrol, is_control = TRUE, suffix = "control")
    
  }
  
  return(FullRegionData)
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

DESeq2Wrap <- function(chicdiff.settings, RU, FullRegionData, suffix = "", theta = NULL){
  
  Grid = chicdiff.settings[["theta_grid"]]
  rmapfile = chicdiff.settings[["rmapfile"]]
  saveRDS = chicdiff.settings[["saveAuxData"]]
  outprefix = chicdiff.settings[["outprefix"]]
  
  if (is.null(theta) & !is.null(chicdiff.settings[["theta"]])){
    theta = chicdiff.settings[["theta"]]
    #message("Mixing parameter theta set to the user-specified value of ", theta)
  }
  
  norm = chicdiff.settings[["norm"]]
  if(!norm %in% c("standard", "fullmean", "combined")){
    stop("DESeq2Wrap error: Unknown normalisation method.")
  }
  
  if (!is.null(theta)){
    if(theta == 1 & norm!="standard"){
      warning("Mixing parameter theta set to 1, equivalent to norm = \"standard\". The norm method has been reset accordingly.")
      norm = "standard"
    }
    
    if(!theta & norm!="fullmean"){
      warning("Mixing parameter theta set to 0, equivalent to norm = \"fullmean\". The norm method has been reset accordingly.")
      norm = "fullmean"
    }
  }
  
  ##Input data:
  
  fragData <- copy(FullRegionData) ## As otherwise it is altered by reference below and the IHW part breaks
  setkey(fragData, otherEndID)
  
  #message("saved a copy of FullRegionData")
  
  #fragData.impute <- fragData 
  #fragData.impute[,s_j.impute := geoMean(s_j, na.rm = TRUE), by=c("baitID", "otherEndID")]
  #fragData.impute[is.na(s_j), s_j:=s_j.impute]
  #fragData.impute[,FullMean.impute := geoMean(FullMean, na.rm = TRUE), by=c("baitID", "otherEndID")]
  #fragData.impute[is.na(FullMean), FullMean:=FullMean.impute]
  
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
  
  nSamples <- length(unique(regionData$sample))
  
  regionDataMatrix <- matrix(regionData$N,
                             ncol=nSamples,
                             byrow = TRUE)
  rownames(regionDataMatrix) <- unique(regionData$regionID)
  colnames(regionDataMatrix) <- unique(regionData$sample)
  colData <- data.frame(condition = fragData$condition[1:ncol(regionDataMatrix)])
  dds <- DESeqDataSetFromMatrix(countData = regionDataMatrix,
                                colData = colData,
                                design = ~ condition)
  
  dds.nullModel <- estimateSizeFactors(dds)
  nullSizeFactors <- sizeFactors(dds.nullModel)
  
  if(norm != "standard"){
    dds.M3 <- copy(dds.nullModel)
  }
  if(norm == "combined"){
    dds.M5 <- copy(dds.nullModel)
  }
  
  ##model 1) standard DESeq2 model
  if (norm == "standard"){
    dds.nullModel <- estimateDispersions(dds.nullModel)
    dds.nullModel <- nbinomWaldTest(dds.nullModel)
  }
  
  ##model 2) skipped
  
  #browser()
  
  ##model 3) use fullMean scaling factors
  if (norm != "standard"){
    normFactorsFull <- matrix(regionData$FullMean, ncol=nSamples, byrow = TRUE)
    
    normFactorsM3 <- normFactorsFull
    normFactorsM3 <- normFactorsM3 / exp(rowMeans(log(normFactorsM3))) ##<- normalize to get geometric mean = 1
    ##remove NA rows
    selNA <- apply(normFactorsM3, 1, function(x){any(is.na(x))})
    normFactorsM3[selNA,] <- rep(nullSizeFactors, each=sum(selNA))
    
    
    # sf2 <- estimateSizeFactorsForMatrix(counts(dds.M3)/normFactorsM3)
    # sf2mat <- matrix(rep(sf2,nrow(normFactorsM3)), 
    #                  ncol=nSamples, nrow=nrow(normFactorsM3), byrow=T) 
    
  }
  
  if (norm == "fullmean"){
    
    normalizationFactors(dds.M3) <- normFactorsM3
    
    dds.M3 <- estimateDispersions(dds.M3)
    dds.M3 <- nbinomWaldTest(dds.M3)
  }
  
  ## model 4) skipped
  
  ##model 5) use a weighted mean of fullMean or DESeq2 scaling factors
  ##that minimises the overall variance of the counts 
  ##for the whole sample
  
  if (norm=="combined"){
    
    nsf = matrix(rep(nullSizeFactors,nrow(normFactorsM3)), 
                 ncol=nSamples, nrow=nrow(normFactorsM3), byrow=T) 
    
    tt = theta
    
    if(is.null(tt)){
      
      message("Optimising scaling factors...")
      
      glen = length(Grid)
      # grid_normFactorsM5 <- matrix(nrow=nrow(normFactorsM3),ncol=glen*nSamples)
      # varsM5 = matrix(nrow = nrow(normFactorsM3),ncol=glen)
      i = 1
      deviances = vector("numeric", length=glen)
      
      ddsTest <- DESeqDataSetFromMatrix(countData = regionDataMatrix,
                                        colData = colData,
                                        design = ~ 1) # sic!
      
      for(tt in Grid){
        
        sc <- normFactorsM3*(1-tt)+nsf*tt
        # sc <- exp (log(normFactorsM3)*(1-tt)+log(nsf)*tt)
        
        sc <- sc / exp(rowMeans(log(sc)))
        
        
        normalizationFactors(ddsTest) <- sc
        
        ddsTest <- estimateDispersions(ddsTest)
        ddsTest <- nbinomWaldTest(ddsTest)
        
        #deviances[i] = median(mcols(ddsTest)$deviance)
        deviances[i] = sum(mcols(ddsTest)$deviance)#, na.rm = TRUE)
        
        # grid_normFactorsM5[,(1+(i-1)*nSamples):(i*nSamples)] <- sc
        
        # varsM5[,i] <- t(rowVars(counts(dds.M5)/sc))
        
        i=i+1
        
      }
      
      message("Total deviances by theta (Fullmean --> Standard):")
      cat(sprintf("%f", deviances), "\n", file=stderr())
      
      tt = Grid[which(deviances==min(deviances)[1])]
      
    } 
    
    message("Theta=", tt)
    
    sc <- normFactorsM3*(1-tt)+nsf*tt
    # sc <- exp (log(normFactorsM3)*(1-tt)+log(nsf)*tt)
    
    sc <- sc / exp(rowMeans(log(sc)))
    
    normalizationFactors(dds.M5) <- sc
    
    dds.M5 <- estimateDispersions(dds.M5)
    dds.M5 <- nbinomWaldTest(dds.M5)
    
    # varsM5 <- cbind(varsM5, grid_normFactorsM5)
    # 
    # normFactorsM5 <- t(apply(varsM5,1,function(x){
    #   whichminvar <- which(x[1:glen]==min(x[1:glen]))[1]
    #   x[(1+glen+(whichminvar-1)*nSamples):(glen+whichminvar*nSamples)]
    # }))
    # 
    # # Purely for diagnostic purposes    
    # whichMinVar5 <- apply(varsM5, 1, function(x)
    #   which(x[1:glen]==min(x[1:glen]))[1])
    
    # normalizationFactors(dds.M5) <- normFactorsM5
    # 
    # dds.M5 <- estimateDispersions(dds.M5)
    # dds.M5 <- nbinomWaldTest(dds.M5)
  }
  
  
  
  ##get annotation information -----------
  
  message("Processing model output")
  
  ##other ends:
  rmap <- fread(rmapfile)
  colnames(rmap) <- c("OEchr", "OEstart", "OEend", "otherEndID")
  
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
  
  
  if (norm == "standard" ){
    res = results(dds.nullModel)
    message("Standard DESeq2 normalisation: # unweighted interactions with padj<0.05: ",
            nrow(res[res$padj<0.05 & !is.na(res$padj),]))
    if(saveRDS == TRUE){
      message("Saving the DESeq object")
      saveRDS(dds.nullModel, paste0(outprefix, "_DESeqObj", suffix, ".Rds"))
    }
  }
  if (norm == "fullmean"){
    res = results(dds.M3)
    message("Chicago full mean-based normalisation: # unweighted interactions with padj<0.05: ",
            nrow(res[res$padj<0.05 & !is.na(res$padj),]))
    if(saveRDS == TRUE){
      message("Saving the DESeq object")
      saveRDS(dds.M3, paste0(outprefix, "_DESeqObj", suffix, ".Rds"))
    }
  }
  if (norm == "combined"){
    res = results(dds.M5)
    message("combined normalisation: # unweighted interactions with padj<0.05: ",
            nrow(res[res$padj<0.05 & !is.na(res$padj),]))
    # message("Weights of standard [vs Fullmean] scaling factors used for normalisation")
    # w = matrix(c(Grid, table(whichMinVar5)/length(whichMinVar5)), nrow=2, byrow=T, 
    #            dimnames = list(c("weight", "% cases")))
    # print(w)
    if(saveRDS == TRUE){
      message("Saving the final DESeq object")
      saveRDS(dds.M5, paste0(outprefix, "_DESeqObj", suffix, ".Rds"))
    }
  }
  
  results <- as.data.table(as.data.frame(res), keep.rownames = "id")
  results$id <- as.integer(results$id)
  setkey(results, id)
  
  out <- cbind(results,annoData)
  out[,id := NULL]
  
  if(exists("tt"))
    attributes(out)$theta <- tt
  
  out 
  
  # results.null <- as.data.table(as.data.frame(results(dds.nullModel)), keep.rownames = "id")
  # results.null$id <- as.integer(results.null$id)
  # setkey(results.null, id)
  # 
  # out.nullModel <- cbind(
  #   results.null,
  #   annoData
  # )
  # 
  # out.nullModel[,id := NULL]
  # 
  # message("out.nullModel")
  
}

#-----------------------------IHW and plotting functions-------------------------------#

plotDiffBaits <- function(output, countput, baitmapfile, baits = NULL, 
                          n = 3, plotBaitNames = TRUE, 
                          plotBaitIDs = TRUE, pcol="weighted_padj", 
                          plevel1 = 5, plevel2 = 3, xlim=c(-1e6,1e6), bgCol = "black", 
                          lev1Col = "red", lev2Col = "blue", plotBmean = FALSE, plotKey = TRUE,
                          ...){
  
  if(any(class(output) == "character")){
    if(file.exists(output)){
      out <- readRDS(output)
    } else{
      "No output file found at specified location"
    }
  } else{
    out <- copy(output)
  }
  
  if(any(class(countput) == "character")){
    if(file.exists(countput)){
      countput_coord <- readRDS(countput)
    } else{
      "No countput file found at specified location"
    }
  } else{
    countput_coord <- copy(countput)
  }
  
  if(any(class(baitmapfile) == "character")){
    if(file.exists(baitmapfile)){
      bmap <- fread(baitmapfile)
    } else{
      "No baitmap file found at specified location"
    }
  } else{
    bmap <- copy(baitmapfile)
  }
  
  conditions <- unique(countput_coord[,condition])
  
  out[, names(out)[!(names(out) %in% c("OEstart", "OEend","baitID", pcol))] := NULL]
  
  colnames(bmap) <- c("chr", "start", "end", "baitID", "name")
  bmap[, chr := NULL]
  
  plot_list <- vector("list", length = n)
  
  if(is.null(baits)){
    baits <- sample(unique(out$baitID), n)
  } else if(!all(baits %in% unique(out$baitID))){
    stop("One or more specified baits not in final Chicdiff output")
  }
  
  chroms_bait <- bmap[bmap[,baitID %in% baits]]
  chroms_bait[, midpoint := (start+end)/2]
  chroms_bait[, c("start", "end") := NULL]
  
  if(plotKey){  
    scale_data <- data.frame(value = 1L, name = factor(c("> 0.05", "0.005 - 0.05", "0.0005 - 0.005", "< 0.0005"), levels = c("< 0.0005", "0.0005 - 0.005", "0.005 - 0.05", "> 0.05")))
    scalebar <- ggplot(scale_data, aes(x = value, y = name)) +
      geom_tile(aes(fill = name)) +
      scale_fill_manual(name = "p-values", values = c("> 0.05" = "#2b8cbe", "0.005 - 0.05" = "#fee8c8", "0.0005 - 0.005" = "#fdbb84", "< 0.0005" = "#e34a33"))
  }
  
  for (i in 1:nrow(chroms_bait))
  {  
    count_temp <- countput_coord[baitID == chroms_bait[i][,baitID]]
    count_temp[score > plevel1, big_score := "red"]
    count_temp[score > plevel2 & score < plevel1, big_score := "blue"]
    count_temp[is.na(big_score), big_score := "black"]
    count_temp[,big_score := factor(big_score, levels = c("black", "blue", "red"))]
    
    xlimit <- c(chroms_bait[i][,midpoint + xlim[1]], chroms_bait[i][,midpoint + xlim[2]])
    ylim_max <- count_temp[,max(Nav)]
    
    p1 <- ggplot(count_temp[condition == conditions[[1]]], aes(oeID_mid, Nav)) +
      geom_point(aes(colour=big_score, alpha = 0.4, stroke = 0)) + 
      scale_colour_manual(values = c(bgCol, lev2Col, lev1Col)) 
    
    if (plotBmean){
      p1 <- p1 + geom_line(aes(oeID_mid, Bav))
    }
    
    p1 <- p1 + 
      xlim(xlimit[[1]], xlimit[[2]]) +
      ylim(0, ylim_max) +
      theme(panel.background = element_rect(fill = "white", 
                                            colour = NA),
            panel.grid = element_line(colour = "grey87"), 
            panel.grid.major = element_line(size = rel(0.5)), 
            panel.grid.minor = element_line(size = rel(0.25)), 
            legend.position = "none",
            axis.text.x = element_text(size = 8),
            axis.text.y = element_text(size = 8),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 8),
            plot.title = element_text(size = 10, hjust = 0.5)) +
      ylab("N")
    
    if(plotBaitNames & plotBaitIDs){
      p1 <- p1 + ggtitle(paste0(strsplit(strsplit(as.character(chroms_bait[i][,name]), ",")[[1]], "-")[[1]], " (", chroms_bait[i][,baitID], ")"), subtitle = conditions[[1]]) + theme(plot.subtitle = element_text(hjust = 1))
    }else if(plotBaitNames){
      p1 <- p1 + ggtitle(strsplit(strsplit(as.character(chroms_bait[i][,name]), ",")[[1]], "-")[[1]], subtitle = conditions[[1]]) + theme(plot.subtitle = element_text(hjust = 1))
    } else if(plotBaitIDs){
      p1 <- p1 + ggtitle(chroms_bait[i][,baitID], subtitle = conditions[[1]]) + theme(plot.subtitle = element_text(hjust = 1))
    }
    
    p2 <- ggplot(count_temp[condition == conditions[[2]]], aes(oeID_mid, Nav)) +
      geom_point(aes(colour=big_score, alpha = 0.4, stroke = 0)) + 
      scale_colour_manual(values = c(bgCol, lev2Col, lev1Col)) + 
      labs(caption = conditions[[2]]) + 
      theme(plot.caption = element_text(hjust = 1))
    
    if(plotBmean){
      p2 <- p2 + geom_line(aes(oeID_mid, Bav))
    }
    
    
    p2 <- p2 +       #     xlim(xlim[[1]], xlim[[2]]) +
      scale_y_reverse(limits = c(ylim_max, 0)) +
      theme(
        panel.background = element_rect(fill = "white", 
                                        colour = NA),
        panel.grid = element_line(colour = "grey87"), 
        panel.grid.major = element_line(size = rel(0.5)), 
        panel.grid.minor = element_line(size = rel(0.25)), 
        legend.position = "none", 
        plot.title = element_text(hjust = 0.5, size=10),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8)) + 
      ylab("N") +
      xlab("        ") +
      scale_x_continuous(position = "top", limits = c(xlimit[[1]], xlimit[[2]]))
    
    
    temp_out <- out[baitID == chroms_bait[i][,baitID]]
    ir_chr <- IRanges(start = temp_out[,OEstart], width = (temp_out[, OEend - OEstart] + 1))
    bins_chr <- disjointBins(IRanges(start(ir_chr), end(ir_chr)))
    dat <- as.data.table(ir_chr)[, bin := bins_chr]
    setnames(dat, c("start", "end"), c("OEstart", "OEend"))
    
    merged_dat <- merge(dat, temp_out, by = c("OEstart", "OEend"))
    merged_dat[,width := NULL]
    merged_dat[,minuslogpvalue := -log10(get(pcol))]
    
    merged_dat[minuslogpvalue <= -log(0.05), pvalue_param := "blue"]
    merged_dat[minuslogpvalue <= -log10(0.005) & minuslogpvalue > -log10(0.05), pvalue_param := "light red"]
    merged_dat[minuslogpvalue <= -log10(0.0005) & minuslogpvalue > -log10(0.005), pvalue_param := "mid red"]
    merged_dat[is.na(pvalue_param), pvalue_param := "red"]
    merged_dat[,pvalue_param := factor(pvalue_param, levels = c("blue", "light red", "mid red", "red"))]   
    
    pvalue_plot <- ggplot(merged_dat[OEstart >= xlimit[[1]] & OEstart <= xlimit[[2]]]) + 
      geom_rect(aes(xmin = OEstart, xmax = OEend,
                    ymin = bin, ymax = bin + 0.9,
                    fill = pvalue_param), colour = "white", size = 0.3) +
      scale_fill_manual(values = c("blue" = "#2b8cbe", "light red" = "#fee8c8", "mid red" = "#fdbb84", "red" = "#e34a33")) +
      xlim(xlimit[[1]], xlimit[[2]]) +
      theme_void() +
      theme(legend.position = "none")
    
    suppressWarnings(plot_list[[i]] <- plot_grid(p1, pvalue_plot, p2, align = "v", nrow = 3, rel_heights = c(1/2, 1/8, 1/2)))
  }
  
  suppressWarnings(chicdiffPlots <- plot_grid(plotlist = plot_list, align = "h"))
  if (plotKey){
    suppressWarnings(scalebar_grid <- plot_grid(get_legend(scalebar)))
    suppressWarnings(plot_grid(chicdiffPlots, scalebar_grid, rel_widths = c(10, 1), ...))
  }else{
    suppressWarnings(plot_grid(chicdiffPlots, rel_widths = c(10, 1), ...))
  }
}

#---------------------------------------------------#

IHWcorrection <- function(chicdiff.settings, DESeqOut, FullRegionData, DESeqOutControl, FullControlRegionData,
                          countput, DiagPlot = TRUE, diffbaitPlot = TRUE, suffix = ""){
  
  baitmapfile = chicdiff.settings[["baitmapfile"]]
  device = chicdiff.settings[["device"]]
  outprefix = chicdiff.settings[["outprefix"]]
  
  out <- DESeqOut ##DESeqOut and DESeqOutControl
  RU.recast <- FullRegionData ##FullRegionData and FullControlRegionData
  RU.distances <- RU.recast[,list(avDist = mean(distSign)),by="regionID"]
  
  out$avDist <- RU.distances$avDist
  
  ##Comparison against genuinely uniform p-vals
  out$uniform <- runif(nrow(out))
  out$shuff <- sample(out$pvalue)
  
  message("Comparison against p-vals for out")
  
  
  setDT(out)  
  out.control <- DESeqOutControl ##DESeqOut and DESeqOutControl
  RU.recastControl <- FullControlRegionData ##FullRegionData and FullControlRegionData
  RU.distancesControl <- RU.recastControl[,list(avDist = mean(distSign)),by="regionID"]
  
  out.control$avDist <- RU.distancesControl$avDist
  
  ##Comparison against genuinely uniform p-vals
  out.control$uniform <- runif(nrow(out.control))
  out.control$shuff <- sample(out.control$pvalue)
  
  message("Comparison against p-vals for outcontrol")
  
  
  ##Convert
  as.data.frame(out.control) #This line is necessary even if the one for out may not be
  
  ##Train weights on the control sample
  ihwRes <- ihw(pvalue ~ abs(avDist),  data = out.control, alpha = 0.05)
  
  message("Trained weights on the control sample")
  
  
  if(DiagPlot == TRUE){
    plot(ihwRes)
    ggsave(paste0(outprefix, "_IHWweightPlot.png"), device = "png", path = "./")
    plot(ihwRes, what = "decisionboundary")
    ggsave(paste0(outprefix, "_IHWdecisionBoundaryPlot.png"), device = "png", path = "./")
  }
  
  
  
  ##---------------------------------------------------------------------
  ##Learn distance dependency
  
  test <- ihwRes@df
  setDT(test)
  test[,log(covariate), by = "group"]
  
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
  
  message("Learned distance dependency")
  
  
  ##--------------------------------------------------------------
  ##Apply to test data
  
  out[,avgLogDist := log(abs(avDist))]
  breaks <- (c(distLookup$minLogDist, Inf) + c(0, distLookup$maxLogDist))/2
  out$group <- as.integer(cut(log(abs(out$avDist)), breaks))
  
  setkey(distLookup, group)
  setkey(out, group)
  
  out <- merge(out, distLookup[,c("group", "avWeights"),with=FALSE], all.x=TRUE, by="group")
  out$weight <- out$avWeights/mean(out$avWeights) ##renormalize
  out[,weighted_pvalue := pvalue/weight]
  
  out[, weighted_padj := p.adjust(weighted_pvalue, method="BH")]
  
  message("applied to test data")
  
  
  if (diffbaitPlot == TRUE){
    sel <- order(out$weighted_padj)
    baits <- sample(head(unique(out[sel]$baitID), 100), 4)
    plotDiffBaits(output = out, countput = countput, baitmapfile = baitmapfile, baits = baits)
    cowplot::ggsave2(paste0(outprefix, "_diffbaitPlot",".",device), device = device, path = "./") 
    dev.off()
  }
  
  saveRDS(out, paste0(outprefix, "_results", suffix, ".Rds"))
  
  return(out)
}


getCandidateInteractions <- function(output, peakFiles, 
                                     chicdiff.settings, 
                                     pcol = "weighted_padj",
                                     method = c("min", "hmp")[1], 
                                     minDeltaAsinhScore=1,                                      
                                     pvcut = 0.05){
  
  chicagoData = chicdiff.settings[["chicagoData"]]
  score <- chicdiff.settings[["score"]]
  targetColumns = chicdiff.settings[["targetColumns"]]
  peakFiles <- fread(peakFiles)
  
  sel <- which(colnames(peakFiles) %in% targetColumns)
  peakFiles <- peakFiles[,c(1:11, sel), with=FALSE] # assumes the standard peakFiles format where condition-specific Chicago scores start from col 12
  sel <- rep(FALSE, nrow(peakFiles))
  for(cl in targetColumns)
  {
     sel <- sel | peakFiles[,get(cl) > score & !is.na(get(cl))] ##Get any rows where at least one score > 5.
  }
  peakFiles <- peakFiles[sel,]
    
    
    
  if(is.character(output)){
    output <- readRDS(output)
  }
  
  if (!method %in% c("min", "hmp")){
    stop ("getCandidateInteractions error: Unknown method to combine p-values (should be 'min' or 'hmp')")
  }
  setkey(output, baitID, minOE, maxOE)
  
  peakFiles <- peakFiles[,oeID1:=oeID]
  
  cond1names <- character(0)
  if(length(names(chicagoData[[1]]))){ 
    cond1names <- names(chicagoData[[1]])
  }else{
    cond1names <- names(chicagoData)[1]
  }
  
  cond2names <- character(0)
  if(length(names(chicagoData[[2]]))){
    cond2names <- names(chicagoData[[2]])
  }else{
    cond2names <- names(chicagoData)[2]
  } 
  
  #res <- output[get(pcol) < pvcut]
  
  if(!is.null(names(chicagoData[[1]]))){ # replicate-level peakFile
    peakFiles[,cond1mean:= asinh(rowMeans(.SD)), .SDcols=cond1names]
    peakFiles[,cond2mean:= asinh(rowMeans(.SD)), .SDcols=cond2names]
    peakFiles[,delta:=abs(cond1mean-cond2mean)]  
    peakFiles[,cond1mean:=NULL]
    peakFiles[,cond2mean:=NULL]
  }
  else{ # merged peakFile
    peakFiles[,delta:=abs(get(names(chicagoData)[2])-get(names(chicagoData)[1]))] 
  }
  
  outpeak <- foverlaps(peakFiles, output, by.x=c("baitID", "oeID", "oeID1"), 
                       by.y=c("baitID", "minOE", "maxOE"), nomatch =0, mult = "all")
  
  pcol_out = ifelse(method=="min", paste0("min_", pcol), 
                    paste0("hm_", pcol))
  
  if(method=="hmp"){
    outpeak[is.na(get(pcol)) | get(pcol)>1, c(pcol):=1]
  }
  
  # This is to include the columns whose names are listed in names(chicagoData[[n]]) and are unknown a priori
  expr = paste0('list( 
                baitChr = baitChr[1], baitstart = baitstart[1], 
                baitend = baitend[1], baitName = baitName[1],',
                paste0(cond1names, "=", cond1names, "[1]", collapse=","), ',',
                paste0(cond2names, "=", cond2names, "[1]", collapse=","), ',',
                ifelse(method=="min", paste0(pcol_out, "=min(", pcol, ")"), 
                       paste0(pcol_out, "=p.hmp(", pcol, ")")),      
                ', 
                deltaAsinhScore = delta[1],
                regionIDs = paste(regionID, collapse=","),
                log2FoldChanges = paste(log2FoldChange, collapse = ","), 
                ',
                pcol, ' = paste(', pcol, ', collapse = ","),  
                OEranges = paste(OEstart, OEend, sep = "-", collapse= ","))
                ')
  
  #print (parse(text=expr) )
  
  setkey(outpeak, baitID, oeID)
  final <- outpeak[, eval(parse(text=expr)),  by=c("baitID", "oeID")]
  
  final[get(pcol_out) <= pvcut & deltaAsinhScore >= minDeltaAsinhScore]
  
}
