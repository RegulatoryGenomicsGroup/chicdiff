## Settings functions -----------------

defaultchicSettings <- function()
 {
   list(
     inputfiles= NA,
     peakfiles= NA,
     targetRDSorRDAs= NA,
     targetChs= NA,
     rmapfile= NA,
     targetColumns= NA,
     baitmapfile= NA,
     RUexpand= 5L, #number of fragments to expand by in either direction
     score= 5, #Give more informative name (see getRegionUniverse())
     saveRDS = TRUE,
     parallel = FALSE,
     printMemory = TRUE
   )
 }

### This is where we now set the defaults
### Order of priority:
### settings override settings from settingsFile
### both override defchic.settings

setchicExperiment = function(designDir="", targetRDSorRDAs = NA, targetChs = NA, peakfiles = NA, settings=list(), settingsFile=NULL, inputfiles = NA,
                             defchic.settings=defaultchicSettings())
{
  
  if(designDir == "" & identical(settings, list()) & is.null(settingsFile) & identical(defchic.settings, defaultchicSettings()))
  {
    stop("Design not specified. Please specify a design directory or design files.")
  }
  defchic.settings = .updatechicSettings(designDir, targetRDSorRDAs, targetChs, peakfiles, 
                                         settings, settingsFile, inputfiles, defchic.settings)
  
}

.updatechicSettings = function(designDir, targetRDSorRDAs, targetChs, peakfiles, 
                               settings, settingsFile, inputfiles, defchic.settings, updateDesign=FALSE){
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
    defchic.settings[[s]] = modSettings[[s]]
  }
  
  defchic.settings[["targetRDSorRDAs"]] = targetRDSorRDAs
  defchic.settings[["targetChs"]] = targetChs
  
  if(is.na(targetRDSorRDAs) | is.na(targetChs)){
    if(is.na(defchic.settings[["inputfiles"]])){
      defchic.settings[["inputfiles"]] = inputfiles
    }else{
      if (!file.exists(defchic.settings[["inputfiles"]])){
        stop(paste("No config file found at the specified location", defchic.settings[["inputfiles"]]))
      }
    }
  } else{
    if(!(all.equal(names(targetRDSorRDAs), names(targetChs)))){
      stop("Conditions for the RDS/RDA files and chinputs must be the same")
    }
    temp <- unlist(targetRDSorRDAs)
    temp_2 <- unlist(targetChs)
    if(!(length(temp) == length(temp_2))){
      stop("Must provide the same number of RDS/RDA files as chinputs")
    }
  }
  
  if(!is.na(inputfiles)){
    input_lists <- .makeTargetFilesList(inputfiles)
    targetRDSorRDAs <- input_lists[[1]]
    targetChs <- input_lists[[2]]
  }

  defchic.settings[["peakfiles"]] = peakfiles

  if(is.na(defchic.settings[["peakfiles"]])){
    stop("No peak files provided")
  } else if(!file.exists(defchic.settings[["peakfiles"]])){
      stop(paste("No peak files found at the specified location", defchic.settings[["peakfiles"]]))
  }
  
  targetColumns <- .getTargetColumns(targetRDSorRDAs)
  defchic.settings[["targetColumns"]] = targetColumns
  peakfile_columns <- colnames(.multimerge(peakfiles))
  
  if(!all(targetColumns %in% peakfile_columns)){
    stop("All specified columns must be present in the peak files")  
  }

  if(is.na(defchic.settings[["baitmapfile"]]) | (updateDesign & !is.null(designDir))){
    defchic.settings[["baitmapfile"]] = .locateFile("<baitmapfile>.baitmap", designDir, "\\.baitmap$")
  }else{
    if (!file.exists(defchic.settings[["baitmapfile"]])){
      stop(paste("No baitmap file found at the specified location", defchic.settings[["baitmapfile"]]))
    }
  }
  
  if(is.na(defchic.settings[["rmapfile"]]) | (updateDesign & !is.null(designDir))){
    defchic.settings[["rmapfile"]] = .locateFile("<rmapfile>.rmap", designDir, "\\.rmap$")
  }else{
    if (!file.exists(defchic.settings[["rmapfile"]])){
      stop(paste("No rmap file found at the specified location", defchic.settings[["rmapfile"]]))
    }
  }
  
  

  message("Checking the design files...")

  baitmapfile = Chicago:::.readBaitmap(defchic.settings)
  if(ncol(baitmapfile)<max(c(defchic.settings[["baitmapFragIDcol"]], defchic.settings[["baitmapGeneIDcol"]]))){
    stop("There are fewer columns in the baitmapfile than expected. Check that this file lists the genomic coordinates as well as both the IDs and names for each baited fragment,
          and that the corresponding columns are specified in baitmapFragIDcol and baitmapGeneIDcol, respectively.")
  }
  rmapfile = Chicago:::.readRmap(defchic.settings)
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

  defchic.settings
  }

# ----------

.getTargetColumns <- function(targetRDSorRDAs){
  
  conditions <- names(targetRDSorRDAs)
  targetColumns <- vector("list", length(conditions))
  for(i in seq_along(conditions)){
    temp <- conditions[[i]]
    targetColumns[[i]] <- names(targetRDSorRDAs[[temp]])
  }
  targetColumns <- unlist(targetColumns)
  if(is.null(targetColumns)){
    targetColumns <- conditions
  } 
  okay = (length(targetColumns) == length(conditions) | length(targetColumns) == length(unlist(targetRDSorRDAs)))
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

.multimerge <- function(peakFiles)
{
  templist <- lapply(peakFiles, fread)
  peakMatrix <- Reduce(merge, templist)
  peakMatrix
}

readAndFilterPeakMatrix <- function(peakFiles, targetColumns, conditions, score){
## Essentially reads in the peak matrix, taking the target columns if specified (else just takes everything) and filters for rows where at least
##one score is > 5 (or whatever is specified) and filters out trans interactions. 
      if(length(peakFiles > 1)){
        x <- .multimerge(peakFiles)
      } else {
        x <- fread(peakFiles)
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
          cols <- colnames(x)[colnames(x) %in% names(targetRDSorRDAs[[temp]])]
          
          sel2 <- sel2 & rowSums(x[,!is.na(.SD[,cols, with = FALSE])]) >= 2 ##Remove any rows where there isn't at least 2 replicates for each condition with non-NA values
        } 
        
        x <- x[sel2]
      }
  
  x <- x[!is.na(dist),] ## <- FILTER OUT TRANS INTERACTIONS - Assumed in (*)
  
  filtered_baits <- all_baits[which(!(all_baits %in% unique(x$baitID)))]
  filtered_baits <- list(filtered_baits)
  fwrite(filtered_baits, "./filtered_baits.txt")
  
  x
}

.printMemoryFunction <- function(printMemory){
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

chicdiffPipeline <- function(defchic.settings, outprefix=NULL, printMemory=TRUE)
{
  
  message("\n*** Running getRegionUniverse\n")
  RU <- getRegionUniverse(defchic.settings)  
  
  message("\n*** Running getControlRegionUniverse\n")
  RUcontrol <- getControlRegionUniverse3(defchic.settings, RU)

  message("\n*** Running getFullRegionData\n")
  FullRegionData <- getFullRegionData(defchic.settings, RU, RUcontrol, suffix = "")

  message("\n*** Running DESeq2Wrap for FullRegion\n")
  DESeqOut <- DESeq2Wrap(defchic.settings, RU, FullRegionData[[1]])

  message("\n*** Running DESeq2Wrap for FullControlRegion\n")
  DESeqOutControl <- DESeq2Wrap(defchic.settings, RUcontrol, FullRegionData[[2]], suffix = "Control") 

  message("\n*** Running IHWcorrection\n")
  output <- IHWcorrection(defchic.settings, DESeqOut, FullRegionData[[1]], DESeqOutControl, FullRegionData[[2]], countput = FullRegionData[[3]]) 

  print(defchic.settings)
  print(sessionInfo())
  
  output

}

#-------------------------------Pipeline functions----------------------------------#

####------------------------------- getRegionUniverse-------------------------------#

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

getRegionUniverse <- function(defchic.settings, suffix = ""){
 
  RUexpand = defchic.settings[["RUexpand"]]
  rmapfile = defchic.settings[["rmapfile"]]
  peakFiles = defchic.settings[["peakfiles"]]
  targetColumns = defchic.settings[["targetColumns"]]
  targetRDSorRDAs = defchic.settings[["targetRDSorRDAs"]]
  score = defchic.settings[["score"]]
  saveRDS = defchic.settings[["saveRDS"]]
  
  conditions <- names(targetRDSorRDAs)
  
  x <- readAndFilterPeakMatrix(peakFiles = peakFiles, score = score, targetColumns = targetColumns, conditions = conditions)
  
  ## Expand the "point universe" to get "region universe" by "window" mode------------------
 

    ##Just expand calls by s in each direction
    
    baits <- unique(x$baitID)
    #baits <- head(baits, 5000) ##Reduce to small subset for now:
    x <- x[baitID %in% baits,]
    
    RU.DT <- x[,c("baitID", "oeID"),with=FALSE]
    RU.DT$regionID <- 1:nrow(RU.DT)
    RU.DT <- RU.DT[,
                   list(otherEndID = expandAvoidBait(baitID, oeID, RUexpand)), 
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
    saveRDS(RU.DT, paste0("./RegionUniverse", suffix, ".Rds"))
  }
  RU.DT
}

#------------------------------- getControlRegionUniverse----------------------------------#

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
                  if(length(which(temp))==1){
                    unlist(sample(list(vIDs[temp])))
                  } else {
                    sample(vIDs[temp], size=1)
                  }
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

getControlRegionUniverse <- function(defchic.settings, RU, suffix = "")
{
  
  rmapfile = defchic.settings[["rmapfile"]]
  baitmapfile = defchic.settings[["baitmapfile"]]
  saveRDS = defchic.settings[["saveRDS"]]
  
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

#-----------------------Functions needed for getControlRegionUniverse2 and 3 ----------------------------#

giveRandomBait <- function(chrom, bmap, regionID){
  vIDs <- bmap[chr == chrom, ID]
  unlist(lapply(regionID, function(x, vIDs){sample(vIDs, 1)}, vIDs))
}

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


#-------------------------------------------getControlRegionUniverse2----------------------------------------#

getControlRegionUniverse2 <- function(defchic.settings, RU){
  
  RUexpand <- defchic.settings[["RUexpand"]]
  saveRDS <- defchic.settings[["saveRDS"]]
  bmap <- fread(defchic.settings[["baitmapfile"]])
  rmap <- fread(defchic.settings[["rmapfile"]])
  colnames(rmap) <- c("chr", "start", "end", "ID")
  colnames(bmap) <- c("chr", "start", "end", "ID", "baitname")
  
  RUshorter <- merge(RU, rmap[,c("chr", "ID")], by.x = "baitID", by.y = "ID")
  max_contacts <- RUshorter[, list(max_contact = max(abs(baitID - otherEndID))),by = "chr"]
  RUshorter <- unique(RUshorter[,c("chr", "baitID", "regionID")])
  
  
  for(ch in unique(RUshorter[,chr])){
    
    std <- max_contacts[chr == ch, max_contact]/3
    max <- rmap[chr == ch, max(ID)]
    min <- rmap[chr == ch, min(ID)]
    
    RUshorter[chr == ch, rand_bait := giveRandomBait(ch, bmap, regionID)]
    RUshorter[chr == ch, seed := giveManySeeds(rand_bait, min = min, max = max, std = std)]
  }
  
  RUcontrol <- RUshorter[,c("rand_bait", "seed", "regionID")]
  colnames(RUcontrol) <- c("baitID", "oeID", "regionID")
  
  RUcontrol <- RUcontrol[,
                         list(otherEndID = expandAvoidBait(baitID, oeID, RUexpand)), 
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
    saveRDS(RUcontrol, "ControlRegionUniverse.Rds")
  }
  
  return(RUcontrol)
}

#-------------------------------------------getControlRegionUniverse3----------------------------------------#

getControlRegionUniverse3 <- function(defchic.settings, RU){
  
  RUexpand <- defchic.settings[["RUexpand"]]
  saveRDS <- defchic.settings[["saveRDS"]]
  bmap <- fread(defchic.settings[["baitmapfile"]])
  rmap <- fread(defchic.settings[["rmapfile"]])
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
                         list(otherEndID = expandAvoidBait(baitID, oeID, RUexpand)), 
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
    saveRDS(RUcontrol, "ControlRegionUniverse.Rds")
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
      message("Cannot read file. Must be an RDA or RDA file.")
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

chicestimateDistFun <- function (x, settings=defaultSettings()) {
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

getFullRegionData1 <- function(defchic.settings, RU, is_control = FALSE, suffix = ""){
  
  targetChs = defchic.settings[["targetChs"]]
  targetRDSorRDAs = defchic.settings[["targetRDSorRDAs"]]
  rmapfile = defchic.settings[["rmapfile"]]
  saveRDS = defchic.settings[["saveRDS"]]
  
  targetRDSorRDAFiles <- unlist(targetRDSorRDAs)
  targetChFiles <- unlist(targetChs)
  
  ## 03getInteractionParameters
  
  RU_IntParams <- copy(RU)
  RU_IntParams[,regionID:=NULL]
  RU_IntParams <- unique(RU_IntParams)
  
  ##prep output
  outData <- vector("list", length(targetRDSorRDAFiles))
  dispersions <- numeric(length(targetRDSorRDAFiles))
  tempForCounts <- vector("list", length(targetRDSorRDAFiles))
  
  if(!is_control){
    
    conditions <- rep(names(targetRDSorRDAs), sapply(targetRDSorRDAs, length))
    countput <- lapply(unique(conditions), assign, 
                       value = vector("list", length(targetRDSorRDAFiles)))
    names(countput) <- unique(conditions)
  }
  
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
    
    message("Finished reading")
    print(mem_used())
    
    ##1b) collect the Bmean, Tmean information that is present
    
    setkey(x, baitID, otherEndID)
    setkey(RU_IntParams, baitID, otherEndID)
    out <- x[RU_IntParams, c("baitID", "otherEndID","distSign","Bmean", "Tmean", "score"), with=FALSE]
    ##(note that score=NA occurs when N=0 - hence replace score=0)
    
    message("Collected Bmean and Tmean")
    print(mem_used())
    
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
    
    message("Recalculated distSign")
    print(mem_used())
    
    ##1d) bait annotation: collect s_js, tblb
    ##Note that s_j can be NA (for baits that were filtered out!)
    temp <- x[,list(s_j=s_j[1], tblb=tblb[1]),by="baitID"]
    setkey(temp, baitID)
    setkey(out, baitID)
    out <- merge(out, temp, all.x=TRUE)
    
    message("Collected s_js and tblb")
    print(mem_used())
    
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
    
    message("Collected s_is and tlb")
    print(mem_used())
    
    ##1f) reconstruct Tmean information
    setkey(x, tlb, tblb)
    setkey(out, tlb, tblb)
    temp <- x[,Tmean[1],by=c("tblb", "tlb")]
    setnames(temp, "V1", "Tmean")
    out[,Tmean:=NULL]
    out <- merge(out, temp, all.x=TRUE)
    
    message("Found Tmean")
    print(mem_used())
    
    ##1g) impute some of the missing Tmeans
    ##    missing tlb suggests that the other end was rarely observed, so assume lowest
    ##    Tmean corresponding to the appropriate tblb.
    templowest <- temp[,c(Tmean = min(Tmean)),by="tblb"]
    tempDictionary <- templowest$V1
    names(tempDictionary) <- templowest$tblb
    out[is.na(tlb) & !is.na(tblb), Tmean := tempDictionary[tblb]]
    
    message("Imputed missing Tmean")
    print(mem_used())
    
    ##2) Reconstruct the distance function from distbin info
    distFunParams <- chicestimateDistFun(x)
    
    ##3) Reconstruct Bmean information (But!! When s_j = NA,
    ##   ensure that Bmean comes out as NA too.)
    out <- Chicago:::.estimateBMean(out, distFunParams)
    out[is.na(s_j),Bmean:=NA]
    setkey(out, baitID, otherEndID)
    
    outData[[i]] <- out
    
    message("Estimated Bmean")
    print(mem_used())
    
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
      
      message("Assigning conditions")
      for(n in 1:length(countput)){
        if(conditions[[i]] == names(countput)[[n]]){
          countput[[n]][[i]] <- x
        }
      }
      
      
    }
    
    #This just keeps all of the loaded RDAs if we still need them for the count data because we don't have chinputs
    if(is.null(targetChs)){  #my lame place holder conditional for checking whether chinputs have been provided. Easily changed. 
      x <- x[, c("baitID", "otherEndID", "N"), with = FALSE]
      columnNames <- paste0("N.", names(targetRDSorRDAFiles))
      setnames(x, old="N", new=columnNames[i])
      setkey(x, baitID, otherEndID)
      tempForCounts[[i]] <- x
      rm(x)
    }
    
    message("Saved counts")
    print(mem_used())
    
  }
  
  if(!is_control){
    for(i in 1:length(countput)){
      countput[[i]] <- rbindlist(countput[[i]])
      gc()
      countput[[i]] <- countput[[i]][,list(Nav = mean(N), Bav = mean(Bmean), score = max(score), midpoint = midpoint[1]), by = c("baitID", "otherEndID")]
      gc()
      countput[[i]][, condition := rep(names(countput)[[i]], nrow(countput[[i]]))]
    }
    
    countput <- rbindlist(countput)
    setnames(countput, "midpoint", "oeID_mid")
    
    saveRDS(countput, "countput.Rds")
  }
  
  message("between (1)")
  print(mem_used())
  
  
  saveRDS(outData, paste0("03interactionsParameters", suffix, ".Rds"))   ##REMEMBER TO REMOVE THESE TWO SAVES
  saveRDS(dispersions, paste0("03dispersionParamters", suffix, ".Rds"))
  
  message("out of iterative reading function")
  print(mem_used())
  
  baits <- sort(unique(RU$baitID)) ##baits to get information from
  
  if(is.null(targetChs)){
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
      setkey(countData[[i]], baitID)
      
      message("countData info collected")
      print(mem_used())
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
      
      message("countData info merged")
      print(mem_used())
    }
    
    message("no chinput case finished")
    print(mem_used())
  }
  
  #Recall that outData and dispersions are our useful outputs from this part above 
  
  ## 02getCounts.R
  
  ##Given a region universe, collect the number of reads in each region from the chinputs.
  
  ##Inputs and parameters:
  
  ##requires targetChs
  
  ##2) Collect all of the read count information ----------------
  if(!is.null(targetChs)){
    
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
      
      message("Finished reading of chinputs")
      print(mem_used())
      gc()
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
      
      message("count the number of reads in each replicate for CountOut")
      print(mem_used())
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
  
  message("collected distances for CountOut")
  print(mem_used())
  
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
    
    message("count + interaction parameter data together")
    print(mem_used())
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
  
  message("recast")
  print(mem_used())
  
  if(saveRDS == TRUE){
    saveRDS(recast, paste0("./FullRegionData", suffix, ".Rds"))
    
    message("saved FullRegionData")
    print(mem_used())
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


#---------------------------First the version which acts on (non)control in parallel-------------#

getFullRegionData2 <- function(defchic.settings, RU, RUcontrol, suffix = ""){
  
  targetChs = defchic.settings[["targetChs"]]
  targetRDSorRDAs = defchic.settings[["targetRDSorRDAs"]]
  targetColumns = defchic.settings[["targetColumns"]]
  rmapfile = defchic.settings[["rmapfile"]]
  saveRDS = defchic.settings[["saveRDS"]]
     
  targetRDSorRDAFiles <- unlist(targetRDSorRDAs)
  targetChFiles <- unlist(targetChs)
  
  ## a.k.a. 03getInteractionParameters.R  
    
  ##Given a region universe, collect Bmean + Tmean  
  ##Inputs and parameters:
  ##requires targetRDAs

  RU_IntParams <- copy(RU) 
  RU_IntParams[,regionID:=NULL]
  RU_IntParams <- unique(RU_IntParams)

  RUcontrol_IntParams <- copy(RUcontrol)
  RUcontrol_IntParams[,regionID:=NULL]
  RUcontrol_IntParams <- unique(RUcontrol_IntParams)
  #Past this point, any mention of RU or RUcontrol in the 03 section of this function is a mistake. Should be RU_IntParams etc. (except for getCount)
    
  ##prep output
  outData <- vector("list", length(targetRDSorRDAFiles))
  outDataControl <- vector("list", length(targetRDSorRDAFiles)) #New to collect control data
  dispersions <- numeric(length(targetRDSorRDAFiles)) #only need one of these
  tempForCounts <- vector("list", length(targetRDSorRDAFiles))
  
  conditions <- rep(names(targetRDSorRDAs), sapply(targetRDSorRDAs, length))
  countput <- lapply(unique(conditions), assign, 
                             value = vector("list", length(targetRDSorRDAFiles)))
  names(countput) <- unique(conditions)
    
  ##1) Read in each file...
  for(i in 1:length(targetRDSorRDAFiles))
  { 
    message("File ",i, " of ", length(targetRDSorRDAFiles))
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

    message("Finished reading")
    print(mem_used())

    ##1b) collect the Bmean, Tmean information that is present
    setkey(x, baitID, otherEndID)
    setkey(RU_IntParams, baitID, otherEndID)
    setkey(RUcontrol_IntParams, baitID, otherEndID)
    out <- x[RU_IntParams, c("baitID", "otherEndID","distSign","Bmean", "Tmean", "score"), with=FALSE]
    outControl <- x[RUcontrol_IntParams, c("baitID", "otherEndID","distSign","Bmean", "Tmean", "score"), with=FALSE]

    message("Collected Bnmean and Tmean")
    print(mem_used())
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

      message("recalculated distSign for out")
      print(mem_used())
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

      message("recalculated distSign for outControl")
      print(mem_used())
    }

    ##1d) bait annotation: collect s_js, tblb
    ##Note that s_j can be NA (for baits that were filtered out!)
    temp <- x[,list(s_j=s_j[1], tblb=tblb[1]),by="baitID"]
    setkey(temp, baitID)
    setkey(out, baitID)
    setkey(outControl, baitID)
    out <- merge(out, temp, all.x=TRUE)
    outControl <- merge(outControl, temp, all.x = TRUE)

    message("recalculated distSign for outControl")
    print(mem_used())

    ##1e) OE annotation: collect s_is, tlb
    ##s_i not allowed to be NA - just assume 1.
    ##(Could be an issue if too many OEs were filtered out => s_i should be 0)
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

    message("collected s_is and tlb")
    print(mem_used())

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

    message("collected s_is and tlb")
    print(mem_used())

    ##1g) impute some of the missing Tmeans
    ##    missing tlb suggests that the other end was rarely observed, so assume lowest
    ##    Tmean corresponding to the appropriate tblb.
    templowest <- temp[,c(Tmean = min(Tmean)),by="tblb"]
    tempDictionary <- templowest$V1
    names(tempDictionary) <- templowest$tblb
    out[is.na(tlb) & !is.na(tblb), Tmean := tempDictionary[tblb]]
    outControl[is.na(tlb) & !is.na(tblb), Tmean := tempDictionary[tblb]]

    message("imputed some of the missing Tmeans")
    print(mem_used())
    
    ##2) Reconstruct the distance function from distbin info
    distFunParams <- chicestimateDistFun(x)

    message("Reconstructed the distance function")
    print(mem_used())
    
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

    message("Reconstructed Bmean")
    print(mem_used())
    
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
    
    message("Assigning conditions")
    for(n in 1:length(countput)){
      if(conditions[[i]] == names(countput)[[n]]){
        countput[[n]][[i]] <- x
      }
    }
    
    if(is.null(targetChs)){ 
      x <- x[, c("baitID", "otherEndID", "N"), with = FALSE]
      columnNames <- paste0("N.", names(targetRDSorRDAFiles))
      setnames(x, old="N", new=columnNames[i])
      setkey(x, baitID, otherEndID)
      tempForCounts[[i]] <- x
      rm(x)

      message("saved counts")
      print(mem_used())
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
  
  saveRDS(countput, "countput.Rds")
  
  message("between (1)")
  print(mem_used())

  # outData, outDataControl and dispersions are our useful outputs from this part above  

  saveRDS(outData, paste0("03interactionsParameters", suffix, ".Rds"))   ##REMEMBER TO REMOVE THESE TWO SAVES
  saveRDS(dispersions, paste0("03dispersionParamters", suffix, ".Rds"))

  message("out of iterative reading function")
  print(mem_used())

  baits <- sort(unique(RU$baitID)) 
  baitsControl <- sort(unique(RUcontrol$baitID))

  ## 02getCounts.R  
  ## Given a region universe, collect the number of reads in each region

  ##2-Rds) If Chinputs are not given: Collect all of the read count information from Rds

  if(is.null(targetChs)){
    gc()
    mergedFiles <- Reduce(merge, tempForCounts)
    countData <- vector("list", length(targetRDSorRDAFiles))
    countDataControl <- vector("list", length(targetRDSorRDAFiles))

    for(i in seq_along(countData)){
      countData[[i]] <- mergedFiles[, c("baitID", "otherEndID", columnNames[i]), with = FALSE]
      setkey(countData[[i]], baitID)
      countData[[i]] <- countData[[i]][J(baits),] #issue that rds may not contain all baitid otherendid pairs which are in RU
      setkey(countData[[i]], baitID)

      message("countData info collected")
      print(mem_used())
    }

    for(i in seq_along(countDataControl)){
      countDataControl[[i]] <- mergedFiles[, c("baitID", "otherEndID", columnNames[i]), with = FALSE]
      setkey(countDataControl[[i]], baitID)
      countDataControl[[i]] <- countDataControl[[i]][J(baitsControl),] #issue that rds may not contain all baitid otherendid pairs which are in RU
      setkey(countDataControl[[i]], baitID)

      message("countDatacontrol info collected")
      print(mem_used())
    }
    
    CountOut <- RU
    CountOutControl <- RUcontrol
    setkey(CountOut, baitID, otherEndID)
    setkey(CountOutControl, baitID, otherEndID)

    message("between (2)")
    print(mem_used())
    
    for(i in seq_along(countData))
    {
      temp <- countData[[i]] #[,c("baitID","otherEndID","N"), with=FALSE]
      setnames(temp, old = columnNames[i], new = "N")
      setkey(temp, baitID, otherEndID)
      CountOut <- merge(CountOut, temp, all.x=TRUE)
      CountOut[is.na(N), N:=0]
      setnames(CountOut, old="N", new=columnNames[i])

      message("countData info merged")
      print(mem_used())
    }

    for(i in seq_along(countDataControl))
    {
      temp <- countDataControl[[i]] #[,c("baitID","otherEndID","N"), with=FALSE]
      setnames(temp, old = columnNames[i], new = "N")
      setkey(temp, baitID, otherEndID)
      CountOutControl <- merge(CountOutControl, temp, all.x=TRUE)
      CountOutControl[is.na(N), N:=0]
      setnames(CountOutControl, old="N", new=columnNames[i])

      message("countDatacontrol info merged")
      print(mem_used())
    }

    message("no chinput case finished")
    print(mem_used())
  }

  ##2-Chi) If Chinputs are given: Collect all of the read count information ----------------
  if(!is.null(targetChs)){
    
    countData <- vector("list", length(targetChFiles))
    countDataControl <- vector("list", length(targetChFiles))
    gc()

    for(i in 1:length(countData))
  {
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

    message("Finished reading of chinputs")
    print(mem_used())
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

    message("count the number of reads in each replicate for CountOut")
    print(mem_used())
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

    message("count the number of reads in each replicate for CountOutControl")
    print(mem_used())
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
  
  CountOut <- merge(CountOut, rmap, by.x="otherEndID", by.y="fragID")
  CountOut <- merge(CountOut, rmap, by.x="baitID", by.y="fragID")
  setkey(CountOut, baitID, otherEndID)

  CountOut[, distSign :=
             ifelse(chr.x == chr.y,
                    midpoint.x - midpoint.y,
                    NA)]
  CountOut[, c("chr.x","midpoint.x","chr.y","midpoint.y"):=NULL]

  message("collected distances for CountOut")
  print(mem_used())
  
  #Just repeating the above code for the control
  CountOutControl <- merge(CountOutControl, rmap, by.x="otherEndID", by.y="fragID")  
  CountOutControl <- merge(CountOutControl, rmap, by.x="baitID", by.y="fragID")
  setkey(CountOutControl, baitID, otherEndID)
  
  CountOutControl[, distSign :=
                    ifelse(chr.x == chr.y,
                           midpoint.x - midpoint.y,
                           NA)]
  CountOutControl[, c("chr.x","midpoint.x","chr.y","midpoint.y"):=NULL]

  message("collected distances for CountOutcontrol")
  print(mem_used())

  # 04mergeData.R
  
  ##Synthesise the count + interaction parameter data together
  
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

    message("count + interaction parameter data together")
    print(mem_used())
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
    rep(names(targetRDSorRDAs), sapply(targetRDSorRDAs, length))
    , each=nrow(CountOut))

  setkey(recast, regionID)

  message("recast")
  print(mem_used())
  
  # Again - just a copy and paste for the control here 
  recastControl = melt(CountOutControl, measure = sels, value.name = myCols,
                       variable.name="sample")
  recastControl$sample <- names(targetRDSorRDAFiles)[recastControl$sample]
  
  recastControl$condition <- rep(
    rep(names(targetRDSorRDAs), sapply(targetRDSorRDAs, length))
    , each=nrow(CountOutControl))
  
  setkey(recastControl, regionID)

  message("recastControl")
  print(mem_used())
  
  if(saveRDS == TRUE){
    saveRDS(recast, paste0("./FullRegionData", suffix, ".Rds"))
    saveRDS(recastControl, paste0("./FullControlRegionData", suffix, ".Rds"))

    message("saved FullRegionData and FullControlRegionData")
    print(mem_used())
  }

  message("reported FullRegionDatalist")
  print(mem_used())
  return(list(recast, recastControl, countput))
}

#-------------------------------The full getFullRegionDataWrapper----------------------------------------#

getFullRegionData <- function(defchic.settings, RU, RUcontrol, suffix = ""){
  
  parallel = defchic.settings[["parallel"]]
  
  if(parallel){
    
    FullRegionData <- getFullRegionData2(defchic.settings = defchic.settings, RU = RU, RUcontrol = RUcontrol, suffix = "")
    
  } else{
    
    FullRegionData <- getFullRegionData1(defchic.settings = defchic.settings, RU = RU, is_control = FALSE, suffix = "")
    FullRegionData[[2]] <- getFullRegionData1(defchic.settings = defchic.settings, RU = RUcontrol, is_control = TRUE, suffix = "")
    
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

DESeq2Wrap <- function(defchic.settings, RU, FullRegionData, suffix = ""){
  ##By default, DESeq2 uses sample-specific scaling factors.
  
  rmapfile = defchic.settings[["rmapfile"]]
  saveRDS = defchic.settings[["saveRDS"]]
  
  ##Input data:
  
  fragData <- copy(FullRegionData) ## As otherwise it is altered by reference below and the IHW part breaks
  setkey(fragData, otherEndID)

  message("saved a copy of FullRegionData")
  print(mem_used())
  
  fragData.impute <- copy(fragData)
  fragData.impute[,s_j.impute := geoMean(s_j, na.rm = TRUE), by=c("baitID", "otherEndID")]
  fragData.impute[is.na(s_j), s_j:=s_j.impute]
  fragData.impute[,FullMean.impute := geoMean(FullMean, na.rm = TRUE), by=c("baitID", "otherEndID")]
  fragData.impute[is.na(FullMean), FullMean:=FullMean.impute]

  ##This is calculated but not used? So I have commented it out
  #fragData[is.na(s_j),N := 0]
  #fragData[is.na(s_j),s_j := 0]

  ##NOTE: Standard normalization factors from DESeq2 used for certain sites later.
  
  ##construct an appropriate DESeq object
  
  # regionData <- fragData[,list(
  #   N=sum(N),
  #   avDist=(min(distSign)+max(distSign))/2,
  #   Bmean=sum(Bmean),
  #   Tmean=sum(Tmean),
  #   FullMean=sum(FullMean),
  #   s_j=s_j[1]
  # ),by=c("baitID", "regionID", "sample")]
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
  rownames(regionDataMatrix) <- unique(regionData.impute$regionID)
  colnames(regionDataMatrix) <- unique(regionData.impute$sample)
  colData <- data.frame(condition = fragData$condition[1:ncol(regionDataMatrix)])
  dds <- DESeqDataSetFromMatrix(countData = regionDataMatrix,
                                colData = colData,
                                design = ~ condition)
  
  dds.nullModel <- estimateSizeFactors(dds)
  dds.M3 <- copy(dds.nullModel)
  
  nullSizeFactors <- sizeFactors(dds.nullModel)

  message("construct an appropriate DESeq object")
  print(mem_used())
  
  ##model 1) standard DESeq2 model
  dds.nullModel <- estimateDispersions(dds.nullModel)
  dds.nullModel <- nbinomWaldTest(dds.nullModel)

  ##model 2 is defunct
  ##model 3) use fullMean scaling factors
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
  
  results <- as.data.table(as.data.frame(results(dds.M3)), keep.rownames = "id")
  results$id <- as.integer(results$id)
  setkey(results, id)
  
  out <- cbind(
    results,
    annoData
  )
  
  out[,id := NULL]

  message("outModel")
  print(mem_used())
  
  results.null <- as.data.table(as.data.frame(results(dds.nullModel)), keep.rownames = "id")
  results.null$id <- as.integer(results.null$id)
  setkey(results.null, id)
  
  out.nullModel <- cbind(
    results.null,
    annoData
  )
  
  out.nullModel[,id := NULL]

  message("out.nullModel")
  print(mem_used())

  
  if(saveRDS == TRUE){
    saveRDS(out, paste0("./05outPilot", suffix, ".Rds"))
    saveRDS(out.nullModel, paste0("./05outStandardDESeq", suffix, ".Rds"))
  }

  message("savedRDS")
  print(mem_used())
  out 
}

#-----------------------------IHW and plotting functions-------------------------------#

plotdiffBaits <- function(output, countput, baitmapfile, n = 3, baits = NULL, plotBaitNames = TRUE, plotBaitIDs = TRUE, plevel1 = 5, plevel2 = 3, xlim=c(-1e6,1e6), bgCol = "black", lev1Col = "red", lev2Col = "blue", ...){
  
  if(class(output) == "character"){
    if(file.exists(output)){
      out <- readRDS(output)
    } else{
      "No output file found at specified location"
    }
  } else{
    out <- copy(output)
  }
  
  if(class(countput) == "character"){
    if(file.exists(countput)){
      countput_coord <- readRDS(countput)
    } else{
      "No countput file found at specified location"
    }
  } else{
    countput_coord <- copy(countput)
  }
  
  if(class(baitmapfile) == "character"){
    if(file.exists(baitmapfile)){
      bmap <- fread(baitmapfile)
    } else{
      "No baitmap file found at specified location"
    }
  } else{
    bmap <- copy(baitmapfile)
  }
  
  conditions <- unique(countput_coord[,condition])
  
  out[, names(out)[!(names(out) %in% c("OEstart", "OEend","baitID", "weighted_pvalue"))] := NULL]
  
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
  
  scale_data <- data.frame(value = 1L, name = factor(c("> 0.05", "0.05 - 0.005", "0.005 - 0.0005", "< 0.0005"), levels = c("< 0.0005", "0.005 - 0.0005", "0.05 - 0.005", "> 0.05")))
  scalebar <- ggplot(scale_data, aes(x = value, y = name)) +
    geom_tile(aes(fill = name)) +
    scale_fill_manual(name = "p-values", values = c("> 0.05" = "#2b8cbe", "0.05 - 0.005" = "#fee8c8", "0.005 - 0.0005" = "#fdbb84", "< 0.0005" = "#e34a33"))
  
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
      geom_point(aes(colour=big_score, alpha = 0.4)) + 
      scale_colour_manual(values = c(bgCol, lev2Col, lev1Col)) +
      geom_line(aes(oeID_mid, Bav)) +
      xlim(xlimit[[1]], xlimit[[2]]) +
      ylim(0, ylim_max) +
      theme(legend.position = "none",
            axis.title.x = element_blank()) +
      ylab("N")
    
    if(plotBaitNames & plotBaitIDs){
      p1 <- p1 + ggtitle(paste0(chroms_bait[i][,name], " (", chroms_bait[i][,baitID], ")"))
    }else if(plotBaitNames){
      p1 <- p1 + ggtitle(chroms_bait[i][,name])
    } else if(plotBaitIDs){
      p1 <- p1 + ggtitle(chroms_bait[i][,baitID])
    }
    
    p2 <- ggplot(count_temp[condition == conditions[[2]]], aes(oeID_mid, Nav)) +
      geom_point(aes(colour=big_score, alpha = 0.4)) + 
      scale_colour_manual(values = c(bgCol, lev2Col, lev1Col)) +
      geom_line(aes(oeID_mid, Bav)) +
      #     xlim(xlim[[1]], xlim[[2]]) +
      scale_y_reverse(limits = c(ylim_max, 0)) +
      theme(legend.position = "none") + 
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
    merged_dat[,minuslogpvalue := -log10(weighted_pvalue)]
    
    merged_dat[minuslogpvalue <= -log(0.05), pvalue_param := "blue"]
    merged_dat[minuslogpvalue <= -log10(0.005) & minuslogpvalue > -log10(0.05), pvalue_param := "light red"]
    merged_dat[minuslogpvalue <= -log10(0.0005) & minuslogpvalue > -log10(0.005), pvalue_param := "mid red"]
    merged_dat[is.na(pvalue_param), pvalue_param := "red"]
    merged_dat[,pvalue_param := factor(pvalue_param, levels = c("blue", "light red", "mid red", "red"))]
    
    #     merged_dat[minuslogpvalue < 0, minuslogpvalue := 0]
    #     merged_dat[minuslogpvalue > 3, minuslogpvalue := 3]
    #     
    
    pvalue_plot <- ggplot(merged_dat[OEstart >= xlimit[[1]] & OEstart <= xlimit[[2]]]) + 
      geom_rect(aes(xmin = OEstart, xmax = OEend,
                    ymin = bin, ymax = bin + 0.9,
                    fill = pvalue_param), colour = "white", size = 0.3) +
      scale_fill_manual(values = c("blue" = "#2b8cbe", "light red" = "#fee8c8", "mid red" = "#fdbb84", "red" = "#e34a33")) +
      xlim(xlimit[[1]], xlimit[[2]]) +
      theme_void() +
      theme(legend.position = "none")
    
    plot_list[[i]] <- plot_grid(p1, pvalue_plot, p2, align = "v", nrow = 3, rel_heights = c(1/2, 1/8, 1/2))
  }
  
  chicdiffPlots <- plot_grid(plotlist = plot_list, align = "h")
  scalebar_grid <- plot_grid(get_legend(scalebar))
  plot_grid(chicdiffPlots, scalebar_grid, rel_widths = c(10, 1), ...)
}

#---------------------------------------------------#

IHWcorrection <- function(defchic.settings, DESeqOut, FullRegionData, DESeqOutControl, FullControlRegionData,
                          countput, DiagPlot = TRUE, diffbaitPlot = TRUE, suffix = ""){
  
  saveRDS = defchic.settings[["saveRDS"]]
  baitmapfile = defchic.settings[["baitmapfile"]]
  
  out <- DESeqOut ##DESeqOut and DESeqOutControl
  RU.recast <- FullRegionData ##FullRegionData and FullControlRegionData
  RU.distances <- RU.recast[,list(avDist = mean(distSign)),by="regionID"]
  
  out$avDist <- RU.distances$avDist
  
  ##Comparison against genuinely uniform p-vals
  out$uniform <- runif(nrow(out))
  out$shuff <- sample(out$pvalue)

  message("Comparison against p-vals for out")
  print(mem_used())
  
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

  message("Comparison against p-vals for outcontrol")
  print(mem_used())
  
  ##Convert
  as.data.frame(out.control) #This line is necessary even if the one for out may not be
  
  ##Train weights on the control sample
  ihwRes <- ihw(pvalue ~ abs(avDist),  data = out.control, alpha = 0.05)

  message("Trained weights on the control sample")
  print(mem_used())
  
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
  print(mem_used())
  
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

  message("applied to test data")
  print(mem_used())
  
  if (diffbaitPlot == TRUE){
    sel <- order(out$weighted_pvalue)
    baits <- sample(head(unique(out[sel]$baitID), 100), 4)
    plotdiffBaits(output = out, countput = countput, baitmapfile = baitmapfile, baits = baits)
    ggsave("diffbaitPlot.pdf", device = "pdf", path = "./")
  }
  
  if(saveRDS == TRUE){
    saveRDS(out, paste0("06Weighted", suffix, ".Rds"))
  }
  return(out)
}