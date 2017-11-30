##Read in the Region Universe, then scramble it
RU <- readRDS(file.path(sub("/control", "", filePath), "01regionUniverse.Rds"))

##Return to short format
RUshort <- RU[,list(START=min(otherEndID), STOP=max(otherEndID)),by=c("baitID", "regionID")]

##Scramble enhancer locations
##get chromosome lengths
rmap <- Chicago:::.readRmap(list(rmapfile="Resources/digest/GRCh37_HindIII.rmap"))
bmap <- Chicago:::.readBaitmap(list(baitmapfile="Resources/digest/GRCh37_HindIII.baitmap"))
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

#RUshortBackup <- copy(RUshort)
RUshort <- shuffle(RUshort, bmap=bmap, rmap=rmap)

##Return to long format
RUshuffled <- RUshort[,c("baitIDshuff", "STARTshuff", "STOPshuff"),with=FALSE]
RUshuffled$regionID <- 1:nrow(RUshuffled)
names(RUshuffled) <- c("baitID","start","end","regionID")
RU <- RUshuffled[,list(otherEndID = start:end),by=c("baitID", "regionID")]

##Output
saveRDS(RU, file = file.path(filePath, "01regionUniverse.Rds"))
