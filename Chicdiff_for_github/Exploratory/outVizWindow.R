##Take a look at the top DESeq hits
##1) work out which regions to plot
##2) go through the Rds files, plotting those regions

source("Pipeline/00ConfigMonVsTCD4ActWindow.R")

library(ggplot2)
library(data.table)
library(gridExtra)
library(lattice)
library(Hmisc)
library(Chicago)

geoMean <- function(x) {exp(sum(log(x)))}

RU <- readRDS("Resources/monVsTCD4ActWindow/01regionUniverse.Rds")
dispersions <- readRDS("Resources/monVsTCD4ActWindow/03dispersionParameters.Rds")
outRecast <- readRDS("Resources/monVsTCD4ActWindow/04recast.Rds")
outPilot <- readRDS("Resources/monVsTCD4ActWindow/05outPilot.Rds")
outNull <- readRDS("Resources/monVsTCD4ActWindow/05outStandardDESeq.Rds")

figFile <- "Resources/monVsTCD4ActWindow/06Figures.Rds" ##to be created if it doesn't exist

rmap <- fread("Resources/digest/GRCh37_HindIII.rmap")
names(rmap) = c("chr","start","end","otherEndID")

regionData <- outRecast[,list(
  N=sum(N),
  nfrags=.N,
  avDist=(min(distSign)+max(distSign))/2,
  Bmean=sum(Bmean),
  Tmean=sum(Tmean),
  FullMean=sum(FullMean),
  s_j=s_j[1]
),by=c("baitID", "regionID", "sample")]
setkey(regionData, regionID)

regionSummary <- regionData[,list(
  avN=mean(N),
  nfrags=nfrags[1],
  avDist=avDist[1],
  avBmean=mean(Bmean),
  avTmean=mean(Tmean),
  avFullMean=mean(FullMean),
  avs_j=mean(s_j)
), by=c("baitID", "regionID")]

setkey(regionSummary, regionID)

##Order of calls
sel <- order(outPilot$pvalue)
sel.Null <- order(outNull$pvalue)

##Take a look at things called by one method but not the other
sel.sig <- sel[1:sum(outPilot$padj < 0.05, na.rm=TRUE)]
sel.Null.sig <- sel.Null[1:sum(outNull$padj < 0.05, na.rm=TRUE)]
sel.excl <- sel.sig[which(!sel.sig %in% sel.Null.sig)]
sel.Null.excl <- sel.Null.sig[which(!sel.Null.sig %in% sel.sig)]

##Collect the IDs
#custom <- 79615
custom <- c(147534,75977,50084,90564,87127,171096,42931,45638,155834,115642)
IDsToPlot <- unlist(
    lapply(list(sel.sig, sel.Null.sig, sel.excl, sel.Null.excl), head)
  )
IDsToPlot <- c(IDsToPlot, custom)
mynames <- paste0(rep(
  c("Sig_NF_", "Sig_Null_", "Sig_NF_not_Null_", "Sig_Null_not_NF_"), each=6),
  rep(1:6, 4))
if(!is.null(custom)) {mynames <- c(mynames, paste0("Custom_", 1:length(custom)))}
names(IDsToPlot) <- mynames

##Work out what distance range to use for each interaction
maxDval <- baitIDval <- rep(0, length(IDsToPlot))
for(i in seq_along(IDsToPlot))
{
  rID <- IDsToPlot[i]
  maxDval[i] <- outPilot[regionID == rID, max(abs(c(baitstart - OEend, baitend - OEstart)))]
  baitIDval[i] <- outPilot[regionID == rID, baitID]
}
maxDval <- 2*maxDval

##Plotting function
plotBaitsStep1 <- function(cd, pcol="score", Ncol="N", n=16, baits=NULL,
    plotBaitNames=TRUE, plotBprof=FALSE,plevel1 = 5, plevel2 = 3,
    outfile=NULL, removeBait2bait=TRUE, width=20, height=20,
    maxD=1e6, bgCol="black", lev2Col="blue", lev1Col="red",
    bgPch=1, lev1Pch=20, lev2Pch=20, disp=NULL, ...)
{
  if(plotBaitNames){
    baitmap = Chicago:::.readBaitmap(cd@settings)
  }
  if (is.null(baits)){
    baits = sample(unique(cd@x$baitID),n)
  }
  else{
    n = length(baits)
  }
  
  ncols = ceiling(n/4)
  nrows = ceiling(n/ncols)
  par(mfrow=c(nrows, ncols))
  
  setkey(cd@x, baitID)
  
  myArgs = list(...)
  xlimInArgs = ("xlim" %in% names(myArgs))
  
  these <- vector("list", n)
  titles <- vector("list", n)
  for(i in 1:n){

    myminD = NULL
    if(xlimInArgs){
      myminD = myArgs$xlim[1]
      mymaxD = myArgs$xlim[2]
    }else{
      if(!is.null(maxD)){
        mymaxD = maxD[i]
        myminD = -maxD[i]
      }
    }
    
    this = cd@x[J(baits[i])]
    this = this[is.na(distSign)==FALSE]
    
    if (!is.null(mymaxD)){
      this = this[distSign<=mymaxD]
    }
    if(!is.null(myminD)){
      this = this[distSign>=myminD]
    }
    
    if (removeBait2bait){
      this = this[isBait2bait==FALSE]
    }
    
    setDF(this)
    this = this[order(this$distSign),]
    
    cols <- rep(bgCol, nrow(this))
    pchs <- rep(bgPch, nrow(this))
    sel1 <- this[,pcol] >=plevel1
    sel2 <- this[,pcol] >=plevel2
    cols[sel2] <- lev2Col ##less stringent first
    cols[sel1] <- lev1Col
    pchs[sel2] <- lev2Pch
    pchs[sel1] <- lev2Pch
    this$pchs <- pchs
    this$cols <- cols
    
    title = paste(baits[i], sep="")
    if(plotBaitNames){
      baitName = baitmap[baitmap$V4==baits[i]][, cd@settings$baitmapGeneIDcol, with=FALSE]
      if (length(grep(",",baitName))){
        baitName = gsub("(\\S+,).+","\\1", baitName)
        baitName = paste0(baitName, "...")
      }
      title = paste0(baitName, " (", title, ")")
    }
    these[[i]] <- this
    titles[[i]] <- title
  }
  list(these, titles)
}

plotBaitsStep2 <- function(this, title, pcol="score", Ncol="N", n=16, baits=NULL,
                           plotBaitNames=TRUE, plotBprof=FALSE,plevel1 = 5, plevel2 = 3,
                           outfile=NULL, removeBait2bait=TRUE, width=20, height=20,
                           maxD=1e6, bgCol="black", lev2Col="blue", lev1Col="red",
                           bgPch=1, lev1Pch=20, lev2Pch=20, disp=NULL, ...)
{
  plot(this$distSign, this[,Ncol], xlab="Distance from viewpoint", ylab=Ncol, main=title, col=this$cols, pch=this$pchs, ...)
  abline(v=0, col="grey", lwd=1)
  
  if(plotBprof){
    lines(this$distSign, this$Bmean, lwd=1, col="darkgrey")
    lines(this$distSign, this$Bmean+1.96*sqrt(this$Bmean+this$Bmean^2/disp), ##need dispersion
          lwd=1, lty=2, col="darkgrey")
  }
}

if(!file.exists(figFile))
{
  ##Systematically work through the Rda files, collecting data to plot.
  outPlots <- vector("list", length(targetRDAFiles))
  names(outPlots) <- targetRDAFiles
  
  for(i in 1:length(targetRDAFiles))
  {
    outPlots[[i]] <- vector("list", length(baitIDval))
    
    #myDisp <- dispersions[i]
    message("File ",i, " of ", length(targetRDAFiles))
    file <- targetRDAFiles[i]
    load(file) ##this produces the object 'x'
  
    cd = setExperiment(designDir="Resources/digest/")
    cd@params$dispersion = attributes(x)$dispersion
    setDT(x)
  
    cd@x = x
    cd@x[, score:=newScore]
    cd@x[, newScore:=NULL]
    for(j in seq_along(baitIDval))
    {
      outPlots[[i]] <- plotBaitsStep1(cd, baits = baitIDval, maxD=maxDval, plotBprof=TRUE)
    }
  }
  
  saveRDS(outPlots, file = figFile)
}

outPlots <- readRDS(figFile)

##part two -----------------------------------------------------------
##plot:
nSites <- length(outPlots[[1]][[2]])
nSamples <- length(outPlots)
pdf("Figures/FullPlotsWindows.pdf", height=14, width=14)
  for(j in 1:nSites)
  {
    ##get region information
    minFrag <- min(RU[regionID == IDsToPlot[j], otherEndID])
    maxFrag <- max(RU[regionID == IDsToPlot[j], otherEndID])
    baitFrag <- unique(RU[regionID == IDsToPlot[j], baitID])
    minX <- rmap[otherEndID == minFrag, start]
    maxX <- rmap[otherEndID == maxFrag, end]
    baitX <- rmap[otherEndID == baitFrag, (start+end)/2]
    myXlim <- c(minX - baitX - 100000, maxX - baitX + 100000)
    
    par(mfrow=c(4,2))
    for(i in 1:nSamples)
    {
      myThis <- copy(outPlots[[i]][[1]][[j]])
      ##paint the relevant restriction sites
      oes <- RU[regionID == IDsToPlot[j],otherEndID]
      if(any(myThis$otherEndID %in% oes))
      {
        myThis[myThis$otherEndID %in% oes,]$cols <- "magenta"
      }
      
      ##plot
      plotBaitsStep2(this=myThis,
                     title=paste0(outPlots[[i]][[2]][[j]],
                                  " ", names(IDsToPlot[j]), " R", IDsToPlot[j]),
                     xlim=c(-maxDval[j], maxDval[j]),
                     plotBprof=TRUE,
                     disp=dispersions[i])
    }
    
    par(mfrow=c(4,2))
    for(i in 1:nSamples)
    {
      myThis <- copy(outPlots[[i]][[1]][[j]])
      ##paint the relevant restriction sites
      oes <- RU[regionID == IDsToPlot[j],otherEndID]
      if(any(myThis$otherEndID %in% oes))
      {
        myThis[myThis$otherEndID %in% oes,]$cols <- "magenta"
      }
      
      ##plot
      plotBaitsStep2(this=myThis,
                     title=paste0(outPlots[[i]][[2]][[j]],
                                  " ", names(IDsToPlot[j]), " R", IDsToPlot[j]),
                     xlim=myXlim,
                     plotBprof=TRUE,
                     disp=dispersions[i])
    }
  }
dev.off()