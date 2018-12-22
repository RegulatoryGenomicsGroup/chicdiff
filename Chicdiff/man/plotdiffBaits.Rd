\name{plotdiffBaits}
\alias{plotdiffBaits}

\title{
Plot differential baits
}
\description{
This function is a part of IHW function. It plots n (n=3 by default) from top 100 differential baits. 
}
\usage{
plotdiffBaits(output, countput, baitmapfile, n = 3, baits = NULL, plotBaitNames = TRUE, plotBaitIDs = TRUE, pcol = "weighted_padj", plevel1 = 5, plevel2 = 3, xlim=c(-1e6,1e6), bgCol = "black", lev1Col = "red", lev2Col = "blue", ...)
}

\arguments{
  \item{output}{final output of IHW function, containing all the information about a bait - region interaction: coordinates, counts, log2FoldChange, DESeq pvalue, weight and weighted pvalue and lots of other data.}
  
  \item{countput}{countput file, containing counts recovered from the .Rda Chicago files, produced in  \code{\link{getFullRegionData1}} or \code{\link{getFullRegionData2}}.}

  \item{baitmapfile}{The location of the bait map file. Recovered from defchic.settings, produced by  \code{\link{setChicExperiment}}.}

  \item{n}{The number of baits to plot (ignored if \code{baits} is specified). The differential landscape of n random unique baits will be displayed.}

  \item{baits}{The IDs of particular baits to plot.}

  \item{plotBaitNames}{plot the bait names at the top of the plots? TRUE by default.}

  \item{plotBaitIDs}{plot the bait IDs at the top of the plots? TRUE by default.} 
  
  \item{pcol}{the name of the p-value column to plot (default: weighted_padj).}

  \item{plevel1}{the hard score threshold (default: 5).}

  \item{plevel2}{the soft score threshold (default: 3).}

  \item{xlim}{the area around the bait to plot. By default it is set up to 1Mb in each direction.}

  \item{bgCol}{The color of background contacts. It is black by default.}

  \item{lev1Col}{The color of contacts that were call significant at least in one of the replicates of a corresponding condition. It is red by default.}

  \item{lev2Col}{The color of contacts that got a subthreshold score at least in one of the replicates of a corresponding condition. It is blue by default.}
  
  \item{\dots}{Additional arguments passed to \code{\link{ggplot}}}
  
}

\value{
  Vector of the baitIDs plotted}


\author{
Jonathan Cairns, William Orchard, Valeriya Malysheva and Mikhail Spivakov
}

\seealso{
  \code{\link{setChicExperiment}}, \code{\link{getFullRegionData1}}, \code{\link{getFullRegionData2}}, \code{\link{IHWcorrection}}
}

\examples{

#load the output and countput
dataPath <- system.file("extdata", package="ChicdiffData")
resultsPath <- file.path(dataPath, "CD4_Mono_results")

mini_output <- readRDS(file.path(resultsPath, "06Weighted.Rds"))
mini_countput <- readRDS(file.path(resultsPath, "countput.Rds"))

testDesignDir <- file.path(dataPath, "designDir")
bmap <- fread(file.path(testDesignDir, "chr19_GRCh37_HindIII.baitmap"))
baits <- c(331172, 320562, 327700)

#plot bait interactions
plotdiffBaits(output = mini_output, countput = mini_countput, baitmapfile = bmap, baits = baits)

}
