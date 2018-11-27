\name{plotdiffBaits}
\alias{plotdiffBaits}

\title{
Plot differential baits
}
\description{
This function is a part of IHW function. It plots n (n=3 by default) from top 100 differential baits. 
}
\usage{
plotdiffBaits(output, countput, baitmapfile, n = 3, baits = NULL, plotBaitNames = TRUE, plotBaitIDs = TRUE, plevel1 = 5, plevel2 = 3, xlim=c(-1e6,1e6), bgCol = "black", lev1Col = "red", lev2Col = "blue", ...)
}

\arguments{
  \item{output}{final output of IHW function, containing all the information about the bait - other end - region: coordinates, counts, log2FoldChange, DESeq pvalue, weight and weighted pvalue and lots of other data.}
  
  \item{countput}{countput file, containing counts recovered from the .Rda Chicago files, produced in  \code{\link{getFullRegionData1}} or \code{\link{getFullRegionData2}}.}

  \item{baitmapfile}{The location of the bait map file. Recovered from defchic.settings, produced by  \code{\link{setchicExperiment}}.}

  \item{n}{The number of baits to plot (ignored if \code{baits} is specified). The differential landscape of n random unique baits will be displayed.}

  \item{baits}{The IDs of particular baits to plot.}

  \item{plotBaitNames}{plot the bait names at the top of the plots? TRUE by default.}

  \item{plotBaitIDs}{plot the bait IDs at the top of the plots? TRUE by default.} 
  
  \item{plevel1}{the hard score threshold = 5}

  \item{plevel2}{the soft score threshold = 3}

  \item{xlim}{the area around the bait to plot. By default it is set up to 1Mb in each direction.}

  \item{bgCol}{The color of background contacts. It is black by default.}

  \item{lev1Col}{The color of contacts that were call significant at least in one of the replicates of a corresponding condition. It is red by default.}

  \item{lev2Col}{The color of contacts that got a subthreshold score at least in one of the replicates of a corresponding condition. It is blue by default.}
  
  \item{\dots}{Additional arguments passed to \code{\link{ggplot}}}
  
}

\value{
  Vector of the baitIDs plotted}


\author{
Jonathan Cairns, Will Orchard, Valeriya Malysheva and Mikhail Spivakov
}

\seealso{
  \code{\link{setchicExperiment}}, \code{\link{getFullRegionData1}}, \code{\link{getFullRegionData2}}, \code{\link{IHWcorrection}}
}

\examples{
#load here the output and countput...and dedault settings...

plotdiffBaits(output = out, countput = countput, n = 4, baitmapfile = baitmapfile)

# if you want to save the plot:
ggsave("diffbaitPlot.pdf", device = "pdf", path = "./")

}
