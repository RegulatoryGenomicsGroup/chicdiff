\name{IHWcorrection}
\alias{IHWcorrection}

\title{
IHW correction
}
\description{
IHW weight estimation and pvalue correction.
}
\usage{
IHWcorrection(defchic.settings, DESeqOut, FullRegionData, DESeqOutControl, FullControlRegionData,
                          countput, DiagPlot = TRUE, diffbaitPlot = TRUE, suffix = "")
}

\arguments{
  \item{defchic.settings}{settings defined by \code{\link{setchicExperiment}}.}
  \item{DESeqOut}{output of the \code{\link{DESeq2Wrap}} function run on the Region Universe.}
  \item{FullRegionData}{Full Region Universe data of the \code{\link{getFullRegionData}}. It is the first element of the FullRegionData vector of three elements.}
  \item{DESeqOutControl}{output of the \code{\link{DESeq2Wrap}} function run on the Control Region Universe}
  \item{FullControlRegionData}{Full Control Region Universe data of the \code{\link{getFullRegionData}}. It is the second element of the FullRegionData vector of three elements.}
  \item{countput}{Countput produced by \code{\link{getFullRegionData}}. It is the third element of the FullRegionData vector of three elements.}
  \item{DiagPlot}{Whether to plot or not the IHW control diagrams: decision boundary and weight plots? TRUE by default.}
  \item{diffbaitPlot}{Whether to plot or not the top differential baits? TRUE by default. The plot shows four random unique baits from the 100 "most differential baits" (baits containing one or more differential contacts with the lowest weighted pvalue).}
  \item{suffix}{Optional. The text to add to the name of the final output file, containig the weighted pvalue for every "bait - other end - region" contact.}
}

\value{

 \item{weighted.Rds }{final output file, containing all the information about the "bait - other end - region" contact: coordinates, counts, log2FoldChange, DESeq pvalue, weight and weighted pvalue and lots of other data.}
 \item{diffbaitplot.pdf }{if diffbaitPlot is TRUE - the plot with four of 100 most differential baits will be saved'.}

}

\author{
Jonathan Cairns, William Orchard, Valeriya Malysheva and Mikhail Spivakov
}

\seealso{
  \code{\link{setChicExperiment}}, \code{\link{getFullRegionData}}, \code{\link{DESeq2Wrap}}, \code{\link{plotdiffBaits}}
}

\examples{
# path to recover test data
dataPath <- system.file("extdata", package="ChicdiffData")

## read in the settings
testDefChicSettings <- readRDS(file.path(dataPath, "defChicSettings.Rds"))

# path to recover pre-processes test data
resultsPath <- file.path(dataPath, "CD4_Mono_results")

# load pre-processed FullRegionData and FullControlRegionData
testFullRegionData <- readRDS(file.path(resultsPath, "FullRegionData.Rds"))
testFullControlRegionData <- readRDS(file.path(resultsPath, "FullControlRegionData.Rds"))
testCountput <- readRDS(file.path(resultsPath, "countput.Rds"))

testDESeqOut <- readRDS(file.path(resultsPath, "outPilot.Rds"))
testDESeqOutControl <- readRDS(file.path(resultsPath, "outPilotControl.Rds"))


# Run IHW correction
output <- IHWcorrection(defchic.settings = testDefChicSettings, DESeqOut = testDESeqOut, FullRegionData = testFullRegionData, 
              DESeqOutControl = testDESeqOutControl, FullControlRegionData = testFullControlRegionData, countput = testCountput)

}
