\name{DESeq2Wrap}
\alias{DESeq2Wrap}

\title{
DESeq2 wrapper
}
\description{
A wrapper around DESeq2 to run on Region Universe or Control Region Universe.
}
\usage{
DESeq2Wrap(defchic.settings, RU, FullRegionData, suffix = "")
}

\arguments{
  \item{defchic.settings}{settings that were set up by running \code{\link{setchicExperiment}}}
  \item{RU}{Region Universe set up by running \code{\link{getRegionUniverse}} or Control Region Universe set by by running \code{\link{getControlRegionUniverse}}}
  \item{FullRegionData}{If running for RU -> FullRegionData[[1]] - Full Region Data for RU. If running for RUcontrol - use FullRegionData[[2]] - Full Region Data for RU control}
  \item{suffix}{Optional: suffix to add to the name of the saved files}
}


\value{
  \item{out}{DESeq results file (in a data.table format, not a DESeq object file!) with normalization factors based on the FullMean}
  \item{outPilot.Rds}{if saveRDS is TRUE - saves DESeq results file (=out) as a .Rds data.table (not a DESeq object file!) with normalization factors based on the FullMean}
  \item{outStandardDESeq.Rds}{if saveRDS is TRUE - saves standard DESeq results file as a .Rds data.table (not a DESeq object file!) with the size factors estimated using the "median ratio method"}

}

\author{
Jonathan Cairns, William Orchard, Valeriya Malysheva and Mikhail Spivakov
}

\seealso{
  \code{\link{setChicExperiment}}, \code{\link{getRegionUniverse}}, \code{\link{getControlRegionUniverse}}, \code{\link{getFullRegionData1}}, \code{\link{getFullRegionData2}}}
  
\examples{
# path to recover test data
dataPath <- system.file("extdata", package="ChicdiffData")

## read in the settings
testDefChicSettings <- readRDS(file.path(dataPath, "defChicSettings.Rds"))
                               
# path to recover pre-processes test data
resultsPath <- file.path(dataPath, "CD4_Mono_results")
                               
## load pre-processed RU and RUcontrol data
testRU <- readRDS(file.path(resultsPath, "RegionUniverse.Rds"))
testRUcontrol <- readRDS(file.path(resultsPath, "ControlRegionUniverse.Rds"))
                               
# load pre-processed FullRegionData and FullControlRegionData
testFullRegionData <- readRDS(file.path(resultsPath, "FullRegionData.Rds"))
testFullControlRegionData <- readRDS(file.path(resultsPath, "FullControlRegionData.Rds"))
                               
# Run DESeq2Wrap
DESeqOut <- DESeq2Wrap(defchic.settings=testDefChicSettings, RU=testRU, FullRegionData=testFullRegionData)
DESeqOutControl <- DESeq2Wrap(defchic.settings=testDefChicSettings, RU=testRUcontrol, FullRegionData=testFullControlRegionData, suffix = "Control")
  
}

