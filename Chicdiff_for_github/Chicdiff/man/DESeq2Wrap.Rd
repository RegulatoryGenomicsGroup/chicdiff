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
  \item{defchic.settings} {settings that were set up by running \code{\link{setchicExperiment}}}
  \item{RU} {Region Universe set up by running \code{\link{getRegionUniverse}} or Control Region Universe set by by running \code{\link{getControlRegionUniverse}}}
  \item{FullRegionData} {If running for RU -> FullRegionData[[1]] - Full Region Data for RU. If running for RUcontrol - use FullRegionData[[2]] - Full Region Data for RU control}
  \item{suffix} {Optional: suffix to add to the name of the saved files}}
}


\value{
  \item{out}{DESeq results file (in a data.table format, not a DESeq object file!) with normalization factors based on the FullMean}
  \item{05outPilot.Rds }{if saveRDS is TRUE - saves DESeq results file (=out) as a .Rds data.table (not a DESeq object file!) with normalization factors based on the FullMean}
  \item{05outStandardDESeq }{standard DESeq results file as a .Rds data.table (not a DESeq object file!) with the size factors estimated using the "median ratio method"}

}

\author{
Jonathan Cairns, Will Orchard, Valeriya Malysheva and Mikhail Spivakov
}

\seealso{
  \code{\link{setchicExperiment}}, \code{\link{getRegionUniverse}}, \code{\link{getControlRegionUniverse}}, \code{\link{getFullRegionData1}}, \code{\link{getFullRegionData2}}}
  
\examples{
defchic.settings <- setchicExperiment(designDir="path/to/DesignDir", settingsFile = "path/to/settingsFile")
  
RU <- getRegionUniverse(defchic.settings)
  
RUcontrol <- getControlRegionUniverse3(defchic.settings, RU)
  
FullRegionData <- getFullRegionData2(defchic.settings = defchic.settings, RU = RU, RUcontrol = RUcontrol, suffix = "")
  
DESeqOut <- DESeq2Wrap(defchic.settings, RU, FullRegionData[[1]])
  
DESeqOutControl <- DESeq2Wrap(defchic.settings, RUcontrol, FullRegionData[[2]], suffix = "Control") 
}

