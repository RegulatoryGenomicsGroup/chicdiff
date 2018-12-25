\name{chicdiffPipeline}
\alias{chicdiffPipeline}

\title{
Chicdiff pipeline function
}
\description{
This function runs data through the Chicdiff pipeline.
}
\usage{
chicdiffPipeline(chicdiff.settings)
}

\arguments{
  \item{chicdiff.settings}{
  chicdiff.settings defined using the \code{\link{setChicdiffExperiment}} function.
}
}

\value{
\code{setChicdiffExperiment} returns a data.table containing the full output of the chicdiff pipeline.

In addition, it saves the following data into files:
  \item{<outprefix>_results.Rds}{Data.table containing all information for contacts that went through the chicdiff comparative analysis. This file includes the contact details, the weight attributed to every contact related region and weighted pvalues (raw and adjusted).}
  \item{<outprefix>_settings.Rds}{File containing all settings used for the run.}
  \item{<outprefix>_countput.Rds}{Contains the counts, Bmean and the max score per contact}
  \item{<outprefix>_diffbaitPlot.pdf}{By default this plot shows random 4 differentially interacting regions.}
  \item{<outprefix>_IHWdecisionBoundaryPlot.png}{IHW Boundary decision plot, produced on the control data.}
  \item{<outprefix>_IHWweightPlot.png}{IHW weight plot, produced on the control data.}

Otherwise, if the saveAuxData setting in setChicdiffExperiment was set to TRUE, Chicdiff returns 4 more objects:
  \item{<outprefix>_RegionUniverse.Rds}{Contains contact - region annotations}
  \item{<outprefix>_ControlRegionUniverse.Rds}{Contains control contact - region annotations}
  \item{<outprefix>_FullRegionData.Rds}{Contains technical Chicago info on contact - region annotations}
  \item{<outprefix>_FullControlRegionData.Rds}{Contains technical Chicago info on control contact - region annotations}
  \item{<outprefix>_DESeqObj.Rds}{The DESeq object from the (final) run on the test data}
  \item{<outprefix>_DESeqObjControl.Rds}{The DESeq object from the run on the control data}
}

\author{
Jonathan Cairns, William Orchard, Valeriya Malysheva and Mikhail Spivakov
}
\details{
This pipeline runs the following functions in order:

- \code{\link{getRegionUniverse}}

- \code{\link{getControlRegionUniverse}}

- \code{\link{getFullRegionData}} 

- \code{\link{DESeq2Wrap}}

- \code{\link{IHWcorrection}}

}
\seealso{
  \code{\link{setChicdiffExperiment}},   \code{\link{plotDiffBaits}},   \code{\link{getCandidateInteractions}}
}

\examples{
# path to recover test data
dataPath <- system.file("extdata", package="ChicdiffData")

## read in the settings
testDefChicSettings <- readRDS(file.path(dataPath, "defChicSettings.Rds"))

# launch the Chicdiff pipeline
output <- chicdiffPipeline(defchic.settings = testDefChicSettings)
}

