\name{chicdiffPipeline}
\alias{chicdiffPipeline}

\title{
Chicdiff pipeline function
}
\description{
This function runs data through the Chicdiff pipeline.
}
\usage{
chicdiffPipeline(defchic.settings)
}

\arguments{
  \item{defchic.settings}{
  defchic.settings set up using setchicExperiment function.
}
}

\value{
If saveRDS argument of setchicExperiment was set to FALSE, then Chicdiff returns 6 objects:
  \item{<test_name>.out}{File containing all settings used for the run, information about the version of all packages used to run the pipeline, comments of the pipeline as passing through every crucial function}
  \item{Weighted.Rds}{Data.table containing all information for contacts that went through the chicdiff comparative analysis. This file includes the contact details, the weight attributed to every contact related region and weighted pvalue.}
  \item{countput.Rds}{Contains the counts, Bmean and the max score per contact}
  \item{diffbaitPlot.pdf}{By default this plot shows random 4 differentially interacting regions.}
  \item{IHWdecisionBoundaryPlot.png}{IHW Boundary decision plot, produced on the control data.}
  \item{IHWweightPlot.png}{IHW weight plot, produced on the control data.}

Otherwise, if saveRDS argument of setchicExperiment was set to TRUE, Chicdiff returns 4 more objects:
  \item{RegionUniverse.Rds}{Contains contact - region annotations}
  \item{ControlRegionUniverse.Rds}{Contains control contact - region annotations}
  \item{FullRegionData.Rds}{Contains technical Chicago info on contact - region annotations}
  \item{FullControlRegionData.Rds}{Contains technical Chicago info on control contact - region annotations}
  \item{DESeq_object}{Standard DESeq object}
}

\author{
Jonathan Cairns, William Orchard, Valeriya Malysheva and Mikhail Spivakov
}
\details{
This pipeline runs the following functions in order:

- \code{\link{getRegionUniverse}}

- \code{\link{getControlRegionUniverse3}}

- \code{\link{getFullRegionData}} 

- \code{\link{DESeq2Wrap}}

- \code{\link{IHWcorrection}}

}
\seealso{
  \code{\link{setChicExperiment}}
}

\examples{
# path to recover test data
dataPath <- system.file("extdata", package="ChicdiffData")

## read in the settings
testDefChicSettings <- readRDS(file.path(dataPath, "defChicSettings.Rds"))

# launch the Chicdiff pipeline
output <- chicdiffPipeline(defchic.settings = testDefChicSettings)
}

