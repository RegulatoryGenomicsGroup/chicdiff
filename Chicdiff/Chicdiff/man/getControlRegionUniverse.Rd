\name{getControlRegionUniverse}
\alias{getControlRegionUniverse}

\title{
Get control Region Universe.
}
\description{
For every contact from Region Universe this function:
  1. randomly selects a bait from the same chromsome, 
  2. attributes to this bait a distance, selected from a normal distribution limited by the size of the chromosome 
  3. recovers the count data for this contact
}
\usage{
getControlRegionUniverse(defchic.settings, RU)
}
\arguments{
  \item{defchic.settings}{defchic.settings defined by running \code{\link{setChicExperiment}}}
  \item{RU}{Region Universe defined by running \code{\link{getRegionUniverse}}}
}

\value{
  \item{RegionUniverse.Rds}{a data.table containing contact - region annotations}
}

\author{
Jonathan Cairns, William Orchard, Valeriya Malysheva and Mikhail Spivakov
}

\seealso{
  \code{\link{setChicExperiment}}
}

\examples{
# path to recover test data
dataPath <- system.file("extdata", package="ChicdiffData")

## read in the test settings
testDefChicSettings <- readRDS(file.path(dataPath, "defChicSettings.Rds")

# path to recover pre-processes test data
resultsPath <- file.path(dataPath, "CD4_Mono_results")

## load pre-processed RU and RUcontrol data
testRU <- file.path(resultsPath, "RegionUniverse.Rds")

RUcontrol <- getControlRegionUniverse(defchic.settings = testDefChicSettings, RU = testRU)
}
