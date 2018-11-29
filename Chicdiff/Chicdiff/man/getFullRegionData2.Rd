\name{getFullRegionData2}
\alias{getFullRegionData2}

\title{
Get Full Region Data and Full Control Region Data(parallel version)
}
\description{
Recovers the full data for every region - contact from targetRDSorRDAs and chinput files. Processes Region Universe and Control Region IUniverse in one go.
}
\usage{
getFullRegionData2(defchic.settings, RU, RUcontrol, suffix = "")
}
\arguments{
  \item{defchic.settings}{settings that were set up by running \code{\link{setChicExperiment}}}
  \item{RU}{Region Universe set up by running \code{\link{getRegionUniverse}}}
  \item{RUcontrol}{Control Region Universe set by by running \code{\link{getControlRegionUniverse}}}
  \item{suffix}{Optional: suffix to add to the name of the saved files}
}

\value{
  \item{FullregionData}{a list of 3 elements: [[1]] - FullRegionData for RU, [[2]] - FullRegionData for RUcontrol, [[3]] - countput, containing data for plotting}
  \item{countput.Rds}{countput, containing data for plotting}
  \item{FullRegionData.Rds}{FullRegionData for RU, if saveRDS is TRUE}
  \item{FullControlRegionData.Rds}{FullControlRegionData for RUcontrol, if saveRDS is TRUE}
}

\author{
Jonathan Cairns, William Orchard, Valeriya Malysheva and Mikhail Spivakov
}

\seealso{
  \code{\link{setChicExperiment}}, \code{\link{getRegionUniverse}}, \code{\link{getControlRegionUniverse}}
}
\examples{
# path to recover test data
dataPath <- system.file("extdata", package="ChicdiffData")

## read in the settings
testDefChicSettings <- readRDS(file.path(dataPath, "defChicSettings.Rds")

# path to recover pre-processes test data
resultsPath <- file.path(dataPath, "CD4_Mono_results")

## load pre-processed data
testRU <- file.path(resultsPath, "RegionUniverse.Rds")
testRUcontrol <- file.path(resultsPath, "ControlRegionUniverse.Rds")

# get FUll Region Data and Full Control Region Data in parallel 
FullRegionData <- getFullRegionData2(defchic.settings = defchic.settings, RU = RU, RUcontrol = RUcontrol, suffix = "")
}

