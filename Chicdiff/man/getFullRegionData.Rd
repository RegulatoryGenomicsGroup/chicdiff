\name{getFullRegionData}
\alias{getFullRegionData}

\title{
Get Full Region Data
}
\description{
Recovers the full data for every region - contact from targetRDSorRDAs and chinput files.
If parallel = TRUE, reads in targetRDSorRDAs files only once for both "control" and "actual" contacts. If FALSE - reads twice: once for actual, and once for control.
Parallel mode is faster, but requires more memory.
}
\usage{
getFullRegionData(defchic.settings, RU, RUcontrol, suffix = "")
}

\arguments{
  \item{defchic.settings}{settings that were set up by running \code{\link{setChicdiffExperiment}}}
  \item{RU}{Region Universe set up by running \code{\link{getRegionUniverse}}}
  \item{RUcontrol}{Control Region Universe set up by running \code{\link{getControlRegionUniverse}}}
  \item{suffix}{Control Region Universe set up by running \code{\link{getControlRegionUniverse}}}
}

\value{
  \item{FullregionData}{a list of 3 elements. [[1]] - FullRegionData for RU, [[2]] - FullRegionData for RUcontrol, [[3]] - countput, containing data for plotting}
  \item{countput.Rds}{countput, containing data for plotting}
  \item{FullRegionData.Rds}{ FullRegionData for RU, if saveRDS is TRUE}
  \item{FullControlRegionData.Rds}{ FullControlRegionData for RUcontrol, if saveRDS is TRUE}
}

\author{
Jonathan Cairns, Willia Orchard, Valeriya Malysheva and Mikhail Spivakov
}


\seealso{
  \code{\link{setChicdiffExperiment}}, \code{\link{getRegionUniverse}}, \code{\link{getControlRegionUniverse}}
}

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

# Get Full Region Data and Full Control Region Data
FullRegionData <- getFullRegionData(defchic.settings = testDefChicSettings, RU = testRU , RUcontrol = testRUcontrol, suffix = "")

}


