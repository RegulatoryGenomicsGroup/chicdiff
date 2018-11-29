\name{getFullRegionData1}
\alias{getFullRegionData1}

\title{
Get Full Region Data (non-parallel version)
}
\description{
Recovers the full data for every region - contact from targetRDSorRDAs and chinput files.
}
\usage{
getFullRegionData1(defchic.settings, RU, is_control = FALSE, suffix = "")
}

\arguments{
  \item{defchic.settings}{settings that were set up by running \code{\link{setChicExperiment}}}
  \item{RU}{Region Universe set up by running \code{\link{getRegionUniverse}} or Control Region Universe set by by running \code{\link{getControlRegionUniverse}}}
  \item{is_control}{is the data provided - RU or RU control? FALSE if RU = RU, TRUE if RU = RUcontrol}
  \item{suffix}{Optional: suffix to add to the name of the saved files}
}

\value{
  \item{FullregionData}{a list of 3 elements: [[1]] - FullRegionData for RU, [[2]] - FullRegionData for RUcontrol, [[3]] - countput, containing data for plotting}
  \item{countput.Rds}{countput, containing data for plotting}
  \item{FullRegionData.Rds}{FullRegionData for RU, if saveRDS is TRUE}
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
testDefChicSettings <- readRDS(file.path(dataPath, "defChicSettings.Rds"))

# path to recover pre-processes test data
resultsPath <- file.path(dataPath, "CD4_Mono_results")

## load pre-processed data
testRU <- readRDS(file.path(resultsPath, "RegionUniverse.Rds"))
testRUcontrol <- readRDS(file.path(resultsPath, "ControlRegionUniverse.Rds"))

#Get Full region Data
FullRegionData <- getFullRegionData1(defchic.settings = testDefChicSettings, RU = testRU, is_control = FALSE, suffix = "")

#Get Full control region Data
FullRegionData[[2]] <- getFullRegionData1(defchic.settings = testDefChicSettings, RU = testRUcontrol, is_control = TRUE, suffix = "")

}

