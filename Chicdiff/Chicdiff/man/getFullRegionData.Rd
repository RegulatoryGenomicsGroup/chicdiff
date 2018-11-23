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
  \item{defchic.settings}{settings that were set up by running \code{\link{setchicExperiment}}}
  \item{RU}{Region Universe set up by running \code{\link{getRegionUniverse}}}
  \item{RUcontrol}{Control Region Universe set up by running \code{\link{getControlRegionUniverse3}}}
  \item{suffix}{Control Region Universe set up by running \code{\link{getControlRegionUniverse3}}}
}

\value{
  \item{FullregionData}{a list of 3 elements. [[1]] - FullRegionData for RU, [[2]] - FullRegionData for RUcontrol, [[3]] - countput, containing data for plotting}
  \item{countput.Rds}{countput, containing data for plotting}
  \item{FullRegionData.Rds}{ FullRegionData for RU, if saveRDS is TRUE}
  \item{FullControlRegionData.Rds}{ FullControlRegionData for RUcontrol, if saveRDS is TRUE}
}

\author{
Jonathan Cairns, Will Orchard, Valeriya Malysheva and Mikhail Spivakov
}


\seealso{
  \code{\link{setchicExperiment}}, \code{\link{getRegionUniverse}}, \code{\link{getControlRegionUniverse3}}
}
\examples{
defchic.settings <- setchicExperiment(designDir="path/to/DesignDir", settingsFile = "path/to/settingsFile")

RU <- getRegionUniverse(defchic.settings)

RUcontrol <- getControlRegionUniverse3(defchic.settings, RU)

FullRegionData <- getFullRegionData(defchic.settings, RU, RUcontrol, suffix = "")
}


