\name{getFullRegionData2}
\alias{getFullRegionData2}

\title{
Get Full Region Data (parallel version)
}
\description{
Recovers the full data for every region - contact from targetRDSorRDAs and chinput files. Processes Region Universe and Control Region IUniverse in one go.
}
\usage{
getFullRegionData2(defchic.settings, RU, RUcontrol, suffix = "")
}
\arguments{
  \item{defchic.settings} {settings that were set up by running \code{\link{setchicExperiment}}}
  \item{RU} {Region Universe set up by running \code{\link{getRegionUniverse}}}
  \item{RUcontrol} {Control Region Universe set by by running \code{\link{getControlRegionUniverse}}}
  \item{suffix} {Optional: suffix to add to the name of the saved files}}
}

\value{
  \item{FullregionData} {a list of 3 elements: [[1]] - FullRegionData for RU, [[2]] - FullRegionData for RUcontrol, [[3]] - countput, containing data for plotting}
  \item{countput.Rds} {countput, containing data for plotting}
  \item{FullRegionData.Rds} {FullRegionData for RU, if saveRDS is TRUE}
  \item{FullControlRegionData.Rds} {FullControlRegionData for RUcontrol, if saveRDS is TRUE}
}

\author{
Jonathan Cairns, Will Orchard, Valeriya Malysheva and Mikhail Spivakov
}

\seealso{
  \code{\link{setchicExperiment}}, \code{\link{getRegionUniverse}}, \code{\link{getControlRegionUniverse}}
}
\examples{
defchic.settings <- setchicExperiment(designDir="path/to/DesignDir", settingsFile = "path/to/settingsFile")
 
RU <- getRegionUniverse(defchic.settings)

RUcontrol <- getControlRegionUniverse3(defchic.settings, RU)
 
FullRegionData <- getFullRegionData2(defchic.settings = defchic.settings, RU = RU, RUcontrol = RUcontrol, suffix = "")
}

