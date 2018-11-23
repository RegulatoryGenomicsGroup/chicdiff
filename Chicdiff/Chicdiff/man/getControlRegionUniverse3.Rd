\name{getControlRegionUniverse3}
\alias{getControlRegionUniverse3}

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
getControlRegionUniverse3(defchic.settings, RU)
}
\arguments{
  \item{defchic.settings}{defchic.settings defined by running \code{\link{setchicExperiment}}}
  \item{RU}{Region Universe defined by running \code{\link{getRegionUniverse}}}
}

\value{
  \item{RegionUniverse.Rds}{a data.table containing contact - region annotations}
}

\author{
Jonathan Cairns, Will Orchard, Valeriya Malysheva and Mikhail Spivakov
}

\seealso{
  \code{\link{setchicExperiment}}
}

\examples{
defchic.settings <- setchicExperiment(designDir="path/to/DesignDir", settingsFile = "path/to/settingsFile")
RU <- getRegionUniverse(defchic.settings)
getControlRegionUniverse3(defchic.settings, RU)
}
