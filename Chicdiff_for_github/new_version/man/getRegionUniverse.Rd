\name{getRegionUniverse}
\alias{getRegionUniverse}

\title{
Get Region Universe.
}
\description{
Filters contacts by score, so that at least one replicate has the score higher than the threshold. It then expands the otherend (seed) by s (RUexpand) in both direction (if possible) to get "region universe". 
}
\usage{
getRegionUniverse(defchic.settings)  
}

\arguments{
  \item{defchic.settings} {defchic.settings defined by \code{\link{setchicExperiment}}
}
}

\value{
  \item{RegionUniverse.Rds} {a data.table contains contact - region annotations}
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
}

