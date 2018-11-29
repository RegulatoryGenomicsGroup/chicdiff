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
  \item{defchic.settings}{defchic.settings defined by \code{\link{setChicExperiment}}}
}

\value{
  \item{RegionUniverse.Rds}{a data.table contains contact - region annotations}
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

# create Region Universe
RU <- getRegionUniverse(defchic.settings=testDefChicSettings) 
}

