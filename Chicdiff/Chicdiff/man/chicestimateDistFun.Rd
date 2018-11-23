\name{chicestimateDistFun}
\alias{chicestimateDistFun}

\title{
Estimate the Distance Function
}
\description{
Estimates the function that models how the expected number of counts decreases with increasing distance. Used in \code{\link{getFullRegionData}}.
}
\usage{
chicestimateDistFun(x)
}

\arguments{
  \item{x}{Chicago .Rds file read with \code{\link{readRDSorRDA}}
}
}

\value{
A list object with the parameters of the distance function.
}

\author{
Jonathan Cairns, Will Orchard, Valeriya Malysheva and Mikhail Spivakov
}

\seealso{
  \code{\link{readRDSorRDA}}
}

\examples{
setchicExperiment(designDir="path/to/DesignDir", settingsFile = "path/to/settingsFile")
targetRDSorRDAs = defchic.settings[["targetRDSorRDAs"]]
file <- targetRDSorRDAFiles[i]
x <- readRDSorRDA(file)
chicestimateDistFun(x)
}
