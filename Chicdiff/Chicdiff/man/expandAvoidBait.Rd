\name{expandAvoidBait}
\alias{expandAvoidBait}

\title{
Expands the region from the other end
}
\description{
Expands the region from the seed by 5 restriction fragments (RUexpand = 5, by default) both directions (if possible). But if the bait is too close to the seed (other end), make the range shorter.
}
\usage{
expandAvoidBait <- function(bait, oe, s)
}

\arguments{
  \item{bait}{baitID}
  \item{oe}{otherEndID, the seed}
  \item{s}{RUexpand}  
}
\details{
This function is a part of \code{\link{getRegionUniverse}} and \code{\link{getControlRegionUniverse3}}
}
\value{
  Function usually returns (oe-s):(oe+s) - range around oe, if the bait and oe are >s+1 apart.
  Otherwise (bait + 2L):(oe + s), (oe - s):(bait - 2L) , depending on a direction of the seed relatively to the bait.
  }

\author{
Jonathan Cairns, Will Orchard, Valeriya Malysheva and Mikhail Spivakov
}

\seealso{
\code{\link{getRegionUniverse}}, \code{\link{getControlRegionUniverse3}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x) 
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
