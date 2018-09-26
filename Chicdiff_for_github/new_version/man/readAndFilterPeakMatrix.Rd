\name{readAndFilterPeakMatrix}
\alias{readAndFilterPeakMatrix}

\title{
Read and filter the peakmatrix
}
\description{
This function reads through the peakfile(s) and applies the defind threshold (score) on Chicago scores recovered from target Columns (targetColumns). It retains the rows (contacts) where at least one score is larger than the specified score and filters out trans interactions. 
}
\usage{
readAndFilterPeakMatrix(peakFiles, targetColumns, score)
}

\arguments{
  \item{peakFiles}{
The location of the Chicago peakfile or several comma separated peakfiles.
}
  \item{targetColumns}{
The target columns containing the Chicago scores of replicates under the study.
}
  \item{score}{
The Chicago score threshold to apply for filtering the contacts.
}
}
\details{
This function is applied within the getRegionUniverse function. Withing that function all arguments are recovered from the final set of settings (defchic.settings). 
}
\value{
Returns the data.table object of filtered contacts.
}

\author{
Jonathan Cairns, Will Orchard, Valeriya Malysheva and Mikhail Spivakov
}

\seealso{
  \code{\link{getRegionUniverse}}
}
\examples{
  designDir <- "~/path/to/DesignDir"
  peakFiles <- "~/path/to/peakfile.txt"
  targetColumns <- c("condition_1", "condition_1")
  x <- readAndFilterPeakMatrix(peakFiles = peakFiles, score = score, targetColumns = targetColumns)
  }


