\name{getCandidateInteractions}
\alias{getCandidateInteractions}
\title{
Get Candidate Interactions
}
\description{
This function will indicate interactions that are most probably causing the difference between samples.
}
\usage{
getCandidateInteractions(x)
}
\arguments{
  \item{defchic.settings}{
The location of the Chicago peakfile or several comma separated peakfiles.
}
  \item{output}{final output of IHW function, containing all the information about a bait - region interaction: coordinates, counts, log2FoldChange, DESeq pvalue, weight and weighted pvalue and lots of other data.
}
  \item{peakFiles}{
The location of the Chicago peakfile or several comma separated peakfiles.
}
  \item{pvcut}{
User-defined confidence cutoff of weighted pvalue.
}
  \item{deltaAsinhScore}{
User-defined minimal difference between the asinh of Chicago scores of studied samples.
}
}

\value{
Data.frame object with full full inofrmation about candidate interactions, including chromosome, coordinates, genes, Chicago scores and weighted pvalues.
}
\author{
Jonathan Cairns, William Orchard, Valeriya Malysheva and Mikhail Spivakov
}

\seealso{
  \code{\link{setChicExperiment}}, \code{\link{chicdiffPipeline}}
}
\examples{
## path to recover test data
dataPath <- system.file("extdata", package="ChicdiffData")

## read in some raw data and settings
peakFiles <- file.path(dataPath,"peakMatrix_cd4_mono_unmerged.txt")
defchic.settings <- readRDS(file.path(dataPath, "defChicSettings.Rds")

#load pre-processed data
resultsPath <- file.path(dataPath, "CD4_Mono_results")
output <- file.path(resultsPath, "weighted.Rds")

#get candidate interactions
out <- getCandidateInteractions(defchic.settings, output = output, peakFiles = peakFiles, pvcut = 0.05, deltaAsinhScore = 1)
}

