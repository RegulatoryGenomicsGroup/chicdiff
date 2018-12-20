\name{defaultchicSettings}
\alias{defaultchicSettings}
\title{
  Default Chicdiff settings
}
\description{
  A function that gives the default settings used for a Chicdiff experiment.
}
\usage{
  defaultchicSettings()
}
\value{
  A list of the following settings:
  \item{inputfiles}{Default: NA. The location of the input file; see the package vignette for a description of the file format.}
  \item{peakfiles}{Default: NA. The location of the Chicago peakfile.}
  \item{targetRDSorRDAs}{Default: NA. A named list of lists of .Rda files, corresponding to replicates of interest; see the package vignette for details.}
  \item{targetChs}{Default: NA. A named list of lists of Chinput files, corresponding to replicates of interest; see the package vignette for details.}
  \item{rmapfile}{Default: NA. The location of the restriction map file; see the package vignette for a description of the file format.}
  \item{targetColumns}{Default: NA. A list of target coloumns containing Chicago scores for corresponding replicates.}
  \item{baitmapfile}{Default: NA. The location of the bait map file; see the package vignette for a description of the file format.}
  \item{RUexpand}{Default: 5L. The number of fragments to expand by in either direction.}
  \item{score}{Default: 5. The filtering score threshold applied to contacts.}
  \item{norm}{Default: "combined" - using a weighted sum of sample-level and interaction-level scaling factors for normalisation. Alternatives are "standard" (sample-level only) and "fullmean" (interaction-level only).}
  \item{theta}{Default: NULL (corresponding to choosing theta automatically if norm="combined"). The mixing parameter for combined normalisation. Note that theta=1 is equivalent to norm="standard", and theta=0 is equivalent to norm="fullmean", with norm settings updated accordingly in these cases.}
  \item{theta_grid}{Default: seq(0,1,0.25). A grid of mixing parameter settings, from which the optimal value will be chosen if norm="combined".}
  \item{saveRDS}{Default: FALSE. Save the output of every step of chicdiff pipeline?}
  \item{parallel}{Default: FALSE. Run in parallel or stepwise mode? Parallel mode would require less time, but more memory.}
  \item{device}{Default: png. The format of the output figures.}
  \item{printMemory}{Default: TRUE. A garbage collection and print of memory usage statistics.}
}
\author{
  Jonathan Cairns, William Orchard, Valeriya Malysheva and Mikhail Spivakov
}
\seealso{
  \code{\link{setChicExperiment}}, \code{\link{.updateChicSettings}}
}
\examples{
DefChicSettings <- defaultchicSettings()
print(DefChicSettings)
}
