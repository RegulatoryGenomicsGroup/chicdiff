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

## setup designDir
testDesignDir <- file.path(dataPath, "designDir")

## read in some raw data
peakFiles <- file.path(dataPath,"peakMatrix_cd4_mono_unmerged.txt")
testDataPath_CD4 <- file.path(dataPath, "CD4")
testDataPath_Mono <- file.path(dataPath, "monocytes")
targetChs <- list(
  CD4 = c(NCD4_22 = file.path(testDataPath_CD4, "unitTest_CD41.chinput"),
          NCD4_23 = file.path(testDataPath_CD4, "unitTest_CD42.chinput")
  ),
  
  Mono = c(Mon_2 = file.path(testDataPath_Mono, "unitTest_Mono2.chinput"),
           Mon_3 = file.path(testDataPath_Mono, "unitTest_Mono3.chinput")
  )
)
testDataPath_rda <- system.file("data", package="ChicdiffData")

targetRDSorRDAs <- list(
  CD4 = c(NCD4_22 = file.path(testDataPath_rda, "unitTest_CD41.RDa"),
          NCD4_23 = file.path(testDataPath_rda, "unitTest_CD42.RDa")
  ),
  
  Mono = c(Mon_2 = file.path(testDataPath_rda, "unitTest_Mono2.RDa"),
           Mon_3 = file.path(testDataPath_rda, "unitTest_Mono3.RDa")
  )
)

settingsFile = file.path(dataPath,"SettingsFile.txt")


#set up Chicdiff settings
defchic.settings <- setChicExperiment(designDir = testDesignDir, targetRDSorRDAs = targetRDSorRDAs, targetChs = targetChs, peakfiles = peakFiles, settingsFile = settingsFile)

#load pre-processed data
resultsPath <- file.path(dataPath, "CD4_Mono_results")
output <- file.path(resultsPath, "weighted.Rds")

#get candidate interactions
out <- getCandidateInteractions(defchic.settings, output = output, peakFiles = peakFiles, pvcut = 0.05, deltaAsinhScore = 1)
}

