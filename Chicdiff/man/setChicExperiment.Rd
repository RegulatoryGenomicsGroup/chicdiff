\name{setchicExperiment}
\alias{setchicExperiment}
\title{
Set experiment
}
\description{
Creates the full set of settings for Chicdiff
}
\usage{
setchicExperiment(designDir="", targetRDSorRDAs = NA, targetChs = NA, peakfiles = NA, targetColumns = NA, settings=list(), settingsFile=NULL, inputfiles = NA, defchic.settings=defaultchicSettings())
}

\arguments{
  \item{designDir}{
  The location of the design directory, e.g "~/resources/path". (Should not end with a slash.)
}
  \item{targetRDSorRDAs}{
  A named list of lists, containing paths to Chicago .Rds or .Rda files, corresponding to the biological replicates of the study, e.g. "targetRDSorRDAs = list(
  condition_1 = c("~path/rep1.Rda","~path/rep2.Rda"), 
  condition_2 = c("~path/rep1.Rda","~path/rep2.Rda")
  )".
}
  \item{targetChs}{
  A named list of lists, containing paths to Chicago chinput files, corresponding to the biological replicates of the study, e.g. "targetChs = list(
  condition_1 = c("~path/rep1.txt","~path/rep2.txt"), 
  condition_2 = c("~path/rep1.txt","~path/rep2.txt")
  )".
}
  \item{peakfiles}{
  The location of the Chicago peakfile or several comma separated peakfiles.
}
  \item{targetColumns}{
  A list of target Columns to read from peakfile(s). These columns should contain Chicago scores.
}
  \item{inputfiles}{
  This is an alternative to described above way of providing .Rda/.Rds and chinput files. It is a three coloumn .txt file of tab or space separated condition, .Rda and chinput files. For more info check the Vignette.
}
  \item{settings}{
  A named list containing settings to apply. This overrides anything specified in \code{settingsFile}, or in \code{def.settings}.
}
  \item{settingsFile}{
  The location of a file containing settings. Each row should contain the name of a setting, followed by whitespace, followed by the value of that setting. Overrides anything specified in \code{def.settings}.
}
  \item{def.settings}{
  These are the default settings.
}
}

\value{
A list updated settings for a chicdiff experiment
}

\author{
  Jonathan Cairns, William Orchard, Valeriya Malysheva and Mikhail Spivakov
}
\seealso{
  \code{\link{defaultchicSettings}}
}
\examples{

# path to recover test data
dataPath <- system.file("extdata", package="ChicdiffData")

# path to DesignDir
testDesignDir <- file.path(dataPath, "designDir")

# recover paths to test data
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
testDefChicSettings <- setChicExperiment(designDir = testDesignDir, targetRDSorRDAs = targetRDSorRDAs, targetChs = targetChs, peakfiles = peakFiles, settingsFile = settingsFile)

}

