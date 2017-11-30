#Broadly what the pipeline function might look like

# chicdiffPipeline <- function(, ){
#   RU <- getRegionUniverse(, )
#   RUcontrol <- getControlRegionUniverse(RU, )
#   FullRegionData <- FullRegionData(, )
#   FullControlRegionData <- getFullRegionData(, ) #This could be combined into one which produces both control and non-control region data (see chicdiffCombined)
#   DESeqOut <- DESeq2Wrap(, )
#   DESeqOutControl <- DESeq2Wrap(, ) #Good idea to just do DESeq2 twice - once for control and once for not, but could just put this into a tiny wrapper which does both
#   final_output <- IHWcorrection(, )
#   final_output #This is the 06Weighted.Rds output from the original scripts
# }


# Job 2584667 (chicdiff_with_chinputs) Complete
# User             = orchardw
# Queue            = all.q@compute-0-6.local
# Host             = compute-0-6.local
# Start Time       = 09/07/2017 19:35:26
# End Time         = 09/07/2017 22:14:34
# User Time        = 02:44:03
# System Time      = 00:07:39
# Wallclock Time   = 02:39:08
# CPU              = 02:51:42
# Max vmem         = 52.342G
# Exit Status      = 0



source("/bi/group/sysgen/chicdiff/Will/CollectionOfFunctions.R") #Mimics doing library("CHiCdiff")
source("/bi/group/sysgen/chicdiff/Pipeline/00ConfigMonVsTCD4ActWindowPkg.R") #Mimics using setExperiment() 
printMemory = TRUE
printMemoryFunction <- function(printMemory){
  if(printMemory == TRUE){
    print(gc(reset=TRUE))
  }
}
#Pipeline - example workflow - many of these arguments have defaults, but I have assigned them explicitly here for 
#clarity. Also, the pipline funcion should just take the settings as default which contains all this info: 

#'targetRDSorRDAFiles' was a stupendous placeholder name by the way

RU <- getRegionUniverse(RUmergeMode = "window", RUexpand = 5L, rmapfile = rmapfile, peakMatrix = peakMatrix, targetColumns = targetColumns, score = 5, saveRDS = TRUE)
  printMemoryFunction(printMemory)

FullRegionData <- getFullRegionData(RU, targetChFiles = targetChFiles, targetRDSorRDAFiles = targetRDSorRDAFiles, targetRDSorRDAs = targetRDSorRDAs, rmapfile = rmapfile, saveRDS = TRUE)
  printMemoryFunction(printMemory)

DESeqOut <- DESeq2Wrap(RU, FullRegionData, rmapfile = rmapfile, saveRDS = TRUE)
  printMemoryFunction(printMemory)

RUcontrol <- getControlRegionUniverse(RU, rmapfile = rmapfile, baitmapfile = baitmapfile, saveRDS = TRUE)
  printMemoryFunction(printMemory)

FullControlRegionData <- getFullRegionData(RUcontrol, targetChFiles = targetChFiles, targetRDSorRDAFiles = targetRDSorRDAFiles, targetRDSorRDAs = targetRDSorRDAs, rmapfile = rmapfile, saveRDS = TRUE, suffix = "Control")
  printMemoryFunction(printMemory)

DESeqOutControl <- DESeq2Wrap(RUcontrol, FullControlRegionData, rmapfile = rmapfile, saveRDS = TRUE, suffix = "Control")
  printMemoryFunction(printMemory)

source("/bi/group/sysgen/chicdiff/Pipeline/00ConfigMonVsTCD4ActWindowPkg.R")
output <- IHWcorrection(DESeqOut, FullRegionData, DESeqOutControl, FullControlRegionData) #The only one which, by default, outputs an RDS (to working directory)
  printMemoryFunction(printMemory)
#Included suffixes for my own convenience here 

## Settings functions ----------------- #Just copies from Chicago at the moment, so not integrated.#------------------ ##

#Commented out means, from Chicago but perhaps applicable? 
# defaultSettings <- function()
# {
#   list(
#     rmapfile= NA,
#     baitmapfile= NA,
#     peakMatrix= NA,
#     targetColumns= NA,
#     targetRDSorRDAFiles= NA,
#     targetChFiles= NA,
#     RUmergeMode= "window",
#     RUexpand= 5L,
#     score= 5, #Give more informative name (see getRegionUniverse())
#     #Ncol = "N",
#     #baitIDcol = "baitID",
#     #otherEndIDcol = "otherEndID",
#     #otherEndLencol = "otherEndLen", 
#     #distcol = "distSign",
# 
#   )
# }