#! /usr/bin/env Rscript
message("\n***runChicdiff.R\n\n")

### Parsing the command line ###

library(argparser)

if (packageVersion("argparser") < 0.3) {
  stop("argparser version (", packageVersion("argparser"), ") is out of date - 0.3 or later is required. Please open R and run install.packages('argparser') to update.")
}


message("\n")

defaultchicSettings <- function()
 {
   list(
     inputfiles= NA,
     peakfiles= NA,
     targetRDSorRDAs= NA,
     targetChs= NA,
     targetColumns= NA,
     rmapfile= NA,
     baitmapfile= NA,
     RUexpand= 5L, 
     score= 5, 
     saveRDS = TRUE,
     parallel = FALSE,
     printMemory = TRUE
   )
 }

args = commandArgs(trailingOnly=TRUE)
spec = matrix(c("<input-files>", "Full path to the configuration file, containing info about input file", 
                "<output-prefix>", "Output file names prefix (cannot contain folders)"),
                "<peak-files>", "Full path to Chicago peak file (or comma-separated list of files)",  byrow=TRUE, ncol=3)

p = arg_parser("Run Chicdiff from input files", name="Rscript runChicdiff.R")

p = add_argument(p, arg=spec[,1], help=spec[,2])

p = add_argument(p, arg="--design-dir", help = "Directory with restriction and bait maps (note the settings file has priority over these)", default = "")

p = add_argument(p, arg="--settings-file", help = "Full path to Chicdiff settings file", default = NA)

p = add_argument(p, arg="--target-coloumns", help = "Comma-separated list of target coloumns to take from peakFiles", default = NA)

p = add_argument(p, arg="--RUexpand", help = "Number of fragments to expand by in either direction", default=5L)

p = add_argument(p, arg="--score", help = "Score cutoff for filtering contacts in peakmatrix", default = 5)

p = add_argument(p, arg="--examples-prox-dist", help = "The distance limit for plotting \"proximal\" examples", default=1e6)

p = add_argument(p, arg="--saveRDS", help = "Save the intermediate files or not", default=TRUE)

p = add_argument(p, arg="--output-dir", help = "The name of the output directory (can be a full path)", default="<output-prefix>")

p = add_argument(p, arg="--run-mode", help = "Running in parallel (requires more memory, but is faster) or sequential mode (less memory, but slower)", default="<output-prefix>")

p = add_argument(p, arg="--print-memory", help = "Should chicdiffPipeline print out memory use?", flag=TRUE) # Possibly remove

opts = parse_args(p, args)


### Auxilliary functions ###
dir.create.ifNotThere = function(path, ...){
  if(file.exists(path)){
    if(!file.info(path)$isdir){
      message("Cannot create directory ", path, " as it coincides with an existing file\n")
      return(0)
    }else{
      message("\nWarning: directory ", path, " exists and will be reused.")
      return(1)
    }
  }
  else{
    dir.create(path, ...)
  }
}

moveToFolder = function(pattern, where){
  files = list.files(".", pattern)
  if (!length(files)){
    message("No files found matching the pattern ", pattern)
    return(0)
  }
  res = vector("numeric", length(files))
  i=1
  for(f in files){
    res[i] = file.rename(f, file.path(where, f))
    if (!res[i]){
      message("Warning: did not succeed in moving ", f, " into ", where)
    }
    i = i+1
  }
  min(res)
}

message("Loading the Chicdiff package and dependencies...\n")
library(Chicdiff)

### argparser read-in ###
if(packageVersion("argparser") < 0.4)
{
  names(opts) <- gsub("-", "_", names(opts))
}


inputFiles = opts[["<input_files>"]]
outPrefix_rel = opts[["<output_prefix>"]]
peakFiles = strsplit(opts[["<peak_files>"]], "\\,")[[1]] 

settingsFile = opts[["settings_file"]]
designDir = opts[["design_dir"]] 
tarcols = strsplit(opts[["target-coloumns"]], "\\,")

printMemory = opts[["print_memory"]] # Possibly remove

score = opts[["score"]]

proxLim = opts[["examples_prox_dist"]]

mode = opts[["run-mode"]]

outDir = ifelse(opts[["output_dir"]]=="<output-prefix>", outPrefix_rel, opts[["output_dir"]])

outPrefix = file.path(outDir, outPrefix_rel)

### main code ###
  
if(is.na(settingsFile)){
  settingsFile = NULL
}

message("Setting the experiment...\n")
setchic = setExperiment(designDir = designDir, settingsFile = settingsFile)


# if(length(inputFiles)>1){
#   message("\nReading config input files...\n")
#   cd = readAndMerge(inputFiles, cd)
# }else{
#   message("\n")
#   cd = readSample(inputFiles, cd)
# }
  
if(!dir.create.ifNotThere(outDir, recursive = TRUE)){
  stop(paste("Couldn't create folder", outDir, "\n"))
}



message("\n\nStarting chicdiffPipeline...\n")
output = chicdiffPipeline(setchic, inputfiles, targetChFiles, targetRDSorRDAFiles, printMemory = printMemory)

# if(saveDESeq){
#   saveRDS(DESeqOut, paste0(outPrefix, ".Rds"))
# }

saveRDS(output, paste0(outPrefix, ".Rds"))

logfile = file(paste0(outPrefix, "_params.txt"), "a")
sink(logfile, append=TRUE)
print(setchic)
print(sessionInfo())
sink(NULL)
close(logfile)

# FIX ME!!!!!!!!!!!!!!!!!!!!!!
# message("\n\nPlotting examples...\n")
# baits=plotBaits(, outfile=paste0(outPrefix, "_proxExamples.pdf"), maxD = proxLim)  


message("\n\nSorting output files into folders...\n")

curDir = getwd()
setwd(outDir)


if(!dir.create.ifNotThere("data")){
  stop("Couldn't create folder data\n")  
}


if(!dir.create.ifNotThere("diag_plots")){
  stop("Couldn't create folder diag_plots\n")
}
if(!dir.create.ifNotThere("examples")){
  stop("Couldn't create folder examples\n")      
}

  
if (!moveToFolder("\\.txt", "data")){
  stop("Couldn't move txt files to data folder\n")
}

if(!moveToFolder("\\.Rds", "data")){
    stop("Couldn't move the Rds file to data folder\n")
}

if(!moveToFolder("xamples\\.pdf", "examples")){
  stop("Couldn't move pdf file(s) to example folder\n")
}

if(!moveToFolder("\\.pdf", "diag_plots")){
  stop("Couldn't move pdf files to diag_plots folder\n")
}

setwd(curDir)

message("All done!\n")
