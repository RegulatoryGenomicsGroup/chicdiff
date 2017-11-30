#source("http://bioconductor.org/biocLite.R")
#biocLite("IRanges")
#biocLite("Chicago")
#biocLite("data.table")
#biocLite("argparser")
#biocLite("IHW")
#biocLite("ggplot2")
#biocLite("hexbin")

source("/bi/group/sysgen/malysheva/chicdiff/chicdiff_pkg/Chicdiff_v2.R")

library(data.table)
library(argparser)
library(IRanges)
library(Chicago)

targetRDSorRDAs <- list(
   CD4 = c("/bi/group/fraser/GWAS_Processing/Samples/Samples/Batch3/NCD4/CHiCAGO/res_Step2_NCD4_22_chicago2_reweight/data/Step2_NCD4_22.RDa",
                      "/bi/group/fraser/GWAS_Processing/Samples/Samples/Batch3/NCD4/CHiCAGO/res_Step2_NCD4_23_chicago2_reweight/data/Step2_NCD4_23.RDa",
                       "/bi/group/fraser/GWAS_Processing/Samples/Samples/Batch3/NCD4/CHiCAGO/res_Step2_NCD4_24_chicago2_reweight/data/Step2_NCD4_24.RDa",
                       "/bi/group/fraser/GWAS_Processing/Samples/Samples/Batch3/NCD4/CHiCAGO/res_Step2_NCD4_64_chicago2_reweight/data/Step2_NCD4_64.RDa"),
   Mono = c("/bi/group/fraser/GWAS_Processing/Samples/Samples/Mon/CHiCAGO/res_Mon_1.1.Mon_1.2.Mon_1.3_Step2_chicago2_reweight/data/Mon_1.1.Mon_1.2.Mon_1.3_Step2.RDa",
                           "/bi/group/fraser/GWAS_Processing/Samples/Samples/Mon/CHiCAGO/res_Mon_2.2.Mon_2.3.Mon_2.4_Step2_chicago2_reweight/data/Mon_2.2.Mon_2.3.Mon_2.4_Step2.RDa",            
            "/bi/group/fraser/GWAS_Processing/Samples/Samples/Mon/CHiCAGO/res_Mon_3.1.Mon_3.2.Mon_3.3_Step2_chicago2_reweight/data/Mon_3.1.Mon_3.2.Mon_3.3_Step2.RDa")
   )

targetRDSorRDAFiles <- unlist(targetRDSorRDAs)

designDir <- "/bi/group/sysgen/malysheva/chicdiff/chicdiff_pkg/DesignDir"
peakMatrix <- "/bi/group/fraser/GWAS_Processing/Samples/Samples/Rda_reweighted/full/individual_samples_full.txt"
outprefix = "chicdiff_v2"

targetChFiles = NULL
settingsFile = NULL

defchic.settings = setchicExperiment(designDir = designDir, settingsFile = settingsFile)
output = chicdiffPipeline(defchic.settings, peakMatrix, targetChFiles, targetRDSorRDAFiles, outprefix = outPrefix, printMemory = printMemory)
