##-----ChicagoDiff pipeline----

library(data.table)
library(Chicago)
library(IRanges)

#x <- fread(input = "/bi/group/fraser/GWAS_Processing/Samples/Rda_reweighted/full/individual_samples_full.txt")
#colnames(x)
resourceDir <- "monVsTCD4Act"
filePath <- file.path("Resources/", resourceDir)
if(!dir.exists(filePath)) dir.create(filePath)

#targetColumns <- c("Mon_1","Mon_2","Mon_3","TCD4Act_40","TCD4Act_42","TCD4Act_44")
targetColumns <- c("Mon_1","Mon_2","Mon_3","NCD4_22", "NCD4_23", "NCD4_24", "NCD4_64")

targetRDAs <- list(
  CD4 = c("/bi/group/fraser/GWAS_Processing/Samples/Samples/Batch3/NCD4/CHiCAGO/res_Step2_NCD4_22_chicago2_reweight/data/Step2_NCD4_22.RDa",
          "/bi/group/fraser/GWAS_Processing/Samples/Samples/Batch3/NCD4/CHiCAGO/res_Step2_NCD4_23_chicago2_reweight/data/Step2_NCD4_23.RDa",
          "/bi/group/fraser/GWAS_Processing/Samples/Samples/Batch3/NCD4/CHiCAGO/res_Step2_NCD4_24_chicago2_reweight/data/Step2_NCD4_24.RDa",
          "/bi/group/fraser/GWAS_Processing/Samples/Samples/Batch3/NCD4/CHiCAGO/res_Step2_NCD4_64_chicago2_reweight/data/Step2_NCD4_64.RDa"),
  
  Mono = c("/bi/group/fraser/GWAS_Processing/Samples/Samples/Mon/CHiCAGO/res_Mon_1.1.Mon_1.2.Mon_1.3_Step2_chicago2_reweight/data/Mon_1.1.Mon_1.2.Mon_1.3_Step2.RDa",
           "/bi/group/fraser/GWAS_Processing/Samples/Samples/Mon/CHiCAGO/res_Mon_2.2.Mon_2.3.Mon_2.4_Step2_chicago2_reweight/data/Mon_2.2.Mon_2.3.Mon_2.4_Step2.RDa",
           "/bi/group/fraser/GWAS_Processing/Samples/Samples/Mon/CHiCAGO/res_Mon_3.1.Mon_3.2.Mon_3.3_Step2_chicago2_reweight/data/Mon_3.1.Mon_3.2.Mon_3.3_Step2.RDa")
)
targetRDAFiles <- unlist(targetRDAs)


targetChs <- list(
  CD4 = c("Data/Step1_NCD4_22_bait_otherEnd_N_len_distSign.txt",
          "Data/Step1_NCD4_23_bait_otherEnd_N_len_distSign.txt",
          "Data/Step1_NCD4_24_bait_otherEnd_N_len_distSign.txt",
          "Data/Step1_NCD4_64_bait_otherEnd_N_len_distSign.txt"
  ),
  Mono = c("Data/Mon_1.1.Mon_1.2.Mon_1.3_Step1_bait_otherEnd_N_len_distSign.txt",
           "Data/Mon_2.2.Mon_2.3.Mon_2.4_Step1_bait_otherEnd_N_len_distSign.txt",
           "Data/Mon_3.1.Mon_3.2.Mon_3.3_Step1_bait_otherEnd_N_len_distSign.txt")
)
targetChFiles <- unlist(targetChs)

peakMatrix <- "/bi/group/fraser/GWAS_Processing/Samples/Samples/Rda_reweighted/full/individual_samples_full.txt"

rmapfile <- "Resources/digest/GRCh37_HindIII.rmap"

#manner to expand in
RUmergeMode <- "union"
#number of fragments to expand by in either direction
RUexpand = 5L