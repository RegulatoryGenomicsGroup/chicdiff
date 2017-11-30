##-----ChicagoDiff pipeline----

library(data.table)
library(Chicago)
library(IRanges)

resourceDir <- "macrophage"
filePath <- file.path("Resources/", resourceDir)
if(!dir.exists(filePath)) dir.create(filePath)

#x <- fread(input = "/bi/group/fraser/GWAS_Processing/Samples/Rda_reweighted/full/individual_samples_full.txt")
#colnames(x)
targetColumns <- c("Mac0_13","Mac0_14","Mac0_15","Mac1_16","Mac1_18","Mac1_63","Mac2_19","Mac2_20","Mac2_21")

targetRDAs <- list(
  Mac = c("/bi/group/fraser/GWAS_Processing/Samples/Samples/Mac0/CHiCAGO/res_Mac0_13_Step2_chicago2_reweight/data/Mac0_13_Step2.RDa",
          "/bi/group/fraser/GWAS_Processing/Samples/Samples/Mac0/CHiCAGO/res_Mac0_14_Step2_chicago2_reweight/data/Mac0_14_Step2.RDa",
          "/bi/group/fraser/GWAS_Processing/Samples/Samples/Mac0/CHiCAGO/res_Mac0_15_Step2_chicago2_reweight/data/Mac0_15_Step2.RDa",
          "/bi/group/fraser/GWAS_Processing/Samples/Samples/Batch4/Mac1/CHiCAGO/res_Mac1_16_Step2_chicago2_reweight/data/Mac1_16_Step2.RDa",
          "/bi/group/fraser/GWAS_Processing/Samples/Samples/Batch4/Mac1/CHiCAGO/res_Mac1_18_Step2_chicago2_reweight/data/Mac1_18_Step2.RDa",
          "/bi/group/fraser/GWAS_Processing/Samples/Samples/Batch4/Mac1/CHiCAGO/res_Mac1_63_Step2_chicago2_reweight/data/Mac1_63_Step2.RDa",
          "/bi/group/fraser/GWAS_Processing/Samples/Samples/Batch2/Mac2/CHiCAGO/res_Step2_Mac2_19_chicago2_reweight/data/Step2_Mac2_19.RDa",
          "/bi/group/fraser/GWAS_Processing/Samples/Samples/Batch2/Mac2/CHiCAGO/res_Step2_Mac2_20_chicago2_reweight/data/Step2_Mac2_20.RDa",
          "/bi/group/fraser/GWAS_Processing/Samples/Samples/Batch2/Mac2/CHiCAGO/res_Step2_Mac2_21_chicago2_reweight/data/Step2_Mac2_21.RDa"
          )
)
targetRDAFiles <- unlist(targetRDAs)


targetChs <- list(
  Mac = c("Data/Mac0_13_Step1_bait_otherEnd_N_len_distSign.txt",
          "Data/Mac0_14_Step1_bait_otherEnd_N_len_distSign.txt",
          "Data/Mac0_15_Step1_bait_otherEnd_N_len_distSign.txt",
          "Data/Mac1_16_Step1_bait_otherEnd_N_len_distSign.txt",
          "Data/Mac1_18_Step1_bait_otherEnd_N_len_distSign.txt",
          "Data/Mac1_63_Step1_bait_otherEnd_N_len_distSign.txt",
          "Data/Step1_Mac2_19_bait_otherEnd_N_len_distSign.txt",
          "Data/Step1_Mac2_20_bait_otherEnd_N_len_distSign.txt",
          "Data/Step1_Mac2_21_bait_otherEnd_N_len_distSign.txt"
  )
)
targetChFiles <- unlist(targetChs)

peakMatrix <- "/bi/group/fraser/GWAS_Processing/Samples/Samples/Rda_reweighted/full/individual_samples_full.txt"

rmapfile <- "Resources/digest/GRCh37_HindIII.rmap"

#manner to expand in
RUmergeMode <- "union"
#number of fragments to expand by in either direction
RUexpand = 5L