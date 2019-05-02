# Community detection

source("code/contingencyTable.R")

# 1. Load datesets ----
insects <- read.table("datasets/arthropods_clean.txt")
treats  <- read.table("datasets/treatments_clean.txt")
plants  <- read.table("datasets/plants_clean.txt")
size_dat <-read.table("datasets/size_dat_bio.txt")

source("code/weighted-modularity-LPAwbPLUS/code/R/LPA_wb_plus.R")
source("code/weighted-modularity-LPAwbPLUS/code/R/MODULARPLOT.R") #read in plotting function

MAT <- subinsct

# Modularity plots - seems like 
for (i in 1:36){
  MAT <- gardnets[[i]]
  if(dim(MAT)[1] == 1){next}
  # MOD1 = LPA_wb_plus(MAT) # find labels and weighted modularity using LPAwb+
  MOD2 = DIRT_LPA_wb_plus(MAT) # find labels and weighted modularity using DIRTLPAwb+
  # MOD3 = DIRT_LPA_wb_plus(MAT>0, 2, 20) 
  patch <- paste("figs/", "gard", i, sep="")
  pdf(patch, width = 5, height=5) 
  MODULARPLOT(MAT,MOD2)
  dev.off()
  # show the modular network configuration found in MOD1. Row and column numbering indicates the ordering of rows and columns in MAT. Modules are highlighted in red rectangles.
}
MAT <- gardnets[[19]]
MAT <- subinsct
MOD1 = LPA_wb_plus(MAT) # find labels and weighted modularity using LPAwb+
MOD2 = DIRT_LPA_wb_plus(MAT) # find labels and weighted modularity using DIRTLPAwb+
MOD3 = DIRT_LPA_wb_plus(MAT>0, 2, 20) 
MODULARPLOT(MAT,MOD2) # show the modular network configuration found in MOD1. Row and column numbering indicates the ordering of rows and columns in MAT. Modules are highlighted in red rectangles.

#use with R library 'bipartite'
library("bipartite")
source("code/weighted-modularity-LPAwbPLUS/code/R/convert2moduleWeb.R") # read in conversion function 

MOD1modWeb = convert2moduleWeb(MAT,MOD2) # converts the configuration found in MOD1 to a moduleWeb object
plotModuleWeb(MOD1modWeb) # plot the corresponding moduleWeb object using plotModuleWeb in library bipartite


source("code/weighted-modularity-LPAwbPLUS/code/R/GetModularInformation.R") #read in function for finding additional information

MOD1information = GetModularInformation(MAT,MOD1)
print(MOD1information$normalised_modularity)  # normalised modularity score for configuration found by MOD1 for MAT
print(MOD1information$realized_modularity)  # realized modularity score for configuration found by MOD1 for MAT

print(MOD1information$RowNodesInModules)  # Shows row nodes per module
print(MOD1information$ColNodesInModules)  # Shows column nodes per module
