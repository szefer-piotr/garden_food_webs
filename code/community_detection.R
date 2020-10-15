# Community detection

source("code/contingencyTable.R")
# source("code/bio_log_ratio.R")
source("code/data_processing_code.R")
source("code/weighted-modularity-LPAwbPLUS/code/R/LPA_wb_plus.R")
source("code/weighted-modularity-LPAwbPLUS/code/R/MODULARPLOT.R") #read in plotting function
source("code/weighted-modularity-LPAwbPLUS/code/R/GetModularInformation.R") 
library("bipartite")
source("code/weighted-modularity-LPAwbPLUS/code/R/convert2moduleWeb.R") # read in conversion 

# 1. Load datesets ----
insects <- read.table("datasets/arthropods_clean.txt")
treats  <- read.table("datasets/treatments_clean.txt")
plants  <- read.table("datasets/plants_clean.txt")
size_dat <-read.table("datasets/size_dat_bio.txt")

# From the bio log ratio
# MAT <- subinsct
# gardnetsorig <- gardnets
# gardnets  <- gardnetsfam
# gardnets <- gardnetsorig



# Analysis below should be performed only for C,  I, P plots.
# par(mfrow = c(3,12))
# Modularity plots - seems like 
# for (i in 1:36){
#   MAT <- gardnets[[i]]
#   if(dim(MAT)[1] == 1){next}
#   # MOD1 = LPA_wb_plus(MAT) # find labels and weighted modularity using LPAwb+
#   MOD2 = DIRT_LPA_wb_plus(MAT) # find labels and weighted modularity using DIRTLPAwb+
#   # MOD3 = DIRT_LPA_wb_plus(MAT>0, 2, 20) 
#   patch <- paste("figs/", "gard", i, sep="")
#   tiff(patch) 
#   MODULARPLOT(MAT,MOD2)
#   dev.off()
#   # show the modular network configuration found in MOD1. Row and column numbering indicates the ordering of rows and columns in MAT. Modules are highlighted in red rectangles.
# }

# Before that we need to get rid of intermediate predators
pcplots <- treats[treats$treat %in% c("CONTROL",
                                       "PREDATOR"), ]$codes

pcnets_noip <- gardnets[pcplots]
length(pcnets_noip)

# Remove intermediate predators
for(i in pcplots){
  sn <- pcnets_noip[[i]]
  herbnames <- grep("aran|mant", colnames(sn))
  colnames(sn[,-herbnames])
  pcnets_noip[[i]] <- sn[,-herbnames]
}

# See if it worked
for(gard in pcplots){
  print(colnames(pcnets_noip[[gard]]))
}


for (i in pcplots){
  MAT <- pcnets_noip[[i]]
  if(dim(MAT)[1] <= 1){next}
  # find labels and weighted modularity using DIRTLPAwb+
  MOD2 = DIRT_LPA_wb_plus(MAT)
  
  patch <- paste("figs/", "gard", i, sep="")
  
  treatment <- as.character(treats[treats$codes == i, "treat"])
  
  png(patch, width = 1000,height = 350) 
  # MODULARPLOT(MAT,MOD2)
  # converts the configuration found in MOD1 to a moduleWeb object
  MOD1modWeb = convert2moduleWeb(MAT,MOD2) 
  # plot the corresponding moduleWeb object using plotModuleWeb in library bipartite
  plotModuleWeb(MOD1modWeb, labsize = 0.8)
  text(5,2,treatment)
  dev.off()
}

# Highlight modules
library(RColorBrewer)

mod <- unique(MOD2$Col_labels)
colors <- brewer.pal(length(mod),"BrBG")
barplot(c(1,1,1,1),col = colors)

upperCol <- c(rep("black", dim(MAT)[2]))
lowerCol <- c(rep("black", dim(MAT)[1]))

i = 2
# Find out which are in different modules and color them
modrows <- rownames(MAT)[which(MOD2$Row_labels == mod[i])]
modcols <- colnames(MAT)[which(MOD2$Col_labels == mod[i])]
upperCol[colnames(MAT) %in% modcols]  <- colors[1]
lowerCol[rownames(MAT) %in% modrows]  <- colors[1]

plotweb(MAT, 
        col.high = upperCol, bor.col.high = upperCol,
        col.low = lowerCol, bor.col.low = lowerCol)

#use with R library 'bipartite'




#read in function for finding additional information

MOD1information = GetModularInformation(MAT,MOD2)
# normalised modularity score for configuration found by MOD1 for MAT
print(MOD1information$normalised_modularity)  

# realized modularity score for configuration found by MOD1 for MAT
print(MOD1information$realized_modularity)  

print(MOD1information$RowNodesInModules)  # Shows row nodes per module
print(MOD1information$ColNodesInModules)  # Shows column nodes per module

