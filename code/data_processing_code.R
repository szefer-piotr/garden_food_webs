# Script that produces datasets used in the paper:
# 1. Measurments for each morphotype
# 2. Abundance for each family in each garden
# 3. Plant biomass in each plot
# 4. Treatment assignments to the plot codes

# 0. Contingency table function ----

# This function creates a contingency table for a given row category
# column categoory and sums values.
contingencyTable2 <- function(dataset, ROW, COL, VALUE,rm.null=TRUE){
  # Get rid of the empty factors
  dataset[, colnames(dataset) == ROW] <- as.character(dataset[, colnames(dataset) == ROW])
  dataset[, colnames(dataset) == COL] <- as.character(dataset[, colnames(dataset) == COL])
  # Make a table, get rid of the empty rows and columns
  plants <- table(dataset[, colnames(dataset) == ROW], dataset[, colnames(dataset) == COL])
  
  if(rm.null){
    plants <- plants[rowSums(plants) != 0, colSums(plants) != 0]
  }
  
  if(is.null(dim(plants))){plants <- t(plants)}
  # See where to insert values
  allSpecCodes <- colnames(plants)
  allPlotCodes <- rownames(plants)
  entries <- which(plants != 0, arr.ind = TRUE)
  # Loop through the entries and insert values
  for (entry in 1:dim(entries)[1]){
    plot <- entries[entry,1]
    plant <- entries[entry,2]
    specCode <- allSpecCodes[plant]
    plotCode <- allPlotCodes[plot]
    #res <- dataset[dataset$ROW == plotCode & dataset$COL == specCode,VALUE]
    res <- dataset[dataset[,ROW] == plotCode & dataset[,COL] == specCode,VALUE]
    # print(sum(res))
    plants[plot,plant] <- sum(res, na.rm = TRUE)
  }
  plants[is.na(plants)] <- 0
  
  # Change the table to a matrix or data.frame
  
  mat_a <- matrix(0, nrow = dim(plants)[1], ncol = dim(plants)[2])
  colnames(mat_a) <- colnames(plants)
  rownames(mat_a) <- rownames(plants)
  for (row in 1:dim(plants)[1]){
    for (col in 1:dim(plants)[2])
      mat_a[row,col] <- plants[row,col]
  }
  
  #return(plants)
  return(mat_a)
}

# 1. Load, fix, and save datasets ----

# Measurements: insect body length ----
measur <- read.csv("datasets/csv_measurments_all.csv")
# Count data: incidence of insects on individual plant species ----
arthro <- read.csv("datasets/csv_wng_all.csv")

# Extract family abbreviation
arthro$family <- substr(arthro$morphotype, 1, 4)
measur$family <- substr(measur$morphotype, 1, 4)

# table(arthro$family) # Number of interactions within each family
# table(measur$family) # Number of measured individuals from each family

# Fixing some of the entries
measur$Size <- as.numeric(measur$Size)
names(measur) <- c("morphotype","no","size","scale","scl","notes", "group")
measur$morphotype <- as.character(measur$morphotype)

# Check and unify plant names
longnames <- which(sapply(as.character(measur$morphotype), nchar)>7)
chlist <- strsplit(as.character(measur[longnames, ]$morphotype), split=",")

for (i in 1:length(longnames)){
  row <- longnames[i]
  measur[row, 1] <- as.character(chlist[[i]][1])
  measur[row, 2] <- as.numeric(chlist[[i]][2])
  measur[row, 3] <- as.numeric(chlist[[i]][3])
  measur[row, 5] <- as.numeric(chlist[[i]][6])
}

# measur <- measur[,c(1,2,3,5,6)]
# NOTE: Scale codes are as follows
# for MANT, ARAN, HEMI, HOMO: 1 cm = 10 [mm], 0.5 cm = 20 [mm], cm = 1
# for COLE, ORTH and LEPI check the data again!!!
# For cole scale 2 = 20, 1 = 10, 

# Transform the scale factors
ent <- as.character(unique(measur$scale))
measur$nscl <- 0
measur[measur$scale %in% ent[c(2,5)], ]$nscl <- 10
measur[(measur$scale %in% ent[c(1,3,4,7)] & measur$scl == 0.5),]$nscl <- 20
measur[(measur$scale %in% ent[c(1,3,6,7)] & measur$scl == 1),]$nscl <- 1
measur[(measur$scale %in% ent[c(1,3,6,7)] & measur$scl == 2),]$nscl <- 2

# Some corrections for missing scale factors
measur[measur$morphotype == "cole065", ]$nscl <- 20
measur[2108, ]$nscl <- 10

# Is everything ok with coleoptera?
measur[measur$group == "cole",] #instead of 2 and 1 there are 1.0 and 2.0

# Real scale [cm] insect sizes
measur$rsize <- measur$size/measur$nscl

# Write a clean measurments table
# write.table(measur, "datasets/wng_measurements.txt")


# 2. Arthropod dataset cleaning ----

sort(unique(arthro$tree))

# Remove Sida rhombifolia - this is not a tree
arthro <- arthro[arthro$tree != "sida",]

# Change the names of trees
tochange <- data.frame(a = sort(unique(arthro$tree)),
           b=sort(unique(arthro$tree)))

# himibi??? "mimodi"
# arthro[arthro$tree == "premna",]
arthro[arthro$tree == "costsp",]
# arthro[arthro$plot == "w1g1p2",]
# 
# premna
tochange$b[1] <- "breyce"
tochange$b[5] <- "cordte"
tochange$b[10] <- "ficuco"
tochange$b[12] <- "mimodi"
tochange$b[14] <- "homano"
tochange$b[15] <- "homano"
tochange$b[25] <- "pipeum"
tochange$b[29] <- "prems1"
tochange$b[32] <- "solatu"
tochange$b[36] <- "tremor"
tochange$b[41] <- "viteco"

arthro$tree <- as.character(arthro$tree)

# Then use the corrected names to change names in the dataset
for(name in tochange$a){
  print(name)
  arthro[arthro$tree == name, ]$tree <- as.character(tochange[tochange$a == name, ]$b)
}

arthro[arthro$plot == "wg3p6",]$plot <- "w1g3p6"
arthro$plot <- as.character(arthro$plot)

# Save the corrected dataset (clean one)
# write.table(arthro, "datasets/wng_arthro_clean.txt")

# 3. Biomass and treatments ----
main <- read.table("datasets/wng_main_clean.txt", header = T)
main$TREAT <- as.character(main$TREAT)
treats <- as.data.frame(tapply(main$TREAT, main$CODE, unique))
treats$code <- rownames(treats)
names(treats) <- c("treat", "codes")
treats$codes <- gsub("W","W1", treats$codes)
treats$codes <- tolower(treats$codes)

# write.table(treats, "datasets/treats_clean.txt")

main$CODE <- gsub("W","W1", main$CODE)
main$CODE <- tolower(main$CODE)

main_biomass <- main[,c("CODE","PLOT","BLOCK","TREAT","SPEC","SP_CODE","LIFE.FORM","BASAL_A","HEIGHT_M", "LEAVES","TRUNK","WEIGHT")]

# write.table(main_biomass, "datasets/wng_main_bio.txt")

# Food webs

# source("code/contingencyTable.R")
library("bipartite")
library("igraph")

# 1. Load datesets ----

insects <- read.table("datasets/arthropods_clean.txt")
treats  <- read.table("datasets/treatments_clean.txt")
plants  <- read.table("datasets/plants_clean.txt")
size_dat <-read.table("datasets/size_dat_bio.txt")
# sizes   <- read.table("datasets/sizes_clean.txt")

# Attach biomass measurments to the main insects dataset
rownames(size_dat) <- size_dat$morph
arthbio  <- size_dat[as.character(insects$morphotype), ]
ins_bio <- cbind(insects, arthbio[, c("morph", "bio")])
ins_bio$totbio <- ins_bio$amount * ins_bio$bio
# ins_bio[1341,]

# Missing stuff, try to measure myself
# ins_bio[!complete.cases(ins_bio),]

# Insect abundances
# tapply(ins_bio$amount,ins_bio$plot, sum)
# ins_bio[ins_bio$plot == pcode,]

# 1.1 Biomass based networks ----

# Biomasses of insects per plot per plant
biofulldf <- data.frame()
gardnets <- list()
gardnetsfam <- list()

# Plots with biomass
for(pcode in as.character(treats$codes)){
  print(pcode)
  subinsdat <- ins_bio[ins_bio$plot == pcode,] # get insect biomass
  plantcodes <- unique(subinsdat$tree) # see which plants have to be extracted
  plantbio <- plants[(plants$CODE == pcode & plants$SP_CODE %in% plantcodes), c("SP_CODE","WEIGHT")] # in KG
  
  plantbio$SP_CODE <- as.character(plantbio$SP_CODE)
  cumWeight <- tapply(plantbio$WEIGHT,
                      plantbio$SP_CODE, 
                      sum)
  
  pbio <- data.frame(SP_CODE = rownames(cumWeight),
                     WEIGHT = cumWeight)
  
  rownames(pbio) <- pbio[,1]
  plantb <- pbio[,2]
  names(plantb) <- pbio[,1]
  subinsct <- contingencyTable2(subinsdat, "tree", "family", "totbio",FALSE)
  listnet <- list(subinsct)
  names(listnet) <- pcode
  gardnetsfam <- append(gardnetsfam, listnet) 
  
  # Add to the list
  # By family or by species
  subinsctsp <- contingencyTable2(subinsdat, "tree", "morph", "totbio",FALSE)
  listnet <- list(subinsctsp)
  names(listnet) <- pcode
  gardnets <- append(gardnets, listnet) 
  
  subdf <- data.frame()
  # Collect data for a given plant within a plot
  for(row in 1:nrow(subinsct)){
    plnm <- rownames(subinsct)[row]
    nms <- rownames(as.matrix(subinsct[row,]))
    if(is.null(nms)){nms <- colnames(subinsct)}
    bio <- as.matrix(subinsct[row,])
    trt <- as.character(treats[treats$codes == pcode, "treat"])
    plbio <- pbio[pbio$SP_CODE == plnm, "WEIGHT"]
    ssdf <- data.frame(plot=pcode,bio=bio, nms=nms, plnm=plnm,
                       trt=trt,plbio=plbio)
    
    subdf <- rbind(subdf, ssdf)
  }
  biofulldf <- rbind(biofulldf, subdf)
}
# pcode <- "w1g5p1" # assign plot name
# biofulldf   # dataframe
# gardnets    # morphotype based networks for each garden
# gardnetsfam # family aagregated networks for all garden


# 1.2 Abundance (number of individuals) based networks ----
# Abundances of insects per plot per plant
abufulldf <- data.frame()
abugardnets <- list()
abugardnetsfam <- list()
# Plots with biomass
for(pcode in as.character(treats$codes)){
  print(pcode)
  subinsdat <- ins_bio[ins_bio$plot == pcode,] # get insect biomass
  plantcodes <- unique(subinsdat$tree) # see which plants have to be extracted
  plantbio <- plants[(plants$CODE == pcode & plants$SP_CODE %in% plantcodes), c("SP_CODE","WEIGHT")] # in KG
  
  plantbio$SP_CODE <- as.character(plantbio$SP_CODE)
  cumWeight <- tapply(plantbio$WEIGHT,
                      plantbio$SP_CODE, 
                      sum)
  
  pbio <- data.frame(SP_CODE = rownames(cumWeight),
                     WEIGHT = cumWeight)
  
  rownames(pbio) <- pbio[,1]
  plantb <- pbio[,2]
  names(plantb) <- pbio[,1]
  subinsct <- contingencyTable2(subinsdat, "tree", "family", "amount",FALSE)
  listnet <- list(subinsct)
  names(listnet) <- pcode
  abugardnetsfam <- append(abugardnetsfam, listnet) 
  
  # Add to the list
  # By family or by species
  subinsctsp <- contingencyTable2(subinsdat, "tree", "morph", "amount",FALSE)
  listnet <- list(subinsctsp)
  names(listnet) <- pcode
  abugardnets <- append(abugardnets, listnet) 
  
  subdf <- data.frame()
  # Collect data for a given plant within a plot
  for(row in 1:nrow(subinsct)){
    plnm <- rownames(subinsct)[row]
    nms <- rownames(as.matrix(subinsct[row,]))
    if(is.null(nms)){nms <- colnames(subinsct)}
    bio <- as.matrix(subinsct[row,])
    trt <- as.character(treats[treats$codes == pcode, "treat"])
    plbio <- pbio[pbio$SP_CODE == plnm, "WEIGHT"]
    ssdf <- data.frame(plot=pcode,bio=bio, nms=nms, plnm=plnm,
                       trt=trt,plbio=plbio)
    
    subdf <- rbind(subdf, ssdf)
  }
  abufulldf <- rbind(abufulldf, subdf)
}


#### bio_log_ratio CUTTED line 22
# Dataset containing biomasses for the log ratio comparisons between predator exclosures and control plots
biollcp <- biofulldf[biofulldf$trt %in% c("CONTROL", "PREDATOR"),]
biollcp$plot <- as.character(biollcp$plot)
biollcp$plnm <- as.character(biollcp$plnm)
biollcp$trt <- as.character(biollcp$trt)
biollcp$gard <- substr(biollcp$plot, 3,4)
# see which species are present in both treatment plots
# table(biollcp$trt, biollcp$plnm) 
