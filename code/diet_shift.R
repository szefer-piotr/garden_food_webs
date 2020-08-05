# diet shift ordination

# rm(list = ls())
# insects <- read.table("datasets/arthropods_clean.txt")
# treats  <- read.table("datasets/treatments_clean.txt")
# plants  <- read.table("datasets/plants_clean.txt")
# size_dat <-read.table("datasets/size_dat_bio.txt")

library("bipartite")
library("vegan")
# source("code/bio_log_ratio.R")
source("code/data_processing_code.R")


# This shift maybe should be beased on abundance rather than like in the plot
ins_bio$morphplot <- paste(ins_bio$family, ins_bio$plot, sep="")
ins_bio$tree <- as.character(ins_bio$tree)
abuinsmat <- contingencyTable2(ins_bio, "tree", "morphplot", "amount" , FALSE)
#dist(x, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)

# For the whole dataset
# abunmds <- metaMDS(abuinsmat, distance = "bray")
# plot(abunmds, display = "species")
# text(abunmds, display = "species")
# abuvd <- as.matrix(vegdist(t(abuinsmat), method = "bray"))

# Maybe this is not the best representation
tp <- as.character(treats[treats$treat %in% c("PREDATOR"), ]$codes)
tc <- as.character(treats[treats$treat %in% c("CONTROL") , ]$codes)

# abuvd[substr(rownames(abuvd), 5,10) %in% tp, 
#       substr(colnames(abuvd), 5,10) %in% tc]

# DIET SWITCHING ----
# Herbivores present in P and C treatments
predsites <- treats[treats$treat == "PREDATOR",]$codes
contsites <- treats[treats$treat == "CONTROL",]$codes

ips <- grep("aran|mant", ins_bio$morphotype)

ins_bioOrig <- ins_bioins_bio <- ins_bio[-ips, ]
ins_bio <- ins_bio[-ips, ]

# Contingency tables for predator and control sites
pabumat <- contingencyTable2(ins_bio[(ins_bio$plot %in% predsites), ],
                             "plot","morphotype","amount")
pbiomat <- contingencyTable2(ins_bio[ins_bio$plot %in% predsites, ],
                             "plot","morphotype","totbio")

cabumat <- contingencyTable2(ins_bio[ins_bio$plot %in% contsites, ],
                             "plot","morphotype","amount")
cbiomat <- contingencyTable2(ins_bio[ins_bio$plot %in% contsites, ],
                             "plot","morphotype","totbio")

# Species present in both predator and control sites
comparable <- colnames(cbiomat)[colnames(cbiomat) %in% colnames(pbiomat)]

# Remove intermediate predators
# comparable

# Food plants for comparable herbivores
# treats_trimmed <- treats[!(treats$treat %in% c("FUNGICIDE", "INSECTICIDE",
#                                                "WEEVIL25",  "WEEVIL125")), ]
# treats_trimmed$codes <- as.character(treats_trimmed$codes)
# 
# cp_treats <- treats_trimmed[treats_trimmed$treat %in% c("CONTROL","PREDATOR"),]$codes
# csites <- treats_trimmed[treats_trimmed$treat %in% c("CONTROL"),]$codes
# psites <- treats_trimmed[treats_trimmed$treat %in% c("PREDATOR"),]$codes

csites <- as.character(contsites)
psites <- as.character(predsites)

# Filter only control and predator sites
ins_bio_cp <- ins_bio[ins_bio$plot %in% c(csites, psites), ]
ins_bio_cp_comparable <- ins_bio_cp[ins_bio_cp$morphotype %in% comparable,]
ibc <- ins_bio_cp_comparable
ibc <- ibc[complete.cases(ibc),]

ccompFood <- contingencyTable2(ibc[ibc$plot %in% csites, ],
                               "tree",
                               "morphotype",
                               "totbio")

pcompFood <- contingencyTable2(ibc[ibc$plot %in% psites, ],
                               "tree",
                               "morphotype",
                               "totbio")

dim(pcompFood)
dim(ccompFood)

# Combine dataset 
rownames(pcompFood) <- paste("p", rownames(pcompFood), sep="_")
rownames(ccompFood) <- paste("c", rownames(ccompFood), sep="_")

compFood <- rbind(pcompFood,ccompFood)

# ALL GOOD

# compvec <- compFood[,'cole080']
compvec <- compFood[, 'orth019']# Shift of the diet histograms
compare_row <- function(row){
  compvec <- compFood[,row]
  pred <- grep("p_", names(compvec))
  cont <- grep("c_", names(compvec))
  # Predator treatment will always have more woody plants
  # However, I need to fix the names.
  predvec <- compvec[pred]
  contvec <- compvec[cont]
  # Remove indicator
  names(predvec) <- substr(names(predvec), 3, 8)
  names(contvec) <- substr(names(contvec), 3, 8)
  
  allnames <- c(names(predvec),names(contvec))
  # Merge two vectors but need to take into account plant names 
  # that are unique for BOTH
  compDat <- data.frame(names = unique(allnames))
  
  ##### ISSUE 1
  compDat$pred <- 0
  compDat$cont <- 0
  
  for(i in compDat$names){
    print(i)
    if(i %in% names(predvec)){
      compDat[compDat$names == i,  ]$pred <- predvec[i]
    }
    if(i %in% names(contvec)){
      compDat[compDat$names == i,  ]$cont <- contvec[i]
    }
  }
  # compDat$pred <- predvec[compDat$names]
  # compDat$cont <- contvec[compDat$names]
  
  # Here is another problem. What about palnts thet were, not
  # present in the control plot but were in the predator exclosure.
  # If they were not present, but favored would that count as a shift?
  # I will initially make these values zeros.
  compDat[is.na(compDat)] <- 0
  
  # Caluculate proportions for the diets
  compDatStand <- as.data.frame(apply(compDat[,-1], 2, function(x){x/sum(x)}))
  # I calculate differences between proportions in utilizing different
  # plant species. I use absolute values. This way values would vary between
  # 0 (no shift) and 2 - compelete shift. Therefore I could divide by 2 to 
  # standardize
  compDatStand$diff <- abs(compDatStand$pred - compDatStand$cont)
  shift <- sum(compDatStand$diff)/2
  return(shift)
}

shiftDf <- data.frame()
for (row in 1:dim(compFood)[2]){
  print(row)
  shiftSpec <- data.frame(species = colnames(compFood)[row],
                          shift = compare_row(row))
  shiftDf <- rbind(shiftDf, shiftSpec)
}

hist(shiftDf$shift, breaks = 50)

# Solving issue 1

# See if it does what it says
# compFood[, 'orth063'] # yes it does
# Where NaNs came up these are zeros...needs to be removed from the dataset
# compFood[, 'homo003'] 
# This species has a log ratio, but no shif in the diet... why?
# pts <- as.character(insects[insects$morphotype == "cole080",]$plot)
# treats[treats$codes %in% pts, ]
# The insects dataset is ok
# insects[insects$morphotype == "cole080",]
# what about ins_bio
# ins_bio[ins_bio$morphotype == "cole080",]
# Filtered data is also ok
# ins_bio_cp[ins_bio_cp$morphotype == "cole080",]
# ibc[ibc$morphotype == "cole080",]
# ccompFood[, 'cole080']
# pcompFood[, 'cole080']
# # Also ok
# compFood[, 'cole080']

