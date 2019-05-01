# Data analysis ----

# Direct and indirect effects on plants can be calculated using the log ratio. ln(Vp+/Vp-) where V's are community variables (herbivore abundance and plant biomass in the presence and absence of herbivores). These effect magnitudes can be plotted in relation to each other on an x-y plane and in relation to a 45 deg reference line tha represents the equivalence in strength of direct and indirect effect of carnivores.

# The log ratio effect of carnivores on herbivores should always be negative if carnivores are limiting herbivores regardles of the of the way herbivores are resource limited.
# page 36 Resloving ecosystem complexity

# I need biomass of herbivores, arthropod predators and plant for each plot.

# Inspect that and see if there are any changes in the food web structure. What can i expect?

source("code/contingencyTable.R")
library("bipartite")
library("igraph")

# 1. Load datesets ----

insects <- read.table("datasets/arthropods_clean.txt")
treats  <- read.table("datasets/treatments_clean.txt")
plants  <- read.table("datasets/plants_clean.txt")
size_dat <-read.table("datasets/size_dat_bio.txt")
# sizes   <- read.table("datasets/sizes_clean.txt")

## 2. Log response ratios ----

# Attach biomass measurments to othe main insects dataset
rownames(size_dat) <- size_dat$morph
arthbio  <- size_dat[as.character(insects$morphotype), ]
ins_bio <- cbind(insects, arthbio[, c("morph", "bio")])
ins_bio$totbio <- ins_bio$amount * ins_bio$bio
# ins_bio[1341,]

# Missing stuff, try to measure myself
# ins_bio[!complete.cases(ins_bio),]


# Plots with biomass
pcode <- "w1g5p1" # assign plot name
subinsdat <- ins_bio[ins_bio$plot == pcode,] # get insect biomass
plantcodes <- unique(subinsdat$tree) # see which plants have to be extracted
plantbio <- plants[(plants$CODE == pcode & plants$SP_CODE %in% plantcodes), c("SP_CODE","WEIGHT")] # in KG

rownames(plantbio) <- plantbio[,1]
plantb <- plantbio[,2]
names(plantb) <- plantbio[,1]
subinsct <- contingencyTable2(subinsdat, "tree", "family", "totbio")

subdf <- data.frame()
# Collect data for a given plant within a plot
for(row in 1:nrow(subinsct)){
  plnm <- rownames(subinsct)[row]
  nms <- rownames(as.matrix(subinsct[row,]))
  bio <- as.matrix(subinsct[row,])
  trt <- as.character(treats[treats$codes == pcode, "treat"])
  plbio <- plantbio[plantbio$SP_CODE == plnm, "WEIGHT"]
  ssdf <- data.frame(bio=bio, nms=nms, plnm=plnm,
                      trt=trt,plbio=plbio)
  
  subdf <- rbind(subdf, ssdf)
}




bipartite::plotweb(subinsct,low.abun = plantb,
                   high.abun = colSums(subinsct))



# Vulnerability etc... for fully resolved networks
networklevel(subinsct)

# Motifs
testGraph = barabasi.game(10, 
                          m = 5,
                          power = 2, 
                          out.pref = TRUE,
                          zero.appeal = 0.5,
                          directed = TRUE)

graph.motifs(testGraph, 
               size = 3)
