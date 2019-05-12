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

# Insect abundances
# tapply(ins_bio$amount,ins_bio$plot, sum)
# ins_bio[ins_bio$plot == pcode,]

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

# Dataset containing biomasses for the log ratio comparisons
biollcp <- biofulldf[biofulldf$trt %in% c("CONTROL", "PREDATOR"),]
biollcp$plot <- as.character(biollcp$plot)
biollcp$plnm <- as.character(biollcp$plnm)
biollcp$trt <- as.character(biollcp$trt)
table(biollcp$trt, biollcp$plnm)

# Herbivore log-ratio for each plant species and each group of insects
# within garden

# Assume that each plant hosts unique community of insects.

# Log ratio analyses

biollcp$gard <- substr(biollcp$plot, 3,4)
# 
# block = "g1"
# plnt = unique(subbl$plnm)[1]
# fam = unique(subblpl$nms)[1]
# 
logratiodf <- data.frame()
# Go through each block
for(block in unique(biollcp$gard)){
  subbl <- biollcp[biollcp$gard == block, ]
  # within each block obtain community for each plant species
  for(plnt in unique(subbl$plnm)){
    subblpl <- subbl[subbl$plnm == plnt, ]
    subblpl$nms <- as.character(subblpl$nms)
    subpllr <- tapply(subblpl$plbio, subblpl$trt, mean)
    pltlr <- log(subpllr["CONTROL"]/subpllr["PREDATOR"])
    if(is.na(pltlr)){next}
    arthrodf <- data.frame()
    print(plnt)
    for(fam in unique(subblpl$nms)){
      famsub <- subblpl[subblpl$nms == fam, c("bio","trt")]
      cont <- famsub[famsub$trt == "CONTROL",]$bio
      pred <- famsub[famsub$trt == "PREDATOR",]$bio
      famlr <- log(cont/pred)
      if(length(famlr) == 0){next}
      print(fam)
      print(famlr)
      arthrodf<- rbind(arthrodf, data.frame(plnt=plnt,fam=fam, lr=famlr,
                                            pltlr=pltlr, gard=block))
    }
  }
  logratiodf <- rbind(logratiodf, arthrodf)
}

logratiodf
logratiodf <- logratiodf[!(is.nan(logratiodf$lr) | is.infinite(logratiodf$lr)), ]

# Biomasses (from the control plots) for each fam and tree
logratiodf$gard <- as.character(logratiodf$gard)
logratiodf$plnt <- as.character(logratiodf$plnt)
logratiodf$fam <- as.character(logratiodf$fam)

biosize <- data.frame()
for(row in 1:dim(logratiodf)[1]){
  loc <- c(logratiodf[row,]$gard,
           logratiodf[row,]$plnt,
           logratiodf[row,]$fam)
  bcpc <- biollcp[biollcp$trt == "CONTROL",]
  bio <- bcpc[(bcpc$gard == loc[1] & bcpc$plnm == loc[2] & 
               bcpc$nms == loc[3]),]$bio
  biosize <- rbind(biosize, data.frame(garrd = loc[1],
                                       plnt = loc[2],
                                       fam = loc[3], 
                                       bio = bio))
}

biosize
logratiodf$size <- biosize$bio

library(ggplot2)
p <- ggplot(logratiodf, 
            aes(x = lr, y = pltlr, label = fam, color=plnt))
p + coord_cartesian(xlim=c(-5,5), ylim = c(-5,5)) +
  geom_jitter(size=logratiodf$size, width=0.5) +  
  geom_text() + 
  geom_abline(slope = 1, intercept = 0, linetype="dashed") +
  geom_abline(slope = -1, intercept = 0,linetype="dashed") +
  xlab("Direct effect of predator removal on insects") + 
  ylab("Indirect effect of predator removal on plants")
  
# I would like to also indicate the abundance in general of these insects on the plot to show their relative importance for the general result.



# plot(pltlr~lr, data=logratiodf, col=logratiodf$fam, pch=19)
# abline(0,1, lty=2)
# abline(0,-1, lty=2)
# abline(h=0, lty=1)
# abline(v=0, lty=1)

# 
# bipartite::plotweb(subinsct,low.abun = plantb,
#                    high.abun = colSums(subinsct))



# Vulnerability etc... for fully resolved networks
# networklevel(subinsct)

# Motifs
# testGraph = barabasi.game(10, 
#                           m = 5,
#                           power = 2, 
#                           out.pref = TRUE,
#                           zero.appeal = 0.5,
#                           directed = TRUE)
# 
# graph.motifs(testGraph, 
#                size = 3)
