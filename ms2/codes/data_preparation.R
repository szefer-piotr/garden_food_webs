# Data preparation script

# 1. Clear the environment, load libraries and datasets
rm(list=ls())

# Arthropod data, no fungicide, no predator
arthropods <- read.table("data/ins_bioOrig.txt")   # arthropods
# plantData <- read.table("data/main_biomass.txt")  # plants
treatments <- read.table("data/treats_clean.txt")  # treatments
arthropodSizes <- read.table("data/size_dat_bio.txt")  # treatments
leafAreaOrig <- read.table("data/wng_main_clean.txt")

library(dplyr)
library(vegan) # for invsimpson
library(data.table)

# 1.1Density

# 1.1.1 Caluclate leaf area per individual tree
leafArea <- leafAreaOrig %>%
  filter(LIFE.FORM %in% c("shrub", "tree")) %>%
  dplyr::select(CODE, PLOT, BLOCK, TREAT, SP_CODE,
                LDMC,WATER,HERB,LEAVES,SLA)%>%
  mutate(area.m2 = ((LEAVES*1000)*SLA)/10000)

# 1.1.2 Calculate leaf area per individual plot
leafAreaPlot <- leafArea %>%
  group_by(CODE) %>%
  summarise(total.leaf = sum(area.m2, na.rm = T))

# 1.1.3 Append area to each row in the arthropod dataset
arthropods <- arthropods %>%
  mutate(ch.plot = as.character(plot),
         ch.tree = as.character(tree),
         leaf.area.m2 = NA)


# 1.2 General descriptors data, grouped by treatment and 
genDescriptors <- arthropods %>%
  group_by(plot, guild, treat) %>%
  summarise(bio = sum(totbio, na.rm = T),
            abu = sum(amount,na.rm = T),
            rich = length(unique(morphotype)),
            ric = n(), # what is this?
            diva = diversity(amount[!is.na(amount)], 
                             index = "invsimpson"),
            divb = diversity(totbio[!is.na(totbio)], 
                             index = "shannon"))


genDescriptorsPlot <- rbind(genDescriptors[1:3],
                            genDescriptors[1:3],
                            genDescriptors[1:3],
                            genDescriptors[1:3],
                            genDescriptors[1:3],
                            genDescriptors[1:3]
                            )

genDescriptorsPlot$ind <- stack(genDescriptors)[,2]
genDescriptorsPlot$val <- stack(genDescriptors)[,1]

genDescriptorsPlot$treat <- factor(genDescriptorsPlot$treat,
                                   levels = c("insecticide",
                                              "control",
                                              "weevil25",
                                              "weevil125"))

genDescriptorsPlot <- genDescriptorsPlot %>%
  mutate(block = substr(plot, 3,4)) %>%
  mutate(treat = as.character(treat))

# genDescriptorsPlot <- as.data.frame(genDescriptorsPlot)

# 3. Data on log ratio between W125/C, W25/C, C/I
# logRatioOrders <- arthropods %>% 
#   summarise()


# Notes ------

# This adds guild identity to the raw arthropods data
# arthropods$guild <- "herbivore"
# arthropods[arthropods$family %in% c("mant", "aran"), ]$guild <- "art_pred"

# This adds treatment identity to arthropod data
# rownames(treatments) <- treatments$codes
# arthropods$treat <- tolower(treatments[as.character(arthropods$plot),]$treat)

# Remove predator and fungicide plots
# treats_to_remove <- c("PREDATOR", "FUNGICIDE")
# arthropods <- arthropods %>%
#   filter(plot %in% treatments[!(treatments$treat %in% treats_to_remove) , ]$codes)

# This adds the plant biomass colum in 'arthropods' datasets
# # Log ratios of plants and 
# plantBio <- function(species, plot){
#   plbiomass <- main_biomass[main_biomass$CODE == plot, ]
#   plsitebio <- plbiomass[plbiomass$SP_CODE == toupper(species), ]
#   bio <- (plsitebio$WEIGHT)
#   if(length(bio) != 0){
#     print(paste("Weight of:",species, "at the plot", plot, "is", bio))
#     return(sum(bio))
#   }else{
#     print(paste("There is no",species, "at the plot", plot))
#     return(0)
#   }
# }
# 
# arthropods$plant_bio <- 0
# head(arthropods)
# 
# for(row in 1:dim(arthropods)[1]){
#   arthropods[row, ]$plant_bio <- plantBio(arthropods[row, ]$tree,
#                                           arthropods[row, ]$plot)
# }
# 
# write.table(arthropods, "data/ins_bioOrig.txt")

# Add density
leafAreaPlot <- leafAreaPlot %>%
  mutate(CODE = tolower(CODE))

abu <- genDescriptorsPlot %>%
  filter(ind == "abu")

# go through the 'abu' plot column and extract the leaf area for that plot from the leafAreaPlot to calculate densities.
# rw = 1

densDf <- data.frame()

for (rw in 1:nrow(abu)){
  subrw <- abu[rw,]
  
  ###### HERE NEEDS TO WORK
  areal = tolower(leafAreaPlot$CODE) %like% 
    substr(as.character(subrw$plot), 3,6) # extract the area
  
  area = leafAreaPlot[areal, ]$total.leaf
  
  gDnR <- data.frame(
    plot = subrw$plot,
    guild = subrw$guild,
    treat = subrw$treat,
    ind = "dens",
    val = subrw$val/area, #ind/m2
    block = subrw$block)
  
  densDf <- rbind(densDf, gDnR)
  
}

genDescriptorsPlotDens <- rbind(as.data.frame(genDescriptorsPlot),
                                densDf)

# Clean-up
genDescriptorsPlotDens$treat <- as.factor(genDescriptorsPlotDens$treat)
genDescriptorsPlotDens$treat <- factor(genDescriptorsPlotDens$treat, 
                                   levels = c("insecticide", 
                                              "control", 
                                              "weevil25",
                                              "weevil125"),
                                   labels = c("I","C","H1","H2"))
rm(list = ls()[!(ls() %in% c("genDescriptorsPlotDens"))])
