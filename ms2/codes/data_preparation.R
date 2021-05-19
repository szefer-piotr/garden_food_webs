# Data preparation script

# 1. Clear the environment, load libraries and datasets
rm(list=ls())

# Arthropod data, no fungicide, no predator
arthropods <- read.table("data/ins_bioOrig.txt")   # arthropods
# plantData <- read.table("data/main_biomass.txt")  # plants
treatments <- read.table("data/treats_clean.txt")  # treatments
arthropodSizes <- read.table("data/treats_clean.txt")  # treatments

library(dplyr)
library(vegan) # for invsimpson



# 1. General descriptors data, grouped by treatment and 
genDescriptors <- arthropods %>%
  group_by(plot, guild, treat) %>%
  summarise(bio = sum(totbio, na.rm = T),
            abu = sum(amount,na.rm = T),
            rich = length(unique(morphotype)),
            ric = n(),
            diva = diversity(amount[!is.na(amount)], 
                             index = "invsimpson"),
            divb = diversity(totbio[!is.na(totbio)], 
                             index = "invsimpson"))


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


# Clean-up

rm(list = ls()[!(ls() %in% c("genDescriptorsPlot"))])
