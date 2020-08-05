# How different are sites.
rm(list=ls())
source("code/data_processing_code.R")

rankAbuData <- function(data, block, treatment, woody = T){
  # Speciefic function which plots rank abundance 
  # plot for given garden and treatment
  
  print(paste(block, treatment, sep = " "))
  
  ca <- data$BLOCK == block
  cb <- data$TREAT == treatment
  cc <- data$LIFE.FORM %in% c("tree", "shrub")

  if(woody){conditions <- (ca&cb&cc)
  } else {conditions <- (ca&cb)}
  
  sub_data <- data[conditions,
                   c("SP_CODE", "WEIGHT")]
  
  return(sub_data)
  
}

data <- main_biomass
blocks <- unique(data$BLOCK)
treatments <- unique(data$TREAT)

rankAbuData(data, blocks[1], treatments[4], woody=T)

# Dissimilarity matrix for woody plants only
trt <- treatments[4] 
main_sub_treat <- main_biomass[main_biomass$TREAT == trt, ]
main_sub_treat_woody <- main_sub_treat[main_sub_treat$LIFE.FORM %in% c("shrub","tree"),]
ct_trt <- contingencyTable2(main_sub_treat_woody, "CODE", "SP_CODE", "WEIGHT")
vegdist(ct_trt)

#w1g1p1 and w1g4p6

treats[treats$codes == "w1g1p1", ]

rankAbuData(data, 'WG1', "FUNGICIDE")
rankAbuData(data, 'WG4', "FUNGICIDE")

# Variability between blocks is high.

# Between treatment and control within the block is low

compTreatToControl <- function(trt){
  dist_vals <- data.frame()
  ct_comparisons <- list()
  for(bck in blocks){
    
    # Rank abundance for two treatments within given block
    radtrt <- rankAbuData(data, bck, trt)
    # For a given block which
    radctr <- rankAbuData(data, bck, "CONTROL")
    
    # Merge some repeating species
    t_vals <- tapply(radtrt$WEIGHT, radtrt$SP_CODE, sum, na.rm=T)
    c_vals <- tapply(radctr$WEIGHT, radctr$SP_CODE, sum, na.rm=T)
    corrt <- data.frame(WEIGHT = t_vals,
                        SP_CODE = names(t_vals))
    corrc <- data.frame(WEIGHT = c_vals,
                        SP_CODE = names(c_vals))
    corrt <- corrt[complete.cases(corrt),]
    corrc <- corrc[complete.cases(corrc),]
    
    
    mdset <- merge(corrt,corrc, by="SP_CODE", all =  T)
    distset <- mdset[,-1]
    rownames(distset) <- mdset[,1]
    dv <- as.numeric(vegdist(t(distset),na.rm = T))
    dist_vals <- rbind(dist_vals, 
                       data.frame(TREAT = trt,
                                  BLOCK = bck,
                                  DIST = dv))
    ct_comparisons[[bck]] <- mdset
  }
  return(list(distances = dist_vals, 
              barplots = ct_comparisons))
}

# Dissimilarity between treatment and control plots within block
allTrtComp <- data.frame()
for(trt in treatments){
  allTrtComp <- rbind(allTrtComp,
                      compTreatToControl(trt)$distances)
}

# Very low value of dissimilarity is caused by a strong dominance of 
# Pipturus argenteum
# INSECTICIDE   WG3 0.04967525


td <- rankAbuData(data, "WG3", "PREDATOR")[,2]
names(td) <- rankAbuData(data, "WG3", "PREDATOR")[,1]

# Barplots
sel_treat <- "PREDATOR"

X11(45,10)
par(mfrow = c(3,2))
for(block_no in 1:6){
  comparisons <- compTreatToControl(sel_treat)
  block_names <- names(comparisons$barplots)
  my.df <- comparisons$barplots[[block_no]]
  
  names(my.df)<- c("species", "treatment", "control")
  
  barplot(t(as.matrix(log(my.df[, 2:3]))), 
          beside = TRUE,
          names.arg = my.df$species,
          legend.text = FALSE,
          ylab = "Biomass [kg]",
          xlab = "",
          las=2,
          main = paste(sel_treat, block_names[block_no], sep = " "))
}

main_biomass[main_biomass$TREAT == sel_treat & main_biomass$BLOCK == block_names[block_no], ]
main_biomass[main_biomass$SP_CODE == "CALOMU", ]
