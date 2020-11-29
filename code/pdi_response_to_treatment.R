# Individual species PDI in response to treatments

rm(list=ls())
source("code/data_processing_code.R")
source("code/contingencyTable.R")
library(ggplot2)
library(emmeans)
library(multcomp)
library(vegan)
library(bipartite)
library(lme4)
library(lmerTest)


# gardnets 

# Random networks for the relative comparison
# vaznull(10, abugardnets[[3]])

# within garden 

treats$block <- substr(treats$codes, 3,4)

# Initiate an empty data frame
pdi_changedf <- data.frame()

comps <- list(toupper(c("control", "predator")),
              toupper(c("control", "weevil25")),
              toupper(c("control", "weevil125")),
              toupper(c("control", "insecticide")))

comps <- list(toupper(c("control", "predator")))

for(num in 1:length(comps)){
  comparison <- comps[[num]]
  print(comparison)
  for(garden in unique(treats$block)){
    
    print(garden)
    
    ctr_plot <- treats[treats$block == garden & treats$treat == comparison[1], ]$code
    trt_plot <- treats[treats$block == garden & treats$treat == comparison[2], ]$code
    
    ctr_net <- abugardnets[[ctr_plot]]
    trt_net <- abugardnets[[trt_plot]]
    
    # print("ASS")
    c1 <- is.null(dim(ctr_net))
    c2 <- is.null(dim(trt_net))
    
    if(c1 | c2){
      print("One of two is empty")
      next
    }
    
    # Remove one row matrices
    if (dim(ctr_net)[1] == 1 | dim(trt_net)[1] == 1){
      print("One or two matrices have only one row")
      next
    }
    
    print("All good")
    
    ctr_net <- ctr_net[,-grep(c("mant|aran"),
                              colnames(ctr_net))]
    trt_net <- trt_net[,-grep(c("mant|aran"),
                              colnames(trt_net))]
    
    # Get species in both treatments
    names_in_both <- colnames(ctr_net)[colnames(ctr_net) %in% colnames(trt_net)]
    ctr_net_f <- ctr_net[,names_in_both]
    trt_net_f <- trt_net[,names_in_both]
    
    # Calculate PDI's
    
    pdictr <- specieslevel(ctr_net_f, index = "PDI")$`higher level`
    pditrt <- specieslevel(trt_net_f, index = "PDI")$`higher level`
    
    log_ratio <- c(log(pditrt/pdictr)$PDI)
    
    spec_abu_ctr <- colSums(ctr_net_f)
    spec_abu_trt <- colSums(trt_net_f)
      
    spec_abu_both <- colSums(rbind(ctr_net_f,
                                   trt_net_f))
      
    gdf <- data.frame(species = rownames(pdictr),
                      vals = c(pdictr$PDI,pditrt$PDI),
                      trt = rep(c(comparison[1],
                                  comparison[2]),
                                each=length(pdictr$PDI)),
                      dupl_lr = rep(log_ratio,2),
                      abu_cont = spec_abu_ctr,
                      abu_tret = spec_abu_trt,
                      abu = rep(spec_abu_both,2),
                      garden = garden,
                      comp = paste(comparison, collapse = " vs. "))
    
    pdi_changedf <- rbind(pdi_changedf, gdf)
    
  }
  
}

pdi_changedf$trt <- as.character(pdi_changedf$trt)

# Remove zero values
pdi_change_nozero <- pdi_changedf[pdi_changedf$vals != 0, ]
pdi_nocole <- pdi_change_nozero[!(pdi_change_nozero$species %in% "cole001"), ]
pdi_only_cole <- pdi_change_nozero[(pdi_change_nozero$species %in% "cole001"), ]

ggplot(pdi_only_cole , aes(x=trt, y = vals, group = garden))+
  geom_jitter(width = 0.1)+
  geom_line(lty = 2, lwd=0.1)+
  facet_wrap(~comp, scales = "free")+
  ylab("Paired Distance Index (specialization")+
  xlab("Treatment")

# Compare log ratios and test wether their mean value is different from 0.

# Select a treatment
pdi_change_nozero

# Filter species which abundance is above some treshold
treshold <- 5 # number of individuals

cond1 <- pdi_change_nozero$abu_cont >= 5
cond2 <- pdi_change_nozero$abu_tret >= 5
pdi_filtered <- pdi_change_nozero[cond1 & cond2, ]

pdi_filtered$sp_gard <- paste(pdi_filtered$species, pdi_filtered$garden, sep = "_")
pdi_filtered$fam <- substr(pdi_filtered$species,1,4)

# There is some relationship between plant sp richness in C vs I comparison.
ggplot(pdi_filtered , aes(x=trt, y = vals, 
                          group = sp_gard, 
                          colour = fam))+
  geom_jitter(width = 0.05, 
              size = log(pdi_filtered$abu), 
              alpha = 0.4)+
  geom_line(lty = 2, lwd=0.9,
            alpha = 0.4)+
  facet_wrap(~comp, scales = "free")+
  ylab("Paired Distance Index (specialization")+
  xlab("Treatment")

logit <- function(x){log(x/(1-x))}
ggplot(pdi_filtered , aes(x=trt, y = logit(vals), 
                          group = sp_gard, 
                          colour = fam))+
  geom_jitter(width = 0.05, 
              size = log(pdi_filtered$abu), 
              alpha = 0.4)+
  geom_line(lty = 2, lwd=0.9,
            alpha = 0.4)+
  facet_wrap(~fam, scales = "free")+
  ylab("Paired Distance Index (specialization")+
  xlab("Treatment")+
  theme(legend.position = "none")

# Test for individual orders paired in morpho-species

# See if change is spread randomly around zero or is there an evidence for more positive or negative shifts

# I think I need to break this dataset into families manually.
mod_dat <- pdi_filtered[pdi_filtered$comp == unique(pdi_filtered$comp)[1], ]
# summary(lmer(dupl_lr~1+(1|garden), data=mod_dat)) # Not different  from zero
summary(lmer(dupl_lr~0+fam+(1|garden), data=mod_dat)) # Not different  from zero

# Break tests into families to help understand these results
mdf <- mod_dat[mod_dat$fam == "cole", ]
summary(lmer(dupl_lr~1+(1|garden), data=mdf))

mod_dat <- pdi_filtered[pdi_filtered$comp == unique(pdi_filtered$comp)[2], ]
# summary(lmer(dupl_lr~1+(1|garden), data=mod_dat)) # Not different  from zero
summary(lmer(dupl_lr~ 0 + fam + (1|garden), data=mod_dat)) # Lepidoptera and orthoptera marginally

mod_dat <- pdi_filtered[pdi_filtered$comp == unique(pdi_filtered$comp)[3], ]
mod_dat <- mod_dat[is.finite(mod_dat$dupl_lr), ]
# summary(lmer(dupl_lr~1+(1|garden), data=mod_dat)) # Not different  from zero
summary(lmer(dupl_lr~0+fam+(1|garden), data=mod_dat)) # Not different

mod_dat <- pdi_filtered[pdi_filtered$comp == unique(pdi_filtered$comp)[4], ]
summary(lmer(dupl_lr~1+(1|garden), data=mod_dat)) # Marginally significant...
summary(lmer(dupl_lr~0+fam+(1|garden), data=mod_dat)) # orthoptera and marginally: cole, hemi, homo.

hist(mod_dat$dupl_lr)

# Diet change of selected species ----

# remove cole and ip
insects_no_cole <- insects[insects$morphotype != "cole001", ]
insdat <- insects_no_cole[-grep("aran|mant", insects_no_cole$morphotype), ]

# most prevalent species
which(table(insdat$morphotype) == max(table(insdat$morphotype)))
"cole002"

getBioFromGarden <- function(species,plot){
  sp <- main_biomass$CODE == plot
  plt <- main_biomass$SP_CODE == toupper(species)
  sb <- sum(main_biomass[sp & plt, ]$WEIGHT)
  return(sb)
}

treatSpecDietVec <- function(species, plot){
  mps <- insects[insects$morphotype %in% c(species),]
  rownames(treats) <- treats$codes
  mps$treatment <- treats[mps$plot, "treat"]
  df_w25 <- mps[mps$plot == plot, ]
  df_w25$tree <- as.character(df_w25$tree)
  net_w25 <- tapply(df_w25$amount, df_w25$tree, sum, na.rm = T)
  retdf <- data.frame()
  for(nm in names(net_w25)){
    print(nm)
    row
  }
}

plot <- "w1g3p1"
mbt <- main_biomass[main_biomass$LIFE.FORM %in% c("tree","shrub"),]
mbt$SP_CODE <- tolower(mbt$SP_CODE)

plotPlantComposition <- function(plot){
  return(mbt[mbt$CODE == plot, c("SP_CODE","WEIGHT")])
}

# 1. Get plots for a comparison
comparison <- c('control', 'weevil25')
cplots <- treats[treats$treat %in% toupper(comparison[1]), ]$codes
tplots <- treats[treats$treat %in% toupper(comparison[2]), ]$codes

# 2. PLots within a block 
bl = "g1"
plotsFromaBlock <- c(as.character(cplots[grep(bl, cplots)]),
                     as.character(tplots[grep(bl, tplots)]))

# 3. Get species names

treshold <- 10

cnet <- abugardnets[[plotsFromaBlock[1]]]
tnet <- abugardnets[[plotsFromaBlock[2]]]
comp_sp <- colnames(cnet)[colnames(cnet) %in% colnames(tnet)]
comp_sp_noip_nocole <- comp_sp[-grep("aran|mant|cole001", comp_sp)]
ins_dat <- insects[((insects$plot %in% c(plotsFromaBlock[1],
                              plotsFromaBlock[2])) & (insects$morphotype %in% comp_sp_noip_nocole)), ]
ins_dat$morphotype <- as.character(ins_dat$morphotype)
ins_abu_vec <- tapply(ins_dat$amount, ins_dat$morphotype, sum)
ins_abu_vec_tr <- ins_abu_vec[ins_abu_vec >= treshold]
selected.species <- names(ins_abu_vec_tr)

# Get background plant composition
cplants <- plotPlantComposition(plotsFromaBlock[1])
tplants <- plotPlantComposition(plotsFromaBlock[2])

# For a given species get diet
hsp <- selected.species[5]
hsp
sp_cont <- treatSpecDietVec(hsp, plotsFromaBlock[1])
sp_tret <- treatSpecDietVec(hsp, plotsFromaBlock[2])

# Calculate diet dissimilarity Sorensen Index
cnet[,hsp]
tnet[,hsp]
cplants
tplants

networklevel(cnet, index = "weighted connectance")
networklevel(cnet, index = "connectance")
