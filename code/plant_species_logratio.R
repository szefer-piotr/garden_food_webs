rm(list=ls())
source("code/data_processing_code.R")
source("code/pdi.R")
source("code/diet_shift.R")
# source("code/bio_log_ratio.R")
### !!!! diet shift should be a function for any two treatments

library(dplyr)
library(ggplot2)
# Select plant species and treatment

# plant species 
species <- c("piptar") # if there is no species selected use all of them
species <- c()
# logratio = 1:(treatment) / 2(control) [addition]
# logratio = 1:(control) / 2(treatment) [removal, exclosure]

treatments <- c("weevil125", "predator")
treatments <- c("weevil25", "predator")
treatments <- c("weevil125", "control")
treatments <- c("weevil25", "control")
treatments <- c("control", "predator")

# Add block to treats
treats$block <- substr(treats$codes, 3,4)

# Species selection
main_tree <- main_biomass[main_biomass$LIFE.FORM %in% c("tree","shrub"), ]
main_tree$SP_CODE <- as.character(main_tree$SP_CODE)
main_tree$CODE <- as.character(main_tree$CODE)
pltab <- table(main_tree$SP_CODE, main_tree$CODE)
species <- rownames(pltab[rowSums(pltab)> 10, ])
# species <- c("piptar", "melamu")
# bl = "g5"
# bl = "g2"

unique(ins_bioOrig$family)


makeLogratioData <- function(species, treatments, abundance = TRUE){
  
  treatments <- toupper(treatments)
  gen_log_ratio_dat  <- data.frame()
  species_log_ratio <- data.frame()
  
  # Filter species from the dataset
  if(!is.null(species)){
    species_filtered_ds <- ins_bio[ins_bio$tree %in% species, ]
    pl_bio_dat <- main_biomass[main_biomass$SP_CODE %in% toupper(species), ]
  }else{ # empty set of species uses all of them
    species_filtered_ds <- ins_bio
    pl_bio_dat <- main_biomass
  }
  
  species_filtered_ds$block <- substr(species_filtered_ds$plot,3,4)
  
  # Numerator site codes
  num_sites <- treats[treats$treat == treatments[1],]$codes 
  num_sites <- as.character(num_sites)
  # Denominator site codes
  den_sites <- treats[treats$treat == treatments[2],]$codes 
  den_sites <- as.character(den_sites)
  
  # For each block calculate log ratio based on abundance or biomass (these are probably the same)
  for(bl in unique(species_filtered_ds$block)){
    print(bl)
    
    bl_subdat <- species_filtered_ds[species_filtered_ds$block == bl, ]
    
    num_site <- num_sites[grep(bl, num_sites)]
    den_site <- den_sites[grep(bl, num_sites)]
    
    # Numerator and denominator sites
    num_subdat <- bl_subdat[bl_subdat$plot == num_site,]
    den_subdat <- bl_subdat[bl_subdat$plot == den_site,]
    
    # General log ratio, sum all species froma a given plot
    if(abundance){
      gen_num <- sum(num_subdat$amount)
      gen_den <- sum(den_subdat$amount)
    }else{ # there migth be more smaller species so there will be a difference between abundance based and biomass based log ratio
      gen_num <- sum(num_subdat$totbio)
      gen_den <- sum(den_subdat$totbio)
    }
    
    # Add how plant biomass changes
    num_pl_bio <- sum(pl_bio_dat[pl_bio_dat$CODE == num_site, ]$WEIGHT)
    den_pl_bio <- sum(pl_bio_dat[pl_bio_dat$CODE == den_site, ]$WEIGHT)
    
    glrrow <- data.frame(block = bl, 
                         lratio=log(gen_num/gen_den),
                         treat = paste(treatments[1],"/",
                                       treatments[2]),
                         spec = paste(species, collapse = ","),
                         plant_abu_logratio = log(num_pl_bio/den_pl_bio))
    gen_log_ratio_dat <- rbind(gen_log_ratio_dat, glrrow)
    
    # Individual herbivore species lratio
    for(spec in species){
      # print(spec)
      
      num_sub_spec <- num_subdat[num_subdat$tree == spec, ]
      den_sub_spec <- den_subdat[den_subdat$tree == spec, ]
      
      if(dim(num_sub_spec)[1] == 0 | dim(den_sub_spec)[1] == 0){
        next
      }
      
      numchars <- as.character(num_sub_spec$morphotype)
      denchars <- as.character(den_sub_spec$morphotype)
      boolherb <- numchars[numchars %in% denchars]
      # herbnms <- den_sub_spec[boolstr, ]$morphotype
      herbnms <- as.character(boolherb)
      herb_num_subdat <- num_sub_spec[num_sub_spec$morphotype %in% herbnms, ]
      herb_den_subdat <- den_sub_spec[den_sub_spec$morphotype %in% herbnms, ]
      # there should be no difference here
      if(abundance){
        gen_num <- herb_num_subdat[, c("morphotype","amount")]
        gen_den <- herb_den_subdat[, c("morphotype","amount")]
      }else{
        gen_num <- herb_num_subdat[, c("morphotype","totbio")]
        gen_den <- herb_den_subdat[, c("morphotype","totbio")]
      }
      
      # DEBUG
      print(gen_num)
      
      mergeddf <- merge(gen_num, gen_den, by="morphotype")
      
      #DEBUG
      print(mergeddf)
      
      ind_spec_lratio <- log(mergeddf[,2]/mergeddf[,3])
      
      # DEBUG
      print(ind_spec_lratio)
      
      names(ind_spec_lratio) <- as.character(mergeddf[,1])
      
      # Plant biomass change
      num_bool<-(pl_bio_dat$CODE==num_site & pl_bio_dat$SP_CODE==toupper(spec))
      den_bool<-(pl_bio_dat$CODE == den_site & pl_bio_dat$SP_CODE==toupper(spec))
      num_pl_bio <- sum(pl_bio_dat[num_bool, ]$WEIGHT)
      den_pl_bio <- sum(pl_bio_dat[den_bool, ]$WEIGHT)
      
      # Put the row together
      slr_row <- data.frame(species = names(ind_spec_lratio),
                            lratio = ind_spec_lratio,
                            tree = spec, 
                            block = bl,
                            plant_lr = log(num_pl_bio/den_pl_bio),
                            treat = paste(treatments[1],"/",
                                          treatments[2]))
      species_log_ratio <- rbind(species_log_ratio, slr_row)
    }
  }
  return(list(species_log_ratio, gen_log_ratio_dat))
}

# No difference
cplratio <- makeLogratioData(tolower(species), treatments, abundance = T)

# Biomass based
cplratio <- makeLogratioData(tolower(species), treatments, abundance = FALSE)



# cplratio <- makeLogratioData(c("piptar", "melamu"), treatments, abundance = T)
cplratio_spec <- cplratio[[1]]
cplratio_gen <- cplratio[[2]]

cplratio_spec$pdi <- diet_breadth_ab[as.character(cplratio_spec$species)]
cplratio_spec$weight <- pdiss[as.character(cplratio_spec$species)]

cplratio_spec$order <- substr(cplratio_spec$species,1,4)


# ggplot(cplratio_spec, aes(y=lratio, x=pdi, color = tree))+
#   geom_jitter(size = log(cplratio_spec$weight))
# 
# ggplot(cplratio_spec, aes(y=lratio, x=pdi, color = order))+
#   geom_jitter(size = log(cplratio_spec$weight))

cplratio_spec$pshape <- as.numeric(cplratio_spec$tree)+20

cbf_2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
           "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

fulltree <- c("Breynia cernua", 
              "Melanolepis multiglandulosa",
              "Pipturus argenteus", 
              "Trichospermum pleiostigma",
             "Trema orientalis")
names(fulltree) <- c("breyce", "melamu","piptar", "tricpl",
                    "tremor")

p1 <- ggplot(cplratio_spec, aes(y=lratio, x=pdi))+
  geom_point(aes(fill= tree, shape = tree))+
  # scale_fill_manual(values=cbf_2[4:8])+
  # scale_color_manual(values=cbf_2[4:8])+
  scale_shape_manual(values=c(21,22,23,24,25))+
  facet_grid(vars(tree),
             labeller = labeller(tree = fulltree,label_wrap_gen()))
  # stat_smooth(mapping=aes(weight=log(weight)),
  #             method = "lm", formula = y~x,
  #             lty = 2, se=F, col="grey40", lwd=0.5)

p1
fullord <- c("Coleoptera", "Hemiptera","Homoptera", "Lepidoptera",
             "Orthoptera")
names(fullord) <- c("cole", "hemi","homo", "lepi",
             "orth")
  
p2 <- ggplot(cplratio_spec, aes(y=lratio, x=pdi))+
  geom_point(aes(fill= order, shape = order))+
  # scale_fill_manual(values=cbf_2[4:8])+
  # scale_color_manual(values=cbf_2[4:8])+
  scale_shape_manual(values=c(21,22,23,24,25))+
  facet_grid(vars(order),
             labeller = labeller(order = fullord))+
  stat_smooth(data = subset(cplratio_spec, order == "homo"),
              mapping=aes(weight=log(weight)),
              method = "lm", formula = y~x,
              lty = 1, se=T, col="grey40", lwd=0.5)
  # scale_linetype_manual(values=c(cole=2,
  #                                hemi=2,
  #                                homo=1,
  #                                lepi=2,
  #                                orth=2))

library(ggpubr)

ggarrange(p1+theme(legend.position = "none")+
            xlab("Specialization (PDI)")+
            ylab("Predator indirect effect on herbivore abundance"),
          p2+theme(legend.position = "none")+
            xlab("Specialization (PDI)")+
            ylab("Predator direct effect on herbivore abundance"))
# I need a genearl pattern and is there an interaction in different orders and plant species

# ggplot(cplratio_spec, aes(y=lratio, x=pdi, color = order))+
#   geom_point(aes(size = log(weight), 
#                  shape = order, fill=order))+
#   # scale_fill_manual(values=cbf_2[4:8])+
#   # scale_color_manual(values=cbf_2[4:8])+
#   scale_shape_manual(values=c(21,22,23,24,25))+
#   stat_smooth(mapping=aes(weight=log(weight)), 
#               method = "lm", formula = y~x+I(x^2),
#               lty = 1, se=F, col="grey80")
# 
# ggplot(cplratio_spec, aes(y=lratio, x=pdi, color = order))+
#   geom_jitter(size = log(cplratio_spec$weight))+
#   stat_smooth(method = "lm", lty = 1, se=T)

# General test

modpip1 <- nlme::lme(lratio~pdi*order, random = ~1|block, data=cplratio_spec,
                     weights = ~ as.vector(weight))

modpip1 <- nlme::lme(lratio~pdi, random = ~0+pdi|order, data=cplratio_spec,
                     weights = ~ as.vector(weight))

summary(modpip1)
modpip1 <- nlme::lme(lratio~pdi*tree, random = ~1|block, data=cplratio_spec,
                  weights = ~ as.vector(weight))

modpip1 <- nlme::lme(lratio~pdi, random = ~1+pdi|tree, data=cplratio_spec,
                     weights = ~ as.vector(weight))

modpip1 <- nlme::lme(lratio~pdi, random = ~1+pdi|tree, 
                     data=cplratio_spec[cplratio_spec$tree == "tremor",],
                     weights = ~ as.vector(weight))

modpip1 <- nlme::lme(lratio~pdi, random = ~1|block, 
                     data=cplratio_spec[cplratio_spec$tree == "tremor",])
modpip1 <- nlme::lme(lratio~pdi, random = ~1|block, data=cplratio_spec)

# modpip2 <- lm(lratio~pdi,data=cplratio_spec)
summary(modpip1)
summary(modpip2)


mod1 <- nlme::lme(lratio~pdi*tree, random = ~1|block, data=cplratio_spec,
                  weights = ~ as.vector(weight))
mod2 <- lm(lratio~pdi*tree,data=cplratio_spec)
summary(mod1)
summary(mod2)

# This just looks the same... positive response from plants in correlated with popsitive response from herbivores
plot(lratio~plant_abu_logratio, data=cplratio_gen)

# Each plant each o0rder: nothing significant
for(plt in unique(cplratio_spec$tree)[-1]){
  print(paste(plt,"_________________"))
  subdat <- cplratio_spec[cplratio_spec$tree == plt, ]
  mod <- nlme::lme(lratio ~ pdi, random = ~1|block, 
                   data = subdat, weights = ~ as.vector(weight))
  print(summary(mod))
}

for(plt in unique(cplratio_spec$order)){
  print(paste(plt,"_________________"))
  subdat <- cplratio_spec[cplratio_spec$order == plt, ]
  mod <- nlme::lme(lratio ~ pdi, random = ~1|block, 
                   data = subdat, weights = ~ as.vector(weight))
  print(summary(mod))
}

# Here homoptera is significant


# Each plant each order plant interacton nothing significant

for(plt in unique(cplratio_spec$tree)){
  for(ord in unique(cplratio_spec$order)){
    
    print(paste(plt, ord))
    
    c1 <- cplratio_spec$tree == plt
    c2 <- cplratio_spec$order == ord
    subdat <- cplratio_spec[c1&c2, ]
    
    isenough <- TRUE

    huh <- tryCatch({
      mod <- nlme::lme(lratio ~ pdi,
                       random = ~1|block,
                       data = subdat,
                       weights = ~ as.vector(weight))
      }, error = function(e) {
        print(paste(title, "not enough obs"))
        isenough <- FALSE
        return(isenough)
      })
    
    print(summary(huh))
    
    # if(huh){
    #   mod <- nlme::lme(lratio ~ pdi, random = ~1|block,
    #                    data = subdat, weights = ~ as.vector(weight))
    #   print(summary(mod))
    # }
    
  }
  
}

ggplot(cplratio_spec, aes(y = lratio,x = pdi))+
  geom_point()+
  geom_smooth(type = "lm")


inter <- ggplot(cplratio_spec, aes(y=lratio, x=pdi))+
  geom_point(aes(fill= order, shape = order))+
  scale_shape_manual(values=c(21,22,23,24,25))+
  facet_grid(vars(order), vars(tree))+
  stat_smooth(mapping=aes(weight=log(weight)),
              method = "lm", formula = y~x,
              lty = 1, se=T, col="grey40", lwd=0.5)
# scale_linetype_manual(values=c(cole=2,

# for(plt in unique(cplratio_spec$order)){
#   print(paste(plt,"_________________"))
#   subdat <- cplratio_spec[cplratio_spec$order == plt, ]
#   mod <- nlme::lme(lratio ~ I(pdi^2), random = ~1|block, 
#                    data = subdat, weights = ~ as.vector(weight))
#   print(summary(mod))
# }
