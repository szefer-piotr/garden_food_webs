# Individual species PDI in response to treatments

rm(list=ls())
source("code/data_processing_code.R")
source("code/contingencyTable.R")
library(ggplot2)
library(emmeans)
library(multcomp)
library(vegan)
library(bipartite)


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

ggplot(pdi_filtered , aes(x=trt, y = vals, group = species))+
  geom_jitter(width = 0.1, size = sqrt(abu))+
  geom_line(lty = 2, lwd=0.1)+
  facet_wrap(~comp, scales = "free")+
  ylab("Paired Distance Index (specialization")+
  xlab("Treatment")

mod_dat <- pdi_filtered[pdi_filtered$comp == unique(pdi_filtered$comp)[1], ]
summary(lmer(dupl_lr~1+(1|garden), data=mod_dat)) # Not different  from zero

mod_dat <- pdi_filtered[pdi_filtered$comp == unique(pdi_filtered$comp)[2], ]
summary(lmer(dupl_lr~1+(1|garden), data=mod_dat)) # Not different  from zero

mod_dat <- pdi_filtered[pdi_filtered$comp == unique(pdi_filtered$comp)[3], ]
mod_dat <- mod_dat[is.finite(mod_dat$dupl_lr), ]
summary(lmer(dupl_lr~1+(1|garden), data=mod_dat)) # Not different  from zero

mod_dat <- pdi_filtered[pdi_filtered$comp == unique(pdi_filtered$comp)[4], ]
summary(lmer(dupl_lr~1+(1|garden), data=mod_dat)) # Not different  from zero
