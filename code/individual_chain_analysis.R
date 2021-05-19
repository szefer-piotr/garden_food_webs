# Each possible chain.

rm(list = ls())

source("code/data_processing_code.R")
source("code/pdi.R")
source("code/diet_shift.R")

library(dplyr)
library(ggplot2)

csites <- treats[treats$treat %in% c("CONTROL"), ]$codes
psites <- treats[treats$treat %in% c("PREDATOR"), ]$codes
csites <- as.character(csites)
psites <- as.character(psites)

# Identify possible chains

main <- main[main$LIFE.FORM %in% c("shrub", "tree"), ]

chains <- data.frame()
for(garden in unique(main$BLOCK)){
  print(garden)
  subdat <- main[main$BLOCK == garden, ]
  ptrees <- subdat[subdat$CODE %in% psites, ]$SP_CODE
  ctrees <- subdat[subdat$CODE %in% csites, ]$SP_CODE
  chains <- rbind(chains,
                  data.frame(tree = as.character(ctrees[ctrees %in% ptrees]),
                             garden = garden,
                             psite = unique(subdat[subdat$CODE %in% psites,]$CODE),
                             csite = unique(subdat[subdat$CODE %in% csites,]$CODE)))
}

# Get plant biomass for each 
# For each plot and plant get IAP and herbivore community (abu and bio) + 
# plant species biomas
chains$garden <- substr(chains$psite,3,4)

chain_vals <- data.frame()
for(row in 1:dim(chains)[1]){
  tsp <- as.character(chains[row, ]$tree)
  gard <- chains[row, ]$garden
  psite <- as.character(chains[row,]$psite)
  csite <- as.character(chains[row,]$csite)
  
  c1 <- grepl(gard, main$BLOCK, ignore.case = T)
  c2 <- main$SP_CODE == tsp
  c3 <- main$CODE %in% c(psite, csite)
  
  vars <- c("CODE","SP_CODE", "TREAT","WEIGHT", "LEAVES","SLA", "LDMC","WATER")
  plsubset <- main[c1 & c2 & c3, vars]
  
  cplbio <- plsubset[plsubset$CODE %in% csites, 4:8]
  pplbio <- plsubset[plsubset$CODE %in% psites, 4:8]
  colnames(cplbio) <- paste(tolower(colnames(cplbio)), "c", sep = "_")
  colnames(pplbio) <- paste(tolower(colnames(pplbio)), "p", sep = "_")
  
  chain_vals <- rbind(chain_vals,data.frame(cbind(cplbio, pplbio)))
}

chainsdf <- cbind(chains, chain_vals)

#
abu_iap_list <- list()
abu_herb_list <- list()

abu_list <- list()

# Abundance based IAP/H ratio

chainsdf$habu_p <- NA
chainsdf$iapabu_p <- NA
chainsdf$habu_c <- NA
chainsdf$iapabu_c <- NA
chainsdf$pdicwm <- NA

for(row in 1:dim(chains)[1]){
  tsp <- as.character(chains[row, ]$tree)
  gard <- chains[row, ]$garden
  psite <- as.character(chains[row,]$psite)
  csite <- as.character(chains[row,]$csite)
  
  insdat <- ins_bioOrig %>%
    filter(plot == psite & tree == tolower(tsp)) %>%
    mutate(guild = ifelse(grepl("aran|mant", family),
                          "predator",
                          "herbivore")) %>%
    group_by(guild) %>%
    summarize(sum = sum(amount, na.rm=T))
  
  # abu_list[[paste(tsp_)]] <- insdat
  
  if(dim(insdat)[1] == 0){
    print("empty")
    next
  } else {
    hval <- insdat[insdat$guild == "herbivore", ]$sum
    pval <- insdat[insdat$guild == "predator", ]$sum
    
    # herbivore diet breadth
    insdatPDI <- ins_bioOrig %>%
      filter(plot == psite & tree == tolower(tsp)) %>%
      mutate(guild = ifelse(grepl("aran|mant", family),
                            "predator",
                            "herbivore"))%>%
      select(morph, amount) %>%
      mutate(pdi = diet_breadth[morph]) %>%
      summarise(cwm = weighted.mean(pdi, amount, na.rm = T))
    
    chainsdf[row, ]$pdicwm <- insdatPDI$cwm
    chainsdf[row, ]$habu_p <- hval
    chainsdf[row, ]$iapabu_p <- pval
  }
}


for(row in 1:dim(chains)[1]){
  tsp <- as.character(chains[row, ]$tree)
  gard <- chains[row, ]$garden
  psite <- as.character(chains[row,]$psite)
  csite <- as.character(chains[row,]$csite)
  
  insdat <- ins_bioOrig %>%
    filter(plot == csite & tree == tolower(tsp)) %>%
    mutate(guild = ifelse(grepl("aran|mant", family),"predator","herbivore")) %>%
    group_by(guild) %>%
    summarize(sum = sum(amount, na.rm=T))
  
  # abu_list[[paste(tsp_)]] <- insdat
  
  if(dim(insdat)[1] == 0){
    print("empty")
    next
  } else {
    
    hval <- insdat[insdat$guild == "herbivore", ]$sum
    pval <- insdat[insdat$guild == "predator", ]$sum
    
    if(length(hval)==1){
      chainsdf[row, ]$habu_c <- hval
    }
    if(length(pval)==1){
      chainsdf[row, ]$iapabu_c <- pval
    }
  }
}

chainsdf.1 <- chainsdf %>%
  mutate(plantlrr = log(weight_c/weight_p),
         iap_h_ratio = iapabu_c/habu_c,
         iap_h_ratio_p = iapabu_p/habu_p,
         iaplrr = log(iapabu_c/iapabu_p),
         hlrr = log(habu_c/habu_p))

# CWM PDI does not affect cascading effects
plot(plantlrr~pdicwm,data = chainsdf.1)
summary(lm(plantlrr~pdicwm,data = chainsdf.1))
abline(lm(plantlrr~pdicwm,data = chainsdf.1))

chdf.pl <- chainsdf %>%
  group_by(tree) %>%
  summarise(hp = sum(habu_p, na.rm = T),
            hc = sum(habu_c, na.rm = T),
            app = sum(iapabu_p, na.rm = T),
            apc = sum(iapabu_c, na.rm = T),
            pp = sum(weight_p, na.rm = T),
              pc = sum(weight_c, na.rm = T)) %>%
  mutate(aphr.exb = app/hp,
         aphr.ctb = apc/hc,
         plrr = pp/pc)

chdf.bl <- chainsdf %>%
  group_by(garden) %>%
  summarise(hp = sum(habu_p, na.rm = T),
            hc = sum(habu_c, na.rm = T),
            app = sum(iapabu_p, na.rm = T),
            apc = sum(iapabu_c, na.rm = T),
            pp = sum(weight_p, na.rm = T),
            pc = sum(weight_c, na.rm = T)) %>%
  mutate(aphr.exb = app/hp,
         aphr.ctb = apc/hc,
         plrr = pp/pc)


# ggplot(chdf.bl, aes(x = aphr.exb, y = plrr))+
#   geom_point()+
#   geom_smooth(method = "lm",col= "grey30",
#               se = F) + xlab("AP to herbivore ratio")+ylab("LRR for plants")

# ABSOLUTE VALS
# 1. Plot based on control ratio
# cr.a <- ggplot(chainsdf.1, aes(x = iap_h_ratio, y = abs(plantlrr)))+
#   geom_point()+
#   geom_smooth(method = "lm",col= "grey30",
#               se = F) + xlab("AP to herbivore ratio")+ylab("LRR for plants")
# 
# # 1.1 without one extreme ratio
# cr.nex.a <- ggplot(chainsdf.1[-which(chainsdf.1$iap_h_ratio == max(chainsdf.1$iap_h_ratio, na.rm = T)),], aes(x = iap_h_ratio, y = abs(plantlrr)))+
#   geom_point()+
#   geom_smooth(method = "lm",
#               col= "grey30",
#               se = F)+ xlab("AP to herbivore ratio")+ylab("LRR for plants")
# 
# # 2. Plot based on predator ratio
# pr.a <- ggplot(chainsdf.1, aes(x = iap_h_ratio_p, y = abs(plantlrr)))+
#   geom_point()+
#   geom_smooth(method = "lm",col= "grey30",
#               se = F)+ xlab("AP to herbivore ratio")+ylab("LRR for plants")
# 
# # 2.2 withouth extreme observation
# pr.nex.a <- ggplot(chainsdf.1[-which(chainsdf.1$iap_h_ratio_p == max(chainsdf.1$iap_h_ratio_p, na.rm = T)),], aes(x = iap_h_ratio_p, y = abs(plantlrr)))+
#   geom_point()+
#   geom_smooth(method = "lm",col= "grey30",
#               se = T)+ xlab("AP to herbivore ratio")+ylab("LRR for plants")

# RAW VALUES
# 1. Plot based on control ratio
cr <- ggplot(chainsdf.1, aes(x = log(iap_h_ratio), y = plantlrr))+
  geom_point(size = 5, alpha= 0.5)+
  # geom_smooth(method = "lm",col= "grey30",
  #             se = F) + 
  xlab("AP to herbivore ratio")+ylab("LRR for plants")

cr
# # 1.1 without one extreme ratio
# cr.nex <- ggplot(chainsdf.1[-which(chainsdf.1$iap_h_ratio == max(chainsdf.1$iap_h_ratio, na.rm = T)),], aes(x = iap_h_ratio, y = plantlrr))+
#   geom_point()+
#   geom_smooth(method = "lm",
#               col= "grey30",
#               se = F)+ xlab("AP to herbivore ratio")+ylab("LRR for plants")

# 2. Plot based on control ratio
# pr <- ggplot(chainsdf.1, aes(x = iap_h_ratio_p, y = plantlrr))+
#   geom_point()+
#   geom_smooth(method = "lm",col= "grey30",
#               lty=2, se = F)+
#   xlab("AP to herbivore ratio")+ylab("LRR for plants")

# 2.2 withouth extreme observation
# pr.nex <- ggplot(chainsdf.1[-which(chainsdf.1$iap_h_ratio_p == max(chainsdf.1$iap_h_ratio_p, na.rm = T)),], aes(x = iap_h_ratio_p, y = plantlrr))+
#   geom_point()+
#   geom_smooth(method = "lm",col= "grey30",
#               se = F)+ xlab("AP to herbivore ratio")+ylab("LRR for plants")
# 
# cr
# cr.nex
# pr
# pr.nex
# 
# # ABSOLUTE
# # Based on control ratio
lm1.cr <- lm(abs(plantlrr)~iap_h_ratio,
            data = chainsdf.1)
summary(lm1.cr)
# 
# # Withouth extreme obs
# lm1.cr.nex <- lm(abs(plantlrr)~iap_h_ratio,
#    data = chainsdf.1[-which(chainsdf.1$iap_h_ratio == max(chainsdf.1$iap_h_ratio, na.rm = T)),])
# summary(lm1.cr.nex)
# 
# # Exclosure based AP/H ratio
# lm1.pr <- lm(abs(plantlrr)~iap_h_ratio, 
#             data = chainsdf.1)
# summary(lm1.pr)
# # Withouth and extreme obs
# lm1.pr.nex <- lm(abs(plantlrr)~iap_h_ratio_p, 
#             data = chainsdf.1[-which(chainsdf.1$iap_h_ratio_p == max(chainsdf.1$iap_h_ratio_p, na.rm = T)),])
# summary(lm1.pr.nex)
# 
# # RAW
# 
# # Based on control ratio
lm1.cr <- lm(plantlrr~log(iap_h_ratio),
             data = chainsdf.1)
summary(lm1.cr)
# 
# # Withouth extreme obs
# lm1.cr.nex <- lm(plantlrr~iap_h_ratio, 
#                  data = chainsdf.1[-which(chainsdf.1$iap_h_ratio == max(chainsdf.1$iap_h_ratio, na.rm = T)),])
# summary(lm1.cr.nex)
# 
# # Exclosure based AP/H ratio
lm1.pr <- lm(plantlrr~iap_h_ratio,
             data = chainsdf.1)
summary(lm1.pr)
# # Withouth and extreme obs
# lm1.pr.nex <- lm(plantlrr~iap_h_ratio_p, 
#                  data = chainsdf.1[-which(chainsdf.1$iap_h_ratio_p == max(chainsdf.1$iap_h_ratio_p, na.rm = T)),])
# summary(lm1.pr.nex)
# 
# ggpubr::ggarrange(cr, cr.nex, pr, pr.nex, nrow = 2, ncol = 2,
#                   labels = c("A","B","C","D"))
# 
# ggpubr::ggarrange(cr.a, cr.nex.a, pr.a, pr.nex.a, nrow = 2, ncol = 2,
#                   labels = c("A","B","C","D"))
