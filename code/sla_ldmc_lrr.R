# SLA and LDMC relationships with LRR

rm(list=ls())
source("code/data_processing_code.R")
# source("code/pdi.R")
# source("code/diet_shift.R")

# I need evaluation of LRR for each plant species
# biollcp

library(dplyr)
library(ggplot2)
library(nlme)
library(MASS)

# Within each garden check which plant species are comparable
biollcp <- biollcp %>%
  mutate(plnat_gard = paste(plnm, gard,trt, sep="_"))

bdat <- biollcp
mb <- main

traitdf <- data.frame()

for(ga in unique(bdat$gard)){
  print(ga)
  cd <- bdat[bdat$gard == ga & bdat$trt == "CONTROL", ]
  pd <- bdat[bdat$gard == ga & bdat$trt == "PREDATOR", ]
  comppl <- unique(cd$plnm)[unique(cd$plnm) %in% unique(pd$plnm)]
  
  print(comppl)
  
  # Conditions
  gcond <- grepl(ga,mb$BLOCK, ignore.case = T)
  trtcond <- mb$TREAT %in% "PREDATOR" # should this be the base for our analysis
  pltcond <- mb$SP_CODE %in% toupper(comppl)
  
  slaval <- mb[gcond & trtcond & pltcond, ]
  
  traitdf <- rbind(traitdf, data.frame(
    garden = ga,
    plant = tolower(slaval$SP_CODE),
    sla = slaval$SLA,
    ldmc = slaval$LDMC,
    water = slaval$WATER/slaval$WET..g., #This is perfectlycorrelated with LDMC
    dam = slaval$HERB
  ))
}

traitdf$pllrr <- NA
arthlrr <- data.frame()

for (row in 1:dim(traitdf)[1]){
  gd <- traitdf[row, ]$garden
  plant <- traitdf[row, ]$plant
  subdat <- bdat %>%
    filter(gard == gd, plnm == plant)
  
  # print(gd)
  # print(plant)

  trtindc <- subdat$trt == "CONTROL"
  trtindp <- subdat$trt == "PREDATOR"
  cpl <- unique(subdat[trtindc, ]$plbio)
  ppl <- unique(subdat[trtindp, ]$plbio)

  plantlrr <- log(cpl/ppl)
  
  traitdf[row, ]$pllrr <- plantlrr
  
  # subdatAP <- bdat %>%
  #   filter(gard == gd, plnm == plant,
  #          nms %in% c("aran","mant"))
  # subdatH <- bdat %>%
  #   filter(gard == gd, plnm == plant,
  #          !(nms %in% c("aran","mant")))
  
  # print(subdat)
  
  for (order in unique(subdat$nms)){
    
    ordind <- subdat$nms %in% order
    trtindc <- subdat$trt == "CONTROL"
    trtindp <- subdat$trt == "PREDATOR"
      
    cbio <- subdat[ordind & trtindc, ]$bio
    pbio <- subdat[ordind & trtindp, ]$bio
    
    if (length(cbio) == 0){cbio <- NA}
    if (length(pbio) == 0){pbio <- NA}
    
    # print(subdat[subdat$trt == "CONTROL", ])
    # print(subdat[subdat$trt == "PREDATOR", ])
    # 
    arthlrr <- rbind(arthlrr, data.frame(
      gard = gd, plnm = plant, ord = order,
      lrr = log(cbio/pbio)
    ))
    
  }
}

# Traits vs LRR of plants

# 1. General model for traits and plant LRR ----
summary(gls(pllrr ~ sla + water, 
            data=traitdf[complete.cases(traitdf), ]))

multi_traitdf <- traitdf[rep(1:14, each = length(unique(bdat$nms))), ]

tlrrdf <- cbind(multi_traitdf, arthlrr)

# 2. Damage vs traits ----
lmdam1 <- gls(dam~sla+water, method = "ML",
              data = traitdf[complete.cases(traitdf), ]) 
summary(lmdam1) 

# Stepwise selection
step.lm <-stepAIC(lmdam1, direction = "both")
summary(step.lm) # Water content is marginally significant predictor of leaf damage.

# 3. Traits vs plant LRR ----

test_lrr <- tlrrdf %>%
  filter(lrr != Inf & lrr != -Inf & !(is.na(lrr)))

# Replace g3 melamu missing val with an average
# herb_lrr <- test_lrr[!(grepl("aran|mant", test_lrr$ord)), ]
# mv <- mean(main[main$SP_CODE == "MELAMU",]$WATER/main[main$SP_CODE == "MELAMU",]$WET..g.,na.rm=T)
# herb_lrr[is.na(herb_lrr$water), ]$water <-  mv
# 
# # Basic model
# B1 <- gls(lrr~1+sla+water+sla*ord+water*ord, 
#           method = "REML",
#           data = herb_lrr)
# 
# B2 <- lme(lrr~1+sla+water+sla*ord+water*ord,
#           random = ~1|gard,
#           method = "REML",
#           data = herb_lrr)
# 
# AIC(B1,B2) # Random effect of garden seems unnecessary
# summary(B1)
# 
# # Refit B1 with ML
# B1 <- gls(lrr~1+sla+water+sla*ord+water*ord, 
#           method = "ML",
#           data = herb_lrr)
# 
# drop1(B1)
# 
# B3 <- gls(lrr~1+sla+water, 
#           method = "ML",
#           data = herb_lrr)
# 
# drop1(B3)
# anova(B1,B3)
# 
# B4 <- gls(lrr~1+sla, 
#           method = "ML",
#           data = herb_lrr)
# 
# drop1(B4)
# anova(B3,B4)
# 
# B5 <- gls(lrr~1, 
#           method = "ML",
#           data = herb_lrr)
# 
# drop1(B5)
# anova(B4,B5)

# 4. Interaction of PDI with quality ----
lrr.pdi <- read.table("datasets/cplratio_spec.txt")
lrr.pdi.1 <- lrr.pdi[,c("species","lratio","tree","block",
                        "fam","plant_lr", "pdi", "weight")]

getQuality <- function (plant, 
                        garden, 
                        treatment, 
                        trait){
  print(paste(plant, garden, treatment, trait))
  subdat <- main %>% 
    filter(SP_CODE == toupper(plant),
           grepl(garden, BLOCK, ignore.case = T),
           TREAT == toupper(treatment))
  if(trait == "water"){
    return(subdat$WATER/subdat$WET..g.)
  }
  if(trait == "ldmc"){
    return(subdat$LDMC)
  }
  if(trait == "sla"){
    return(subdat$SLA)
  }
}

# getQuality(plant = "melamu",
#            garden = "g3",
#            treatment = "predator",
#            trait = "ldmc")

lrr.pdi.1$water <- NA
lrr.pdi.1$sla <- NA
lrr.pdi.1$ldmc <- NA
for (row in 1:dim(lrr.pdi.1)[1]){
  lrr.pdi.1[row, ]$water <- getQuality(lrr.pdi.1[row,]$tree,
             lrr.pdi.1[row,]$block,
             "predator",
             "water")
  lrr.pdi.1[row, ]$sla <- getQuality(lrr.pdi.1[row,]$tree,
                                       lrr.pdi.1[row,]$block,
                                       "predator",
                                       "sla")
  lrr.pdi.1[row, ]$ldmc <- getQuality(lrr.pdi.1[row,]$tree,
                                       lrr.pdi.1[row,]$block,
                                       "predator",
                                       "ldmc")
}

lrr.pdi.2 <- lrr.pdi.1 %>%
  filter(!grepl("aran|mant",fam))

lrr.pdi.2[is.na(lrr.pdi.2$water), ]$water <- mv

# ggplot(lrr.pdi.2, aes(x = pdi, y = lratio))+
#   geom_point()+
#   stat_smooth(method = "lm")+
#   facet_wrap(~water)
# library(boot)
# ggplot(lrr.pdi.2, aes(x = pdi, y = lratio))+
#   geom_point(size = log(lrr.pdi.2$weight))+
#   stat_smooth(method = "lm")+
#   facet_wrap(~inv.logit(ldmc))

# Model
lrr.pdi.test <- lrr.pdi.2 %>%
  filter(!is.na(lratio)) %>%
  filter(lratio != 0) %>%
  filter(!is.na(water))# delete this in case

L1 <- gls(lratio~pdi+sla+water+sla*pdi+water*pdi, 
          data = lrr.pdi.test, method = "ML")
L2 <- lme(lratio~pdi+sla+water+sla*pdi+water*pdi,
          random = ~1|block,
          data = lrr.pdi.test, method = "ML")

AIC(L1,L2) # Block random effect improves model
drop1(L2)
summary(L2)

L3 <- lme(lratio~pdi+sla+water+sla*pdi,
          random = ~1|block,
          data = lrr.pdi.test, method = "ML")

AIC(L2,L3)
drop1(L3)

L4 <- lme(lratio~pdi+sla+water,
          random = ~1|block,
          data = lrr.pdi.test, method = "ML")

AIC(L3,L4) # Mode with sla seems better
drop1(L4) 

L5 <- lme(lratio~sla+water,
          random = ~1|block,
          data = lrr.pdi.test, method = "ML")
AIC(L4,L5) # Mode with sla seems better
drop1(L5) 

L6 <- lme(lratio~water,
          random = ~1|block,
          data = lrr.pdi.test, method = "ML")

AIC(L5,L6) # Mode with sla seems better

anova(L5,L6) 

L7 <- lme(lratio~1,
          random = ~1|block,
          data = lrr.pdi.test, method = "ML")

AIC(L6,L7)
anova(L6,L7) 
summary(L6) # water content predicts well responses of herbivores to predators

# Plot of water content
lrr.pdi.test$fam <- recode(lrr.pdi.test$fam, cole = "Coleoptera",
         hemi = "Heteroptera",
         homo = "Homoptera",
         lepi = "Lepidoptera",
         orth = "Orthoptera")
lrr.pdi.test.upd <- rename(lrr.pdi.test, "Order" = fam)

colorBlindGrey8   <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
scales::show_col(colorBlindGrey8)

ggplot(lrr.pdi.test.upd,
       aes(x = water, y = lratio))+
  geom_point(aes(color = Order),size=5, alpha = 0.8)+
  geom_smooth(method = "lm", col = "grey40",
              lwd = 1, lty = 1, se=T)+
  ylab("Log Response Ratio")+
  xlab("Water content [%]")+
  scale_colour_manual(values = colorBlindGrey8[2:6])+
  theme(legend.position = "bottom")







# Final model would be th L3 model
# summary(L3) # NEvertheless no significant effects can be found
# plot(L3)
# lrr.pdi.test$preds <- predict(L3)

# ggplot(lrr.pdi.test, aes(x = pdi, y = lratio))+
#   geom_point()+
#   stat_smooth(method = "lm")+
#   facet_wrap(~ldmc)
# 
# ggplot(lrr.pdi.test, aes(x = pdi, y = preds))+
#   geom_point()+
#   stat_smooth(method = "lm")+
#   facet_wrap(~sla)
# 
# ggplot(lrr.pdi.test, aes(x = pdi, y = preds))+
#   geom_point()+
#   stat_smooth(method = "lm")+
#   facet_wrap(~water)


# L3.w <- lme(lratio~pdi+sla+water+ldmc*pdi+water*pdi,
#             random = ~1|block, weights = ~ I(1/weight),
#             data = lrr.pdi.test, method = "ML")
# summary(L3.w) # weighted regresssion supports lack of effects

# Simplification of the model based on assumptions
# Interaction is not significant
# L5 <- lme(lratio~water*pdi,
#           random = ~1|block, weights = ~ I(1/log(weight)),
#           data = lrr.pdi.test, method = "ML")
# L5.1 <- lme(lratio~pdi*ldmc,
#           random = ~1|block, weights = ~ I(1/log(weight)),
#           data = lrr.pdi.test, method = "ML")
# L5.2 <- lme(lratio~sla*pdi,
#           random = ~1|block, weights = ~ I(1/log(weight)),
#           data = lrr.pdi.test, method = "ML")
# summary(L5)
# summary(L5.1)
# summary(L5.2)

# Repeated selection with weights
L1.w <- gls(lratio~pdi+sla+water+sla*pdi+water*pdi, 
          weights = ~ I(1/log(weight)),
          data = lrr.pdi.test, method = "ML")
L2.w <- lme(lratio~pdi+sla+water+sla*pdi+water*pdi,
          random = ~1|block,
          weights = ~ I(1/log(weight)),
          data = lrr.pdi.test, method = "ML")
L2.w.reml <- lme(lratio~pdi+sla+water+sla*pdi+ldmc*pdi+water*pdi,
          random = ~1|block,
          weights = ~ I(1/log(weight)),
          data = lrr.pdi.test, method = "REML")

AIC(L1.w,L2.w) # Block random effect improves the model
drop1(L2.w) # No improvemet - keep the full model

L3.w <- lme(lratio~pdi+sla+water+sla*pdi,
            random = ~1|block,
            weights = ~ I(1/log(weight)),
            data = lrr.pdi.test, method = "ML")

AIC(L2.w,L3.w) # Block random effect improves the model
drop1(L3.w)

L4.w <- lme(lratio~pdi+sla+water,
            random = ~1|block,
            weights = ~ I(1/log(weight)),
            data = lrr.pdi.test, method = "ML")

AIC(L3.w,L4.w) # Block random effect improves the model
drop1(L4.w)

L5.w <- lme(lratio~sla+water,
            random = ~1|block,
            weights = ~ I(1/log(weight)),
            data = lrr.pdi.test, method = "ML")

AIC(L4.w,L5.w) # Block random effect improves the model
drop1(L5.w)

L6.w <- lme(lratio~water,
            random = ~1|block,
            weights = ~ I(1/log(weight)),
            data = lrr.pdi.test, method = "ML")

AIC(L5.w, L6.w)
drop1(L6.w) # no more improvement

summary(L6.w)
R2(L6.w)
# Prediction
# lrr.pdi.test$preds <- predict(L2.w)
# ggplot(lrr.pdi.test, aes(x = pdi, y = preds))+
#   geom_point()+
#   stat_smooth(method = "lm")+
#   facet_wrap(~ldmc)

# After removing water content we loose all significance
# L3.1 <- lme(lratio~pdi+sla+sla*pdi+ldmc*pdi,
#             random = ~1|block,
#             weights = ~ I(1/log(weight)),
#             data = lrr.pdi.test, method = "ML")

# summary(L3.1) 


# library(piecewiseSEM)
# rsquared(L2.w)
# rsquared(L2)
# rsquared(L3)
# however it is difficult to interpret effects of weighted corelation in effects which not include relation of lrr with pdi.

# Higher level interactions are not significna
# L2.max <- lme(lratio~sla*water*ldmc*pdi+sla*water*pdi+sla*ldmc*pdi+water*ldmc*pdi,
#             random = ~1|block,
#             weights = ~ I(1/log(weight)),
#             data = lrr.pdi.test, method = "ML")
# summary(L2.max)
# drop1(L2.max)
# 
# L3.max <- lme(lratio~sla*water*pdi+sla*ldmc*pdi+water*ldmc*pdi,
#               random = ~1|block,
#               weights = ~ I(1/log(weight)),
#               data = lrr.pdi.test, method = "ML")
# drop1(L3.max)
# 
# pairs(lrr.pdi.test[,c(9,10,11)])


# Discrete categories and test for interactions
library(vegan)
quality <- rda(lrr.pdi.test[,c(9,10)])
biplot(quality, type = "text")
text(quality, display = "bp")
qty <- summary(quality)
lrr.pdi.test$qty <- qty$sites[,1]

lrr.pdi.test.int <-lrr.pdi.test %>%
  mutate(disqual = ifelse(qty > 0, "Quality -", "Quality +"),
         dispdi = ifelse(pdi < median(pdi), "Generalists", "Specialists"))

ggplot(lrr.pdi.test.int, aes(x = dispdi, y = lratio))+
  geom_jitter(width = 0.05)+
  stat_summary(method = "pointrange", 
               col = "grey20",
               lwd = 1)+
  facet_wrap(~disqual)

ggplot(lrr.pdi.test.int, aes(x = dispdi, y = lratio))+
  geom_jitter(width = 0.1, alpha=0.3)+
  stat_summary(fun.data=mean_cl_boot, 
               geom="pointrange",
               col = "grey30",
               lwd = 1)+
  facet_wrap(~disqual) + xlab("")+ylab("")

# Test
QL1 <- gls(lratio ~ dispdi*disqual, 
           data = lrr.pdi.test.int,
           method = "ML")
QL1.rand <- lme(lratio ~ dispdi*disqual, 
                random = ~1|block, 
                data = lrr.pdi.test.int,
                method = "ML")
AIC(QL1, QL1.rand) # Random effect improves the model

# QUality influenced log ration. Higher quality was related to average positive log response ratio. However interaction with specialization was not significant
summary(QL1.rand)
drop1(QL1.rand)
