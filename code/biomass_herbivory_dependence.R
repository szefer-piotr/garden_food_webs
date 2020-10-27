rm(list=ls())
source("code/data_processing_code.R")
source("code/contingencyTable.R")
library(ggplot2)
library(emmeans)
library(multcomp)
library(vegan)

# remove cole001
treats_to_plot <- as.character(unique(treats$treat))[c(6,3,4,5,2)]
treats_to_remove_cole <- treats_to_plot
sites_to_remove_cole <- treats[treats$treat %in% treats_to_remove_cole, ]$codes
strc <- as.character(sites_to_remove_cole)

gardnets_nocole <- gardnets
abugardnets_nocole <- abugardnets


# Removing cole001 from all networks
for (net in names(gardnets_nocole)){
  if(net %in% strc){
    subnet <- gardnets_nocole[[net]]
    gardnets_nocole[[net]] <- subnet[, colnames(subnet) != "cole001"]
    
    subnet <- abugardnets_nocole[[net]]
    abugardnets_nocole[[net]] <- subnet[, colnames(subnet) != "cole001"]
  }
}

gardnets <- gardnets_nocole
gardnets[["w1g2p4"]] <- t(gardnets[["w1g2p4"]])
gardnets <- gardnets[names(gardnets) != "w1g1p6"]
gardnets <- gardnets[names(gardnets) != "w1g2p5"]

abugardnets <- abugardnets_nocole
abugardnets[["w1g2p4"]] <- t(abugardnets[["w1g2p4"]])
abugardnets <- abugardnets[names(abugardnets) != "w1g1p6"]
abugardnets <- abugardnets[names(abugardnets) != "w1g2p5"]

# Filter either IPs or herbivores
what_to_filter_out <- "ips"

selected_gardnets <- gardnets
abuselected_gardnets <- abugardnets

# Removing cole001 from all networks
for (net in names(selected_gardnets)){
  if(net %in% strc){
    
    subnet <- selected_gardnets[[net]]
    if(what_to_filter_out == "ips"){
      selected_gardnets[[net]] <- subnet[, -grep("aran|mant", colnames(subnet))]
    }
    if(what_to_filter_out == "herbivores"){
      selected_gardnets[[net]] <- subnet[, grep("aran|mant", colnames(subnet))]
    }
    
    subnet <- abuselected_gardnets[[net]]
    if(what_to_filter_out == "ips"){
      abuselected_gardnets[[net]] <- subnet[, -grep("aran|mant", colnames(subnet))]
    }
    if(what_to_filter_out == "herbivores"){
      abuselected_gardnets[[net]] <- subnet[, grep("aran|mant", colnames(subnet))]
    }
    
  }
}

gardnets <- selected_gardnets
abugardnets <- abuselected_gardnets

# Analyses
getBioFromGarden <- function(species,plot){
  sp <- main_biomass$CODE == plot
  plt <- main_biomass$SP_CODE == toupper(species)
  sb <- sum(main_biomass[sp & plt, ]$WEIGHT)
  return(sb)
}

all_treatments <- c("weevil125", "weevil25", 
               "control", "insecticide",
               "predator")

all_treat_bio_herb_data <- data.frame()

for (treatment in all_treatments) {
  print(treatment)
  
  treatment_plots <- treats[treats$treat == toupper(treatment), ]$code
  
  treat_bio_herb_data <- data.frame()
  
  for(plt in treatment_plots){
    for(spc in rownames(gardnets[[plt]])){
      
      print(paste(plt,spc))
      print(getBioFromGarden(spc, plt))
      

      # Community descriptors
      biomass <- sum(gardnets[[plt]][spc, ])
      abundance <- sum(abugardnets[[plt]][spc, ])
      sp_rich <- sum(gardnets[[plt]][spc, ]>0)
      div_b <- vegan::diversity(gardnets[[plt]][spc, ], 
                         index = "shannon")
      div_a <- vegan::diversity(abugardnets[[plt]][spc, ], 
                         index = "shannon")

      
      # Add summaries for herbivorous communities
      treat_bio_herb_data <- rbind(treat_bio_herb_data, 
                             data.frame(plot = plt,
                                        species = spc,
                                        biomass = getBioFromGarden(spc, plt),
                                        bio = biomass,
                                        abu = abundance,
                                        rich = sp_rich,
                                        div_a = div_a, # abund based
                                        div_b = div_b, # bio based
                                        treatment = treatment))
      
    }
  }
  all_treat_bio_herb_data <- rbind(all_treat_bio_herb_data,
                                   treat_bio_herb_data)
}


# Biomass of plants vs species richness of arthropods
# only species with more than three observations within the treatment
atbhd <- all_treat_bio_herb_data

# ggplot(atbhd, aes(x = log(biomass), 
#                           y=rich))+
#   geom_point()+
#   geom_smooth(method = "lm")+
#   facet_wrap(~treatment)

more_than_x <- tapply(atbhd$species, atbhd$treatment, function(x){which(table(x)>=3)} )

filtered_data <- data.frame()
for(nm in names(more_than_x)){
  c1 <- atbhd$treatment == nm
  c2 <- atbhd$species %in% names(more_than_x[[nm]])
  filtered_data <- rbind(filtered_data,
                         atbhd[c1 & c2, ])
}

# ggplot(filtered_data, aes(x = log(biomass), 
#                           y=rich, 
#                           colour = species))+
#   geom_point()+
#   geom_smooth(method = "lm")+
#   facet_wrap(~treatment)

# ggplot(filtered_data, aes(x = log(biomass), 
#                           y=log(val), 
#                           colour = species))+
#   geom_point()+
#   geom_smooth(method = "lm")+
#   facet_wrap(~treatment)

# That needs to be broken into treatments.

# For each species and for each treatment I can get an estimate of the slope parameter
filtered_slope <- filtered_data[filtered_data$species %in% c("melamu", "piptar"),]

# slope_estims <- data.frame()
# # Forr abundnace 
# for(spc in c("melamu", "piptar")){
#   trt_dat <- data.frame()
#   for(nm in unique(filtered_slope$treatment)){
#     subdat <- filtered_slope[filtered_slope$treatment == nm, ]
#     subdat_pl <- subdat[subdat$species == spc, ]
#     if(dim(subdat_pl)[1] == 0){
#       next
#     }
#     # Model for log-log linear relationship plant bio and invertebrate abundance
#     mod <- lm(log(bio)~log(biomass), 
#               data=subdat_pl)
#     results <- data.frame(t(summary(mod)$coefficients[2, ]),
#                           species = spc,
#                           treatment = nm,
#                           size = dim(subdat_pl)[1])
#     trt_dat <- rbind(trt_dat, results)
#   }
#   slope_estims <- rbind(slope_estims, trt_dat)
# }

# p<- ggplot(slope_estims, aes(x=treatment, 
#                              y=Estimate, 
#                              group=species, 
#                              color=species,
#                              label = size)) + 
#   geom_line() +
#   geom_point()+
#   geom_text(hjust = 3, nudge_x = 0.1,color = "black")+
#   geom_errorbar(aes(ymin=Estimate-Std..Error, 
#                     ymax=Estimate+Std..Error), width=.2,
#                 position=position_dodge(0.05))
# p


# slope_estims <- data.frame()
# # For diversity 
# for(spc in c("melamu", "piptar")){
#   trt_dat <- data.frame()
#   for(nm in unique(filtered_slope$treatment)){
#     subdat <- filtered_slope[filtered_slope$treatment == nm, ]
#     subdat_pl <- subdat[subdat$species == spc, ]
#     if(dim(subdat_pl)[1] == 0){
#       next
#     }
#     # Model for log-log linear relationship plant bio and invertebrate abundance
#     mod <- lm(rich~log(biomass), 
#               data=subdat_pl)
#     results <- data.frame(t(summary(mod)$coefficients[2, ]),
#                           species = spc,
#                           treatment = nm,
#                           size = dim(subdat_pl)[1])
#     trt_dat <- rbind(trt_dat, results)
#   }
#   slope_estims <- rbind(slope_estims, trt_dat)
# }

# p<- ggplot(slope_estims, aes(x=treatment, 
#                              y=Estimate, 
#                              group=species, 
#                              color=species,
#                              label = size)) + 
#   geom_line() +
#   geom_point()+
#   geom_text(hjust = 3, nudge_x = 0.1,color = "black")+
#   geom_errorbar(aes(ymin=Estimate-Std..Error, 
#                     ymax=Estimate+Std..Error), width=.2,
#                 position=position_dodge(0.05))
# p

# Species identity
filtered_slope$block <- substr(filtered_slope$plot,3,4)
filtered_data$block <- substr(filtered_data$plot,3,4)

library(lme4)
library(lmerTest)

# Invertebrate abundance vs species identity
# glmer1 <- glmer.nb(val~log(biomass)+species+(1|block), data=filtered_slope)
# glmer2 <- lmer(log(val)~log(biomass)+species+(1|block), data=filtered_slope)
# glmer3a <- lm(log(val)~log(biomass)+species, data=filtered_slope)

# lm3 <- lm(log(val)~log(biomass), 
#            data=filtered_slope)
# lm3a <- lm(log(val)~log(biomass)+species, 
#            data=filtered_slope)

# For all possible species
# Biomass ----
lmer_random <- lmer(log(bio)~log(biomass)+species+(1|block), data=filtered_data)
lmer_norand <- lm(log(bio)~log(biomass)+species, data=filtered_data)

AIC(lmer_random, lmer_norand)
# Random factor is not necessarry

# Lets see wether biomass predicts herbivore biomass better
lm_bio <- lm(log(bio)~log(biomass), 
                data=filtered_data)
lm_spc <- lm(log(bio)~species, 
             data=filtered_data)

lm_bio_spc <- update(lm_bio, .~.+species)
anova(lm_bio, lm_bio_spc, test = "Chisq")

lm_spc_bio <- update(lm_spc, .~.+log(biomass))
anova(lm_spc, lm_spc_bio, test = "Chisq")


# By the way, this same approach is proposed by N. W. Galwey in 'Introduction to Mixed Modelling: Beyond Regression and Analysis of Variance' on pages 213-214.

# Species identity is an important part!
# Biomass itsef is not explaining all the variability

# Do species differ in their slopes?
# Is there an interaction?
lm_int_bio_sp <-update(lm_bio_spc, .~. + log(biomass)*species)
anova(lm_bio_spc, lm_int_bio_sp)
# There is no signigficant interaction there.

# What if the slope is random? Unconstrained slopes
random_slope <- lmer(log(bio)~biomass+species+(0+biomass|species), data = filtered_data)
summary(random_slope)

# Form the graph it deasn't look like it.
# Variance doesn't seem to be different than zero
lm_int_bio_sp_trt <-update(lm_int_bio_sp, .~. + treatment + treatment*species*log(biomass))
anova(lm_int_bio_sp, lm_int_bio_sp_trt)


# Lets see whether treatement can affect these accumulation curves
ggplot(filtered_data, aes(x = log(biomass), y=log(bio), color = species))+
  geom_jitter()+
  geom_smooth(method = "lm", se=FALSE)+
  facet_wrap(~treatment)
  


lm_int <- lm(log(val)~log(biomass)+species*treatment, 
             data=filtered_slope)

anova(lm_no_int, lm_int, test="Chisq")

lm_no_int <- lm(log(bio)~log(biomass)+species+treatment, 
           data=filtered_slope)
lm_int <- lm(log(val)~log(biomass)+species*treatment, 
              data=filtered_slope)

anova(lm_no_int, lm_int, test="Chisq")
# Interaction term is not significant!




# filtered_slope$trtspec <- paste(filtered_slope$treatment, filtered_slope$species)
# lm3c <- lm(log(val)~log(biomass)+trtspec, 
#               data=filtered_slope)
# 
# summary(glmer3b)
# 
# inter.test <- emmeans(glmer3c, "trtspec")
# pairwise <- cld(inter.test, Letter="abcdefghijklm")


# Abundance vs species identity ----
glmer1 <- glm.nb(abu~log(biomass)+species, data=filtered_data)
# nb model with random effect does not converge

glmer2 <- lmer(log(abu)~log(biomass)+species+(1|block), data=filtered_data)
glmer3a <- lm(log(abu)~log(biomass)+species, data=filtered_data)

AIC(glmer1, glmer2)
AIC(glmer2, glmer3a)
anova(glmer2, glmer3a)
# Negative binomial doesnt look good
# Model withouth random effect of block performs better.

summary(glmer3a)
glmer3b <- lm(log(abu)~log(biomass)*species, data=filtered_data)

anova(glmer3b, glmer3a)
# There is a marginal signifdicant interaction term in case of abundnace

ggplot(filtered_data, aes(x = log(biomass), y=log(abu), color = species))+
  geom_jitter()+
  geom_smooth(method = "lm", se=FALSE)+
  facet_wrap(~treatment)

# Richness ----

# should be Poisson?
lmer_rich <- lmer(rich~log(biomass)+species+(1|block), data=filtered_data)


lmer_rich_norand <- lm(rich~log(biomass)+species, data=filtered_data)

AIC(lmer_rich, lmer_rich_norand) # Here random is better
anova(lmer_rich, lmer_rich_norand) # And here is not :/

lm_rich_norand_int <- lm(rich~log(biomass)*species, data=filtered_data)

anova(lmer_rich_norand, lm_rich_norand_int) # no significant improvement

lm_rich_norand_trt <- lm(rich~log(biomass)+species+treatment, data=filtered_data)
anova(lmer_rich_norand, lm_rich_norand_trt) # teatment not significant
lm_rich_norand_trt_int <- lm(rich~log(biomass)*treatment+species,
                             data=filtered_data)
anova(lmer_rich_norand, lm_rich_norand_trt_int) # interaction not significant

lm_rich_norand_bio <- lm(rich~log(biomass), data=filtered_data)
lm_rich_norand_bio_sp <- update(lm_rich_norand_bio, .~.+species)
anova(lm_rich_norand_bio, lm_rich_norand_bio_sp) # species is significant but only for piptar

# final simplest model approximately 70% of variation
summary(lm_rich_norand_bio_sp)

ggplot(filtered_data, aes(x = log(biomass), y=rich, color = species))+
  geom_jitter()+
  geom_smooth(method = "lm", se=FALSE)

## glmer poisson

lmer_rich <- glmer(rich~log(biomass)+species+(1|block),family = "poisson", data=filtered_data)

lmer_rich_norand <- glm(rich~log(biomass)+species, 
                        family="poisson", 
                        data=filtered_data)

AIC(lmer_rich, lmer_rich_norand) 
anova(lmer_rich, lmer_rich_norand) #Here random is better

# interaction with random effect is to complex for simple dataset.
lm_rich_norand_int <- update(lmer_rich_norand, .~ log(biomass)*species)

anova(lmer_rich_norand, lm_rich_norand_int, test="Chisq") # marginally significant improvement

lm_rich_norand_trt <- lm(rich~log(biomass)+species+treatment, data=filtered_data)
anova(lmer_rich_norand, lm_rich_norand_trt) # teatment not significant
lm_rich_norand_trt_int <- lm(rich~log(biomass)*treatment+species,
                             data=filtered_data)
anova(lmer_rich_norand, lm_rich_norand_trt_int) # interaction not significant

lm_rich_norand_bio <- lm(rich~log(biomass), data=filtered_data)
lm_rich_norand_bio_sp <- update(lm_rich_norand_bio, .~.+species)
anova(lm_rich_norand_bio, lm_rich_norand_bio_sp) # species is significant but only for piptar

# final simplest model approximately 70% of variation
summary(lm_rich_norand_bio_sp)

ggplot(filtered_data, aes(x = log(biomass), y=rich, color = species))+
  geom_jitter()+
  geom_smooth(method = "lm", se=FALSE)

# Diversity ----

lmer_div <- lmer(rich~log(biomass)+species+(1|block), data=filtered_data)
lmer_rich_norand <- lm(rich~log(biomass)+species, data=filtered_data)

AIC(lmer_rich, lmer_rich_norand) # Here random is better
anova(lmer_rich, lmer_rich_norand) # And here is not :/

lm_rich_norand_int <- lm(rich~log(biomass)*species, data=filtered_data)

anova(lmer_rich_norand, lm_rich_norand_int) # no significant improvement

lm_rich_norand_trt <- lm(rich~log(biomass)+species+treatment, data=filtered_data)
anova(lmer_rich_norand, lm_rich_norand_trt) # teatment not significant
lm_rich_norand_trt_int <- lm(rich~log(biomass)*treatment+species,
                             data=filtered_data)
anova(lmer_rich_norand, lm_rich_norand_trt_int) # interaction not significant

lm_rich_norand_bio <- lm(rich~log(biomass), data=filtered_data)
lm_rich_norand_bio_sp <- update(lm_rich_norand_bio, .~.+species)
anova(lm_rich_norand_bio, lm_rich_norand_bio_sp) # species is significant but only for piptar

# final simplest model approximately 70% of variation
summary(lm_rich_norand_bio_sp)

ggplot(filtered_data, aes(x = log(biomass), y=rich, color = species))+
  geom_jitter()+
  geom_smooth(method = "lm", se=FALSE)

# Notes ----
# glmer3b <- lm(log(val)~log(biomass)+species+treatment+species*treatment, 
#               data=filtered_slope)
# lm3d <- lm(rich~log(biomass)+trtspec, 
#               data=filtered_slope)
# 
# inter.test <- emmeans(glmer3c, "trtspec")
# pairwise <- cld(inter.test, Letter="abcdefghijklm")
# 
# plot(pairwise)
 
# ltrs <- data.frame(pt = pairwise$treat,
#                    pg = pairwise$.group)

plot(glmer1)
plot(glmer2)

# Log abundance model is better, no need for random effec neither
anova(glmer1, glmer2)
AIC(glmer1, glmer2, glmer3a, glmer3b,glmer3c)
)

library(insight)
library(MuMIn)
library(rsq)

rsq(glmer3a)
rsq(glmer3b)

rsq.partial(glmer3a)
rsq.partial(glmer3b)

# Both significant but still 

r.squaredGLMM(glmer1)
get_variance_fixed(glmer1)
get_variance_distribution(glmer1)
R2(glmer1)

ggplot(filtered_slope, aes(x = species, y=log(val)))+
  geom_jitter()+
  facet_wrap(~treatment)
