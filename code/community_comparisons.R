# Community descriptors of intermediate predators, herbivores and plants.
rm(list=ls())
source("code/data_processing_code.R")



# [1] "FUNGICIDE"   [2] "WEEVIL125"   
# [3] "CONTROL"     [4] "PREDATOR"   
# [5] "WEEVIL25"    [6] "INSECTICIDE"

# treats_to_plot <- as.character(unique(treats$treat))[c(6,3,4,5,2)]
treats_to_plot <- as.character(unique(treats$treat))[c(3,4)]


# Indicate insect order here ----
# order <- c("cole")

# ins_bio <- ins_bio[ins_bio$family %in% order, ]

# Herbivores and IPS diversity, abumdnace, richness in different treatments
ins_bio$group <- "Herbivore"
ins_bio$group[grep("aran|mant", ins_bio$family)] <- "Intermediate predator"

# Test
# unique(ins_bio[ins_bio$group == "Intermediate predator", ]$morphotype)

# Abundnaces
herbAbundance <- tapply(ins_bio[ins_bio$group == "Herbivore", ]$amount, 
                        ins_bio[ins_bio$group == "Herbivore", ]$plot, 
                        sum, na.rm=T)
ipsAbundance <- tapply(ins_bio[ins_bio$group == "Intermediate predator", ]$amount, 
                    ins_bio[ins_bio$group == "Intermediate predator", ]$plot, 
                    sum, na.rm=T)

# Biomass

#Why there are so many NAs in herbivore biomass?
ins_bio[is.na(ins_bio$totbio), ]$totbio <- 0
# because there are no estimates of biomass for these species :/

herbBio <- tapply(ins_bio[ins_bio$group == "Herbivore", ]$totbio, 
                  ins_bio[ins_bio$group == "Herbivore", ]$plot, 
                  sum, na.rm=T)

ipsBio <- tapply(ins_bio[ins_bio$group == "Intermediate predator", ]$totbio, 
                  ins_bio[ins_bio$group == "Intermediate predator", ]$plot, 
                  sum, na.rm=T)

# Diversity
herbDivBio <- tapply(ins_bio[ins_bio$group == "Herbivore", ]$totbio, 
                  ins_bio[ins_bio$group == "Herbivore", ]$plot, 
                  vegan::diversity)

ipsDivBio <- tapply(ins_bio[ins_bio$group == "Intermediate predator", ]$totbio, 
                 ins_bio[ins_bio$group == "Intermediate predator", ]$plot, 
                 vegan::diversity)

# Richness [test this!!!]
herbRich <- tapply(ins_bio[ins_bio$group == "Herbivore", ]$morphotype, 
                        ins_bio[ins_bio$group == "Herbivore", ]$plot, 
                        length)
ipsRich <- tapply(ins_bio[ins_bio$group == "Intermediate predator", ]$morphotype, 
                       ins_bio[ins_bio$group == "Intermediate predator", ]$plot, 
                       length)

descriptorDf <- data.frame(plot = names(herbAbundance), 
                           habund = herbAbundance,
                           ipabund = ipsAbundance,
                           hbio = herbBio,
                           ipsbio = ipsBio,
                           hrich = herbRich,
                           iprich = ipsRich,
                           hdivbio = herbDivBio,
                           ipsdivbio = ipsDivBio)
rownames(treats) <- tolower(rownames(treats))
descriptorDf$treat <- treats[descriptorDf$plot, ]$treat

# ggplot dataset

panelData <- data.frame(value = c(herbAbundance,
                                  ipsAbundance,
                                  herbBio,
                                  ipsBio,
                                  herbRich,
                                  ipsRich,
                                  herbDivBio,
                                  ipsDivBio),
                        group = rep(c("Herbivore","IPs",
                                    "Herbivore","IPs",
                                    "Herbivore","IPs",
                                    "Herbivore","IPs"),each = 36),
                        descriptor = rep(c("Abundance",
                                           "Biomass",
                                           "Richness",
                                           "Diversity"), each = 2*36),
                        type = rep(c("Herbivore abundance",
                                 "IPs abundance",
                                 "Herbivore biomass",
                                 "IPs biomass",
                                 "Herbivore richness",
                                 "IPs richness",
                                 "Herbivore diversity (biomass based)",
                                 "IPs diversity (biomass bases)"), each = 36),
                        treat = rep(treats[descriptorDf$plot, ]$treat, 8),
                        plots = rep(treats[descriptorDf$plot, ]$codes, 8))
panelData$block <- substr(panelData$plots,3,4)
panelData$tukey <- "" # set up column for labels from Tukey test

# panelData

# 1. Ready ----
library(lme4)
library(MASS)
library(lmerTest)
library(emmeans)
library(multcomp)
library(AER)


# glmer

### Biomass ----

# Group
herb <- panelData$group == "Herbivore"
ips <- panelData$group == "IPs"

# Descriptors
bio <- panelData$descriptor == "Biomass"
abu <- panelData$descriptor == "Abundance"
rich <- panelData$descriptor == "Richness"
div <- panelData$descriptor == "Diversity"

# Selected treatments
trtsel <- panelData$treat %in% treats_to_plot

# Function that performs the test and appends results into the panelData
# analyzeAndAppend <- function(Conditions,
#                              FUN,
#                              Form = formula(value~treat+(1|block)),
#                              Fam = gaussian(link="log"), ...){
# 
#   #Call the function
#   model <- FUN(formula = Form, family = Fam,
#                         data = panelData[Conditions, ], ...)
# 
#   # Get model summary and pairwise comparisons
#   summary(model)
#   inter.test <- emmeans(model, "treat")
#   pairwise <- cld(inter.test, Letter="abcdefghijklm")
# 
#   # Extract letters and append to the dataset
#   ltrs <- data.frame(pt = pairwise$treat,
#                      pg = pairwise$.group)
#   rownames(ltrs) <- ltrs$pt
#   labs <- ltrs[as.character(panelData[condition, "treat"]),]$pg
#   panelData[condition, ]$tukey <<- as.character(labs)
# 
#   # Print model results
#   print(summary(model))
# }

# FUN = glmer
# formula = value~treat+(1|block)
# family = gaussian(link="log")
# Conditions = (herb & bio & trtsel)
# model <- FUN(data = panelData[Conditions, ], formula=formula, family=family)  


analyzeAndAppend2 <- function(Conditions, 
                             FUN,
                             dispTest = F,
                             ...){
  
  # Specific function that analyses the chunk of my data and
  # appends results of a Tukey post hoc test to panelData
  
  #Call the function
  model <- FUN(..., data = panelData[Conditions, ])
  
  print("model calculated")
  if(dispTest){
    print(dispersiontest(model,trafo=1))
  }
  
  # Get model summary and pairwise comparisons
  inter.test <- emmeans(model, "treat", 
                        data = panelData[Conditions, ])
  print("inter.test passed")
  pairwise <- cld(inter.test, Letter="abcdefghijklm")
  
  print("pairwise computed")
  
  # Extract letters and append to the dataset
  ltrs <- data.frame(pt = pairwise$treat,
                     pg = pairwise$.group)
  rownames(ltrs) <- ltrs$pt
  labs <- ltrs[as.character(panelData[Conditions, "treat"]),]$pg
  panelData[Conditions, ]$tukey <<- as.character(labs)
  
  print("data appended")
  
  # Print model results
  print(summary(model))
}


# Analyses

## 1A. Herbivore - biomass ----
analyzeAndAppend2(Conditions = (herb & bio & trtsel),
                 FUN = glmer,
                 formula = value~treat+(1|block),
                 family = gaussian(link="log"),
                 control=glmerControl(optimizer="bobyqa",
                                      optCtrl=list(maxfun=2e5)))

## 1B. IPs - biomass ----
# analyzeAndAppend2( 
#   Conditions = (ips & bio & trtsel),
#   FUN = glmer,
#   formula = value~treat+(1|block),
#   family = gaussian(link="log"))

analyzeAndAppend2( 
  Conditions = (ips & bio & trtsel),
  FUN = lmer,
  formula = value~treat+(1|block))

## 2A. Herbivore - Abundance ----
# no random effects because I was getting singular var-cov matrix
analyzeAndAppend2(
  Conditions = (herb & abu & trtsel),
  FUN = glm.nb,
  formula = value~treat)

## 2B. IPs - Abundnace ----
# no random effects because I was getting singular var-cov matrix
# Barr DJ, Levy R, Scheepers C, Tily HJ. Random effects structure for confirmatory hypothesis testing: Keep it maximal. Journal of Memory and Language, 68(3):255–278, April 2013.
analyzeAndAppend2(
  Conditions = (ips & abu & trtsel),
  FUN = glm.nb,
  formula = value~treat)

## 3A. Herbivore - Richness ----
analyzeAndAppend2(
  Conditions = (herb & rich & trtsel),
  FUN = glmmPQL,
  family = quasipoisson(),
  random = ~ 1|block,
  fixed = value~treat,
  dispTest = F)

## 3B. IPs - Richness ----
# overdispersed random effect models
analyzeAndAppend2(
  Conditions = (ips & rich & trtsel),
  FUN = glmmPQL,
  family = quasipoisson(),
  random = ~ 1|block,
  fixed = value~treat)

## 4A. Herbivore - diversity SW ----
analyzeAndAppend2(Conditions = (herb & div & trtsel),
                  FUN = lmer,
                  formula = value~treat+(1|block))

library(nlme)
analyzeAndAppend2(Conditions = (herb & div & trtsel),
                  FUN = lme,
                  fixed = value~treat,
                  random= ~1|block)


## 4B. IPs - diversity SW ----
analyzeAndAppend2(Conditions = (ips & div & trtsel),
                  FUN = lm,
                  formula = value~treat)

# I should check also wether these models fit data well.
# Check for overdispersion
# data(RecreationDemand)
# rd <- glm(trips ~ ., data = RecreationDemand, family = poisson)
# dispersiontest(rd,trafo=1)

# ### Individual analyses ----
# c1 <- panelData$group == "IPs"
# 
# condition <- (c1 & c2 & c3)
# ips_bio_glmer <- glmer(value~treat+(1|block), family = gaussian(link="log"),
#                        data = panelData[condition, ])
# summary(ips_bio_glmer)
# inter.test <- emmeans(ips_bio_lmer, "treat")
# pairwise <- cld(inter.test, Letter="abcdefghijklm")
# 
# ltrs <- data.frame(pt = pairwise$treat,
#                    pg = pairwise$.group)
# rownames(ltrs) <- ltrs$pt
# labs <- ltrs[as.character(panelData[condition, "treat"]),]$pg
# panelData[condition, ]$tukey <- as.character(labs)
# 
# # Abundnace
# c1 <- panelData$group == "Herbivore"
# c2 <- panelData$descriptor == "Abundance"
# c3 <- panelData$treat %in% treats_to_plot
# 
# # Herbivore
# condition <- (c1 & c2 & c3)
# herb_abu_glmer <- glmer(value~treat+(1|block), family = poisson(),
#                         data = panelData[condition, ])
# herb_abu_quasi <- glmmPQL(value~treat, random = ~ 1|block, family = quasipoisson(),
#                         data = panelData[condition, ])
# summary(herb_abu_glmer)
# summary(herb_abu_quasi)
# 
# # Which one is better?
# # This should be 1: 2053.9/14
# 
# inter.test <- emmeans(herb_abu_quasi, "treat")
# pairwise <- cld(inter.test, Letter="abcdefghijklm")
# 
# ltrs <- data.frame(pt = pairwise$treat,
#                    pg = pairwise$.group)
# rownames(ltrs) <- ltrs$pt
# labs <- ltrs[as.character(panelData[condition, "treat"]),]$pg
# panelData[condition, ]$tukey <- as.character(labs)
# 
# # IPs
# c1 <- panelData$group == "IPs"
# condition <- (c1 & c2 & c3)
# ips_abu_glmer <- glmer(value~treat+(1|block), family = poisson(),
#                         data = panelData[condition, ])
# ips_abu_quasi <- glmmPQL(value~treat, random = ~ 1|block, family = quasipoisson(),
#                           data = panelData[condition, ])
# summary(ips_abu_glmer)
# summary(ips_abu_quasi)
# 
# # Which one is better?
# # This should be 1: 328.1/13
# 
# inter.test <- emmeans(ips_abu_quasi, "treat")
# pairwise <- cld(inter.test, Letter="abcdefghijklm")
# 
# ltrs <- data.frame(pt = pairwise$treat,
#                    pg = pairwise$.group)
# rownames(ltrs) <- ltrs$pt
# labs <- ltrs[as.character(panelData[condition, "treat"]),]$pg
# panelData[condition, ]$tukey <- as.character(labs)
# 
# # Diversity
# 
# # Herbivores
# c1 <- panelData$group == "Herbivore"
# c2 <- panelData$descriptor == "Diversity"
# c3 <- panelData$treat %in% treats_to_plot
# condition <- (c1 & c2 & c3)
# herb_div_lmer <- lmer(value~treat+(1|block), 
#                        data = panelData[condition, ])
# summary(herb_div_lmer)
# inter.test <- emmeans(herb_div_lmer, "treat")
# pairwise <- cld(inter.test, Letter="abcdefghijklm")
# 
# ltrs <- data.frame(pt = pairwise$treat,
#                    pg = pairwise$.group)
# rownames(ltrs) <- ltrs$pt
# labs <- ltrs[as.character(panelData[condition, "treat"]),]$pg
# panelData[condition, ]$tukey <- as.character(labs)
# 
# # IPs
# c1 <- panelData$group == "IPs"
# condition <- (c1 & c2 & c3)
# ips_div_lm <- lm(value~treat, 
#                       data = panelData[condition, ])
# summary(ips_div_lm)
# inter.test <- emmeans(ips_div_lm, "treat")
# pairwise <- cld(inter.test, Letter="abcdefghijklm")
# 
# ltrs <- data.frame(pt = pairwise$treat,
#                    pg = pairwise$.group)
# rownames(ltrs) <- ltrs$pt
# labs <- ltrs[as.character(panelData[condition, "treat"]),]$pg
# panelData[condition, ]$tukey <- as.character(labs)
# 
# # Richness
# 
# # Herbivore
# c1 <- panelData$group == "Herbivore"
# c2 <- panelData$descriptor == "Richness"
# c3 <- panelData$treat %in% treats_to_plot
# condition <- (c1 & c2 & c3)
# herb_rich_glmer <- glmer(value~treat+(1|block), family = poisson(),
#                         data = panelData[condition, ])
# herb_rich_quasi <- glmmPQL(value~treat, random = ~ 1|block, family = quasipoisson(),
#                           data = panelData[condition, ])
# summary(herb_rich_glmer)
# summary(herb_rich_quasi)
# inter.test <- emmeans(herb_rich_quasi, "treat")
# pairwise <- cld(inter.test, Letter="abcdefghijklm")
# 
# ltrs <- data.frame(pt = pairwise$treat,
#                    pg = pairwise$.group)
# rownames(ltrs) <- ltrs$pt
# labs <- ltrs[as.character(panelData[condition, "treat"]),]$pg
# panelData[condition, ]$tukey <- as.character(labs)
# 
# # IPs
# c1 <- panelData$group == "IPs"
# condition <- (c1 & c2 & c3)
# ips_rich_glmer <- glmer(value~treat+(1|block), family = poisson(),
#                          data = panelData[condition, ])
# ips_rich_quasi <- glmmPQL(value~treat, random = ~ 1|block, family = quasipoisson(),
#                            data = panelData[condition, ])
# summary(ips_rich_glmer)
# summary(ips_rich_quasi)
# inter.test <- emmeans(ips_rich_quasi, "treat")
# pairwise <- cld(inter.test, Letter="abcdefghijklm")



# 2. Plots ----

treats_to_plot
library(ggplot2)
panelData$treat <- factor(panelData$treat,
                          levels = treats_to_plot)

ann_text = 

panelData[panelData$treat %in% treats_to_plot,] -> pD
ggplot(pD, 
       aes(x = treat, 
           y = value,
           group = group,
           label = tukey))+
  geom_jitter(width = 0.1, col = rgb(10,10,10,80,maxColorValue = 255))+
  facet_wrap(~type, ncol=4, scales = "free") + 
  stat_summary(fun.data=mean_cl_boot, 
               geom="pointrange", lwd=0.5) +
  # Significance of herbivore biomass
  stat_summary(data = pD[pD$type == "Herbivore biomass",], 
               fun.data=mean_cl_boot, 
               geom="pointrange", 
               lwd=0.8, col = c("black", "red")) +
  stat_summary(data = pD[pD$type == "Herbivore diversity (biomass based)",], 
               fun.data=mean_cl_boot, 
               geom="pointrange", 
               lwd=0.8, col = c("black", "gold"))+
  theme_bw()
  # stat_summary(data = pD[pD$type == "Herbivore biomass",],
  #              fun=mean, geom="text", hjust = 2,
  #              vjust = 1)
  
  scale_color_manual(values = c("black", "black",
                               "black", "red",
                               "black", "gold",
                               "black", "black",
                               "black", "black",
                               "black", "black",
                               "black", "black",
                               "black", "black"))

                    