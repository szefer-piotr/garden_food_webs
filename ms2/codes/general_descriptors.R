# general_descriptors

source("codes/data_preparation.R")
source("codes/functions.R")

library(ggplot2)
library(dplyr)
library(lme4)
library(emmeans)
library(multcomp)
library(lmerTest)

# Plot
gd <- ggplot(genDescriptorsPlot, aes(x = treat, y= log(val), 
                                     group = guild))+
  geom_jitter(width = 0.1, aes(colour = block))+
  facet_grid(cols = vars(guild), 
             rows = vars(ind),
             scales = "free")

gd

# 1. Tests ----

gdp <- genDescriptorsPlot
gdp$tukey <- "none"

## 1.1 Herbivore bio
Conditions <- gdp$guild == "herbivore" & gdp$ind == "bio"

# Log link gaussian.
glmer1 <- glmer(val/1000~treat+(1|block),family = gaussian(link = "log"),
               data = gdp[Conditions, ],
               control=glmerControl(optimizer="bobyqa",
                                    optCtrl=list(maxfun=2e5)))
summary(glmer1)

tg <- TukeyGroups(glmer1, gdp[Conditions, ])
gdp[Conditions, ] <- AppendGroups(tg, gdp[Conditions, ])

## 1.2 AP bio
Conditions <- gdp$guild == "art_pred" & gdp$ind == "bio"

# lmer on logged values quick and dirty.
# lm1 <- lm(log(val)~treat,
#               data = gdp[Conditions, ])

glmer2 <- glmer(val/1000~treat+(1|block),family = gaussian(link = "log"),
                data = gdp[Conditions, ],
                control=glmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=2e5)))

# lmer2 <- lmer(log(val)~treat+(1|block),
#                 data = gdp[Conditions, ],
#                 control=lmerControl(optimizer="bobyqa",
#                                      optCtrl=list(maxfun=2e5)))
summary(glmer2)

tg <- TukeyGroups(lmer2, gdp[Conditions, ])
gdp[Conditions, ] <- AppendGroups(tg, gdp[Conditions, ])
