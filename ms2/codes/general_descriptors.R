# general_descriptors

source("codes/data_preparation.R")
source("codes/functions.R")

library(ggplot2)
library(dplyr)
library(lme4)
library(emmeans)
library(multcomp)
library(lmerTest)
library(scales)

# Plot
# gd <- ggplot(genDescriptorsPlotDens, aes(x = treat, y= val))+
#   geom_jitter(width = 0.05, alpha=0.1)+
#   stat_summary(geom="pointrange", fun.data = "mean_cl_boot")+
#   facet_grid(rows = vars(guild),
#              cols = vars(ind),
#              scales = "free") +
#   theme_bw()
# 
# gd

# 1. Tests ----

gdp <- genDescriptorsPlotDens

## 1.1 Herbivore bio ----
Conditions <- gdp$guild == "herbivore" & gdp$ind == "bio"

# Log link gaussian.
gdp$treat <- factor(gdp$treat,
                    levels = c("C","I","H1","H2"))

glmer1 <- glmer(val~treat+(1|block),family = gaussian(link = "log"),
               data = gdp[Conditions, ],
               control=glmerControl(optimizer="bobyqa",
                                    optCtrl=list(maxfun=2e5)))
summary(glmer1)

tg <- TukeyGroups(glmer1, gdp[Conditions, ])

# gdp[Conditions, ] <- AppendGroups(tg, gdp[Conditions, ])

yvalue <- aggregate(.~treat, 
                    data=gdp[Conditions, c("treat","val")], 
                    mean)
tdat <- cbind(yvalue, tg[,c(2,5,6,7)])


gdp$treat <- factor(gdp$treat,
                    levels = c("I","C","H1","H2"))

herb_bio <- ggplot(gdp[Conditions,], 
                   aes(x = treat, y = val))+
  scale_y_continuous(trans = "log10",
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))+
  stat_summary(geom="pointrange", 
               fun.data = "mean_cl_boot")+
  geom_text(data =  tdat,
            aes(label = .group),vjust=-3.5,hjust=-.5)+
  geom_jitter(width = 0.05, 
              alpha=0.1)+
  theme_bw()

# herb_bio  


# Ordinal scale
# gdp$ord_trt <- as.numeric(factor(gdp$treat,
#                     levels = c("I","C","H1","H2"),
#                     labels = 1,2,3,4))
# glmer1_ord <- glmer(val~ord_trt+(1|block),family = gaussian(link = "log"),
#                 data = gdp[Conditions, ],
#                 control=glmerControl(optimizer="bobyqa",
#                                      optCtrl=list(maxfun=2e5)))
# summary(glmer1_ord)

# Non-parametric: not significant
# cor.test(gdp[Conditions, ]$ord_trt, gdp[Conditions, ]$val)


## 1.2 AP bio ----
Conditions <- gdp$guild == "art_pred" & gdp$ind == "bio"

gdp$treat <- factor(gdp$treat,
                    levels = c("C","I","H1","H2"))

# apbio1 <- glmer(val~treat+(1|block),family = gaussian(link = "log"),
#                 data = gdp[Conditions, ],
#                 control=glmerControl(optimizer="bobyqa",
#                                      optCtrl=list(maxfun=2e5)))

# Dropped random effect because of convergence problems
apbio2 <- glm(val~treat,
              family = gaussian(link = "log"),
              data = gdp[Conditions, ])

summary(apbio2)

# apbio3 <- nlme::lme(log(val)~treat, random = ~1|block,
#                     data = gdp[Conditions, ])
# 
# summary(apbio3)

tg <- TukeyGroups(apbio2, gdp[Conditions, ])

# gdp[Conditions, ] <- AppendGroups(tg, gdp[Conditions, ])

yvalue <- aggregate(.~treat, 
                    data=gdp[Conditions, c("treat","val")], 
                    mean)
tdat <- cbind(yvalue, 
              tg[order(tg$treat),c(2,5,6,7)])

gdp$treat <- factor(gdp$treat,
                    levels = c("I","C","H1","H2"))

ap_bio <- ggplot(gdp[Conditions,], 
                   aes(x = treat, y = val))+
  scale_y_continuous(trans = "log10",
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))+
  stat_summary(geom="pointrange", 
               fun.data = "mean_cl_boot")+
  geom_text(data =  tdat,
            aes(label = .group),vjust=-3.5,hjust=-.5)+
  geom_jitter(width = 0.05, 
              alpha=0.1)+
  theme_bw()

# ap_bio  

# 1.3 Herbivore abundance ----
Conditions <- gdp$guild == "herbivore" & gdp$ind == "abu"
gdp$treat <- factor(gdp$treat,
                    levels = c("C","I","H1","H2"))

# habu1 <- glmer.nb(val~treat+(1|block),
#                 data = gdp[Conditions, ],
#                 control=glmerControl(optimizer="bobyqa",
#                                      optCtrl=list(maxfun=2e10)))
# convergence issues

habu2 <- glm.nb(val~treat,
                data = gdp[Conditions, ])

summary(habu2)

tg <- TukeyGroups(habu2, gdp[Conditions, ])

# gdp[Conditions, ] <- AppendGroups(tg, gdp[Conditions, ])

yvalue <- aggregate(.~treat, 
                    data=gdp[Conditions, c("treat","val")], 
                    mean)
tdat <- cbind(yvalue, 
              tg[order(tg$treat),c(2,5,6,7)])

gdp$treat <- factor(gdp$treat,
                    levels = c("I","C","H1","H2"))

herb_abu <- ggplot(gdp[Conditions,], 
                 aes(x = treat, y = val))+
  scale_y_continuous(trans = "log10",
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))+
  stat_summary(geom="pointrange", 
               fun.data = "mean_cl_boot")+
  geom_text(data =  tdat,
            aes(label = .group),vjust=-3.5,hjust=-.5)+
  geom_jitter(width = 0.05, 
              alpha=0.1)+
  theme_bw()

# herb_abu  

# 1.4 AP abundance ----
Conditions <- gdp$guild == "art_pred" & gdp$ind == "abu"
gdp$treat <- factor(gdp$treat,
                    levels = c("C","I","H1","H2"))

# apabu1 <- glmer.nb(val~treat+(1|block),
#                 data = gdp[Conditions, ],
#                 control=glmerControl(optimizer="bobyqa",
#                                      optCtrl=list(maxfun=2e10)))
# convergence issues, degenerate Hessian

apabu2 <- glm.nb(val~treat,
                data = gdp[Conditions, ])

summary(apabu2)

tg <- TukeyGroups(apabu2, gdp[Conditions, ])

yvalue <- aggregate(.~treat, 
                    data=gdp[Conditions, c("treat","val")], 
                    mean)
tdat <- cbind(yvalue, 
              tg[order(tg$treat),c(2,5,6,7)])

gdp$treat <- factor(gdp$treat,
                    levels = c("I","C","H1","H2"))

ap_abu <- ggplot(gdp[Conditions,], 
                   aes(x = treat, y = val))+
  scale_y_continuous(trans = "log10",
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))+
  stat_summary(geom="pointrange", 
               fun.data = "mean_cl_boot")+
  geom_text(data =  tdat,
            aes(label = .group),vjust=-3.5,hjust=-.5)+
  geom_jitter(width = 0.05, 
              alpha=0.1)+
  theme_bw()

# ap_abu  

# 1.5 Herbivore richness ----
Conditions <- gdp$guild == "herbivore" & gdp$ind == "rich"

# Make the control as base
gdp$treat <- factor(gdp$treat,
                    levels = c("C","I","H1","H2"))

hrich1 <- lmer(val~treat+(1|block),
                 data = gdp[Conditions, ])

summary(hrich1)

tg <- TukeyGroups(hrich1, gdp[Conditions, ])

yvalue <- aggregate(.~treat, 
                    data=gdp[Conditions, c("treat","val")], 
                    mean)
tdat <- cbind(yvalue, 
              tg[order(tg$treat),c(2,5,6,7)])

gdp$treat <- factor(gdp$treat,
                    levels = c("I","C","H1","H2"))

herb_rich <- ggplot(gdp[Conditions,], 
                 aes(x = treat, y = val))+
  stat_summary(geom="pointrange", 
               fun.data = "mean_cl_boot")+
  geom_text(data =  tdat,
            aes(label = .group),vjust=-3.5,hjust=-.5)+
  geom_jitter(width = 0.05, 
              alpha=0.1)+
  theme_bw()

# herb_rich  


# 1.5 AP richness ----
Conditions <- gdp$guild == "art_pred" & gdp$ind == "rich"

# Make the control as base
gdp$treat <- factor(gdp$treat,
                    levels = c("C","I","H1","H2"))

# Well approximated by normal distribtion
# hrich1 <- nlme::lme(val~treat, random = ~1|block,
#                data = gdp[Conditions, ])

# Negative binomial withouth the block works well.
hrich1 <- glm.nb(val~treat,
                    data = gdp[Conditions, ])

# library(brms)
# brm(val~treat+(1|block),
#     data = gdp[Conditions, ])

summary(hrich1)

# Predict values
# require(nlme)
# simulate.lme(hrich1, 
#          nsim = 1)
# 
# gdp[Conditions, ] -> d
# n = 50
# d[rep(seq_len(nrow(d)),n),] -> drep
# drep$predset <- stack(simulate(hrich1, 
#                                nsim = n))[,1]
# ggplot(drep, aes(y = predset,x=treat))+
#   geom_jitter(alpha = 0.3)+
#   geom_jitter(data=d, aes(y = val,x=treat, col = "red"), alpha = 0.9)

tg <- TukeyGroups(hrich1, gdp[Conditions, ])

yvalue <- aggregate(.~treat, 
                    data=gdp[Conditions, c("treat","val")], 
                    mean)
tdat <- cbind(yvalue, 
              tg[order(tg$treat),c(2,5,6,7)])

gdp$treat <- factor(gdp$treat,
                    levels = c("I","C","H1","H2"))

ap_rich <- ggplot(gdp[Conditions,], 
                    aes(x = treat, y = val))+
  stat_summary(geom="pointrange", 
               fun.data = "mean_cl_boot")+
  geom_text(data =  tdat,
            aes(label = .group),vjust=-3.5,hjust=-.5)+
  geom_jitter(width = 0.05, 
              alpha=0.1)+
  theme_bw()

# ap_rich  

# 1.6 Herbivore diversity Simpson ----
Conditions <- gdp$guild == "herbivore" & gdp$ind == "diva"

# Make the control as base
gdp$treat <- factor(gdp$treat,
                    levels = c("C","I","H1","H2"))

# Model
hSimp1 <- lmer(val~treat+(1|block),
                 data = gdp[Conditions, ])

summary(hSimp1)
# plot(hSimp1)

# Predict values
# gdp[Conditions, ] -> d
# n = 100
# d[rep(seq_len(nrow(d)),n),] -> drep
# drep$predset <- stack(simulate(hSimp1,
#                                nsim = n))[,1]
# ggplot(drep, aes(y = predset,x=treat))+
#   geom_jitter(alpha = 0.3)+
#   geom_jitter(data=d, aes(y = val,
#                           x=treat,
#                           col = "red"), alpha = 0.9)

tg <- TukeyGroups(hSimp1, gdp[Conditions, ])

yvalue <- aggregate(.~treat, 
                    data=gdp[Conditions, c("treat","val")], 
                    mean)
tdat <- cbind(yvalue, 
              tg[order(tg$treat),c(2,5,6,7)])

gdp$treat <- factor(gdp$treat,
                    levels = c("I","C","H1","H2"))

h_simp <- ggplot(gdp[Conditions,], 
                  aes(x = treat, y = val))+
  stat_summary(geom="pointrange", 
               fun.data = "mean_cl_boot")+
  geom_text(data =  tdat,
            aes(label = .group),vjust=-3.5,hjust=-.5)+
  geom_jitter(width = 0.05, 
              alpha=0.1)+
  theme_bw()

# h_simp  

# 1.6 AP diversity Simpson ----
Conditions <- gdp$guild == "art_pred" & gdp$ind == "diva"

# Make the control as base
gdp$treat <- factor(gdp$treat,
                    levels = c("C","I","H1","H2"))

apSimp1 <- lmer(val~treat+(1|block),
                 data = gdp[Conditions, ])

# library(brms)
# brm(val~treat+(1|block),
#     data = gdp[Conditions, ])

summary(apSimp1)

# Predict values
# require(nlme)
# simulate.lme(hrich1, 
#          nsim = 1)
# 
# gdp[Conditions, ] -> d
# n = 50
# d[rep(seq_len(nrow(d)),n),] -> drep
# drep$predset <- stack(simulate(hrich1, 
#                                nsim = n))[,1]
# ggplot(drep, aes(y = predset,x=treat))+
#   geom_jitter(alpha = 0.3)+
#   geom_jitter(data=d, aes(y = val,x=treat, col = "red"), alpha = 0.9)

tg <- TukeyGroups(apSimp1, gdp[Conditions, ])

yvalue <- aggregate(.~treat, 
                    data=gdp[Conditions, c("treat","val")], 
                    mean)
tdat <- cbind(yvalue, 
              tg[order(tg$treat),c(2,5,6,7)])

gdp$treat <- factor(gdp$treat,
                    levels = c("I","C","H1","H2"))

ap_simp <- ggplot(gdp[Conditions,], 
                  aes(x = treat, y = val))+
  stat_summary(geom="pointrange", 
               fun.data = "mean_cl_boot")+
  geom_text(data =  tdat,
            aes(label = .group),vjust=-3.5,hjust=-.5)+
  geom_jitter(width = 0.05, 
              alpha=0.1)+
  theme_bw()

# ap_simp  

# 1.7 Herbivore div SW ----
# pass
# 1.7 AP div SW ----
# pass

# 1.8 Herbivore density ----
Conditions <- gdp$guild == "herbivore" & gdp$ind == "dens"

# Make the control as base
gdp$treat <- factor(gdp$treat,
                    levels = c("C","I","H1","H2"))

# Model
hdens1 <- lmer(val~treat+(1|block),
               data = gdp[Conditions, ])
# hdens2 <- glmer(val~treat+(1|block),
#                 family = gaussian(link = "log"),
#                data = gdp[Conditions, ])

# Linear fit seems better
# anova(hdens1,hdens2)
summary(hdens1)
# summary(hdens2)

# plot(hdens1)
# plot(hdens2)

# Predict values
# gdp[Conditions, ] -> d
# n = 100
# d[rep(seq_len(nrow(d)),n),] -> drep
# drep$predset <- stack(simulate(hdens1,
#                                nsim = n))[,1]
# ggplot(drep, aes(y = predset,x=treat))+
#   geom_jitter(alpha = 0.3)+
#   geom_jitter(data=d, aes(y = val,
#                           x=treat,
#                           col = "red"), alpha = 0.9)

tg <- TukeyGroups(hdens1, gdp[Conditions, ])

yvalue <- aggregate(.~treat, 
                    data=gdp[Conditions, c("treat","val")], 
                    mean)
tdat <- cbind(yvalue, 
              tg[order(tg$treat),c(2,5,6,7)])

gdp$treat <- factor(gdp$treat,
                    levels = c("I","C","H1","H2"))

h_dens <- ggplot(gdp[Conditions,], 
                 aes(x = treat, y = val))+
  stat_summary(geom="pointrange", 
               fun.data = "mean_cl_boot")+
  geom_text(data =  tdat,
            aes(label = .group),vjust=-3.5,hjust=-.5)+
  geom_jitter(width = 0.05, 
              alpha=0.1)+
  theme_bw()

# h_dens 

# 1.8 AP density ----
Conditions <- gdp$guild == "art_pred" & gdp$ind == "dens"

# Make the control as base
gdp$treat <- factor(gdp$treat,
                    levels = c("C","I","H1","H2"))

# Model
hdens1 <- lmer(val~treat+(1|block),
               data = gdp[Conditions, ])

# Linear fit seems better
summary(hdens1)
# summary(hdens2)

# Predict values
# gdp[Conditions, ] -> d
# n = 100
# d[rep(seq_len(nrow(d)),n),] -> drep
# drep$predset <- stack(simulate(hdens1,
#                                nsim = n))[,1]
# ggplot(drep, aes(y = predset,x=treat))+
#   geom_jitter(alpha = 0.3)+
#   geom_jitter(data=d, aes(y = val,
#                           x=treat,
#                           col = "red"), alpha = 0.9)

tg <- TukeyGroups(hdens1, gdp[Conditions, ])

yvalue <- aggregate(.~treat, 
                    data=gdp[Conditions, c("treat","val")], 
                    mean)
tdat <- cbind(yvalue, 
              tg[order(tg$treat),c(2,5,6,7)])

gdp$treat <- factor(gdp$treat,
                    levels = c("I","C","H1","H2"))

ap_dens <- ggplot(gdp[Conditions,], 
                 aes(x = treat, y = val))+
  stat_summary(geom="pointrange", 
               fun.data = "mean_cl_boot")+
  geom_text(data =  tdat,
            aes(label = .group),vjust=-3.5,hjust=-.5)+
  geom_jitter(width = 0.05, 
              alpha=0.1)+
  theme_bw()

# ap_dens

ggpubr::ggarrange(herb_bio, herb_abu, h_simp,h_dens,
                  ap_bio, ap_abu, ap_simp, ap_dens, 
                  labels = c("A","B","C","D",
                             "E","F","G","H"),nrow = 2, ncol = 4)
