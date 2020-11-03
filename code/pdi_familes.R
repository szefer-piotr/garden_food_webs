rm(list=ls())
source("code/data_processing_code.R")
source("code/pdi.R")
source("code/diet_shift.R")

familydf <- data.frame(pdi = diet_breadth_ab,
                       spec = names(diet_breadth_ab),
                       fam = substr(names(diet_breadth_ab),1,4))


library(ggplot2)
library(lme4)
library(emmeans)
library(multcomp)

# Weight for each species

treshold <- 20
filtered_species <- names(pdiss[pdiss>treshold])
plotdf <- familydf[(familydf$spec %in% filtered_species), ]
plotdf$abu <- pdiss[plotdf$spec]

plotdf <- plotdf[!(plotdf$fam %in% c("mant", "aran")), ]



# Models
lm1 <- lm(pdi~fam, data=plotdf)

summary(lm1)

inter.test <- emmeans(lm1, "fam")
pairwise <- cld(inter.test, Letter="abcdefghijklm")

# ggplot(plotdf, aes(x = fam, y = pdi))+
#   geom_boxplot(outlier.colour = "grey50")
plotdf$fam
labeldat <- data.frame(labels = c("a","ab","ab","ab","b"),
                       fam = as.character(unique(plotdf$fam)),
                       ypos = tapply(plotdf$pdi,
                                     as.character(plotdf$fam),
                                     max)+0.03)

ggplot(plotdf, aes(x = fam, y = pdi)) + 
  geom_boxplot(outlier.shape = 1)+
  geom_text(data = labeldat, aes(label = labels, y = ypos),
            position = position_dodge(width = .75),
            show.legend = FALSE)


# Classification model

library(tidyverse)
library(caret)
library(nnet)

# there are only three observations with PDI = 1
plotdf[plotdf$pdi == 1,]

# for the rest of the data:
brmod <- betareg::betareg(pdi~fam, data = plotdf[plotdf$pdi != 1,])
summary(brmod)
inter.test <- emmeans(brmod, "fam")
pairwise <- cld(inter.test, Letter="abcdefghijklm")
plot(pairwise)

# Predict 
y_new <- simulate(lm1) # terrible model, predicts values higher than 1
y_new <- simulate(brmod)
y_new$fam <- substr(rownames(y_new),1,4)

ggplot(y_new, aes(x = fam, y = sim_1)) + 
  geom_boxplot(outlier.shape = 1)


# Pdi for species selected in ordination analysis
# from pRDA for herbiviores
selSp <- c("cole003",
           "cole059",
           "cole061",
           "homo004",
           "lepi011",
           "lepi019",
           "orth006",
           "orth013",
           "orth014",
           "orth016",
           "orth018",
           "orth020")

selIns <- insects[insects$morphotype %in% selSp,]
rownames(treats) <- treats$codes
selIns$treat <- treats[selIns$plot, ]$treat
selIns$pdi <- diet_breadth_ab[selIns$morphotype]

ggplot(selIns, aes(x = treat, y=log(amount+1),
                   color = morphotype))+
  geom_jitter(width=0.05)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_wrap(~morphotype, scales = "free")
  
selIns[selIns$morphotype == "lepi011",]$pdi
selIns[selIns$morphotype == "orth020",]$pdi
