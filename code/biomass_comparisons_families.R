# Community comparisons
rm(list=ls())
source("code/data_processing_code.R")
# source("code/pdi.R")
# source("code/diet_shift.R")

library(ggplot2)
library(dplyr)

head(ins_bio)
csites <- treats[treats$treat %in% c("CONTROL"), ]$codes
psites <- treats[treats$treat %in% c("PREDATOR"), ]$codes
csites <- as.character(csites)
psites <- as.character(psites)

cins <-ins_bio[ins_bio$plot %in% csites, ]
pins <-ins_bio[ins_bio$plot %in% psites, ]

# Biomass 

fam_bio <- ins_bio %>%
  group_by(family,plot) %>%
  filter(plot %in% c(psites,csites)) %>%
  summarise(bio = sum(totbio, na.rm=T)) %>%
  mutate(treatment = if_else(plot %in% psites, "Exclosure", "Control")) %>%
  mutate(block = substr(plot, 3,4))

bio_mod <- nlme::lme(log(bio)~family, random = ~1|block, data = fam_bio)
tukey <- emmeans::emmeans(bio_mod, "family")
grouping <- as.data.frame(multcomp::cld(tukey, Letter="abcdefghijklm"))

ord.labs <- c("Orthoptera",
                    "Aranea",
                    "Homoptera",
                    "Heteroptera",
                    "Mantodea", 
                    "Coleoptera","Lepidoptera")

levels(fam_bio$family) <- levels(grouping$family) <- sort(ord.labs)
bp <- ggplot(fam_bio)+
  geom_boxplot(aes(family, log(bio)))+
  geom_text(data = grouping,
            aes(family, sort(tapply(log(fam_bio$bio),fam_bio$family, max)+0.5), 
                label = .group))+
  xlab("")+ylab("Log(biomass)") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Abundance

fam_abu <- ins_bio %>%
  group_by(family,plot) %>%
  filter(plot %in% c(psites,csites)) %>%
  summarise(abu = sum(amount, na.rm=T)) %>%
  mutate(treatment = if_else(plot %in% psites, "Exclosure", "Control")) %>%
  mutate(block = substr(plot, 3,4))

library(lme)
library(lme4)

abu_mod <- MASS::glm.nb(abu~family, data = fam_abu)
abu_mod_rand <- lme4::glmer.nb(abu~family+(1|block), data = fam_abu)
tukey <- emmeans::emmeans(abu_mod_rand, "family")
grouping <- as.data.frame(multcomp::cld(tukey, Letter="abcdefghijklm"))

levels(fam_abu$family) <- levels(grouping$family) <- sort(ord.labs)

ap <- ggplot(fam_abu)+
  geom_boxplot(aes(family, log(abu)))+
  geom_text(data = grouping,
            aes(family, sort(tapply(log(fam_abu$abu),fam_abu$family, max)+1), 
                label = .group))+
  xlab("")+ylab("Log(abundance)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

fam_abu %>% group_by(family) %>% summarise(ABU = sum(abu))

# Species number
fam_sp <- ins_bio %>%
  group_by(family,plot) %>%
  filter(plot %in% c(psites,csites)) %>%
  summarise(spno = length(unique(morph))) %>%
  mutate(treatment = if_else(plot %in% psites, "Exclosure", "Control")) %>%
  mutate(block = substr(plot, 3,4))

sp_mod <- glm(spno~family, family = "quasipoisson", data = fam_sp)
sp_mod_rand <- glmer.nb(spno~family+(1|block), 
                        data = fam_sp)

summary(sp_mod_rand)
# AER::dispersiontest(sp_mod)
tukey <- emmeans::emmeans(sp_mod_rand, "family")
grouping <- as.data.frame(multcomp::cld(tukey, Letter="abcdefghijklm"))

levels(fam_sp$family) <- levels(grouping$family) <- sort(ord.labs)

spa <- ggplot(fam_sp)+
  geom_boxplot(aes(family, spno))+
  geom_text(data = grouping,
            aes(family, sort(tapply(fam_sp$spno,fam_sp$family, max)+3), 
                label = .group))+
  xlab("")+ylab("Number of species") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggpubr::ggarrange(bp,ap,spa, labels = c("A","B","C"), ncol = 3)
