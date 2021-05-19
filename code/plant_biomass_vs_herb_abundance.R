rm(list=ls())
source("code/data_processing_code.R")

library(dplyr)
library(ggplot2)

csites <- as.character(treats[treats$treat %in% c("CONTROL"),]$codes)
psites <- as.character(treats[treats$treat %in% c("PREDATOR"),]$codes)

main_sp_bio <- main %>%
  filter(CODE %in% c(psites, csites) & LIFE.FORM %in% c("shrub", "tree")) %>%
  group_by(CODE, SP_CODE) %>%
  summarise(tot_bio = sum(WEIGHT),
            leaf_weight = sum(LEAVES)) %>%
  mutate(SP_CODE = tolower(SP_CODE)) %>%
  mutate(treat = ifelse(CODE %in% csites, "control", "exclosure"))
  
ins_sp_abubio <- ins_bio %>%
  filter(plot %in% c(psites, csites)) %>%
  group_by(plot, tree, family) %>%
  summarise(ins_totbio = sum(totbio),
            ins_totabu = sum(amount))

fulldf <- left_join(ins_sp_abubio, main_sp_bio, 
          by = c("plot" = "CODE",
                 "tree" = "SP_CODE"))

# ggplot(fulldf, aes(x = tot_bio, y = log(ins_totbio)))+
#   geom_point()+
#   geom_smooth()


# ggplot(fulldf, aes(x = tot_bio*1000, y = ins_totabu))+
#   geom_point(aes(colour = tree))+
#   geom_smooth(method = "lm", aes(colour=treat))+
#   scale_x_continuous(trans = 'log10') + 
#   scale_y_continuous(trans = 'log10')+
#   facet_wrap(~family, scales = "free")

ord.labs <- c("Orthoptera",
              "Homoptera",
              "Heteroptera",
              "Aranea",
              "Mantodea", 
              "Coleoptera",
              "Lepidoptera")

levels(fulldf$family) <- sort(ord.labs)
names(fulldf)[8] <- "Treatment"
fulldf$Treatment <- as.factor(fulldf$Treatment)
levels(fulldf$Treatment) <- c("Control", "Exclosure")

ibrary(ggpmisc)

my.formula <- y~x

ggplot(fulldf, aes(x = leaf_weight*1000, y = ins_totabu))+
  geom_point()+
  geom_smooth(method = "lm",aes(colour=Treatment))+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +
  scale_x_continuous(trans = 'log10') + 
  scale_y_continuous(trans = 'log10')+
  facet_wrap(~family, scales = "free")+
  ylab("Abundance")+xlab("Leaf weight [g]")+
  theme_bw()+
  theme(legend.position="bottom")

fam_df <- data.frame()
for(fam in unique(fulldf$family)){
  ss <- summary(lm(log10(ins_totabu)~ log10(leaf_weight*1000),
                   data = fulldf[fulldf$family == fam, ]))
  row_df <- data.frame(family = fam,
                       rsq = ss$r.squared,
                       pval = round(ss$coefficients[2,4],10),
                       cat = "leaf")
  fam_df <- rbind(fam_df, row_df)
  print(ss$coefficients)
  
}
for(fam in unique(fulldf$family)){
  ss <- summary(lm(log10(ins_totabu)~ log10(tot_bio*1000),
                   data = fulldf[fulldf$family == fam, ]))
  row_df <- data.frame(family = fam,
                       rsq = ss$r.squared,
                       pval = round(ss$coefficients[2,4],10),
                       cat = "total")
  fam_df <- rbind(fam_df, row_df)
  print(ss$coefficients)
  
}


# Total plant biomass

main_bio <- main %>%
  filter(CODE %in% c(psites, csites) & LIFE.FORM %in% c("shrub", "tree")) %>%
  group_by(CODE) %>%
  summarise(tot_bio = sum(WEIGHT),
            leaf_weight = sum(LEAVES)) %>%
  mutate(treat = ifelse(CODE %in% csites, "control", "exclosure"))

ins_abubio <- ins_bio %>%
  filter(plot %in% c(psites, csites)) %>%
  group_by(plot) %>%
  summarise(ins_totbio = sum(totbio, na.rm = T),
            ins_totabu = sum(amount,na.rm = T))%>%
  mutate(treat = ifelse(plot %in% csites, "control", "exclosure"))

cumdf <- left_join(ins_abubio, main_bio, 
                    by = c("plot" = "CODE"))

ggplot(cumdf, aes(x = leaf_weight*1000, y = ins_totabu))+
  geom_point()+
  geom_smooth(method = "lm", aes(colour = treat.x))+
  scale_x_continuous(trans = 'log10') + 
  scale_y_continuous(trans = 'log10')

lmcum <- lm(log10(ins_totabu)~log10(leaf_weight*1000)*treat.x, data = cumdf)
summary(lmcum)
lmcum1 <- lm(log10(ins_totabu)~log10(tot_bio*1000), data = cumdf)
summary(lmcum1)


ggplot(cumdf, aes(x = tot_bio*1000, y = ins_totabu))+
  geom_point()+
  geom_smooth(method = "lm", aes(colour = treat.x))+
  scale_x_continuous(trans = 'log10') + 
  scale_y_continuous(trans = 'log10')

lmcum <- lm(log10(ins_totabu)~log10(tot_bio*1000)*treat.x, 
            data = cumdf)
summary(lmcum)

