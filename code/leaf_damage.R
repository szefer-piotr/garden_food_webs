# Leaf damage analysis
rm(list=ls())
source("code/data_processing_code.R")

library(dplyr)
library(lme4)
library(lmerTest)
library(ggplot2)

ddat <- main %>%
  filter(LIFE.FORM %in% c("shrub", "tree")) %>%
  filter(TREAT %in% toupper(c("control", "predator"))) %>%
  select(toupper(c("code", "block", "treat","leaves","weight",
                   "sp_code","no_stems", "sla","herb"))) %>%
  mutate(TREAT = ifelse(TREAT %in% "CONTROL", "Co", "Ex"),
         sp_block = paste(SP_CODE, BLOCK),
         area = SLA * LEAVES*1000)

colnames(ddat) <- c("code", "block", "treat","leaves","weight",
           "sp_code","no_stems", "sla","herb", "sp_block", "area")
# Analysis
dlmer1 <- lmer(herb ~ treat+(1|sp_code), data = ddat)
dlmer2 <- lmer(herb ~ treat+sla+(1|sp_code), data = ddat)
dlmer3 <- lmer(herb ~ treat+sla+(1|sp_code)+(1|block), data = ddat)
dlmer4 <- lmer(herb ~ treat + (1|sp_code)+(1|block), data = ddat)
dlmer5 <- lmer(herb ~ treat*sla+(1|sp_code)+(1|block), data = ddat)
dlmer6 <- lmer(herb ~ sla+(1|sp_code)+(1|block), data = ddat)
dlmer7 <- lmer(herb ~ sla+(1|sp_code), data = ddat)

anova(dlmer1, dlmer2) #2 is better
anova(dlmer2, dlmer3) #2 still most parsimonious
anova(dlmer2, dlmer4) #2 still the best
anova(dlmer2, dlmer5) # 2nd still
anova(dlmer2, dlmer6)
anova(dlmer2, dlmer7)
anova(dlmer6, dlmer7)

summary(dlmer2)
summary(dlmer5)
summary(dlmer7)

# Weighted analysis - not significant
wlmer1 <- lmer(herb ~ treat + (1|sp_code)+(1|block), 
               data = ddat, weights = weight)

summary(wlmer1)
wlmer2 <- lmer(sla ~ treat + (1|sp_code)+(1|block), 
               data = ddat, weights = weight)

summary(wlmer2)

# Community Weighted Mean differences

# Sla predicts herbivore damage
sla1 <-  lmer(sla ~ treat+(1|sp_code)+(1|block), data = ddat)
summary(sla1)

area1 <-  lmer(area ~ treat+(1|block), 
                data = ddat)
summary(area1)

ggplot(ddat, aes(x = treat, y = log(area)))+
  geom_jitter()

# expr <- expression(paste("Density"," [individuals/", m^{2}, " of foliage]"))

hd <- ggplot(ddat, aes(x = treat, 
                        y = herb))+
  geom_jitter(width = 0.02, size = 5, alpha=0.1)+
  geom_line(aes(group = sp_block), lty = 2, lwd = 1, alpha=0.1)+
  stat_summary(geom="pointrange", fun.data = mean_cl_boot,
               col = "black",  lwd = 1)+
  theme_bw() + 
  xlab("") + ylab("logit(Leaf damage)")
  # scale_x_discrete(labels = c("C","Ex"))+
  # scale_y_log10(breaks = seq(20,150, 15))
# hd

hsla <- ggplot(ddat, aes(x = treat, 
                       y = sla))+
  geom_jitter(width = 0.02, size = 5, alpha=0.1)+
  geom_line(aes(group = sp_block), lty = 2, lwd = 1, alpha=0.1)+
  stat_summary(geom="pointrange", fun.data = mean_cl_boot,
               col = c("black"),  lwd = 1)+
  theme_bw() + 
  xlab("") + ylab("Specific Leaf Area")
# scale_x_discrete(labels = c("C","Ex"))+
# scale_y_log10(breaks = seq(20,150, 15))
# hsla

damsla <- ggplot(ddat, aes(y = sla, 
                           x = herb, col = treat))+
  geom_point(size = 5, alpha=0.1)+
  stat_smooth(method = "lm", se  = T)+
  theme_bw() + 
  ylab("Specific Leaf Area") + xlab("logit(leaf damage)")

ggpubr::ggarrange(hd, hsla, damsla, ncol=3, nrow = 1, 
                  labels = c("A","B","C"))

# vals <- seq(from =0.01, to=0.99, by = 0.05)
# lvals <- log(vals/(1 - vals))

# Do plants that compensate have higher relative abundance?
data_full <- data.frame()
for(gd in unique(ddat$block)){
  print(gd)
  subsC <- ddat %>% 
    filter(block == gd & treat == "Co")
  subsEx <- ddat %>% 
    filter(block == gd & treat == "Ex")
  print(subsC)
  print(subsEx)
  plant_nm <- subsC$sp_code[subsC$sp_code %in% subsEx$sp_code]
  subs_fin <- cbind(subsC[subsC$sp_code %in% plant_nm, ],
                    subsEx[subsEx$sp_code %in% plant_nm, ])
  data_full <- rbind(data_full, subs_fin)
}
colnames(data_full) <- c("codeC", "blockC", "treatC", "leavesC", "weightC", "sp_codeC", "no_stemsC", "slaC", "herbC", "sp_blockC", "areaC","code","block","treat","leaves","weight","sp_code","no_stems","sla","herb","sp_block","area")

relDominance <- function(bl, plant, tr){
  reldat <- ddat %>%
    filter(block == bl & treat == tr)
  reldom <- reldat[reldat$sp_code == plant, ]$weight/sum(reldat$weight)
  return(reldom)
}

data_full$RD <- NA

# Log ratio vs relative dominance in the control
for(row in 1:nrow(data_full)){
  print(row)
  data_full[row, ]$RD <- relDominance(data_full[row, ]$block,
                                      data_full[row, ]$sp_code,
                                      "Co")
}

data_full$RD <- data_full$RD
data_full <- data_full %>%
  mutate(deltaSLA = log(slaC/sla))

ggplot(data_full, aes(y = RD,x = deltaSLA))+
  geom_point(size = 5, alpha=0.5)

summary(betareg::betareg(RD~deltaSLA, data = data_full))
