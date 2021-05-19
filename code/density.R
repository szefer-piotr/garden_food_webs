# Density

# Density of herbivores in treatment and control plots

rm(list=ls())
source("code/data_processing_code.R")

library(dplyr)

# Within each garden check which plant species are comparable
# Abudance of herbivores and AP in each plot.

cp <- treats[treats$treat %in% "CONTROL", ]$codes
pp <- treats[treats$treat %in% "PREDATOR", ]$codes

ins_bio_cp <- ins_bio %>%
  filter(plot %in% c(as.character(cp),
                     as.character(pp)))

ddfh <- ins_bio_cp %>%
  group_by(plot) %>%
  filter(!grepl("aran|mant",family)) %>% 
  summarise(herbabu = sum(amount))

ddfap <- ins_bio_cp %>%
  group_by(plot) %>%
  filter(grepl("aran|mant",family)) %>% 
  summarise(apabu = sum(amount))

ddf <- cbind(ddfh, ddfap)
ddf <- ddf[,-3]

# SLA is in cm2/g
mb <- main[main$CODE %in% c(as.character(cp),
                            as.character(pp)), ]
mbw <- mb %>% 
  filter(LIFE.FORM %in% c("tree", "shrub")) %>%
  mutate(area = SLA * LEAVES*1000) # in cm2

# ggplot(mbw, aes(y = log(area), x = SLA))+
#   geom_point()+
#   geom_smooth(method = "lm", se = T)

mbws <- mbw %>%
  group_by(CODE) %>%
  summarise(area.cm2 = sum(area, na.rm=T)) %>%
  mutate(area.m2 = area.cm2/10000)

# Area comparison
mbsa <- mbws %>%
  mutate(treat = ifelse(CODE %in% cp, "C", "Ex"),
         block = substr(CODE, 3,4))
library(ggplot2)
ggplot(mbsa, aes(x = treat, y=area.m2, group = block))+
  geom_jitter(width = 0.05, size = 5, color = "grey30")+
  geom_line(lty=2,lwd = 1, alpha=0.3)+
  xlab("")+ylab(expression(paste('Total leaf area [', m^2, ']')))

qqnorm(log(mbsa$area.m2))
qqline(log(mbsa$area.m2))

area.lmer <- nlme::lme(area.m2 ~ treat, random=~1|block,
                       data = mbsa)
plot(area.lmer)
area.glmer <- lme4::glmer(area.m2 ~ treat+(1|block),
                   family = gaussian(link = "log"),
                       data = mbsa)


summary(area.lmer)

# Density
ddfp <- cbind(ddf, mbws)

ddfpd <- ddfp %>%
  mutate(hdens = herbabu/area.m2,
         apdens = apabu/area.m2) %>%
  mutate(trt = ifelse(plot %in% as.character(cp), "control", "exclosure"),
         block = substr(plot, 3,4))

# write.table(ddfpd, "datasets/densities.txt")

library(ggplot2)

# Density plots

expr <- expression(paste("Density"," [individuals/", m^{2}, " of foliage]"))

hp <- ggplot(ddfpd, aes(x = trt, 
                  y = hdens))+
  geom_jitter(width = 0.01, size = 5, alpha=0.3)+
  geom_line(aes(group = block), lty = 2, lwd = 1, alpha=0.3)+
  stat_summary(geom="pointrange", col = "black",  lwd = 1)+
  theme_bw() + 
  xlab("") + ylab(expr) + 
  scale_x_discrete(labels = c("C","Ex"))+
  scale_y_log10(breaks = seq(20,150, 15))


hap <- ggplot(ddfpd, aes(x = trt, 
                  y = apdens))+
  geom_jitter(width = 0.01, size = 5, alpha=0.3)+
  geom_line(aes(group = block), lty = 2, lwd = 1, alpha=0.3)+
  stat_summary(geom="pointrange", col = "red",  lwd = 1)+
  theme_bw() + 
  xlab("") + ylab("") +
  scale_x_discrete(labels = c("C","Ex"))+
  scale_y_log10(breaks = seq(0,30,2))
  
# hp
# hap
  
ggpubr::ggarrange(hp, hap, ncol=2, nrow = 1, 
                  labels = c("H","AP"))

# Tests
library(lmerTest)
hlmer <- nlme::lme(log(hdens) ~ trt , 
                   random = ~1|block, data = ddfpd)

hglmer <- glmer(hdens ~ trt +(1|block), 
                data = ddfpd,
                family = gaussian("log"))

apglmer <- glmer(apdens ~ trt +(1|block), 
                data = ddfpd,
                family = gaussian("log"))
summary(hglmer)
summary(apglmer)

# Simulate from a model

plotdata <- ddfpd
simn = 50
preddat <- expand.grid(treat = rep(c("control", "exclosure"),simn),
                         block = c("g1","g2","g3","g4","g5","g6"))
preddat$predicted <- stack(simulate(hglmer, nsim = simn))[,1]
  
ggplot(plotdata, aes(y = hdens, x = trt))+
  geom_jitter(data = preddat, aes(y = predicted, x = treat), 
                size = 2, colour = "blue")+
  geom_jitter(size = 2,colour = "orange")


summary(hlmer)
summary(hglmer) # significant
summary(apglmer) # marginally

# Variances are equal
car::leveneTest(log(hdens)~trt, data = ddfpd)
car::leveneTest(log(apdens)~trt, data = ddfpd)

aplmer <- nlme::lme(log(apdens) ~ trt , random = ~1|block, data = ddfpd)
summary(aplmer) # not significant

bartlett.test(log(apdens)~trt, data = ddfpd)

# Test for the log ratios.
dflrrp <- ddfpd %>% 
  filter(trt == "exclosure")
dflrrc <- ddfpd %>% 
  filter(trt == "control")

lrrherb <- log(dflrrc$hdens/dflrrp$hdens)
lrrap <- log(dflrrc$apdens/dflrrp$apabu)

summary(lm(lrrherb~1))
summary(lm(lrrap~1))

# Density correlation
ggplot(ddfpd, aes(x = hdens, y = apdens))+
         geom_point(aes(color = trt), size = 5)+
  scale_x_log10()+
  scale_y_log10()+
  ylab("Arthropod predator density")+xlab("Herbivore density")+
  scale_colour_discrete(name = "Treatment", 
                        labels = c("Control", "Exclosure"))+
  theme(legend.position = "bottom")

summary(lm(apdens~hdens+trt, data=ddfpd))
       