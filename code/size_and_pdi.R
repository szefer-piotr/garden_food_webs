# size and specialization... can they predict the predator effect/abundance change?

source("code/data_processing_code.R")
source("code/pdi.R")

diet_breadth
ins_bio
treats
biollcp # i should have used this instead
abufulldf # or this for abundnace
plbio <- tapply(biollcp$plbio, biollcp$plnm, max)

# General model, offset for biomass
psites <- as.character(treats[treats$treat %in% c("PREDATOR"), ]$codes)
csites <- as.character(treats[treats$treat %in% c("CONTROL"), ]$codes)

inssub <- ins_bio[ins_bio$plot %in% c(psites, csites), ]
# bcpsub <- biollcp[biollcp$plot %in% c(psites, csites), ]
# inssub <- bcpsub
inssub$plottreeherb <- paste(inssub$plot, inssub$tree, inssub$family, sep = "") 
biollcp$plottreeherb <- paste(biollcp$plot, biollcp$plnm, biollcp$nms, sep = "") 
rownames(biollcp) <- biollcp$plottreeherb

inssub$plbio <- biollcp[inssub$plottreeherb,]$plbio

# db values
inssub$db <- diet_breadth[as.character(inssub$morphotype)]
cnames <- inssub[inssub$plot %in% csites, ]$morphotype
pnames <- inssub[inssub$plot %in% psites, ]$morphotype
compnames <- as.character(unique(cnames[cnames %in% pnames]))
inssub$treat <- "control"
inssub[inssub$plot %in% psites, ]$treat <- "predator"
inssub$block <- substr(inssub$plot, 3,4)


# Treshold how many individuals has to be in each group
treshold <- 5
finalno <- 0
specdat <- c()
nm <- compnames[22]

for (nm in compnames){
  # print(nm)
  nmsub <- inssub[inssub$morphotype == nm, ]
  
  nmsN <- group_by(nmsub, block, tree) %>%
    tally()
  nmsN %>% tbl_df %>% print(n=40)
  nmsN[nmsN$n == 2, ]
  
  subtb <- table(nmsub$treat)
  if(subtb[1] >= treshold & subtb[2] >= treshold){
    msg <- paste(nm, " ", "is a match", sep = "")
    print(msg)
    finalno <- finalno +1
    specdat <- c(specdat, nm)
  }
}

specdat
spsub <- inssub[inssub$morphotype == specdat[1], ] 

specdat_noip <- specdat[-grep("aran|mant", specdat)]

biosub <- inssub[inssub$morphotype %in% specdat_noip, ]
biosub$block <- substr(biosub$plot, 3,4)

# Fill in later
biosub$ipbio <- 0
biosub$ipabu <- 0

# remove unused columns
bsc <- biosub[,c("morphotype", "amount", "bio", "totbio",
                 "ipbio", "ipabu",
                 "tree", "family", "morph", "plbio", "db",
                 "treat", "block", "plot")]

# get some indication of food quality for plants and instead of using species
# use that indincation
source("code/plant_species_quality.R")
plant_quality

# Plant quality measured as percentage water content
# Percentage water content
nms <- tolower(paste(main$SP_CODE, main$CODE, sep = ""))
wat <- main$WATER/main$WET..g.
speccode <- paste(main$CODE, tolower(main$SP_CODE), sep = "")
waterdf <- data.frame(water = wat,name = tolower(nms))
waterdf <- waterdf[!duplicated(waterdf$name), ]
rownames(waterdf) <- waterdf$name

# There are some missing values for melanolepis multiglandulosa at one site
# Fill these with average for the species
nm  <- "melamu"

for (nm in unique(substr(waterdf[is.na(waterdf$water), ]$name, 1,6))){
  print(nm)
  subwat <- waterdf[grep(nm, waterdf$name),]
  meandat <- subwat[!is.na(subwat$water), ]$water
  meandatnames <- as.character(subwat[is.na(subwat$water), ]$name)
  print(meandatnames)
  if((length(meandat[!is.na(meandat)]))!= 0){
    aver <- mean(meandat)
    waterdf[meandatnames, ]$water <- aver
  }
}

bsc$treesite <- paste(bsc$tree, bsc$plot, sep = "")
bsc$water <- waterdf[bsc$treesite, 1]

bsc$quality <- "high"
bsc[bsc$water <= median(bsc$water, na.rm=T), ]$quality <- "low"

# bsc$quality <- plant_quality[bsc$tree, 1]
# bsc$palat <- plant_quality[bsc$tree, 2]

library(ggplot2)
# Evaluate the response
# qqnorm(log(bsc$totbio), col = alpha("black", 0.1))
# qqline(log(bsc$totbio))
# hist(bsc$totbio, breaks = 100)
summary(log(bsc$totbio))
# it is not hugely different from the normal
# qqnorm(log(bsc[bsc$treat == "predator", ]$totbio), col = alpha("black", 0.1))
# qqline(log(bsc[bsc$treat == "predator", ]$totbio))
# 
# qqnorm(log(bsc[bsc$treat == "control", ]$totbio), col = alpha("black", 0.1))
# qqline(log(bsc[bsc$treat == "control", ]$totbio))
# 
# hist(log(bsc[bsc$treat == "predator", ]$totbio), breaks = 20, col="red")
# hist(log(bsc[bsc$treat == "control", ]$totbio), breaks = 20)
# 
# ggplot(bsc, aes(x = log(totbio), fill = treat)) + geom_density(alpha = 0.2)

library(lme4)
library(lmerTest)

# Hypotheses for models:
# TTI dietary specialist herbivore should be less affected by enemies, doesnt matter whether on poor or high quality resources. Generalists on the other hand should be strongly negatively affecte by predators on low quality food. pd*quality.

# Howoever, if generalist can move (mobility), then these effects should not be as strong, and maybe there will be no differences in mobile groups.

# Schmitz evaluated the effect of size (page 128 in "Resolving... "): small grashoppers have a lower capacity tto digest and assimilate much of the vegetation in the fields oowing to its generally poor quality. Smaller grashoppers must sttpend considerable effort seeking high quality plants resources that tend to be rare. This can highten starvation mortality of individuals in smaller size classes relative to larger individuals owing to greater  predation mortality due to more risk prone foraging behaviour by small individuals... however, there was little or no difference in the effect of grasshopper size class on the final abundance of grasses and herbs in the plant community. Grasshoppers in smaller size classes exhibited higher growth rates. These can be only sustained by higher foraging effert. Lower density is compensated by the effect of greater per capita foraging effort of the surviving individuals.

# Does treatment have an overal effect on herbivore biomass
# LMER1 <- glmer(amount ~ treat+(1|block), data=bsc, family = "poisson")
# LMER1nb <- glmer.nb(amount ~ treat+(0+treat|tree), data=bsc)

# summary(LMER1nb)
library(blmeco)
# dispersion_glmer(LMER1) # should not exceed 1.4

library(emmeans)
# plot(emmeans(LMER1nb,"treat",adjust = "tukey"))
# emmeans(LMER1nb, list(pairwise~treat))

# library(lme4)
# effect.out.lme <- effect(term = "treat", 
#                          mod = LMER1nb)


# boxplot(log(bsc[complete.cases(bsc),]$totbio)~bsc[complete.cases(bsc),]$treat)
# LMER1 <- lmer(log(totbio) ~ treat+offset(plbio)+(1|block), data=bsc)

# LMER1b <- lmer(log(totbio) ~ treat*quality+(1+treat|tree) + offset(log(bio)), data=bsc)
# summary(LMER1b)

# LMER1b <- glmer(amount ~ treat*quality+(1+treat|db), data=bsc, family = "poisson")
# summary(LMER1b)

# Why to uses offset ? - log biomass of species is actually larger on more abundant plants, but this account for only 11% of variation
# plot(log(totbio)~log(plbio), data=bsc, col = alpha("black", 0.1), pch=19)
# abline(lm(log(totbio)~log(plbio), data=bsc)
# )
# summary(lm(log(totbio)~log(plbio), data=bsc))
# I dont think I need to use it at all. 

# unique(bsc$db)

# visualize
# qvals <- unique(bsc$quality)
# 
# ut <- unique(bsc$tree)
# sort(plant_quality[as.character(ut), ])
# form <- plbio~plant_quality[names(plbio),1]
# plot(form)

poor <- "breyce"
high <- "pipeum"
average <- "piptar"

# subbsc <- bsc[bsc$tree == poor, ]
# subbsc <- bsc[bsc$tree == high, ]
# subbsc <- bsc[bsc$tree == average, ]
# 
# # 
# ggplot(subbsc, aes(y = log(totbio), x = treat, 
#                 col = db, 
#                 group = db)) + 
#   stat_summary(fun.data=mean_cl_boot, 
#                geom="pointrange", lwd=0.8) +
#   stat_summary(fun=mean, geom="point",cex = 2) +
#   stat_summary(fun=mean, geom="line",lwd=1, lty=2)
#   # stat_summary(fun.y=mean, geom="text", 
#   #              col = rgb(10,10,10,180,maxColorValue = 255),
#   #              hjust = 1.2,
#   #              vjust = -1.5)
# plot(log(bio)~db, data = bsc, xlim = c(0.8, 1))
# 
# plot(bsc$quality~bsc$db)

# Comparison between plant species

# Rename variables
# names(bsc) <- c("species","habu","hsize","hbio",
#                 "ipbio","ipabu","treesp", "hfam",
#                 "spec", "plbio", "db", "treat",
#                 "block","quality","palat")

# ggplot(subbsc, aes(y = log(amount+1), x = treat, 
#                    col = db, 
#                    group = db)) + 
#   stat_summary(fun.data=mean_cl_boot, 
#                geom="pointrange", lwd=0.8) +
#   stat_summary(fun=mean, geom="point",cex = 2) +
#   stat_summary(fun=mean, geom="line",lwd=1, lty=2)

# Full interaction model

# Check how much data I have for each interaction

library(dplyr)
library(tidyr)

# Example
# loc<-c("city1","city2","city1","city2","city1","city1","city2","city2","city1","city2")
# q1<-c("YES","YES","NO","MAYBE","NO","NO","YES","NO","MAYBE","MAYBE")
# q2<-c("YES","NO","MAYBE","YES","NO","MAYBE","MAYBE","YES","YES","NO")
# q3<-c("NO","NO","NO","NO","YES","YES","MAYBE","MAYBE","NO","MAYBE")
# df<-data.frame(loc,q1,q2,q3)

# dfN <- gather(df, quest, answ, q1:q3) %>%
#   complete(loc, quest, answ) %>%
#   unique()
# bsc$db

#db, quality
mean(bsc$db)
bsc$specialization <- "generalist"
bsc[bsc$db <= mean(bsc$db), ]$specialization <- "specialists"

bsc$pqual <- "poor"
bsc[bsc$quality >= mean(bsc$quality, na.rm = T), ]$pqual <- "high"

bscN <- group_by(bsc, specialization, quality, treat) %>%
  tally()
bscN %>% tbl_df %>% print(n=40)

# Are bigger species more specialized?
sizepdi <- data.frame()
for(mft in unique(bsc$morphotype)){
  ss <- bsc[bsc$morphotype == mft, ]
  dfrow <- data.frame(morph = ss$morphotype[1], 
                      fam = ss$family[1], 
                      db = ss$db[1], 
                      bio = ss$bio[1])
  sizepdi <- rbind(sizepdi, dfrow)
}
# bigger species are not more specialized

# Are different families differ in size?
plot(ins_bio$family, log(ins_bio$bio))

#

summary(lm(db~log(bio), data=sizepdi))
plot(db~log(bio), data=sizepdi)

# Bayesian Interaction Model
library(brms)
LMER1 <- lmer(log(totbio) ~ treat*specialization*quality+offset(bio)+(1|morphotype), data=bsc)
summary(LMER1)

LMER2 <- glmer.nb(amount ~ treat*specialization*quality+(1|morphotype), data=bsc)
summary(LMER2)

LMER2b <- glmer.nb(amount ~ treat*specialization+(1|morphotype), data=bsc)
summary(LMER2b)

LMER2c <- glmer.nb(amount ~ treat*quality+(1|morphotype), data=bsc)
summary(LMER2c)


bm1 <-    brm(log(totbio) ~ treat*specialization*quality+offset(bio)+(1|morphotype), data=bsc)
plot(bm1)
summary(bm1)
# Total biomass of a given species in responce
# Change names such that I can read the results

# LMER1 <- lmer(log(totbio) ~ treat*specialization*pqual+offset(totbio)+(1|morphotype), data=bsc)
# summary(LMER1)
# LMER2 <- glmer.nb(amount ~ treat*specialization+(1|morphotype), 
#                data=bsc)
# summary(LMER2)
# LMER3 <- glmer.nb(amount ~treat + (1+treat|pqual/specialization), 
#                   data= bsc)
# summary(LMER3)
ggplot(bsc, aes(x = treat, y = log(amount),
                group = morphotype))+
  geom_jitter() +
  facet_grid(bsc$specialization~bsc$pqual)+
  geom_line()

# NEW IDEA [9.06.2020] Evalute using bayesian methods the beta parameters for the treatment response and then correlate this to the specialization and quality - impossible

# NEW IDEA [1.06.2020] evaluate the strength of the predator effeect using partial RDA and correlate these measurments with quality and specialisation of herbivores
