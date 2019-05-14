insects <- read.table("datasets/arthropods_clean.txt")
treats  <- read.table("datasets/treatments_clean.txt")
plants  <- read.table("datasets/plants_clean.txt")
size_dat <-read.table("datasets/size_dat_bio.txt")

library("bipartite")
library(lme4)
library(lmerTest)
source("code/bio_log_ratio.R")
source("code/contingencyTable.R")
source("code/bio_log_ratio.R")
source("code/weighted-modularity-LPAwbPLUS/code/R/LPA_wb_plus.R")
source("code/weighted-modularity-LPAwbPLUS/code/R/MODULARPLOT.R") #read in plotting function
# gardnets - networks in individual plots
# biofulldf biomass of all insects and trees

# Is it possible to study average number of plants per herbivore species using raw data?
# That could be done for arthropod species present in both predator and control plots. That would be problematic if there are differences in plant communities between blocks

# Simple vulnerability and generality patterns in networks
gnames <- as.character(treats$codes)
genvuldf <- data.frame()
for(gname in gnames){
  print(gname)
  MAT <- gardnets[[gname]]
  
  tt <- tryCatch(networklevel(MAT, index = "vulnerability"),error=function(e) e, warning=function(w) w)
  ifelse(is(tt,"warning"),next,"OK")
  
  nlres <- networklevel(MAT, index = "vulnerability")
  pdi <- mean(PDI(MAT))
  con <-  networklevel(MAT, index = "connectance")
  
  mod1 = DIRT_LPA_wb_plus(MAT)
  MODinformation = GetModularInformation(MAT,mod1)
  mod = MODinformation$realized_modularity
  ispp = ncol(MAT[, substr(colnames(MAT),1,4) %in% c("aran", "mant")])
  isph = ncol(MAT) - ispp
  
  subgvdf <- data.frame(plot = gname, gen = nlres[1], vul = nlres[2], 
                        pdi = pdi, con = con, mod = mod, isph=isph,
                        ispp = ispp)
  genvuldf <- rbind(genvuldf, subgvdf)
}

rownames(treats) <- treats$codes
genvuldf$trt <- treats[genvuldf$plot, ]$treat
genvuldf$block <- substr(genvuldf$plot, 3,4)

# But generality would also change if the number of plant at the plot is smaller!
library(ggplot2)

gen <- genvuldf[,c("plot","gen","trt","block")]
vul <- genvuldf[,c("plot","vul","trt","block")]
pdi <- genvuldf[,c("plot","pdi","trt","block")]
con <- genvuldf[,c("plot","con","trt","block")]
mod <- genvuldf[,c("plot","mod","trt","block")]
isph <- genvuldf[,c("plot","isph","trt","block")]
ispp <- genvuldf[,c("plot","ispp","trt","block")]

# Change names
colnames(gen) <- colnames(vul) <- colnames(pdi) <- colnames(con) <- colnames(mod) <- colnames(isph) <- colnames(ispp) <- c("plot","ind","trt","block")
plotdat <- rbind(gen, vul, pdi, con, mod, isph, ispp)

plotdat$type <- rep(c("Generality", "Vulnerability", "Specialization PDI", "Connectance", 
      "Modularity", "Herbivore Species", 
      "Predator Species"), each=33)

p <- ggplot(plotdat, aes(x = trt, y = ind))
p + stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                                geom="errorbar", width=0.1, 
                 color = "grey60", lwd=1.5) +
  geom_jitter(width = 0.1) + 
  theme_bw() + 
  facet_wrap(~type, scales="free")

# Add plant species richness to the plotdat
sr <- as.data.frame(tapply(plants$SPEC, plants$CODE, function(x){length(unique(x))}))
plotdat$sr <- sr[plotdat$plot,]

# Tests! - consider transformations
# for(type in unique(plotdat$type)){
#   print(type)
#   subdat <- plotdat[plotdat$type == type, ]
#   sublm <- lmer(ind~trt+sr+(1|block), data=subdat)
#   print(summary(sublm))
# }

nms <- unique(plotdat$type)
type <- nms[1]
print(type)
subdat <- plotdat[plotdat$type == type, ]
sublm <- lmer(ind~trt+sr+(1|block), data=subdat)
print(summary(sublm))
plot(sublm)

##########>>>>>>>>>>>>>>>>>>>>>>>>>>>>> WARNING - something wrong here with treatments
# 


# Dealing with singularity because of the optimizers - use nlme::lme()

table(subdat$block, subdat$trt)
table(treats$treat, treats$codes) #

vullme <- lmer(vul~trt+(1|block), data=genvuldf) 
genlme <- lmer(gen~trt+(1|block), data=genvuldf) 
summary(genlme)
summary(vullme)

# Controling for the plant diversity
sr <- as.data.frame(tapply(plants$SPEC, plants$CODE, function(x){length(unique(x))}))
sr$code <- rownames(sr)
colnames(sr) <- c("sr", "code")
genvuldf$sr <- sr[genvuldf$plot, "sr"]
colnames(genvuldf)
# 
vullmesr <- lmer(vul~trt+sr+(1|block), data=genvuldf) 
genlmesr <- lmer(gen~trt+sr+(1|block), data=genvuldf)

summary(vullmesr)
summary(vullme)

anova(vullme, vullmesr)

summary(genlmesr)
summary(genlme)

anova(genlme, genlmesr)

# What can I say using this data?
