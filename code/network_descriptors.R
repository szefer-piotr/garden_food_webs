# Network descriptors ----
rm(list=ls())

insects <- read.table("datasets/arthropods_clean.txt")
treats  <- read.table("datasets/treatments_clean.txt")
plants  <- read.table("datasets/plants_clean.txt")
size_dat <-read.table("datasets/size_dat_bio.txt")

library(bipartite)
library(lme4)
library(lmerTest)
library(emmeans)
library(multcomp)
library(glmmTMB)
library(ggplot2)

# Modularity sourcer code
# source("code/bio_log_ratio.R")
source("code/data_processing_code.R")
source("code/contingencyTable.R")
source("code/weighted-modularity-LPAwbPLUS/code/R/LPA_wb_plus.R")
source("code/weighted-modularity-LPAwbPLUS/code/R/MODULARPLOT.R")
source("code/weighted-modularity-LPAwbPLUS/code/R/convert2moduleWeb.R")
source("code/weighted-modularity-LPAwbPLUS/code/R/GetModularInformation.R")

# Save original list of networks
gardnetsorrig <- gardnets
# gardnets - networks in individual plots
# biofulldf biomass of all insects and trees

# Switch datasets
# gardnets <- gardnetsorrig
gardnets <- abugardnets # abundance based networks


# Dont run this - it will be taken care of later
# Remove intermediate predators from the networks
# for(garden in names(gardnets)){
#   subnet <- as.matrix(gardnets[[garden]])
#   colstorem <- grep(c("mant|aran"),colnames(subnet))
#   cnms <- colnames(subnet)[!(1:dim(subnet)[2] %in% colstorem)]
# 
#   if(length(colstorem) == 0){
#     subnet <- subnet
#   } else {subnet <- (as.matrix(subnet[,cnms]))}
#   rownames(subnet)
#   gardnets[[garden]] <- subnet
# }
# 
# gardnets[["w1g2p4"]] <- t(gardnets[["w1g2p4"]])

# Remove cole001 from sampled networks! This should be done in order to alalyse effect of weevils on the REST of the community! It might be bad idea to do this if cole001 is some kind of key species in the network
treats_to_plot <- as.character(unique(treats$treat))[c(6,3,4,5,2)]
treats_to_remove_cole <- as.character(unique(treats$treat)[c(2,3,4,5)])
treats_to_remove_cole <- treats_to_plot
sites_to_remove_cole <- treats[treats$treat %in% treats_to_remove_cole, ]$codes
strc <- as.character(sites_to_remove_cole)
gardnets_nocole <- gardnets

# Removing cole001 from all networks
for (net in names(gardnets_nocole)){
  if(net %in% strc){
    subnet <- gardnets_nocole[[net]]
    gardnets_nocole[[net]] <- subnet[, colnames(subnet) != "cole001"]
  }
}

gardnets <- gardnets_nocole
gardnets[["w1g2p4"]] <- t(gardnets[["w1g2p4"]])
gardnets <- gardnets[names(gardnets) != "w1g1p6"]
gardnets <- gardnets[names(gardnets) != "w1g2p5"]

# Dataset with descriptors ----
gnames <- as.character(treats$codes)


# See which network has only one plant species and whether randomization would work on them
# for(grd in names(gardnets)){
#   print(grd)
#   print(dim(gardnets[[grd]]))
# }
# 
# vaznull(10, gardnets[["w1g2p4"]])
# One line netwrok can also be randomized

randomizations <- 10
genvuldf <- data.frame()

# gname <- "w1g2p3"
gname <- names(gardnets)[[2]]

# CALCULATE ----
for(gname in names(gardnets)){
  
  print(gname)
  
  # Create a list of randomized networks
  MAT <- gardnets[[gname]]
  
  # Randomization!
  RandMAT <- vaznull(randomizations, MAT)
  
  randomAverageData <- data.frame()

  # Calculate indices for randomized networks
  for(net in 1:length(RandMAT)){
    rMAT <- RandMAT[[net]]

    # rownames(rMAT) <- rownames(MAT)
    colnames(rMAT) <- colnames(MAT)

    rMAT_noip <- rMAT[,-grep(c("mant|aran"),
                             colnames(rMAT))]

    tt <- tryCatch(networklevel(rMAT,
                                index = "vulnerability"),
                   error=function(e) e,
                   warning=function(w) w)
    ifelse(is(tt,"warning"),
           next,
           print(net))

    if(is.null(gardnets[[gname]])){
      next
    }

    nlres <- networklevel(rMAT_noip, index = "vulnerability")
    pdi <- mean(PDI(rMAT_noip, normalise = F, log=T))
    con <-  networklevel(rMAT_noip, index = "connectance")

    # Randomized
    mod1 = DIRT_LPA_wb_plus(rMAT_noip)
    MODinformation = GetModularInformation(rMAT_noip,mod1)
    mod = MODinformation$normalised_modularity

    ispp = ncol(rMAT[, substr(colnames(rMAT),
                              1,
                              4) %in% c("aran",
                                        "mant")])
    isph = ncol(rMAT) - ispp

    ntd = networklevel(rMAT_noip, index = "nestedness")
    asym = networklevel(rMAT_noip, index = "web asymmetry")

    rsubgvdf <- data.frame(plot = gname,
                          gen = nlres[1],
                          vul = nlres[2],
                          pdi = pdi,
                          con = con,
                          mod = mod,
                          isph=isph,
                          ispp = ispp,
                          nestedness = ntd,
                          asym = asym)
    randomAverageData <- rbind(randomAverageData, rsubgvdf)

  }

  rAD <- apply(randomAverageData[, -1],
                             2,
                             mean)
  
  # Remove ips
  MAT_noip <- MAT[,-grep(c("mant|aran"),colnames(MAT))]
  
  # Some networks are too small, print their names, count them!
  tt <- tryCatch(networklevel(MAT, index = "vulnerability"),error=function(e) e, warning=function(w) w)
  ifelse(is(tt,"warning"),next,print("OK"))
  
  if(is.null(gardnets[[gname]])){
    next
  }
  
  nlres <- networklevel(MAT_noip, index = "vulnerability")
  pdi <- mean(PDI(MAT_noip, normalise = F, log=T))
  con <-  networklevel(MAT_noip, index = "connectance")
  
  mod1 = DIRT_LPA_wb_plus(MAT_noip)
  MODinformation = GetModularInformation(MAT_noip,mod1)
  mod = MODinformation$normalised_modularity
  
  # https://royalsocietypublishing.org/doi/full/10.1098/rsos.140536
  # 1 when all links exists between modules
  
  ispp = ncol(MAT[, substr(colnames(MAT),1,4) %in% c("aran", "mant")])
  isph = ncol(MAT) - ispp
  
  ntd = networklevel(MAT_noip, index = "nestedness")
  asym = networklevel(MAT_noip, index = "web asymmetry")
  
  
  subgvdf <- data.frame(plot = gname, 
                        gen = nlres[1],
                        vul = nlres[2], 
                        pdi = pdi, 
                        con = con, 
                        mod = mod, 
                        isph=isph,
                        ispp = ispp, 
                        nestedness = ntd,
                        asym = asym, # here
                        rgen = rAD[1],
                        rvul = rAD[2], 
                        rpdi = rAD[3], 
                        rcon = rAD[4], 
                        rmod = rAD[5], 
                        risph= rAD[6],
                        rispp =rAD[7], 
                        rnestedness = rAD[8],
                        rasym = rAD[9])
  genvuldf <- rbind(genvuldf, subgvdf)
}

#### ------

rownames(treats) <- treats$codes
genvuldf$trt <- treats[as.character(genvuldf$plot),
                       ]$treat
genvuldf$block <- substr(genvuldf$plot, 3,4)

genvuldf$rmod_val <- with(genvuldf, {
  (mod - rmod)/rmod
})
genvuldf$rpdi_val <- with(genvuldf, {
  (pdi - rpdi)/rpdi
})

genvuldf$mod 

genvuldf$rmod 

# Network size and PDI index ----
nsizedf <- data.frame()
for (plt in genvuldf$plot){
  subnet <- gardnets[[plt]]
  subnet <- subnet[,-grep(c("mant|aran"),
                          colnames(subnet))]
  nsizerow <- data.frame(plants = dim(subnet)[1],
                         invert = dim(subnet)[2],
                         pdi = genvuldf[genvuldf$plot == plt, ]$pdi,
                         treatment = treats[treats$codes == plt, ]$treat )
  nsizedf <- rbind(nsizedf, nsizerow)
}

ggplot(nsizedf, aes(x = pdi, y = plants, color = treatment))+
  geom_point() + 
  geom_smooth(method="lm")

lm1 <- lm(pdi~plants*treatment, data = nsizedf)
summary(lm1)
lm1 <- lm(pdi~plants, data = nsizedf[nsizedf$treatment == "INSECTICIDE",])
summary(lm1)


# PDI doesnt seem to be dependent on the network size

# PLots ----
# But generality would also change if the number of plant at the plot is smaller!


gen <- genvuldf[,c("plot","gen","trt","block")]
vul <- genvuldf[,c("plot","vul","trt","block")]
pdi <- genvuldf[,c("plot","pdi","trt","block")]
con <- genvuldf[,c("plot","con","trt","block")]
mod <- genvuldf[,c("plot","mod","trt","block")]
isph <- genvuldf[,c("plot","isph","trt","block")]
ispp <- genvuldf[,c("plot","ispp","trt","block")]
ntd <- genvuldf[,c("plot","nestedness","trt","block")]


# Change names
colnames(gen) <- colnames(vul) <- colnames(pdi) <- colnames(con) <- colnames(mod) <- colnames(isph) <- colnames(ispp) <- colnames(ntd) <- c("plot","ind","trt","block")
plotdat <- rbind(gen, vul, pdi, con, mod, isph, ispp, ntd)

plotdat$type <- rep(c("Generality", "Vulnerability", "Specialization PDI", "Connectance", 
      "Modularity", "Herbivore Species", 
      "IP Species", "Nestedness"), 
      each=length(as.character(genvuldf$plot)))

# pdf("manuscript/figs/descriptors.pdf", 7, 7)

# clip ds to only C vs P comparison

# Facet plot all descriptors ----
# p <- ggplot(plotdat, aes(x = trt, y = ind))
# p + stat_summary(fun.data=mean_cl_boot, 
#                                 geom="pointrange", width=0.1, 
#                  color = "red", lwd=1) +
#   stat_summary(fun.y=mean, geom="point", color="red", cex = 2) +
#   geom_jitter(width = 0.1, col = rgb(128,128,128, alpha = 100, maxColorValue = 255)) + 
#   theme_bw() + 
#   theme(axis.text.x=element_text(angle=90, size=5, hjust=0.5))+
#   facet_wrap(~type, scales="free")

# dev.off()

# Add plant species richness to the plotdat
sr <- as.data.frame(tapply(plants$SPEC, plants$CODE, function(x){length(unique(x))}))

plotdat$sr <- sr[plotdat$plot,]

# Test for differences

# Filter data
desStatMod <- function(data, descriptor){
  print(descriptor)
  subdat <- data[data$type == descriptor, ]
  sublm <- lmer(ind ~ trt + (1|block), data=subdat)
  print(summary(sublm))
  return(subdat)
}

library(lmerTest)
logit <- function(x){log(x/(1-x))}

clippedpd <- plotdat[plotdat$trt %in% treats_to_plot,]

sigcol <- rgb(255,0,0,150,maxColorValue = 255)
nsigcol <- rgb(211,211,211,150,maxColorValue = 255)

clippedpd$colors <- nsigcol

# 1. Connectance - significant somewhat----
condat <- desStatMod(clippedpd, "Connectance")
condat$trt <- factor(condat$trt, 
                     levels = c("CONTROL",
                                "PREDATOR",
                                "WEEVIL125",
                                "WEEVIL25",
                                "INSECTICIDE"))
# conlmer1 <- lmer(logit(ind)~trt+(1|block), condat)
# conlme1 <- nlme::lme(logit(ind)~trt, random = ~1|block, condat)
# summary(conlme1)

# Beta distribution
brrand <- glmmTMB(ind ~ trt + (1|block), data = condat, 
                  family= beta_family(link = "logit"))
sbrand <- summary(brrand)

# Do it manually
cond1 <- clippedpd$trt %in% c("WEEVIL25","WEEVIL125")
cond2 <- clippedpd$type %in% "Connectance"
clippedpd[cond1 & cond2,]$colors <- sigcol

# Grradient significance
# gradvec <- c(1,2,3,4,5)
# names(gradvec) <- treats_to_plot
# condat$graddat <- gradvec[as.character(condat$trt)]
# 
# brrandgrad <- glmmTMB(ind ~ graddat + (1|block), data = condat, 
#                   family= beta_family(link = "logit"))
# summary(brrandgrad)
# library(mgcv)
# gam1 <- gam(ind~s(graddat, k=5), data=condat)
# summary(gam1)
# plot(gam1,pages=1,residuals=TRUE)  ## show partial residuals
# plot(gam1,pages=1,seWithMean=TRUE) ## `with intercept' CIs
# ## run some basic model checks, including checking
# # ## smoothing basis dimensions...
# gam.check(gam1)
# summary(lm(ind ~ poly(graddat, 3), data=condat))
# # linear fits better... 

# 2. Generality - NOT SIGNIFICANT ----
library(truncreg)
gendat <- desStatMod(clippedpd, "Generality")
trunc <- summary(truncreg(ind~trt, gendat))
genlmer <- lm(ind~trt, gendat)
summary(genlmer)
plot(ind~trt,gendat)

# Herbivore species *
hsdat <- desStatMod(clippedpd, "Herbivore Species")
summary(glmer.nb(ind~trt+(1|block), hsdat))



# Modularity NS
# Can modularity be modified by plant species number?
moddat <- desStatMod(clippedpd, "Modularity")
# summary(lmer(logit(ind)~trt+(1|block), moddat))
# plot(moddat$ind~ moddat$trt)
# summary(betareg::betareg(ind~trt, data=moddat))

library(glmmTMB)
plot(ind~trt,moddat)

# relevel treatments
# moddat$trt <- factor(moddat$trt, 
#                      levels = c("PREDATOR",
#                                 "WEEVIL125",
#                                 "CONTROL",
#                                 "WEEVIL25"))

br <- glmmTMB(ind ~ trt+sr, data = moddat, 
              family= beta_family(link = "logit"))
brrand <- glmmTMB(ind ~ trt+sr+(1|block), 
                  data = moddat, 
                  family= beta_family(link = "logit"))
summary(brrand)
summary(br)

brtest <- emmeans(br, "trt")
brrandtest <- emmeans(brrand, "trt")
pairwise <- cld(brtest, Letter="abcdefghijklm")
pairwise <- cld(brrandtest, Letter="abcdefghijklm")

# Something is wrong with blocks in this dataset! - SOLVED
# ggplot(moddat, aes(y=ind, x =trt, col=block)) + 
#   geom_jitter(width = 0.1, cex = 5)

# gradient not significant for modularity
gradvec <- c(1,2,3,4,5)
names(gradvec) <- treats_to_plot
moddat$graddat <- gradvec[as.character(moddat$trt)]
modgrad <- glmmTMB(ind ~ graddat+sr+(1|block), 
                  data = moddat, 
                  family= beta_family(link = "logit"))
summary(modgrad)


# IP species * 
ipdat <- desStatMod(clippedpd, "IP Species")
summary(glmer.nb(ind~trt+(1|block), ipdat))

# Nestedness NS
nsdat <- desStatMod(clippedpd, "Nestedness")
nslmer <- lmer(ind~trt+(1|block), nsdat)
summary(nslmer)

# Specialization PDI not! 
pddat <- desStatMod(clippedpd, "Specialization PDI")
pdlmer1 <- lmer(ind~trt+(1|block), pddat)
pdlme <- nlme::lme(ind~trt, random=~1|block, pddat)


# pdlmer2 <- lmer(logit(ind)~trt+(1|block), pddat)
summary(pdlmer1)
summary(pdlme)

# summary(pdlmer2)

# See if the gradient is significant
gradvec <- c(1,2,3,4,5)
names(gradvec) <- treats_to_plot
pddat$graddat <- gradvec[as.character(pddat$trt)]
pdlme <- nlme::lme(ind~graddat, random=~1|block, pddat)
summary(pdlme)

# gradient seems to result in a decrease
# Gam
# gam1 <- gam(ind~s(graddat, k=3), data=pddat)
# summary(gam1)
# plot(gam1,pages=1,residuals=TRUE)  ## show partial residuals
# plot(gam1,pages=1,seWithMean=TRUE) ## `with intercept' CIs
# ## run some basic model checks, including checking
# ## smoothing basis dimensions...
# gam.check(gam1)

# Vulnerability * 
vuldat <- desStatMod(clippedpd, "Vulnerability")
summary(truncreg(ind~trt, vuldat))
vullmer <- lmer(ind~trt+(1|block), vuldat)
summary(vullmer)

# sigind <- c("Vulnerability",
#             "Specialization PDI",
#             "Herbivore Species",
#             "IP Species",
#             "Modularity")
# 
# clippedpd <- clippedpd[clippedpd$type %in% sigind, ]

# Clipped plot ----

plotdat_filtered <- plotdatrep(, 
                    each=length(as.character(genvuldf$plot)))

clippedpd$trt <- factor(clippedpd$trt, 
                        levels = treats_to_plot)

indices_to_plot <- c("Generality", "Vulnerability", "Specialization PDI", "Connectance", 
                     "Modularity", "Nestedness")

clippedpd_itp <- clippedpd[clippedpd$type %in% indices_to_plot, ]

sigcolors <- rep(nsigcol, 
                 length(indices_to_plot)*length(treats_to_plot))

# Connectance sig
sigcolors[2] <- "black"
sigcolors[c(4,5)] <- sigcol

p <- ggplot(clippedpd_itp, aes(x = trt, y = ind))
p + stat_summary(fun.data=mean_cl_boot,
                 geom="pointrange", width=0.1,
                 color = sigcolors, lwd=1) +
  stat_summary(fun.y=mean, 
               geom="point", 
               color=sigcolors, cex = 2) +
  geom_jitter(width = 0.1, col = rgb(128,128,128, alpha = 100, maxColorValue = 255)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=0, size=5, hjust=0.5))+
  facet_wrap(~type, scales="free")+
  ylab("")+xlab("")

# stat_summary(fun.data=mean_cl_boot,
# geom="pointrange", color= colors, width=0.2, lwd=1.5) +
#   stat_summary(fun.y=mean, geom="point", color=colors, cex = 5) +
#   theme(axis.text.x=element_text(angle=0, size=20, hjust=0.5)

# Tests! - consider transformations
# for(type in unique(plotdat$type)){
#   print(type)
#   subdat <- plotdat[plotdat$type == type, ]
#   sublm <- lmer(ind~trt+sr+(1|block), data=subdat)
#   print(summary(sublm))
# }
# 
# for(dstp in unique(plotdat$type)){
#   print("_________________________")
#   desStatMod(plotdat, dstp)
# }

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
# vullmesr <- lmer(vul~trt+sr+(1|block), data=genvuldf) 
# genlmesr <- lmer(gen~trt+sr+(1|block), data=genvuldf)
# 
# summary(vullmesr)
# summary(vullme)
# 
# anova(vullme, vullmesr)
# 
# summary(genlmesr)
# summary(genlme)
# 
# anova(genlme, genlmesr)

#  IP/herb ----
# gardnets
# for(gard in names(gardnets)){
#   print(dim(gardnets[[gard]]))
# }
# 
# psites <- as.character(treats[treats$treat %in% c("PREDATOR"), ]$codes)
# csites <- as.character(treats[treats$treat %in% c("CONTROL"), ]$codes)
# 
# pihratio <- data.frame()
# for(gard in psites){
#   print(gard)
#   submat <- gardnets[[gard]]
#   ipcols <- grep("aran|mant", colnames(submat))
#   ipbio <- sum(colSums(submat[, ipcols]))
#   hbio <- sum(colSums(submat[, -ipcols]))
#   iphr <- ipbio/hbio
#   phrow <- data.frame(site = gard, 
#                        trt = "predator",
#                        iphr = iphr)
#   # sum(submat) == ipbio+hbio
#   pihratio <- rbind(pihratio, phrow)
# } 
# 
# cihratio <- data.frame()
# for(gard in csites){
#   print(gard)
#   submat <- gardnets[[gard]]
#   ipcols <- grep("aran|mant", colnames(submat))
#   ipbio <- sum(colSums(submat[, ipcols]))
#   hbio <- sum(colSums(submat[, -ipcols]))
#   iphr <- ipbio/hbio
#   chrow <- data.frame(site = gard, 
#                       trt = "control",
#                       iphr = iphr)
#   # sum(submat) == ipbio+hbio
#   cihratio <- rbind(cihratio, chrow)
# } 
# 
# iphratio <- rbind(cihratio,pihratio)
# iphratio$block <- substr(iphratio$site, 3, 4)
# library(ggplot2)
# # p <- ggplot(iphratio, aes(x =trt, y = log(iphr)))
# # p + geom_jitter()
# 
# plot(log(iphr)~trt, iphratio)
# library(lme4)
# library(lmerTest)
# summary(lmer(log(iphr)~trt+(1|block), iphratio))
# 
