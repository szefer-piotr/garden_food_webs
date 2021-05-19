# Resource switching

rm(list=ls())
source("code/data_processing_code.R")
source("code/pdi.R")
# source("code/bio_log_ratio.R")
source("code/diet_shift.R")

library(lme4)
library(lmerTest)
library(emmeans)
library(multcomp)
library(ggplot2)

# Define transparent gray color
gray <- rgb(0,0,0,50,maxColorValue = 255)

ins_bio$block <- substr(ins_bio$plot, 3,4)
treats$treat <- tolower(treats$treat)
treats$block <- substr(treats$codes, 3,4)

# predator vs control
cpsites <- as.character(treats[treats$treat %in% c("control", 
                                      "predator"),]$codes)
cpfull <- ins_bio[ins_bio$plot %in% cpsites, ]
cpfull$morphotype <-as.character(cpfull$morphotype)

# remove ip, it has already been already removed.
cpfull_noip <- cpfull

# For each block get species which are comparable
# bl <- "g2"



# 1. Stayed/lost/gained analysis ----

# 1. to which 2. basal
# treatments <- c("weevil125", "predator")
treatments <- c("weevil25", "predator")
# treatments <- c("weevil125", "control")
# treatments <- c("weevil25", "control")
# treatments <- c("predator", "control")

# remove ip, it has already been already removed.
slgfull_noip <- ins_bio
slgfulldat <- data.frame()

# bl <- unique(slgfull_noip$block)[1]
# bl <- "g2"

`%notin%` <- Negate(`%in%`)

diet_breadth <- diet_breadth_ab

for(bl in unique(slgfull_noip$block)){
  print(bl)
  subbl <- slgfull_noip[slgfull_noip$block == bl,]
  psite <- treats[treats$block == bl & treats$treat == treatments[1],
                  ]$codes
  csite <- treats[treats$block == bl & treats$treat == treatments[2],
                  ]$codes
  psite <- as.character(psite)
  csite <- as.character(csite)
  pdat <- subbl[subbl$plot == psite, ]
  cdat <- subbl[subbl$plot == csite, ]
  
  pdat <- pdat[complete.cases(pdat), ]
  cdat <- cdat[complete.cases(cdat), ]
  
  pmat <- contingencyTable2(pdat, "tree", "morphotype", "totbio")
  cmat <- contingencyTable2(cdat, "tree", "morphotype", "totbio")
  
  #stayed
  stayed <- colnames(cmat)[colnames(cmat) %in% colnames(pmat)]
  #lost
  lost <- colnames(cmat)[colnames(cmat) %notin% colnames(pmat)]
  #gained
  gained <- colnames(pmat)[colnames(pmat) %notin% colnames(cmat)]
  
  allnames <- length(unique(c(colnames(pmat), colnames(cmat))))
  length(c(stayed, lost, gained)) == allnames # do we have all of them?
  
  # Check for NAs
  dbs <- diet_breadth[stayed][!is.na(diet_breadth[stayed])]
  dbl <- diet_breadth[lost][!is.na(diet_breadth[lost])]
  dbg <- diet_breadth[gained][!is.na(diet_breadth[gained])]
  
  # pdi
  staydf <- data.frame(pdi = dbs,
                       type = "stayed")
  lostdf <- data.frame(pdi = dbl,
                       type = "lost")
  gainedf <- data.frame(pdi = dbg,
                        type = "gained")
  slgdat <- rbind(staydf, lostdf, gainedf)
  slgdat$block <- bl
  slgfulldat <- rbind(slgfulldat, slgdat)
}

slgfulldat
dim(slgfulldat)
slgfullnozero <- slgfulldat[slgfulldat$pdi != 0, ]
dim(slgfullnozero)

# plot

slgfullnozero$fam <- substr(rownames(slgfullnozero), 1, 4)
# slgfullnozero_noip <- slgfullnozero[-grep("aran|mant", 
#                                           slgfullnozero$fam), ]
slgfullnozero_noip <- slgfullnozero
slgfullnozero_noip$morpho <- substr(rownames(slgfullnozero_noip), 1,7)
slgfullnozero_noip$abundnce <- pdiss[slgfullnozero_noip$morpho]

# Exclude species observed only minabu-times
minabu <- 10
slgfullnozero_noip_minabu <- slgfullnozero_noip[slgfullnozero_noip$abundnce >= minabu, ]

# 1a. Plot----
slgfullnozero_noip_minabu$fam
slgfullnozero_noip_minabu -> abudat

abudat$fam <- as.factor(abudat$fam)
levels(abudat$fam) <- c("Coleoptera",
                        "Heteroptera",
                        "Homoptera",
                        "Lepidoptera",
                        "Orthoptera")
levels(abudat$type) <- c("C+Ex", "C","Ex")

# Calculate 

bddat <- slgfullnozero_noip_minabu
bddat[bddat$pdi == 0, ]$pdi <- 0.001
bddat[bddat$pdi == 1, ]$pdi <- 0.999

grouping <- expand.grid(levels(abudat$type),
                        levels(abudat$fam))
grouping$label <- c("a","b","b",
                    "","","",
                    "","","",
                    "a","b","ab",
                    "a","b","b")

pdicols <- c("red","grey20","grey20",
             "grey20","grey20","grey20",
             "grey20","grey20","grey20",
             "red","grey20","grey20",
             "red","grey20","grey20")

library(ggsignif)

p1 <- ggplot(abudat,aes(x = type, y = pdi))+
  geom_jitter(width = 0.1, pch=19, alpha = 0.3) +
  stat_summary(fun = mean, geom = "point", col= "red")+
  stat_summary(fun.data = "mean_cl_boot",
               geom = "pointrange",
               col=pdicols, lwd=1.1)+
  facet_grid(cols = vars(fam)) + xlab("")+ylab("Paired Differences Index")

p2 <- ggplot(abudat[abudat$type %in% c("C","Ex"),],aes(x = type, y = log(abundnce)))+
  geom_jitter(width = 0.1, pch=19, alpha = 0.3) +
  stat_summary(fun = mean, geom = "point", col= "red")+
  stat_summary(fun.data = "mean_cl_boot",
               geom = "pointrange",
               col="grey20", lwd=1.1)+
  facet_grid(cols = vars(fam))+ xlab("")+ylab("log[Abundance]")

ggpubr::ggarrange(p1,p2, nrow = 2, labels = c("A", "B"))

# 1b. Test ----
# Beta distribution - theoretically better model
library(glmmTMB)
bddat <- slgfullnozero_noip_minabu
bddat[bddat$pdi == 0, ]$pdi <- 0.001
bddat[bddat$pdi == 1, ]$pdi <- 0.999

brrand <- glmmTMB(pdi ~ type*fam+(1|block), 
                  data = bddat, 
                  family= beta_family(link = "logit"))
summary(brrand)

emmeans(brrand, "type")
emmip(brrand, type ~ type | fam)
emm_s.t <- emmeans(brrand, pairwise ~ type | fam)
plot(emm_s.t)
as.data.frame(summary(emm_s.t)$contrasts)

# Prediction
# newdata <- expand.grid(block = unique(slgfullnozero_noip_minabu$block),
#                        fam = unique(slgfullnozero_noip_minabu$fam),
#                        type = unique(slgfullnozero_noip_minabu$type))
# 
# newdata$prediction <- predict(brrand, newdata, type = 'response')
# ggplot(newdata,aes(x = type, y = prediction,
#                    color = block))+
#   geom_jitter(width = 0.1, pch=19) +
#   stat_summary(fun = mean, geom = "point", col= "red")+
#   stat_summary(fun.data = "mean_cl_boot",
#                geom = "errorbar",
#                width=0.05, col="red", lwd=1.1)+
#   facet_wrap(vars(fam))


# Abundance plot for the lost-stayed-gained

glmnb1 <- glm.nb(abundnce~type*fam,data = abudat[abudat$type %in% c("C","Ex"), ])
# glmernb1 <- glmer.nb(abundnce~type*fam+(1|block),data = abudat)

summary(glmnb1)
emmeans(glmnb1, "type")
emmip(glmnb1, type ~ type | fam)
emm_s.t <- emmeans(glmnb1, pairwise ~ type | fam)
summary(emm_s.t)
# 2. Weighted lratio of individual species change vs PDI ----
treatments <- c("predator", "weevil125")

logratioComp <- function(treatments){
  compdf <- data.frame()
  den <- treatments[1]
  num <- treatments[2]
  for(bl in unique(treats$block)){
    print(bl)
    # Any species, that was found in pred and control site in any block
    subbl <- cpfull_noip[cpfull_noip$block == bl,]
    amount <- subbl$amount
    psite <- treats[treats$block == bl & treats$treat == den, ]$codes
    csite <- treats[treats$block == bl & treats$treat == num, ]$codes
    
    pn <- gardnets[[psite]]
    cn <- gardnets[[csite]]
    
    mtnms <- colnames(cn)[colnames(cn) %in% colnames(pn)]
    # High abundance in predatry exclosure suggest strong effect of predators

    if(is.null(dim(cn[,mtnms])[1])){
      cScn <- colSums(matrix(cn[,mtnms], 
                             nrow = 1))
    }else{cScn <- colSums(cn[,mtnms])}
    
    if(is.null(dim(pn[,mtnms])[1])){
      cSpn <- colSums(matrix(pn[,mtnms], 
                             nrow = 1))
    }else{cSpn <- colSums(pn[,mtnms])}

    ratios <- cScn/cSpn
    # ratios <- colSums(cn[,mtnms])/colSums(pn[,mtnms])
    dbs <- diet_breadth[mtnms]
    
    bldf <- data.frame(spec = mtnms, 
                       ratio = ratios, 
                       pdi = dbs, 
                       effort = pdiss[mtnms], 
                       block = bl)
    compdf <- rbind(compdf, bldf)
  }
  return(compdf)
}

compdf <- logratioComp(c("predator","control"))
# compdf <- logratioComp(c("predator", "weevil125"))
# compdf <- logratioComp(c("predator", "weevil25"))
# compdf <- logratioComp(c("control", "weevil125"))

compdf_niop <- compdf[-grep("aran|mant", compdf$spec), ]

library(lme4)
library(lmerTest)
library(RColorBrewer)
library(ggplot2)

# weights, I could use abundnace as a weighting factor
ins_bio$morphotype <- as.character(ins_bio$morphotype)
morphbios  <- tapply(ins_bio$totbio, ins_bio$morphotype, sum, na.rm=T)
morphbios  <- tapply(ins_bio$amount, ins_bio$morphotype, sum, na.rm=T)
compdf_niop$bio <- morphbios[compdf_niop$spec] 


dotcols <- brewer.pal(6, "BrBG")
dotcols <- alpha(dotcols, 0.5)

# Adding the diet switch ... capabilities?
#Examine log ratios
# lrnm <- "homo004"
# compdf_niop[compdf_niop$spec == lrnm, ]


compdf_niop$lratio <- log(compdf_niop$ratio)
compdf_niop$alratio <- abs(compdf_niop$lratio)
rownames(shiftDf) <- shiftDf$species
compdf_niop$spec <- as.character(compdf_niop$spec)
compdf_niop$shift <- shiftDf[compdf_niop$spec,]$shift

# Remove single- and doubletones
minabu <- 20
compdf_niop <- compdf_niop[compdf_niop$effort > minabu, ]
# log the effort
compdf_niop$effort <- log(compdf_niop$effort)

# PDI and lratio ----
plot(lratio~pdi, data = compdf_niop, 
     pch = 19, col = rgb(180,180,180,80,maxColorValue = 255), 
     cex = log(compdf_niop$effort))

plot(lratio~pdi, data = compdf_niop, 
     pch = 19, col = compdf_niop$block, 
     cex = log(compdf_niop$effort))


# Random slopes weighted by abundance
attach(compdf_niop)
library(nlme)

# names(compdf_niop) <- c("spec","ratio", "pdi", "block", "w", 
#                         "lratio", "alratio", "shift" )

lme1c <- nlme::lme(lratio~pdi, random = ~ pdi|block, 
               data = compdf_niop,
               weights = ~ as.vector(effort))
summary(lme1c)

# Garden loop
par(mfrow=c(2,3))
for(bl in unique(compdf_niop$block)){
  subdf_niop <- compdf_niop[compdf_niop$block == bl, ]
  sub_lm1d <- lm(lratio~pdi,data = subdf_niop,
                     weights = effort)
  smry <- summary(sub_lm1d)
  pval <- smry$coefficients[2,4]
  tit <- paste(bl, smry$coefficients[2,1])
  linetype <- 1
  if(pval >= 0.05){
    linetype <- 2
  }
  plot(lratio~pdi, data=subdf_niop, main = tit)
  abline(sub_lm1d, lwd = 2, lty = linetype)
}

lme1d <- nlme::lme(lratio~1, random = ~ pdi|block, 
                   data = compdf_niop,
                   weights = ~ as.vector(effort),
                   control = lmeControl(msVerbose = TRUE))
summary(lme1d)

fit <- lmList(lratio ~ pdi | block, data=compdf_niop[, -8])
cov(coef(fit))

# No correlation of random effects/ uncostrained slopes
lme1e <- nlme::lme(lratio~pdi, random = ~ 0+pdi|block, 
                   data = compdf_niop,
                   weights = ~ as.vector(effort),
                   control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))

lme1f <- nlme::lme(lratio~pdi, random = ~ pdi|block, 
                   data = compdf_niop,
                   weights = ~ as.vector(effort),
                   control = lmeControl(opt = 'optim'))
summary(lme1e)
summary(lme1f)


detach(compdf_niop)

coef(lme1c, unconstrained = FALSE)

coef(lme1d, unconstrained = FALSE)

summary(lme1c)
VarCorr(lme1c)
plot(lme1c)

library(ggeffects)
lme_pred <- ggpredict(lme1c, terms = c("pdi"), type = "random")
plot(lme_pred)
lme_pred_block <- ggpredict(lme1c, terms = c("block"), type = "random")
plot(lme_pred_block)

#https://cran.r-project.org/web/packages/lme4/vignettes/lmer.pdf
# Unconstrained slopes!
lmer2 <- lmer(lratio~pdi+(1|block)+(0+pdi|block), 
              data = compdf_niop,
              weights = bio)



anova(lmer1, lmer2) # lmer2 is an improvement
summary(lmer2)

summary(lmer1b)


plot(lmer1, type = c("p", "smooth"))

# Within each block there is another 
# for(bl in unique(compdf_niop$block)){
#   print(bl)
#   subs <- compdf_niop[compdf_niop$block == bl, ]
#   fitbl <- lm(lratio~pdi, data=subs)
#   print(summary(fitbl))
# }

# library(brms)
# brm2 <- brm(lratio|weights(bio)~pdi+(1+pdi|block), 
#               data = compdf_niop, 
#             control = list(adapt_delta = 0.99))
# plot(brm2)

# When weighted by abundance (lookup 'morphbios') we are getting highly significant effect
# summ1 <- summary(lmer1)

library(ggeffects)
# pr <- ggpredict(lmer1, "pdi")
# pr <- ggpredict(lmer1, "pdi", type = "random", condition = c(block = factor("g1")))
# 
# 
# plot(pr)

me <- ggpredict(lmer2, terms = c("pdi"), type = "random")
plot(me)
me <- ggpredict(lmer2, terms = c("pdi", "block"), type = "random")
plot(me)
# Ther is no evidence that herbivores' shift in diet affect their log-ratio response
plot(lratio~shift,data=compdf_niop)
lratio_shift <- lmer(lratio~shift+(1|block), data=compdf_niop)
summary(lratio_shift)



sl <- summary(lmer1)
abline(sl$coefficients[1,1],
       sl$coefficients[2,1], lwd = 3,  col = alpha("black", 0.5))
abline(h=0,col = alpha("black", 0.5), lty = 2)

# Need to add confidence intervals
pdiseq <- seq(min(compdf_niop$pdi), max(compdf_niop$pdi), by=0.001)
blocks <- unique(compdf_niop$block)
newdata = expand.grid(pdi = pdiseq, block = blocks)
newdata_noref = expand.grid(pdi = pdiseq)

library(merTools)
preds <- predictInterval(lmer1, newdata = newdata, n.sims = 999)

# pnd <- predict(lmer1, newdata=newdata, interval='prediction')
newdata$lratio.upr <- preds$upr
newdata$lratio.lwr <- preds$lwr
newdata$lratio.fit <- preds$fit

plot(lratio.upr~pdi, newdata)
par(new=T)
plot(lratio.lwr~pdi, newdata)


# # plot
# library(ggplot2)
# gray <- rgb(0,0,0,50,maxColorValue = 255)
# ggplot(slgfullnozero, aes(x = type, y = pdi, 
#                        color = ))+
#   geom_jitter(width = 0.1, color = gray, pch=19) + 
#   stat_summary(fun = mean, geom = "point", col= "red")+
#   stat_summary(fun.data = "mean_cl_boot", 
#                geom = "errorbar",
#                width=0.05, col="red", lwd=1.1)

# There is no difference in PDI for species lost, gained or stayed in the plot


ins_bio$morphplot <- paste(ins_bio$morphotype, 
                            ins_bio$plot, sep="")
rownames(ins_bio) <- ins_bio$morphotype
slgfullnozero$abundance <- ins_bio

lmer1 <-(lmer(pdi~type+(1|block), slgfullnozero))

inter.test1 <- emmeans(lmer1, "type")
# plot(inter.test1)

# PDI vs ABUNDANCE ----
# Were species with high PD more abundnat?

diet_breadth
diet_breadth["orth024"]
abund <- ins_bio[, c("morphotype","totbio")]
abund$morphotype <- as.character(abund$morphotype)
abund$db <- diet_breadth[abund$morphotype]
abund_noip <- abund[-grep("aran|mant", abund$morphotype), ]

anip <- abund_noip[!is.na(abund_noip$db) & abund_noip$db != 0,]

plot(totbio~db ,
     anip, 
     pch=19, col = gray,
     xlab = "PDI", ylab="Biomass")
library(MASS)

glm1 <- glm(totbio~db, anip, family = gaussian(link="log"))
summary(glm1)

ndat <- data.frame(db = seq(0,3.5,by=0.01))
res <- predict(glm1, newdata = ndat, 
               type = "response",
               se.fit = T)
lines(res$fit~ndat$db, lwd =2, col = "red")

# 2.Calculate distances between diets
pdimat

# 3. See wether species reduced or switched its resources

# 4. How abundance changed depending on plant quality.

