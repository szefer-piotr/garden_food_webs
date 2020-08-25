# Resource switching

rm(list=ls())
source("code/data_processing_code.R")
source("code/pdi.R")
# source("code/bio_log_ratio.R")
source("code/diet_shift.R")

ins_bio$block <- substr(ins_bio$plot, 3,4)
treats$treat <- tolower(treats$treat)
treats$block <- substr(treats$codes, 3,4)
# predator vs control
cpsites <- as.character(treats[treats$treat %in% c("control", 
                                      "predator"),]$codes)
cpfull <- ins_bio[ins_bio$plot %in% cpsites, ]
cpfull$morphotype <-as.character(cpfull$morphotype)

# remove ip
cpfull_noip <- cpfull

# For each block get species which are comparable
# bl <- "g1"

slgfulldat <- data.frame()

# bl <- unique(cpfull_noip$block)[1]

`%notin%` <- Negate(`%in%`)

for(bl in unique(cpfull_noip$block)){
  subbl <- cpfull_noip[cpfull_noip$block == bl,]
  psite <- treats[treats$block == bl & treats$treat == "predator", ]$codes
  csite <- treats[treats$block == bl & treats$treat == "control", ]$codes
  psite <- as.character(psite)
  csite <- as.character(csite)
  pdat <- subbl[subbl$plot == psite, ]
  cdat <- subbl[subbl$plot == csite, ]
  
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


# Weighted lratio of individual species change vs PDI ----

compdf <- data.frame()

for(bl in unique(treats$block)){
  
  # Any species, that was found in pred and control site in any block
  subbl <- cpfull_noip[cpfull_noip$block == bl,]
  amount <- subbl$amount
  psite <- treats[treats$block == bl & treats$treat == "predator", ]$codes
  csite <- treats[treats$block == bl & treats$treat == "control", ]$codes
  
  pn <- gardnets[[psite]]
  cn <- gardnets[[csite]]
  
  mtnms <- colnames(cn)[colnames(cn) %in% colnames(pn)]
  # High abundance in predatry exclosure suggest strong effect of predators
  ratios <- colSums(cn[,mtnms])/colSums(pn[,mtnms])
  dbs <- diet_breadth[mtnms]
  plot(log(ratios)~dbs)
  bldf <- data.frame(spec = mtnms, ratio = ratios, pdi = dbs, block = bl)
  compdf <- rbind(compdf, bldf)
}

compdf_niop <- compdf[-grep("aran|mant", compdf$spec), ]

library(lme4)
library(lmerTest)

# weights
ins_bio$morphotype <- as.character(ins_bio$morphotype)
morphbios  <- tapply(ins_bio$totbio, ins_bio$morphotype, sum, na.rm=T)
morphbios  <- tapply(ins_bio$amount, ins_bio$morphotype, sum, na.rm=T)

compdf_niop$bio <- morphbios[compdf_niop$spec] 
library(RColorBrewer)
library(ggplot2)

dotcols <- brewer.pal(6, "BrBG")
dotcols <- alpha(dotcols, 0.5)

# Adding the diet switch ... capabilities?
#Examine log ratios
lrnm <- "homo004"
compdf_niop[compdf_niop$spec == lrnm, ]


compdf_niop$lratio <- log(compdf_niop$ratio)
compdf_niop$alratio <- abs(compdf_niop$lratio)
rownames(shiftDf) <- shiftDf$species
compdf_niop$spec <- as.character(compdf_niop$spec)
compdf_niop$shift <- shiftDf[compdf_niop$spec,]$shift

# PDI and lratio ----
plot(lratio~pdi, data = compdf_niop, 
     pch = 19, col = rgb(180,180,180,80,maxColorValue = 255), 
     cex = log(compdf_niop$bio*10))

# Why is that?
lm1  <- lm(lratio~pdi, 
           weights = bio,
           data =compdf_niop)
summary(lm1)
lmer1 <- lmer(lratio~pdi+(1|block), 
     data = compdf_niop,
     weights = bio)

# Example
fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)

# Full model
lmer1b <- lmer(lratio~pdi+(pdi|block), 
               data = compdf_niop,
               weights = bio)
as.data.frame(VarCorr(lmer1b))
as.data.frame(VarCorr(lmer2))

summary(lmer1b)
summary(fm1)

#nlme model: should I use intercept 0 ???
lme1a <- nlme::lme(lratio~pdi, random = ~ 1|block, 
                   data = compdf_niop)
summary(lme1a)
lme1b <- nlme::lme(lratio~pdi, random = ~ pdi|block, 
                  data = compdf_niop)
summary(lme1b)
VarCorr(lme1b)

attach(compdf_niop)
names(compdf_niop) <- c("spec","ratio", "pdi", "block", "w", 
                        "lratio", "alratio", "shift" )
lme1c <- nlme::lme(lratio~pdi, random = ~ pdi|block, 
               data = compdf_niop,
               weights = ~ as.vector(w))
detach(compdf_niop)
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


# plot
library(ggplot2)
gray <- rgb(0,0,0,50,maxColorValue = 255)
ggplot(slgfullnozero, aes(x = type, y = pdi, 
                       color = ))+
  geom_jitter(width = 0.1, color = gray, pch=19) + 
  stat_summary(fun = mean, geom = "point", col= "red")+
  stat_summary(fun.data = "mean_cl_boot", 
               geom = "errorbar",
               width=0.05, col="red", lwd=1.1)

# There is no difference in PDI for species lost, gained or stayed in the plot
library(lme4)
library(lmerTest)
library(emmeans)

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

