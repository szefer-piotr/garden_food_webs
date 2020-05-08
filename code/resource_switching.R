# Resource switching

rm(list=ls())
source("code/data_processing_code.R")
source("code/pdi.R")

ins_bio$block <- substr(ins_bio$plot, 3,4)
treats$treat <- tolower(treats$treat)
treats$block <- substr(treats$codes, 3,4)
# predator vs control
cpsites <- as.character(treats[treats$treat %in% c("control", 
                                      "predator"),]$codes)
cpfull <- ins_bio[ins_bio$plot %in% cpsites, ]
  
# remove ip
cpfull_noip <- cpfull[-grep("aran|mant", cpfull$morphotype), ]
cpfull_noip$morphotype <- as.character(cpfull_noip$morphotype)

# For each block get species which are comparable
bl <- "g1"

slgfulldat <- data.frame()
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
  `%notin%` <- Negate(`%in%`)
  lost <- colnames(cmat)[colnames(cmat) %notin% colnames(pmat)]
  #gained
  gained <- colnames(pmat)[colnames(pmat) %notin% colnames(cmat)]
  
  allnames <- length(unique(c(colnames(pmat), colnames(cmat))))
  length(c(stayed, lost, gained)) == allnames # do we have all of them?
  
  # pdi
  staydf <- data.frame(pdi = diet_breadth[stayed],
                       type = "stayed")
  lostdf <- data.frame(pdi = diet_breadth[lost],
                       type = "lost")
  gainedf <- data.frame(pdi = diet_breadth[gained],
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
  subbl <- cpfull_noip[cpfull_noip$block == bl,]
  psite <- treats[treats$block == bl & treats$treat == "predator", ]$codes
  csite <- treats[treats$block == bl & treats$treat == "control", ]$codes
  pn <- gardnets[[psite]]
  cn <- gardnets[[csite]]
  mtnms <- colnames(cn)[colnames(cn) %in% colnames(pn)]
  ratios <- colSums(pn[,mtnms])/colSums(cn[,mtnms])
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
compdf_niop$bio <- morphbios[compdf_niop$spec] 
library(RColorBrewer)
library(ggplot2)

dotcols <- brewer.pal(6, "BrBG")
dotcols <- alpha(dotcols, 0.5)

compdf_niop$lratio <- log(compdf_niop$ratio)
compdf_niop$alratio <- abs(compdf_niop$lratio)

plot(lratio~pdi, data = compdf_niop, 
     pch = 19, col = dotcols, cex = log(compdf_niop$bio*10))
lmer1 <- lmer(lratio~pdi+(1|block), 
     data = compdf_niop,
     weights = bio)
summary(lmer1)
abline(-0.04862,0.24944, lwd = 3,  col = alpha("black", 0.5))
abline(h=0,col = alpha("black", 0.5), lty = 2)

# plot
library(ggplot2)
gray <- rgb(0,0,0,50,maxColorValue = 255)
ggplot(slgfullnozero, aes(x = type, y = log(pdi), 
                       color = ))+
  geom_jitter(width = 0.1, color = gray, pch=19) + 
  stat_summary(fun.y = mean, geom = "point", col= "red")+
  stat_summary(fun.data = "mean_cl_boot", 
               geom = "errorbar",
               width=0.05, col="red", lwd=1.1)

# There is no difference in PDI for species lost, gained or stayed in the plot
library(lme4)
library(lmerTest)
library(emmeans)

lmer1 <-(lmer(log(pdi)~type+(1|block), slgfullnozero))

inter.test1 <- emmeans(lmer1, "type")
plot(inter.test1)

# PDA vs ABUNDANCE ----
# Were species with high PD more abundnat?

diet_breadth
diet_breadth["orth024"]
abund <- ins_bio[, c("morphotype","totbio")]
abund$morphotype <- as.character(abund$morphotype)
abund$db <- diet_breadth[abund$morphotype]
abund_noip <- abund[-grep("aran|mant", abund$morphotype), ]

plot(totbio~db ,abund_noip, pch=19, col = gray,
     xlab = "PDI", ylab="Biomass")
library(MASS)
glm1 <- glm(totbio~db, abund_noip, family = gaussian(link="log"))
summary(glm1)
ndat <- data.frame(db = seq(0,3.5,by=0.01))
res <- predict(glm1, newdata = ndat, 
               type = "response",
               se.fit = T)
lines(res$fit~ndat$db, lwd =2, col = "red")

# 2.Calculate distances between diets

# 3. See wether species reduced or switched its resources

# 4. How abundance changed depending on plant quality.