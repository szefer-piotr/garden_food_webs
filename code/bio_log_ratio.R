# Data analysis ----

# Direct and indirect effects on plants can be calculated using the log ratio. ln(Vp+/Vp-) where V's are community variables (herbivore abundance and plant biomass in the presence and absence of herbivores). These effect magnitudes can be plotted in relation to each other on an x-y plane and in relation to a 45 deg reference line tha represents the equivalence in strength of direct and indirect effect of carnivores.

# The log ratio effect of carnivores on herbivores should always be negative if carnivores are limiting herbivores regardles of the of the way herbivores are resource limited.

# page 36 Resloving ecosystem complexity

# I need biomass of herbivores, arthropod predators and plant for each plot.

# Inspect that and see if there are any changes in the food web structure. What can i expect?

# source("code/data_processing_code.R")

# abufulldf # dataframe
# abugardnets # morphotype based networks for each garden
# abugardnetsfam # family aagregated networks for all garden


## * 2. Log response ratios ----

#### CUTTED
rm(list=ls())
source("code/data_processing_code.R")
source("code/pdi.R")
source("code/diet_shift.R")

library(ggplot2)

# Assume that each plant hosts unique community of insects.

# Log ratio analyses

treats$block <- substr(treats$codes,3,4)
# ins_bio with biomass for plants
plot <- "w1g2p4"
species <- "piptar"

plantBio <- function(species, plot){
  plbiomass <- main_biomass[main_biomass$CODE == plot, ]
  plsitebio <- plbiomass[plbiomass$SP_CODE == toupper(species), ]
  bio <- (plsitebio$WEIGHT)
  if(length(bio) != 0){
    print(paste("Weight of:",species, "at the plot", plot, "is", bio))
    return(bio)
  }else{
    print(paste("There is no",species, "at the plot", plot))
    return(0)
  }
}

# * 2.0 Log ratio for given plant, family and set of treatments ----

families <- c("orth", "hemi")
plants <- c("piptar", "melamu")
treatments <- c("weevil125", "control") #numerator/denominator
cumulative = TRUE
bl <- unique(subdat$block)[1]

lratioFamPlantTreat <- function(ins_bio, 
                                families, 
                                plants, 
                                treatments, 
                                cumulative = F,
                                abundance = T){
  
  treatment_plots <- treats[as.character(treats$treat) %in% toupper(treatments),]$codes
  
  ca <- ins_bio$family %in% families
  cb <- ins_bio$tree %in% plants
  cc <- ins_bio$plot %in% treatment_plots
  
  subdat <- ins_bio[(ca&cb&cc), ]
  
  if(is.null(families)){
    families <- as.character(unique(subdat$family))
  }
  if(is.null(plants)){
    plants <- as.character(unique(subdat$tree))
  }
  
  subdat$block <- substr(subdat$plot,3,4)
  num_sites <- treats[treats$treat %in% toupper(treatments[1]), ]$code
  den_sites <- treats[treats$treat %in% toupper(treatments[2]), ]$code
  
  logData <- data.frame()
  
  # For each garden
  for(bl in unique(subdat$block)){
    print(bl)
    num_site <- as.character(num_sites[grep(bl, num_sites)])
    den_site <- as.character(den_sites[grep(bl, den_sites)])
    
    if(!cumulative){
      for(pl in plants){
        print(pl)
        for(fam in families){
          print(fam)
          plbool <- subdat$tree == pl
          fambool <- subdat$family == fam
          numbool <- subdat$plot == num_site
          denbool <- subdat$plot == den_site
          
          numdat <- subdat[ plbool & fambool & numbool, ]
          dendat <- subdat[ plbool & fambool & denbool, ]
          
          if(abundance){
            numval <- sum(numdat$amount)
            denval <- sum(dendat$amount)
          }else{
            numval <- sum(numdat$totbio)
            denval <- sum(dendat$totbio)
          }
          
          plbiomass <- main_biomass[]
          
          logDataRow <- data.frame(numv = numval,
                                   denv = denval,
                                   numr = dim(numdat)[1],
                                   denr = dim(dendat)[1],
                                   plbn = plantBio(pl, num_site),
                                   plbd = plantBio(pl, den_site),
                                   fams = fam,
                                   plnt = pl,
                                   blck = bl,
                                   nplt = num_site,
                                   dplt = den_site,
                                   ntrt = treatments[1],
                                   dtrt = treatments[2])
          
          logData <- rbind(logData, logDataRow)
        }
      }
    }else{
      # cumulative logratio
      
    }
  }
  return(logData)
}

testDat <- lratioFamPlantTreat(ins_bioOrig, 
                    families=c("aran", "mant"),
                    plants,
                    treatments,
                    cumulative = F)
testDat$ilr <- log(testDat$numv/testDat$denv)

# * 2.1 General log ratio for herbivores, intermediate predators and plants ----

# Biomass of each group
# For each garden bioHp (plus), bioHm (minus), bioIPp, bioIPm, bioPp, bioPm
fam <- c("aran")
treatments <- c("weevil125", "control")
groupLogratio <-function(fam, treatments, abundance=T){
  # if no family is selected then collect them all
  genlratio <- data.frame()
  
  if(abundance){
    biodata <- abufulldf
  }else{
    biodata <- biofulldf
  }
  
  biodata <- biodata[biodata$trt %in% toupper(treatments), ]
  biodata$gard <- substr(biodata$plot, 3,4)
  biodata$group <- "herb"
  biodata[as.character(biodata$nms) %in% c("aran","mant"), ]$group <- "ips"
  
  if(length(fam) != 0){ #all families
    biodata <- biodata[biodata$nms %in% fam, ]
  }
  
  num_sites <- treats[treats$treat %in% toupper(treatments[1]), ]$codes
  den_sites <- treats[treats$treat %in% toupper(treatments[2]), ]$codes
  
  num_sites <- as.character(num_sites)
  den_sites <- as.character(den_sites)
  
  for (block in unique(biodata$gard)){
    
    num_site <- num_sites[grep(block, num_sites)]
    den_site <- den_sites[grep(block, den_sites)]
    
    subbl <- biodata[biodata$gard == block, ]
    
    numDat <- subbl[subbl$plot == num_site,  ]
    denDat <- subbl[subbl$plot == den_site,  ]
  
    numDiv <- tapply(numDat$bio, 
                     numDat$group, sum)
    denDiv <- tapply(denDat$bio, 
                     denDat$group, sum)
    numPl <- sum(tapply(numDat$plbio, 
                    as.character(numDat$plnm), unique))
    denPl <- sum(tapply(denDat$plbio, 
                    as.character(denDat$plnm), unique))
    
    numhR <- as.numeric(numDiv["herb"])
    denhR <- as.numeric(denDiv["herb"])
    
    numipR <- as.numeric(numDiv["ips"])
    denipR <- as.numeric(denDiv["ips"])
    
    genlratioRow <- data.frame(numh = numhR,
                               denh = denhR,
                               numip = numipR,
                               denip = denipR,
                               numpl = numPl,
                               denpl = denPl,
                               block = block,
                               comp = paste(treatments[1],"/", treatments[2]),
                               fams = paste(fam, collapse=" "))
    genlratio <- rbind(genlratio,genlratioRow)
  }
  return(genlratio)
}

groupLogratio(c(), c("weevil125", "control"), abundance=T)

treat_pair_list <- list(c("weevil125", "control"),
                        # c("weevil125", "predator"),
                        c("weevil25", "control"),
                        # c("weevil25", "predator"),
                        c("control", "insecticide"),
                        c("control", "predator"))

treat_pair_list <- list(# c("weevil125", "control"),
                        # c("weevil125", "predator"),
                        # c("weevil25", "control"),
                        # c("weevil25", "predator"),
                        # c("control", "insecticide"),
                        c("control", "predator"))

# Plots for families and general log ratio
plotDataFamilieLogRatio <- data.frame()
for(fams in unique(ins_bioOrig$family)){
  print(fams)
  for(trtnum in 1:length(treat_pair_list)){
    print(treat_pair_list[[trtnum]])
    pdfR <- groupLogratio(fams, treat_pair_list[[trtnum]], abundance=T)
    plotDataFamilieLogRatio <- rbind(plotDataFamilieLogRatio, pdfR)
    pdfRG <- groupLogratio(c(), treat_pair_list[[trtnum]], abundance=T)
  }
  
}

# General
for(trtnum in 1:length(treat_pair_list)){
  print(treat_pair_list[[trtnum]])
  pdfR <- groupLogratio(c(), treat_pair_list[[trtnum]], abundance=T)
  plotDataFamilieLogRatio <- rbind(plotDataFamilieLogRatio, pdfR)
}

# Calculate log(ratios)
plotDFLR <- plotDataFamilieLogRatio
plotDFLR$lratioH <- log(plotDFLR$numh/plotDFLR$denh)
plotDFLR$lratioIP <- log(plotDFLR$numip/plotDFLR$denip)
plotDFLR$lratioPL <- log(plotDFLR$numpl/plotDFLR$denpl)

plotDFLR$fams <- as.character(plotDFLR$fams)
plotDFLR[plotDFLR$fams == "", ]$fams <- "cumulative"

# Herbivore vs plants
HvP <- ggplot(plotDFLR[plotDFLR$fams == "cumulative", ], 
       aes(y = lratioH, x=lratioPL))+
  geom_point()+
  geom_hline(yintercept=0,linetype="dotted", color="grey60", size=0.5)+
  geom_vline(xintercept=0,linetype="dotted", color="grey60", size=0.5)+
  stat_smooth(method="lm", se=T, col = "grey50",
              data = plotDFLR[plotDFLR$comp == "control / insecticide",])+
  stat_smooth(method="lm", se=F, col = "grey50",lty=2,
              data = plotDFLR[plotDFLR$comp == "control / predator",])+
  facet_wrap(~comp, scales = "free")

#regression coefs.
# rcoefdf <- plotDFLR[plotDFLR$fams == "cumulative", ]
# lmw125 <- lm(lratioH~lratioPL, 
#           data = rcoefdf[rcoefdf$comp == unique(rcoefdf$comp)[1],])
# summary(lmw125)
# 
# lmw25 <- lm(lratioH~lratioPL, 
#              data = rcoefdf[rcoefdf$comp == unique(rcoefdf$comp)[2],])
# summary(lmw25)
# 
# lmI <- lm(lratioH~lratioPL, 
#              data = rcoefdf[rcoefdf$comp == unique(rcoefdf$comp)[3],])
# summary(lmI)
# 
lmP <- lm(lratioH~lratioPL,
             data = rcoefdf[rcoefdf$comp == unique(rcoefdf$comp)[4],])
summary(lmP)



# Herbivore vs IPS
HvIAP <- ggplot(plotDFLR[plotDFLR$fams == "cumulative", ], 
       aes(x = lratioH, y=lratioIP))+
  geom_point()+
  geom_hline(yintercept=0,linetype="dotted", color="grey60", size=1)+
  geom_vline(xintercept=0,linetype="dotted", color="grey60", size=1)+
  stat_smooth(method="lm", se=T, col = "grey50",
              data = plotDFLR[plotDFLR$comp == "control / insecticide",])+
  stat_smooth(method="lm", se=T, col = "grey50",
              data = plotDFLR[plotDFLR$comp == "control / predator",])+
  facet_wrap(~comp, scales = "free")

rcoefdf <- plotDFLR[plotDFLR$fams == "cumulative", ]
# iplmw125 <- lm(lratioIP~lratioPL, 
#              data = rcoefdf[rcoefdf$comp == unique(rcoefdf$comp)[1],])
# summary(iplmw125)
# 
# iplmw25 <- lm(lratioIP~lratioPL, 
#             data = rcoefdf[rcoefdf$comp == unique(rcoefdf$comp)[2],])
# summary(iplmw25)
# 
# iplmI <- lm(lratioIP~lratioPL, 
#           data = rcoefdf[rcoefdf$comp == unique(rcoefdf$comp)[3],])
# summary(iplmI)
# 
# iplmP <- lm(lratioIP~lratioPL, 
#           data = rcoefdf[rcoefdf$comp == unique(rcoefdf$comp)[4],])
# summary(iplmP)

# Herbivore vs PLant FAMS

# I could try to figure out maybe the average ind size for a given family at a given plot to see wether there is a pattern with this relationship

plotDFLR$some_quality <- 1

HvPfam <- ggplot(plotDFLR[plotDFLR$fams != "cumulative", ], 
       aes(y = lratioH, x=lratioPL, col = fams))+
  geom_point(size = plotDFLR[plotDFLR$fams != "cumulative",
                             ]$some_quality)+
  geom_hline(yintercept=0,linetype="dotted", color="grey60", size=1)+
  geom_vline(xintercept=0,linetype="dotted", color="grey60", size=1)+
  stat_smooth(method="lm", se=F)+
  facet_wrap(~comp, scales = "free")

IAPvP <- ggplot(plotDFLR[plotDFLR$fams == "cumulative", ], 
                       aes(y = lratioIP, x=lratioPL))+
  geom_point()+
  geom_hline(yintercept=0,linetype="dotted", color="grey60", size=0.5)+
  geom_vline(xintercept=0,linetype="dotted", color="grey60", size=0.5)+
  stat_smooth(method="lm", se=T, col = "grey50",
              data = plotDFLR[plotDFLR$comp == "control / predator",])+
  facet_wrap(~comp, scales = "free")

# lmIP <- lm(lratioIP~lratioPL,
#           data = plotDFLR[plotDFLR$fams == "cumulative", ])
# summary(lmIP)

# Herbivore vs IPS FAMS
HvIAPfam <- ggplot(plotDFLR[plotDFLR$fams != "cumulative", ], 
                   aes(x = lratioPL, y=lratioIP, col = fams))+
  geom_point()+
  geom_hline(yintercept=0,linetype="dotted", color="grey60", size=1)+
  geom_vline(xintercept=0,linetype="dotted", color="grey60", size=1)+
  stat_smooth(method="lm", se=FALSE)+
  facet_wrap(~comp, scales = "free")

# Arrange panel plot ----

library(ggpubr)

ggarrange(HvP,
          HvPfam+theme(legend.position = "none"),
          IAPvP,
          HvIAPfam+theme(legend.position = "none"), 
          labels = c("A", "B", "C","D"),
          ncol = 2, nrow = 2,
          legend = "bottom", common.legend = T)

# PAirwise correlations of IAPs and Herb
pairdat <- plotDFLR[plotDFLR$fams != "cumulative", ]
herb <- pairdat[!is.na(pairdat$lratioH), c("block",  "fams", "lratioH")]
iaps <- pairdat[!is.na(pairdat$lratioIP), c("block",  "fams", "lratioIP")]

names(iaps) <- names(herb) <- c("block", "ord", "lratio")

pairwise_df <- rbind(herb, iaps)

bm <- matrix("NA", nrow=length(unique(pairwise_df$block)),
       ncol=length(unique(pairwise_df$ord)),
       dimnames = list(unique(pairwise_df$block),
                       unique(pairwise_df$ord)))

pairwise_mat <- data.frame(bm)
orders <- unique(pairwise_df$ord)
for(ord in orders){
  # cord <- orders[ord]
  ss <- pairwise_df[pairwise_df$ord == ord,]
  rownames(ss) <- as.character(ss$block)
  vals <- ss[as.character(unique(pairwise_df$block)), ]$lratio
  pairwise_mat[, ord] <- vals
}

pairs(pairwise_mat)
names(pairwise_mat)
library(corrplot)
pairwise_cor <- cor(pairwise_mat, 
                    method = "pearson", 
                    use = "pairwise.complete.obs")
pairwise_sig <- cor.mtest(pairwise_mat, 
                          conf.level = .95,
                          use = "pairwise.complete.obs")
corrplot(pairwise_cor,
         p.mat = pairwise_sig$p, 
         insig = "label_sig",
         sig.level = c(.001, .01, .05), 
         pch.cex = .9, 
         pch.col = "white",
         type = "upper")

library("PerformanceAnalytics")
chart.Correlation(pairwise_mat, histogram=F, pch=19)
library(GGally)

ggpairs(pairwise_mat,upper=list(continuous="smooth"),lower=list(continuous="smooth"))

####

genllratio <- data.frame()
for (block in unique(biollcp$gard)){
  subbl <- biollcp[biollcp$gard == block, ]
  # for individual block
  for(plt in unique(subbl$plot)){
    # Values for control
    print(plt)
    subplot <- subbl[subbl$plot == plt, ]
 
    crow <- data.frame(bl = block, 
                       pt = plt,
                       trt = unique(subplot$trt),
                       bioH= NA, 
                       bioIP= NA,
                       bioPp = NA)
    
    # for each plant evaluate H, IP and P and then sum them together
    biovec <- c(0,0,0)
    for (plant in unique(subplot$plnm)){
      #print(plant)
      subplant <- subplot[subplot$plnm == plant,]
      ip <- sum(subplant[subplant$nms %in% c("aran","mant"), ]$bio)
      h <- sum(subplant[!(subplant$nms %in% c("aran","mant")), ]$bio)
      p <- subplant$plbio[1]
      #print(c(h,ip,p))
      biovec <- biovec + c(h,ip,p)
    }
    crow$bioH <- biovec[1]
    crow$bioIP <- biovec[2]
    crow$bioPp <- biovec[3]
    genllratio <- rbind(genllratio,
                        crow)
  }
}


genllratio

# Correlations ----
cor_df_cum <- plotDFLR[plotDFLR$fams == "cumulative", ]
for(comp in unique(cor_df_cum$comp)){
  subdset <- cor_df_cum[cor_df_cum$comp == comp, ]
  print(comp)
  
  # Model
  # Simple correlation.
  sub_mod <- lm(lratioPL~lratioH, data = subdset)
  
  # Print the result
  print(summary(sub_mod))
}

# * 2.1.1 Effects of predator removal on herbivores, intermediate predators and plants ----

# library(ggplot2)
# par(mfrow=c(1,1))
# p1 <- ggplot(genllratio, aes(x=trt, y = bioH))+
#   stat_summary(fun = mean, geom = "bar")+
#   stat_summary(fun.data = "mean_cl_normal",
#                geom = "errorbar",
#                width=0.3)+
#   ggtitle("Herbivores")
# 
# p2 <- ggplot(genllratio, aes(x=trt, y = bioIP))+
#   stat_summary(fun = mean, geom = "bar")+
#   stat_summary(fun.data = "mean_cl_normal",
#                geom = "errorbar",
#                width=0.3)+
#   ggtitle("Intermediate predators")
# 
# p3 <- ggplot(genllratio, aes(x=trt, y = bioPp))+
#   stat_summary(fun = mean, geom = "bar")+
#   stat_summary(fun.data = "mean_cl_normal",
#                geom = "errorbar",
#                width=0.3)+
#   ggtitle("Plants")
# 
# library(gridExtra)
# grid.arrange(p1,p2,p3,nrow = 1)

#

# Evaluate the log ratios: CONTROL/PREDATOR
generallr <- data.frame()
for(block in unique(genllratio$bl)){
  sbl <- genllratio[genllratio$bl == block,]
  ctrlplot <- sbl[sbl$trt == "CONTROL",]
  predplot <- sbl[sbl$trt == "PREDATOR",]
  ctrlvals <- ctrlplot[, c("bioH","bioIP","bioPp")]
  predvals <- predplot[, c("bioH","bioIP","bioPp")]
  llratiovals <- log(ctrlvals/predvals)
  llratiorow <- data.frame(bl = block,
                           H = llratiovals$bioH,
                           IP = llratiovals$bioIP,
                           P = llratiovals$bioPp)
  generallr <- rbind(generallr, llratiorow)
}

# Correlation plot
# library(psych)
# pairs.panels(generallr[,c(2,3,4)],
#              method = "pearson", # correlation method
#              hist.col = "#00AFBB",
#              density = TRUE,  # show density plots
#              ellipses = FALSE # show correlation ellipses
# )

# par(mfrow=c(1,3))
# with(generallr, 
#      plot(H,P, pch=19, cex=1.5, main=))
# abline(h=0)
# abline(v=0)
# abline(0,1, lty =2)
# abline(0, -1, lty = 2)
# 
# with(generallr, 
#      plot(IP,P, pch=19, cex=1.5, main=))
# abline(h=0)
# abline(v=0)
# abline(0,1, lty =2)
# abline(0, -1, lty = 2)
# 
# with(generallr, 
#      plot(IP,H, pch=19, cex=1.5, main=))
# abline(h=0)
# abline(v=0)
# abline(0,1, lty =2)
# abline(0, -1, lty = 2)


# * 2.1.2 General log ratio for abundances ----

# Dataset containing biomasses for the log ratio comparisons between predator exclosures and control plots
abullcp <- abufulldf[abufulldf$trt %in% c("CONTROL", "PREDATOR"),]
abullcp$plot <- as.character(abullcp$plot)
abullcp$plnm <- as.character(abullcp$plnm)
abullcp$trt <- as.character(abullcp$trt)
abullcp$gard <- substr(abullcp$plot, 3,4)
# see which species are present in both treatment plots
# table(abullcp$trt, abullcp$plnm) 


genllratio <- data.frame()

for (block in unique(abullcp$gard)){
  subbl <- abullcp[abullcp$gard == block, ]
  # for individual block
  for(plt in unique(subbl$plot)){
    # Values for control
    print(plt)
    subplot <- subbl[subbl$plot == plt, ]
    
    crow <- data.frame(bl = block, 
                       pt = plt,
                       trt = unique(subplot$trt),
                       bioH= NA, 
                       bioIP= NA,
                       bioPp = NA)
    
    # for each plant evaluate H, IP and P and then sum them together
    biovec <- c(0,0,0)
    for (plant in unique(subplot$plnm)){
      #print(plant)
      subplant <- subplot[subplot$plnm == plant,]
      ip <- sum(subplant[subplant$nms %in% c("aran","mant"), ]$bio)
      h <- sum(subplant[!(subplant$nms %in% c("aran","mant")), ]$bio)
      p <- subplant$plbio[1]
      #print(c(h,ip,p))
      biovec <- biovec + c(h,ip,p)
    }
    crow$bioH <- biovec[1]
    crow$bioIP <- biovec[2]
    crow$bioPp <- biovec[3]
    genllratio <- rbind(genllratio,
                        crow)
  }
}


genllratio

# ggplot(genllratio, aes(x=trt, y = bioH))+
#   stat_summary(fun = mean, geom = "bar")+
#   stat_summary(fun.data = "mean_cl_normal",
#                geom = "errorbar",
#                width=0.3)+
#   ggtitle("Herbivores")
# 
# ggplot(genllratio, aes(x=trt, y = bioIP))+
#   stat_summary(fun = mean, geom = "bar")+
#   stat_summary(fun.data = "mean_cl_normal",
#                geom = "errorbar",
#                width=0.3)+
#   ggtitle("Intermediate predators")
# 
# ggplot(genllratio, aes(x=trt, y = bioPp))+
#   stat_summary(fun = mean, geom = "bar")+
#   stat_summary(fun.data = "mean_cl_normal",
#                geom = "errorbar",
#                width=0.3)+
#   ggtitle("Plants")

# library(gridExtra)
# grid.arrange(p1,p2,p3,nrow = 1)

# Evaluate the log ratios: CONTROL/PREDATOR
generallr <- data.frame()
for(block in unique(genllratio$bl)){
  sbl <- genllratio[genllratio$bl == block,]
  ctrlplot <- sbl[sbl$trt == "CONTROL",]
  predplot <- sbl[sbl$trt == "PREDATOR",]
  ctrlvals <- ctrlplot[, c("bioH","bioIP","bioPp")]
  predvals <- predplot[, c("bioH","bioIP","bioPp")]
  llratiovals <- log(ctrlvals/predvals)
  llratiorow <- data.frame(bl = block,
                           H = llratiovals$bioH,
                           IP = llratiovals$bioIP,
                           P = llratiovals$bioPp)
  generallr <- rbind(generallr, llratiorow)
}

# Correlation plot
# library(psych)
# pairs.panels(generallr[,c(2,3,4)], 
#              method = "pearson", # correlation method
#              hist.col = "#00AFBB",
#              density = TRUE,  # show density plots
#              ellipses = FALSE # show correlation ellipses
# )

# par(mfrow=c(1,3))
# with(generallr, 
#      plot(H,P, pch=19, cex=1.5, main=))
# abline(h=0)
# abline(v=0)
# abline(0,1, lty =2)
# abline(0, -1, lty = 2)
# 
# with(generallr, 
#      plot(IP,P, pch=19, cex=1.5, main=))
# abline(h=0)
# abline(v=0)
# abline(0,1, lty =2)
# abline(0, -1, lty = 2)
# 
# with(generallr, 
#      plot(IP,H, pch=19, cex=1.5, main=))
# abline(h=0)
# abline(v=0)
# abline(0,1, lty =2)
# abline(0, -1, lty = 2)

# * 2.1.3 Herbivore families log ratio responses -----

# fam <- "aran"
# bl <- "g1"
datfam <- data.frame()
for(fam in unique(biollcp$nms)){
  print(fam)
  fambiocp <- biollcp[biollcp$nms == fam,]
  blockdat <- data.frame()
  for (bl in unique(biollcp$gard)){
    print(bl)
    fcpbl <- fambiocp[fambiocp$gard == bl,]
    print(fcpbl)
    cbio <- sum(fcpbl[fcpbl$trt %in% c("CONTROL"), ]$bio)
    pbio <- sum(fcpbl[fcpbl$trt %in% c("PREDATOR"), ]$bio)
    blockrow <- data.frame(fam = fam,
                           val = c(cbio,pbio),
                           trt = c("control", "predator"), 
                           block = bl)
    blockdat <- rbind(blockdat, blockrow)
  }
  datfam <- rbind(datfam, blockdat)
}


# * 2.1.4 Differences between orders based on the abundance ----
treatments <- c("control","predator","weevil25", "weevil125")
treatments <- c("control","predator")


abufulldf$gard <- substr(abufulldf$plot, 3,4)
abudatfam <- data.frame()

# fam <- unique(abufulldf$nms)[3]
# bl <- "g2"

for(fam in unique(abufulldf$nms)){
  print(fam)
  fambiocp <- abufulldf[abufulldf$nms == fam,]
  fambiocp$trt <- as.character(fambiocp$trt)
  blockdat <- data.frame()
  for (bl in unique(abufulldf$gard)){
    print(bl)
    fcpbl <- fambiocp[fambiocp$gard == bl,]
    # print(fcpbl)
    fcpbl_trt <- fcpbl[fcpbl$trt %in% toupper(treatments),]
    vals <- tapply(fcpbl_trt$bio, fcpbl_trt$trt, sum)
    blockrow <- data.frame(fam = fam,
                           val = vals,
                           trt = tolower(names(vals)), 
                           block = bl)
    blockdat <- rbind(blockdat, blockrow)
  }
  abudatfam <- rbind(abudatfam, blockdat)
}

abudatfam_nozero <- abudatfam[abudatfam$val != 0, ]

# * 2.1.5 Diffs based on biomass ----

biofulldf$gard <- substr(biofulldf$plot, 3,4)
biodatfam <- data.frame()

fam <- unique(abufulldf$nms)[3]
bl <- "g2"

for(fam in unique(biofulldf$nms)){
  print(fam)
  fambiocp <- biofulldf[biofulldf$nms == fam,]
  fambiocp$trt <- as.character(fambiocp$trt)
  blockdat <- data.frame()
  for (bl in unique(biofulldf$gard)){
    print(bl)
    fcpbl <- fambiocp[fambiocp$gard == bl,]
    # print(fcpbl)
    fcpbl_trt <- fcpbl[fcpbl$trt %in% toupper(treatments),]
    vals <- tapply(fcpbl_trt$bio, fcpbl_trt$trt, sum)
    blockrow <- data.frame(fam = fam,
                           val = vals,
                           trt = tolower(names(vals)), 
                           block = bl)
    blockdat <- rbind(blockdat, blockrow)
  }
  biodatfam <- rbind(biodatfam, blockdat)
}

biodatfam_nozero <- biodatfam[biodatfam$val != 0, ]



# FIG: Figure biomass for insect groups ----

# make possible to include other treatments

library("ggplot2")
datfam_nozero <- datfam[datfam$val != 0, ]
# ggplot(datfam_nozero, aes(x = trt, y=log(val))) +
#   geom_jitter(width = 0.1) +
#   stat_summary(fun = mean, geom = "point", col="red")+
#   stat_summary(fun.data = "mean_cl_normal",
#                geom = "errorbar",
#                width=0.3, col="red") +
#   facet_grid(~fam)

abudatfam_nozero
biodatfam_nozero

ggplot(abudatfam_nozero, aes(x = trt, y=log(val))) +
  geom_jitter(width = 0.1) +
  stat_summary(fun = mean, geom = "point", col="red")+
  stat_summary(fun.data = "mean_cl_normal",
               geom = "errorbar",
               width=0.1, col="red") +
  facet_wrap(vars(fam), scales = "free")

# Statistical tests
ad <- abudatfam_nozero
library(lme4)
library(MASS)
ord <- "orth"
mod <- glmer.nb(val~trt+(1|block), data = ad[ad$fam == ord, ])
modnr <- glm.nb(val~trt, data = ad[ad$fam == ord, ])
summary(mod)
summary(modnr)




# EACH ORDER Statistical test ----
ad <- biodatfam_nozero
library(lme4)
library(MASS)

ord <- "orth"
mod <- glmer(val*100~trt+(1|block), 
             data = ad[ad$fam == ord, ], 
             family = gaussian(link=log))


mod <- nlme::lme(log(val*111)~trt,
                 random = ~1|block, 
             data = ad[ad$fam == ord, ])

# coledat <- ad[ad$fam == ord, ]
# coledat <- coledat[-1,]
# mod <- glmer(val~trt+(1|block),
#              data = coledat,
#              family = gaussian(link=log))

modnr <- glm(val~trt, data = ad[ad$fam == ord, ],family = gaussian(link=log))
summary(mod)   # Paired
summary(modnr) # Unpaired

ccol <- rgb(10,10,10,150, maxColorValue = 255)
csig <- rgb(255,0,0,150,maxColorValue = 255)
cmar <-rgb(255,215,0,150,maxColorValue = 255)
sigcols <- c(ccol, csig, # aran
             ccol, cmar, # cole
             ccol, csig, # hemi
             ccol, ccol, # homo
             ccol, ccol, # lepi
             ccol, ccol, # mant
             ccol, ccol
             )

ggplot(biodatfam_nozero, aes(x = trt, y=log(val))) +
  geom_jitter(width = 0.05,col = rgb(111,111,111,100, maxColorValue = 255)) +
  stat_summary(fun = mean, geom = "point", 
               size = 3,
               col=sigcols)+
  stat_summary(fun.data = "mean_cl_normal",
               geom = "errorbar",
               width=0.05, lwd = 1.5, 
               col=sigcols) +
  facet_wrap(vars(fam), scales = "free")+
  geom_line(aes(x = trt, y=log(val),group = block), lty = 2, 
            col = rgb(111,111,111,100, maxColorValue = 255))

# test for biomass
# plot(log(val)~trt, data=datfam_nozero)
famlme_custom <- nlme::lme(log(val)~trt,random=~1|block,
                     data=datfam_nozero[datfam_nozero$fam == "lepi", ])
summary(famlme_custom)

famlme1 <- nlme::lme(val~trt*fam,random=~1|block,data=datfam_nozero)
summary(famlme1)

# test for abundance ---- 
abudatfam_nozero$famtrt <- paste(abudatfam_nozero$fam,
                                 abudatfam_nozero$trt, sep = "_")
library(MASS)
library(lme4)

art_fam <- "mant"
sel_fam <- abudatfam_nozero[abudatfam_nozero$fam == art_fam, ]

famlme_allpair <- glmer.nb(val~trt+(1|block),
                           data=sel_fam)

famlme_allpair <- nlme::lme(log(val)~trt, random=~1|block,
                           data=sel_fam)


summary(famlme_allpair)
plot(famlme_allpair)
ggplot(sel_fam, aes(y=log(val), x=trt, col = block))+
  geom_jitter()
# ggplot(sel_fam, aes(y=val, x=trt)) + geom_jitter()

# Individual families log ratio plots ----
herbfams <- unique(biollcp$nms)[-grep("aran|mant", unique(biollcp$nms))]
herbfams <- as.character(herbfams)
#family <- herbfams[1]

# lratioFam <-function(family){
#   genllratio <- data.frame()
#   fambiollcp <- biollcp[biollcp$nms == family, ]
#   for (block in unique(fambiollcp$gard)){
#     subbl <- fambiollcp[fambiollcp$gard == block, ]
#     # for individual block
#     for(plt in unique(subbl$plot)){
#       # Values for control
#       print(plt)
#       subplot <- subbl[subbl$plot == plt, ]
#       
#       crow <- data.frame(bl = block, 
#                          pt = plt,
#                          trt = unique(subplot$trt),
#                          bioH= NA, 
#                          bioIP= NA,
#                          bioPp = NA)
#       
#       # for each plant evaluate H, IP and P and then sum them together
#       biovec <- c(0,0,0)
#       for (plant in unique(subplot$plnm)){
#         #print(plant)
#         subplant <- subplot[subplot$plnm == plant,]
#         ip <- sum(subplant[subplant$nms %in% c("aran","mant"), ]$bio)
#         h <- sum(subplant[!(subplant$nms %in% c("aran","mant")), ]$bio)
#         p <- subplant$plbio[1]
#         #print(c(h,ip,p))
#         biovec <- biovec + c(h,ip,p)
#       }
#       crow$bioH <- biovec[1]
#       crow$bioIP <- biovec[2]
#       crow$bioPp <- biovec[3]
#       genllratio <- rbind(genllratio,
#                           crow)
#     }
#   }
#   
#   
#   genllratio
#   
#   # Evaluate the log ratios: CONTROL/PREDATOR
#   generallr <- data.frame()
#   for(block in unique(genllratio$bl)){
#     sbl <- genllratio[genllratio$bl == block,]
#     ctrlplot <- sbl[sbl$trt == "CONTROL",]
#     predplot <- sbl[sbl$trt == "PREDATOR",]
#     ctrlvals <- ctrlplot[, c("bioH","bioIP","bioPp")]
#     predvals <- predplot[, c("bioH","bioIP","bioPp")]
#     if(dim(predvals)[1] == 0 | dim(ctrlvals)[1] == 0){
#       print(paste("No insects from this family in block",block))
#       llratiorow <- data.frame(bl = block,
#                                H = NA,
#                                IP = NA,
#                                P = NA)
#     }else{
#       llratiovals <- log(ctrlvals/predvals)
#       llratiorow <- data.frame(bl = block,
#                                H = llratiovals$bioH,
#                                IP = llratiovals$bioIP,
#                                P = llratiovals$bioPp)
#     }
#     
#     generallr <- rbind(generallr, llratiorow)
#   }
#   
#   # Correlation plot
#   # library(psych)
#   # pairs.panels(generallr[,c(2,4)], 
#   #              method = "pearson", # correlation method
#   #              hist.col = "#00AFBB",
#   #              density = TRUE,  # show density plots
#   #              ellipses = FALSE # show correlation ellipses
#   with(generallr, plot(H,P, pch=19, cex=1.5,main=family))
#   abline(h=0)
#   abline(v=0)
#   abline(0,1, lty =2)
#   abline(0, -1, lty = 2)
#   
#   generallr <- generallr[,c("bl","H","P")]
#   colnames(generallr) <- c("bl",
#                            paste(family, "H", sep=""),
#                            paste(family, "Pl", sep=""))
#   return(generallr)
# }


#********************************************************

# For a single family
lratioFamTrt <-function(family, 
                        treatments, 
                        abundance=T){
  genllratio <- data.frame()
  
  if(abundance){
    biollcp <- abufulldf
  }else{
    biollcp <- biofulldf
  }
  
  fambiollcp <- biollcp[biollcp$nms == family, ]
  for (block in unique(fambiollcp$gard)){
    subbl <- fambiollcp[fambiollcp$gard == block, ]
    # for individual block
    for(plt in unique(subbl$plot)){
      # Values for control
      print(plt)
      subplot <- subbl[subbl$plot == plt, ]
      
      crow <- data.frame(bl = block, 
                         pt = plt,
                         trt = unique(subplot$trt),
                         bioH= NA, 
                         bioIP= NA,
                         bioPp = NA)
      
      # for each plant evaluate H, IP and P and then sum them together
      biovec <- c(0,0,0)
      for (plant in unique(subplot$plnm)){
        #print(plant)
        subplant <- subplot[subplot$plnm == plant,]
        ip <- sum(subplant[subplant$nms %in% c("aran","mant"), ]$bio)
        h <- sum(subplant[!(subplant$nms %in% c("aran","mant")), ]$bio)
        p <- subplant$plbio[1]
        #print(c(h,ip,p))
        biovec <- biovec + c(h,ip,p)
      }
      crow$bioH <- biovec[1]
      crow$bioIP <- biovec[2]
      crow$bioPp <- biovec[3]
      genllratio <- rbind(genllratio,
                          crow)
    }
  }
  
  
  genllratio
  
  # Evaluate the log ratios: CONTROL/PREDATOR
  generallr <- data.frame()
  for(block in unique(genllratio$bl)){
    sbl <- genllratio[genllratio$bl == block,]
    ctrlplot <- sbl[sbl$trt == "CONTROL",]
    predplot <- sbl[sbl$trt == "PREDATOR",]
    ctrlvals <- ctrlplot[, c("bioH","bioIP","bioPp")]
    predvals <- predplot[, c("bioH","bioIP","bioPp")]
    if(dim(predvals)[1] == 0 | dim(ctrlvals)[1] == 0){
      print(paste("No insects from this family in block",block))
      llratiorow <- data.frame(bl = block,
                               H = NA,
                               IP = NA,
                               P = NA)
    }else{
      llratiovals <- log(ctrlvals/predvals)
      llratiorow <- data.frame(bl = block,
                               H = llratiovals$bioH,
                               IP = llratiovals$bioIP,
                               P = llratiovals$bioPp)
    }
    
    generallr <- rbind(generallr, llratiorow)
  }
  
  # Correlation plot
  # library(psych)
  # pairs.panels(generallr[,c(2,4)], 
  #              method = "pearson", # correlation method
  #              hist.col = "#00AFBB",
  #              density = TRUE,  # show density plots
  #              ellipses = FALSE # show correlation ellipses
  with(generallr, plot(H,P, pch=19, cex=1.5,main=family))
  abline(h=0)
  abline(v=0)
  abline(0,1, lty =2)
  abline(0, -1, lty = 2)
  
  generallr <- generallr[,c("bl","H","P")]
  colnames(generallr) <- c("bl",
                           paste(family, "H", sep=""),
                           paste(family, "Pl", sep=""))
  return(generallr)
}


#********************************************************



jpeg("fam_logratio.jpg", width = 600, height = 800)
par(mfrow=c(3,2))

colelr <- lratioFam(herbfams[1])
hemilr <- lratioFam(herbfams[2])
homolr <- lratioFam(herbfams[3])
lepilr <- lratioFam(herbfams[4])
orthlr <- lratioFam(herbfams[5])

dev.off()

# famdat <- Reduce(function(x, y) merge(x, y, by = "bl"), 
#        list(colelr,
#             hemilr,
#             homolr,
#             lepilr,
#             orthlr))

# pairs.panels(famdat[,grep("H",colnames(famdat))], 
#              method = "pearson", # correlation method
#              hist.col = "#00AFBB",
#              density = TRUE,  # show density plots
#              ellipses = FALSE # show correlation ellipses
# )

# * 2.2.4 Herbivore families biomass/abundance responces ----



# 2.2 Herbivore log-ratio for each plant species and each group of insects ----
# within garden

# block = "g1"
# plnt = unique(subbl$plnm)[1]
# fam = unique(subblpl$nms)[1]

logratiodf <- data.frame()
# Go through each block
for(block in unique(biollcp$gard)){
  subbl <- biollcp[biollcp$gard == block, ]
  # within each block obtain community for each plant species
  for(plnt in unique(subbl$plnm)){
    subblpl <- subbl[subbl$plnm == plnt, ]
    subblpl$nms <- as.character(subblpl$nms)
    subpllr <- tapply(subblpl$plbio, subblpl$trt, mean)
    pltlr <- log(subpllr["CONTROL"]/subpllr["PREDATOR"])
    if(is.na(pltlr)){next}
    arthrodf <- data.frame()
    print(plnt)
    for(fam in unique(subblpl$nms)){
      famsub <- subblpl[subblpl$nms == fam, c("bio","trt")]
      cont <- famsub[famsub$trt == "CONTROL",]$bio
      pred <- famsub[famsub$trt == "PREDATOR",]$bio
      famlr <- log(cont/pred)
      if(length(famlr) == 0){next}
      print(fam)
      print(famlr)
      arthrodf<- rbind(arthrodf, data.frame(plnt=plnt,fam=fam, lr=famlr,
                                            pltlr=pltlr, gard=block))
    }
  }
  logratiodf <- rbind(logratiodf, arthrodf)
}

logratiodf
logratiodf <- logratiodf[!(is.nan(logratiodf$lr) | is.infinite(logratiodf$lr)), ]

# Biomasses (from the control plots) for each fam and tree
logratiodf$gard <- as.character(logratiodf$gard)
logratiodf$plnt <- as.character(logratiodf$plnt)
logratiodf$fam <- as.character(logratiodf$fam)

biosize <- data.frame()
for(row in 1:dim(logratiodf)[1]){
  loc <- c(logratiodf[row,]$gard,
           logratiodf[row,]$plnt,
           logratiodf[row,]$fam)
  bcpc <- biollcp[biollcp$trt == "CONTROL",]
  bio <- bcpc[(bcpc$gard == loc[1] & bcpc$plnm == loc[2] & 
               bcpc$nms == loc[3]),]$bio
  biosize <- rbind(biosize, data.frame(garrd = loc[1],
                                       plnt = loc[2],
                                       fam = loc[3], 
                                       bio = bio))
}

biosize
logratiodf$size <- biosize$bio

library(ggplot2)
# p <- ggplot(logratiodf, 
#             aes(x = lr, y = pltlr, label = fam, color=plnt))
# p + coord_cartesian(xlim=c(-5,5), ylim = c(-5,5)) +
#   geom_jitter(size=logratiodf$size, width=0.5) +  
#   geom_text() + 
#   geom_abline(slope = 1, intercept = 0, linetype="dashed") +
#   geom_abline(slope = -1, intercept = 0,linetype="dashed") +
#   xlab("Direct effect of predator removal on insects") + 
#   ylab("Indirect effect of predator removal on plants")
  
# I would like to also indicate the abundance in general of these insects on the plot to show their relative importance for the general result.

ib_slice <- ins_bio[(ins_bio$plot == "w1g1p1" & ins_bio$family == "orth"), 
        c("amount", "bio")]

# plot(ib_slice$amount)
# plot(ib_slice$bio)
# plot(ib_slice$amount ~ log(ib_slice$bio)) #seems normally dist. for log(bio)
# Small sized grashoppers have higher growth rates and higher foraging efforts. This should be the case if there is a high risk of small gh to complete their development before the end of the seasson (no clear seassonality in tropics, food is available all year round for chewers?)

# General agregeated results for the log ratios for various insect groups
# comptrt <- "FUNGICIDE"
comptrt <- "PREDATOR"
# comptrt <- "WEEVIL125"
# comptrt <- "WEEVIL25"
# comptrt <- "INSECTICIDE"

tpcodes <- treats[treats$treat %in% c(comptrt,"CONTROL"), ]

tpcodes$codes <- as.character(tpcodes$codes)
tpcodes$gard <- substr(tpcodes$codes, 3, 4)
gard <- tpcodes$gard[1]

fams <- as.character(unique(insects$family))
# fam <- as.character(fams[1])

# x11(6,6)
# par(mfrow = c(4,2))

# 2.3 Families agregated plants ---- 

llrodf <- data.frame()

for(fam in fams){
  llro <- data.frame()
  for(gard in unique(tpcodes$gard)){
    
    pc <- tpcodes[tpcodes$gard == gard & tpcodes$treat == comptrt, ]$code
    cc <- tpcodes[tpcodes$gard == gard & tpcodes$treat == "CONTROL", ]$code
    
    VCpt <- sum(plants[(plants$CODE == cc & plants$LIFE.FORM %in% c("shrub","tree")),]$WEIGHT, na.rm = T)
    VCmt <- sum(plants[(plants$CODE == pc & plants$LIFE.FORM %in% c("shrub","tree")),]$WEIGHT, na.rm=T)
    
    VCph <- sum(plants[plants$CODE == cc & !(plants$LIFE.FORM %in% c("shrub","tree")),]$WEIGHT, na.rm = T)
    VCmh <- sum(plants[plants$CODE == pc & !(plants$LIFE.FORM %in% c("shrub","tree")),]$WEIGHT, na.rm = T)
    
    # totC <- sum(plants[(plants$CODE == cc),]$WEIGHT)
    # totP <- sum(plants[(plants$CODE == pc),]$WEIGHT)
    # VCph + VCpt # control plot biomass
    # VCmh + VCmt # predator plot biomass
    # 
    
    HCp <- sum(ins_bio[(ins_bio$family == fam & ins_bio$plot == cc), ]$bio, 
               na.rm = T)
    HCm <- sum(ins_bio[(ins_bio$family == fam & ins_bio$plot == pc), ]$bio,
               na.rm=T)
    llro <- rbind(llro, data.frame(fam=fam, gard=gard, lVt=log(VCpt/VCmt),
                                   lVh = log(VCph/VCmh), lH = log(HCp/HCm)))
    
  }
  
  llrodf <- rbind(llrodf, llro)
  
}

llrodf
# pdf("manuscript/figs/llratio.pdf", 6,6)
# llp <- ggplot(llrodf, aes(x = lVh, y=lVt))
# llp + geom_point() +
#   geom_text(label = llrodf$gard)+
#   geom_point(aes(x = lVh, y=lVt, col="red")) + 
#   facet_wrap(llrodf$fam) + 
#   theme_bw() +
#   xlab("Direct effect of predator removal on insects") + 
#   ylab("Indirect effect of predator removal on plants")
# dev.off()
# Are responces related to the diversity of the control plot (background diveristy?
sr <- as.data.frame(tapply(plants$SPEC, plants$CODE, function(x){length(unique(x))}))

bgrdiv <- sr[treats[treats$treat == "CONTROL", ]$codes,]
names(bgrdiv) <- c("rich","code")

# sr$code <- rownames(sr)
# colnames(sr) <- c("sr", "code")
# genvuldf$sr <- sr[genvuldf$plot, "sr"]
# colnames(genvuldf)

# plot(pltlr~lr, data=logratiodf, col=logratiodf$fam, pch=19)
# abline(0,1, lty=2)
# abline(0,-1, lty=2)
# abline(h=0, lty=1)
# abline(v=0, lty=1)

# 
# bipartite::plotweb(subinsct,low.abun = plantb,
#                    high.abun = colSums(subinsct))



# Vulnerability etc... for fully resolved networks
# networklevel(subinsct)

# Motifs
# testGraph = barabasi.game(10, 
#                           m = 5,
#                           power = 2, 
#                           out.pref = TRUE,
#                           zero.appeal = 0.5,
#                           directed = TRUE)
# 
# graph.motifs(testGraph, 
#                size = 3)

# X.X Simple biomass plots ----

# See where these values are coming from??
# break it into groups! mobile and the rest

biofulldf
abufulldf
rownames(treats) <- treats$codes
library(ggplot2)
library(lme4)
library(lmerTest)
library(MASS)

# loop for the biomass
dst <- biofulldf
nms <- as.character(unique(biofulldf$nms))

predators <- c("aran", "mant")
mobile <- c("cole", "hemi", "homo", "orth")
lepidoptera <- c("lepi")

nm <- nms[5]
nm <- predators

for(nm in unique(biofulldf$nms)){
  
  print(nm)
  
  #subset the data
  datast <- dst[dst$nms %in% nm, ]
  
  # summarise the data
  sums <- tapply(datast$bio, datast$plot, sum)
  
  # add treatments
  subtreats <- treats[rownames(sums), ]
  sumdf <- data.frame(vals = sums, 
             trt = subtreats$treat,
             code = subtreats$codes)
  
  # Process subdata
  sumdf$gard <- substr(sumdf$code, 3,4)
  
  # Plot
  ggplot(sumdf, aes(y = vals, x = trt)) + geom_jitter(width = 0.1)
  
  # Model
  lmer1 <- glmer(vals~trt+(1|gard), family=gaussian(link="log"), data=sumdf)
  print(summary(lmer1))
  
}


# loop for the abundance

dst <- abufulldf
nm <- nms[2]
library(blmeco)

predators <- c("aran", "mant")
mobile <- c("cole", "hemi", "homo", "orth")
lepidoptera <- c("lepi")

nm <- lepidoptera

for(nm in unique(biofulldf$nms)){
  print(nm)
  
  #subset the data
  datast <- dst[dst$nms %in% nm, ]
  
  # summarise the data
  sums <- tapply(datast$bio, datast$plot, sum)
  # add treatments
  subtreats <- treats[rownames(sums), ]
  sumdf <- data.frame(vals = sums, 
                      trt = subtreats$treat,
                      code = subtreats$codes)
  
  # Process subdata
  sumdf$gard <- substr(sumdf$code, 3,4)
  
  # Plot
  ggplot(sumdf, aes(y = vals, x = trt)) + geom_jitter(width = 0.1)
  
  # Model
  glmer1 <- glmer(vals~trt+(1|gard), family="poisson", data=sumdf)
  glmer2 <- glmer.nb(vals~trt+(1|gard), data=sumdf)
  summary(glmer1)
  summary(glmer2)
  # Here is the test for overdispersion and with Poisson there is an overdispersion
  dispersion_glmer(glmer1) # should not exceed 1.4
  dispersion_glmer(glmer2)
  
}

#################################

# use all species

dst <- abufulldf
nm <- nms[2]
predators <- c("aran", "mant")
mobile <- c("cole", "hemi", "homo", "orth")
lepidoptera <- c("lepi")

nm <- lepidoptera

for(nm in unique(biofulldf$nms)){
  print(nm)
  
  #subset the data
  datast <- dst[dst$nms %in% nm, ]
  
  # No blocks here
  datast$gard <- substr(datast$plot, 3,4)
  
  # Plot
  ggplot(datast, aes(y = bio, x = trt)) + geom_jitter(width = 0.1)
  
  # Model
  
  glmer1 <- glmer(bio~trt+(1|gard), family="poisson", data=datast)
  glmer2 <- glmer.nb(vals~trt+(1|gard), data=sumdf)
  
  # glmer1 <- glmer(vals~trt+(1|gard), family="poisson", data=datast)
  # glmer2 <- glmer.nb(vals~trt+(1|gard), data=sumdf)
  summary(glmer1)
  summary(glmer2)
  # Here is the test for overdispersion and with Poisson there is an overdispersion
  dispersion_glmer(glmer1) # should not exceed 1.4
  dispersion_glmer(glmer2)
  
}



dst <- biofulldf
nms <- as.character(unique(biofulldf$nms))

predators <- c("aran", "mant")
mobile <- c("cole", "hemi", "homo", "orth")
lepidoptera <- c("lepi")
nm <- nms[5]
nm <- predators
for(nm in unique(biofulldf$nms)){
  print(nm)
  #subset the data
  datast <- dst[dst$nms %in% nm, ]
  # Process subdata
  datast$gard <- substr(datast$plot, 3,4)
  datast <- datast[complete.cases(datast), ]
  datast <- datast[datast$bio != 0,]
  datast$trt <- factor(datast$trt, labels = c("CONTROL",
                                            "FUNGICIDE", 
                                            "WEEVIL125",
                                            "PREDATOR", 
                                            "WEEVIL25", 
                                            "INSECTICIDE"))
  # Plot
  ggplot(datast, aes(y = log(bio), x = trt)) + geom_jitter(width = 0.1)
 
  # Model
  lmer1 <- lmer(log(bio)~trt+(1|gard),data=datast)
  print(summary(lmer1))
}

#################################


#################################

# Number of stems
# abu_rbl <- glmer(stems ~ TREATMENT + (1|GARDEN),
#                  data = tree_test, family = "poisson")
# 
# abu_nb_rbl <- glmer.nb(stems ~ TREATMENT + (1|GARDEN),
#                        data = tree_test)
# 
# 
# summary(abu_rbl)
# summary(abu_nb_rbl)
# 
# library(blmeco)
# # Here is the test for overdispersion and with Poisson there is an overdispersion
# dispersion_glmer(abu_rbl) # should not exceed 1.4
# dispersion_glmer(abu_nb_rbl)

##################################





ds <- tapply(biofulldf$bio, biofulldf$plot, sum)
dfbio <- data.frame(bio = ds, bcodes = rownames(ds), 
                    trt = treats$treat, trtcodes = treats$codes)
dfbio$gard <- substr(dfbio$bcodes, 3, 4)
# bioplt <- ggplot(dfbio, aes(x = trt, y = bio))
# bioplt + geom_jitter(width = 0.1)

library(lme4)
library(lmerTest)
mod1 <- lmer(bio~trt + (1|gard), data=dfbio)
mod1 <- lm(bio~trt, data=dfbio)
summary(mod1)

# 3. Resource limitation exploration (Schmitz 2010, p.31)
# biofulldf