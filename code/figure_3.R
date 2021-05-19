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

# Based on biomass

pbio <- data.frame()
for(fams in unique(ins_bioOrig$family)){
  print(fams)
  for(trtnum in 1:length(treat_pair_list)){
    print(treat_pair_list[[trtnum]])
    pdfR <- groupLogratio(fams, treat_pair_list[[trtnum]], abundance=F)
    pbio <- rbind(pbio, pdfR)
    pdfRG <- groupLogratio(c(), treat_pair_list[[trtnum]], abundance=F)
  }
  
}

# General
for(trtnum in 1:length(treat_pair_list)){
  print(treat_pair_list[[trtnum]])
  pdfR <- groupLogratio(c(), treat_pair_list[[trtnum]], abundance=T)
  plotDataFamilieLogRatio <- rbind(plotDataFamilieLogRatio, pdfR)
}

# General
for(trtnum in 1:length(treat_pair_list)){
  print(treat_pair_list[[trtnum]])
  pdfR <- groupLogratio(c(), treat_pair_list[[trtnum]], abundance=F)
  pbio <- rbind(pbio, pdfR)
}

# Calculate log(ratios)
plotDFLR <- plotDataFamilieLogRatio
plotDFLR$lratioH <- log(plotDFLR$numh/plotDFLR$denh)
plotDFLR$lratioIP <- log(plotDFLR$numip/plotDFLR$denip)
plotDFLR$lratioPL <- log(plotDFLR$numpl/plotDFLR$denpl)

plotDFLR$fams <- as.character(plotDFLR$fams)
plotDFLR[plotDFLR$fams == "", ]$fams <- "cumulative"

# Biomass df
plotDFLRbio <- pbio
plotDFLRbio$lratioH <- log(plotDFLRbio$numh/plotDFLRbio$denh)
plotDFLRbio$lratioIP <- log(plotDFLRbio$numip/plotDFLRbio$denip)
plotDFLRbio$lratioPL <- log(plotDFLRbio$numpl/plotDFLRbio$denpl)

plotDFLRbio$fams <- as.character(plotDFLRbio$fams)
plotDFLRbio[plotDFLRbio$fams == "", ]$fams <- "cumulative"

# Figure 3 panels ----
# PANE 1 ----

# Run this to change abundance into biomass
# plotDFLR <- plotDFLRbio

# Create plots
pane1 <- ggplot(plotDFLR[plotDFLR$fams == "cumulative", ], 
                aes(y = lratioH, x=lratioPL))+
  geom_point()+
  geom_hline(yintercept=0,linetype="dotted", color="grey60", size=0.5)+
  geom_vline(xintercept=0,linetype="dotted", color="grey60", size=0.5)+
  stat_smooth(method="lm", se=T, col = "grey50",
              data = plotDFLR[plotDFLR$comp == "control / insecticide",])+
  stat_smooth(method="lm", se=F, col = "grey50",lty=2,
              data = plotDFLR[plotDFLR$comp == "control / predator",])

#regression coefs.
rcoefdf <- plotDFLR[plotDFLR$fams == "cumulative", ]

pane1test <- lm(lratioH~lratioPL,
                data = rcoefdf)
summary(pane1test)



# PANE 2 ----
pane2 <- ggplot(plotDFLR[plotDFLR$fams == "cumulative", ], 
                aes(x = lratioH, y=lratioIP))+
  geom_point()+
  geom_hline(yintercept=0,linetype="dotted", color="grey60", size=1)+
  geom_vline(xintercept=0,linetype="dotted", color="grey60", size=1)+
  stat_smooth(method="lm", se=T, col = "grey50",
              data = plotDFLR[plotDFLR$comp == "control / insecticide",])+
  stat_smooth(method="lm", se=T, col = "grey50",
              data = plotDFLR[plotDFLR$comp == "control / predator",])

rcoefdf <- plotDFLR[plotDFLR$fams == "cumulative", ]

pane2test <- lm(lratioIP~lratioH,
                data = rcoefdf)
summary(pane2test)

# Herbivore vs PLant FAMS

# I could try to figure out maybe the average ind size for a given family at a given plot to see wether there is a pattern with this relationship

plotDFLR$some_quality <- 1

# PANE 3 ----
pane3 <- ggplot(plotDFLR[plotDFLR$fams != "cumulative", ], 
                aes(y = lratioH, x=lratioPL, col = fams))+
  geom_point(size = 4)+
  geom_hline(yintercept=0,linetype="dotted", color="grey60", size=1)+
  geom_vline(xintercept=0,linetype="dotted", color="grey60", size=1)+
  stat_smooth(method="lm", se=F)+
  facet_wrap(~comp, scales = "free")

rcoefdf <- plotDFLR[plotDFLR$fams != "cumulative", ]
rcoefdf <- rcoefdf[-grep("aran|mant",rcoefdf$fams), ]
rcoefdf <- rcoefdf[complete.cases(rcoefdf$lratioH),]

p3dat <- plotDFLR[plotDFLR$fams != "cumulative", ]

p3dat <- p3dat[-grep("aran|mant",p3dat$fams),]

# Random slope test for orders
pane3test <- nlme::lme(lratioH~lratioPL,
                       random = ~0+lratioPL|fams,
                       data = rcoefdf)
sp3 <- summary(pane3test)

annotation <- paste("P=", 
                    round(sp3$tTable[2,5],3))

pane3 <- ggplot(p3dat, 
                aes(y = lratioH, x=lratioPL))+
  geom_point(aes(pch = fams, color = fams), size = 4)+
  geom_hline(yintercept=0,linetype="dotted", color="grey60", size=1)+
  geom_vline(xintercept=0,linetype="dotted", color="grey60", size=1)+
  stat_smooth(data = p3dat, 
              mapping = aes(y = lratioH, x=lratioPL), 
              method="lm", se=T, col = "gray40")

# pane3

dummy_data <- rcoefdf[, c("lratioH","lratioPL","fams")]
library(psych)
fams_dummy <- dummy.code(dummy_data$fams)
dummy_data <- cbind(dummy_data, fams_dummy)

pane3DummyTest <- lm(lratioH~lratioPL+hemi+lepi+orth+cole,
                     data = dummy_data)
summary(pane3DummyTest)
library(lmSupport)
# lm.sumSquares(pane3DummyTest)

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

# PANE 4 ----
pane4df <- plotDFLR[plotDFLR$fams != "cumulative", ]
pane4df_noIp <- pane4df[!(pane4df$fams %in% c("aran","mant")), 
                        c("block","fams","lratioH")]
pane4df_Ip <- pane4df[(pane4df$fams %in% c("aran","mant")), 
                      c("block","fams","lratioIP")]

basedf<- expand.grid(unique(pane4df_noIp$fams), unique(pane4df_Ip$fams), unique(pane4df_noIp$block))
colnames(basedf) <- c("hfam","ipfam","block")
basedf$lrIP <- NA
basedf$lrH <- NA

for(row in 1:dim(basedf)[1]){
  print(row)
  ipfam <- as.character(basedf[row,]$ipfam)
  hfam <- as.character(basedf[row,]$hfam)
  block <- as.character(basedf[row,]$block)
  
  ipval <- pane4df[(pane4df$fams %in% ipfam & pane4df$block %in% block), ]$lratioIP
  hval <- pane4df[(pane4df$fams %in% hfam & pane4df$block %in% block), ]$lratioH
  
  basedf[row,]$lrH <- hval
  basedf[row,]$lrIP <- ipval
  
}

basedf$ordcomp <- paste(basedf$hfam, basedf$ipfam, sep="_")
library(RColorBrewer)
# display.brewer.all(n=NULL, type="all", select=NULL, exact.n=TRUE, 
#                    colorblindFriendly=T)

pane4 <- ggplot(basedf, 
                aes(x = lrH, y=lrIP, col = ordcomp))+
  geom_point()+
  geom_hline(yintercept=0,linetype="dotted", color="grey60", size=1)+
  geom_vline(xintercept=0,linetype="dotted", color="grey60", size=1)+
  stat_smooth(method="lm", se=F)+
  scale_color_manual(values =brewer.pal(10,"Paired"))+
  facet_wrap(~ipfam, scales = "free")

fullfams <- c("Araneae","Mantodea")
names(fullfams) <- c("aran","mant")

pane4 <- ggplot(basedf, 
                aes(x = lrH, y=lrIP))+
  geom_point(aes(pch = hfam, col = hfam), size = 4)+
  geom_hline(yintercept=0,linetype="dotted", color="grey60", size=1)+
  geom_vline(xintercept=0,linetype="dotted", color="grey60", size=1)+
  facet_wrap(~ipfam, scales = "free", labeller = labeller(ipfam = fullfams))+
  stat_smooth(data = basedf, 
              mapping = aes(x = lrH, y=lrIP), 
              method="lm", se=T, col = "gray40", lty = 1)


# pane4
pane4testA <- lm(lrIP~lrH+hfam,data = basedf[basedf$ipfam == "aran", ])
pane4testM <- lm(lrIP~lrH+hfam,data = basedf[basedf$ipfam == "mant", ])
summary(pane4testA)
summary(pane4testM)

# Run for biomass
# 
# pane4 <- ggplot(basedf[basedf$ipfam == "aran",], 
#                 aes(x = lrH, y=lrIP))+
#   geom_point(aes(pch = hfam, col = hfam), size = 4)+
#   geom_hline(yintercept=0,linetype="dotted", color="grey60", size=1)+
#   geom_vline(xintercept=0,linetype="dotted", color="grey60", size=1)+
#   stat_smooth(data = basedf[basedf$ipfam == "aran",], 
#               mapping = aes(x = lrH, y=lrIP), 
#               method="lm", se=T, col = "gray40", lty = 1)
# 
pane5 <- ggplot(basedf[basedf$ipfam == "mant",],
                         aes(x = lrH, y=lrIP))+
  geom_point(aes(pch = hfam, col = hfam), size = 4)+
  geom_hline(yintercept=0,linetype="dotted", color="grey60", size=1)+
  geom_vline(xintercept=0,linetype="dotted", color="grey60", size=1)
# 
# pane4
# pane5

# Arrange panel plot ----

library(ggpubr)
# library(gridExtra)
# grid.arrange(arrangeGrob(pane1,pane2, pane3, ncol=3, nrow=1),
#              arrangeGrob(pane4, ncol=1, nrow=1), heights=c(4,1), widths=c(2,1))


# Panels again
library(dplyr)

p3dat <- as_tibble(p3dat)

p3datP <- p3dat %>% 
  rename(Order = fams)
p3datP$Order <- dplyr::recode(p3datP$Order, 
                orth = "Orthoptera",
                homo = "Homoptera",
                hemi = "Heteroptera",
                cole = "Coleoptera",
                lepi = "Lepidoptera")

basedfP <- basedf  %>% 
  rename(Order = hfam)
basedfP$Order <- dplyr::recode(basedfP$Order, 
                        orth = "Orthoptera",
                        homo = "Homoptera",
                        hemi = "Heteroptera",
                        cole = "Coleoptera",
                        lepi = "Lepidoptera")
pane3 <- ggplot(p3datP, 
                aes(y = lratioH, x=lratioPL))+
  geom_point(aes(pch = Order, color = Order), size = 4)+
  geom_hline(yintercept=0,linetype="dotted", color="grey60", size=1)+
  geom_vline(xintercept=0,linetype="dotted", color="grey60", size=1)+
  stat_smooth(data = p3datP, 
              mapping = aes(y = lratioH, x=lratioPL), 
              method="lm", se=T, col = "gray40")
pane4 <- ggplot(basedfP, 
                aes(x = lrH, y=lrIP))+
  geom_point(aes(pch = Order, col = Order), size = 4)+
  geom_hline(yintercept=0,linetype="dotted", color="grey60", size=1)+
  geom_vline(xintercept=0,linetype="dotted", color="grey60", size=1)+
  facet_wrap(~ipfam, scales = "free", labeller = labeller(ipfam = fullfams))+
  stat_smooth(data = basedfP, 
              mapping = aes(x = lrH, y=lrIP), 
              method="lm", se=T, col = "gray40", lty = 1)

lH <- c("LRR of herbivores")
lP <- c("LRR of plants")
lIP  <- c("LRR of AP")

labsize <- 10
axissize <- 10

pane1 <- pane1 + 
  theme_bw()+
  theme(axis.text=element_text(size=axissize),
        axis.title=element_text(size=labsize))
pane2 <- pane2 + 
  theme_bw()+
  theme(axis.text=element_text(size=axissize),
        axis.title=element_text(size=labsize))
pane3 <- pane3 + 
  theme_bw()+
  theme(axis.text=element_text(size=axissize),
        axis.title=element_text(size=labsize))
pane4 <- pane4 + 
  theme_bw()+
  theme(axis.text=element_text(size=axissize),
        axis.title=element_text(size=labsize))

# Run for biomass
# pane5 <- pane5 +
#   theme_bw()+
#   theme(axis.text=element_text(size=axissize),
#         axis.title=element_text(size=labsize))

#labels for the panels
# summary(pane1test)
# summary(pane2test)
# summary(pane3test)
# summary(pane4testA)
# summary(pane4testM)

# Arrange the plot
# svg("ms1/draft_3/figures/fig2.svg", width = 8, height=8)

tiff(filename='ms1/draft_4/figures/fig2.tif',
     height=5600,
     width=5200,
     units='px',
     res=800,compression='lzw') 
ggarrange(pane1+ylab(lH)+xlab(lP),
          pane2+ylab(lIP)+xlab(lH),
          pane3+theme(legend.position = "none")+ylab(lH)+xlab(lP),
          pane4+theme(legend.position = "none")+ylab(lIP)+xlab(lH),
          # pane5+theme(legend.position = "none")+ylab(lIP)+xlab(lH),
          labels = c("A", "B", "C","D","E"),
          ncol = 2, nrow = 2,
          legend = "bottom", common.legend = T)
dev.off()
