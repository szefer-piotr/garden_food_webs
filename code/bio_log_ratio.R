# Data analysis ----

# Direct and indirect effects on plants can be calculated using the log ratio. ln(Vp+/Vp-) where V's are community variables (herbivore abundance and plant biomass in the presence and absence of herbivores). These effect magnitudes can be plotted in relation to each other on an x-y plane and in relation to a 45 deg reference line tha represents the equivalence in strength of direct and indirect effect of carnivores.

# The log ratio effect of carnivores on herbivores should always be negative if carnivores are limiting herbivores regardles of the of the way herbivores are resource limited.

# page 36 Resloving ecosystem complexity

# I need biomass of herbivores, arthropod predators and plant for each plot.

# Inspect that and see if there are any changes in the food web structure. What can i expect?

source("code/data_processing_code.R")

# abufulldf # dataframe
# abugardnets # morphotype based networks for each garden
# abugardnetsfam # family aagregated networks for all garden


## 2. Log response ratios ----

# Dataset containing biomasses for the log ratio comparisons between predator exclosures and control plots
biollcp <- biofulldf[biofulldf$trt %in% c("CONTROL", "PREDATOR"),]
biollcp$plot <- as.character(biollcp$plot)
biollcp$plnm <- as.character(biollcp$plnm)
biollcp$trt <- as.character(biollcp$trt)
biollcp$gard <- substr(biollcp$plot, 3,4)
# see which species are present in both treatment plots
table(biollcp$trt, biollcp$plnm) 

# Assume that each plant hosts unique community of insects.

# Log ratio analyses

# 2.1 General log ratio for herbivores, intermediate predators and plants ----

# For each garden bioHp (plus), bioHm (minus), bioIPp, bioIPm, bioPp, bioPm
genllratio <- data.frame()

for (block in unique(biollcp$gard)){
  subbl <- biollcp[biollcp$gard == block, ]
  # for individual block
  for(plt in unique(subbl$plot)){
    # Values for control
    print(plot)
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

# * 2.1.1 Effects of predator removal on herbivores, intermediate predators and plants ----

library(ggplot2)
par(mfrow=c(1,3))
p1 <- ggplot(genllratio, aes(x=trt, y = bioH))+
  stat_summary(fun.y = mean, geom = "bar")+
  stat_summary(fun.data = "mean_cl_normal", 
               geom = "errorbar",
               width=0.3)+
  ggtitle("Herbivores")
  
p2 <- ggplot(genllratio, aes(x=trt, y = bioIP))+
  stat_summary(fun.y = mean, geom = "bar")+
  stat_summary(fun.data = "mean_cl_normal", 
               geom = "errorbar",
               width=0.3)+
  ggtitle("Intermediate predators")

p3 <- ggplot(genllratio, aes(x=trt, y = bioPp))+
  stat_summary(fun.y = mean, geom = "bar")+
  stat_summary(fun.data = "mean_cl_normal", 
               geom = "errorbar",
               width=0.3)+
  ggtitle("Plants")

library(gridExtra)
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
library(psych)
# pairs.panels(generallr[,c(2,3,4)], 
#              method = "pearson", # correlation method
#              hist.col = "#00AFBB",
#              density = TRUE,  # show density plots
#              ellipses = FALSE # show correlation ellipses
# )

# with(generallr, plot(H,P, pch=19, cex=1.5))
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
table(abullcp$trt, abullcp$plnm) 


genllratio <- data.frame()

for (block in unique(abullcp$gard)){
  subbl <- abullcp[abullcp$gard == block, ]
  # for individual block
  for(plt in unique(subbl$plot)){
    # Values for control
    print(plot)
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

p1 <- ggplot(genllratio, aes(x=trt, y = bioH))+
  stat_summary(fun.y = mean, geom = "bar")+
  stat_summary(fun.data = "mean_cl_normal", 
               geom = "errorbar",
               width=0.3)+
  ggtitle("Herbivores")

p2 <- ggplot(genllratio, aes(x=trt, y = bioIP))+
  stat_summary(fun.y = mean, geom = "bar")+
  stat_summary(fun.data = "mean_cl_normal", 
               geom = "errorbar",
               width=0.3)+
  ggtitle("Intermediate predators")

p3 <- ggplot(genllratio, aes(x=trt, y = bioPp))+
  stat_summary(fun.y = mean, geom = "bar")+
  stat_summary(fun.data = "mean_cl_normal", 
               geom = "errorbar",
               width=0.3)+
  ggtitle("Plants")

library(gridExtra)
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
library(psych)
pairs.panels(generallr[,c(2,3,4)], 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             ellipses = FALSE # show correlation ellipses
)

# with(generallr, plot(H,P, pch=19, cex=1.5))
# abline(h=0)
# abline(v=0)
# abline(0,1, lty =2)
# abline(0, -1, lty = 2)

# * 2.1.3 Herbivore families log ratio responses -----

fam <- "aran"
bl <- "g1"
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

library("ggplot2")
# datfam_nozero <- datfam[datfam$val != 0, ]
# ggplot(datfam_nozero, aes(x = trt, y=log(val))) + 
#   geom_jitter(width = 0.1) +
#   stat_summary(fun = mean, geom = "point", col="red")+
#   stat_summary(fun.data = "mean_cl_normal",
#                geom = "errorbar",
#                width=0.3, col="red") +
#   facet_grid(~fam)

herbfams <- unique(biollcp$nms)[-grep("aran|mant", unique(biollcp$nms))]
herbfams <- as.character(herbfams)
#family <- herbfams[1]
lratioFam <-function(family){
  genllratio <- data.frame()
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

par(mfrow=c(2,3))

colelr <- lratioFam(herbfams[1])
hemilr <- lratioFam(herbfams[2])
homolr <- lratioFam(herbfams[3])
lepilr <- lratioFam(herbfams[4])
orthlr <- lratioFam(herbfams[5])

famdat <- Reduce(function(x, y) merge(x, y, by = "bl"), 
       list(colelr,
            hemilr,
            homolr,
            lepilr,
            orthlr))

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

plot(ib_slice$amount)
plot(ib_slice$bio)
plot(ib_slice$amount ~ log(ib_slice$bio)) #seems normally dist. for log(bio)
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
llp <- ggplot(llrodf, aes(x = lVh, y=lVt))
llp + geom_point() +
  geom_text(label = llrodf$gard)+
  geom_point(aes(x = lVh, y=lVt, col="red")) + 
  facet_wrap(llrodf$fam) + 
  theme_bw() +
  xlab("Direct effect of predator removal on insects") + 
  ylab("Indirect effect of predator removal on plants")
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
bioplt <- ggplot(dfbio, aes(x = trt, y = bio))
bioplt + geom_jitter(width = 0.1)

library(lme4)
library(lmerTest)
mod1 <- lmer(bio~trt + (1|gard), data=dfbio)
mod1 <- lm(bio~trt, data=dfbio)
summary(mod1)

# 3. Resource limitation exploration (Schmitz 2010, p.31)
biofulldf