rm(list=ls())
source("code/data_processing_code.R")
# source("code/pdi.R")
# source("code/diet_shift.R")


# deps <- tools::package_dependencies("ggplot2", recursive = TRUE)$ggplot2
# for (dep in deps)
#   try(install.packages(dep))



library(ggplot2)
library(dplyr)

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
# treatments <- c("control","predator","weevil25", "weevil125")
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



getLogRatiosCP <- function(dataset){
  
  famlogratio <- data.frame()
  
  predset <- dataset[dataset$trt %in% "predator",]
  contset <- dataset[dataset$trt %in% "control",]
  
  rownames(predset) <- paste(predset$fam,predset$block)
  rownames(contset) <- paste(contset$fam,contset$block)
  
  dsrn <- sort(unique(c(rownames(predset),rownames(contset))))
  
  for(name in dsrn){
    dfrow <- data.frame(name = name, 
                        pvals = predset[name,]$val,
                        cvals = contset[name,]$val)
    famlogratio <- rbind(famlogratio, dfrow)
  }
  
  famlogratio$lratio <- log(famlogratio$cvals/famlogratio$pvals)
  famlogratio$fam <- substr(famlogratio$name,1,4)
  famlogratio$block <- substr(famlogratio$name, 6,7)
  
  return(famlogratio)
}


abulr <- getLogRatiosCP(abudatfam_nozero)
biolr <- getLogRatiosCP(biodatfam_nozero)

abulr$type <- "abundance"
biolr$type <- "biomass"

fulldf <- rbind(abulr,biolr)

ns <- table(fulldf[complete.cases(fulldf), ]$fam, fulldf[complete.cases(fulldf), ]$type)
ns <- ns[,1]

# StatN <- ggproto("StatN", Stat,
#                  required_aes = c("x", "y"), 
#                  compute_group = function(data, scales) {
#                    y <- data$y
#                    y <- y[!is.na(y)]
#                    n <- length(y)
#                    data.frame(x = data$x[1], y = min(y), label = paste0("n=", n))
#                  }
# )
# 
# stat_n <- function(mapping = NULL, data = NULL, geom = "text", 
#                    position = "identity", inherit.aes = TRUE, show.legend = NA, 
#                    na.rm = FALSE, ...) {
#   ggplot2::layer(stat = StatN, mapping = mapping, data = data, geom = geom, 
#                  position = position, inherit.aes = inherit.aes, 
#                  show.legend = show.legend, 
#                  params = list(na.rm = na.rm, ...))
# }

csite <- treats[treats$treat == "CONTROL", ]$codes
psite <- treats[treats$treat == "PREDATOR", ]$codes

abullcp <- ins_bio %>%
  filter(plot %in% c(as.character(csite), 
                     as.character(psite))) %>%
  select(c("family","plot", "amount")) %>%
  mutate(garden = substr(plot,3,4)) %>%
  mutate(name = paste(family, garden)) %>%
  mutate(trt = ifelse(plot %in% csite, "C","Ex"))

# SLA is in cm2/g
mb <- main[main$CODE %in% c(as.character(csite),
                            as.character(psite)), ]
mbw <- mb %>% 
  filter(LIFE.FORM %in% c("tree", "shrub")) %>%
  mutate(area = SLA * LEAVES*1000) # in cm2

mbws <- mbw %>%
  group_by(CODE) %>%
  summarise(area.cm2 = sum(area, na.rm=T)) %>%
  mutate(area.m2 = area.cm2/10000)


tcpp <- treats %>%
  filter(treat %in% c("PREDATOR"))
tcpc <- treats %>%
  filter(treat %in% c("CONTROL"))
# Loop
fulldiv <- data.frame()
descriptors <- c("Diversity", "Richness", "Density")
for(desc in descriptors){
  print(desc)
  for(fam in unique(abullcp$family)){
    print(fam)
    for(gard in unique(abullcp$garden)){
      print(gard)
      subdat <- abullcp %>%
        filter(family == fam & garden == gard)
      
      if(desc == "Diversity"){
        pvals <- vegan::diversity(subdat[subdat$trt == "Ex",]$amount, index = "invsimpson")
        cvals <- vegan::diversity(subdat[subdat$trt == "C",]$amount, index = "invsimpson")
      }
      if(desc == "Richness"){
        pvals <- nrow(subdat[subdat$trt == "Ex",])
        cvals <- nrow(subdat[subdat$trt == "C",])
      }
      
      if(desc == "Density"){
        explot <- tcpp[grepl(gard, tcpp$codes), ]$codes
        cplot <-tcpc[grepl(gard, tcpc$codes), ]$codes
        print("plots assigned")
        carea <- mbws %>% 
          filter(CODE == as.character(cplot)) %>%
          select(area.m2)
        parea <- mbws %>% 
          filter(CODE == as.character(explot)) %>%
          select(area.m2)
        print("area obtained")
        pvals <- sum(subdat[subdat$trt == "Ex",]$amount)/parea
        cvals <- sum(subdat[subdat$trt == "C",]$amount)/carea
        pvals <- as.numeric(pvals)
        cvals <- as.numeric(cvals)
        
        print("vals calculated")
      }
      
      fulldiv <- rbind(fulldiv, data.frame(
        name = subdat$name[1], pvals = pvals, cvals = cvals,
        lratio = log(cvals/pvals), fam = fam, block = gard, 
        type = desc
      ))
    }
  }
}


head(fulldiv)
head(fulldf)
# Connect fulldf with fulldiv

full.desc <- rbind(fulldiv, fulldf)

ord.labs <- c("Orthoptera",
              "Homoptera",
              "Heteroptera",
              "Aranea",
              "Mantodea", 
              "Coleoptera",
              "Lepidoptera")
levels(full.desc$fam) <- ord.labs 

names(full.desc)[7] <- "Descriptor"
levels(full.desc$Descriptor) <- c("Diversity",
                                  "Richness",
                                  "Density",
                                  "Abundance",
                                  "Biomass" )

ggplot(data = full.desc,aes(y=lratio,x=fam,
                            fill = Descriptor))+
  stat_summary(aes(color = Descriptor),
    fun.data="mean_sdl",  fun.args = list(mult=1), 
    geom = "pointrange",  size = 0.4,
    position = position_dodge(0.8)
  )+
  # coord_flip() +
  theme_minimal()+
  geom_hline(yintercept = 0, lty=2, 
             col = rgb(150,150,150,150,maxColorValue = 255))+
  ylab("Log-response ratio")+xlab("")


# Tests
desc.test.df <- data.frame()
for (desc in unique(full.desc$Descriptor)){
  print(desc)
  for (ord in unique(full.desc$fam)){
    print(ord)
    
    testvals <- full.desc %>%
      filter(fam == ord & Descriptor == desc) %>%
      select(lratio)
    
    print(testvals)
    
    testvals <- testvals[complete.cases(testvals),]
    testvals <- testvals[!is.infinite(testvals)]
    
    test <- summary(lm(testvals~1))$coefficients
    
    print(test)
    
    desc.test.df <- rbind(desc.test.df, data.frame(
      Descriptor = desc, Order = ord, p = test
    ))
  }
}
