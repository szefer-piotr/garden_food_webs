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

ggplot(data = fulldf,aes(y=lratio,x=fam,
                         fill = type))+
  stat_summary(
    aes(color = type),
    fun.data="mean_sdl",  fun.args = list(mult=1), 
    geom = "pointrange",  size = 0.4,
    position = position_dodge(0.8)
  )+
  coord_flip() +
  theme_minimal()+
  geom_hline(yintercept = 0, lty=2, 
             col = rgb(150,150,150,150,maxColorValue = 255))


# Tests

for (fms in unique(fulldf$fam)){
  testdf <- fulldf[fulldf$fam == fms, ]
  testlm <- lm(lratio~1+type, data=testdf)
  print(summary(testlm)$coefficients)
}
