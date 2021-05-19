# Community descriptors of intermediate predators, herbivores and plants.
rm(list=ls())
source("code/data_processing_code.R")

# treats_to_plot <- as.character(unique(treats$treat))[c(6,3,4,5,2)]
treats_to_plot <- as.character(unique(treats$treat))[c(3,4)]


# Indicate insect order here ----
# order <- c("cole")

# ins_bio <- ins_bio[ins_bio$family %in% order, ]

# Herbivores and IPS diversity, abumdnace, richness in different treatments
ins_bio$group <- "Herbivore"
ins_bio$group[grep("aran|mant", ins_bio$family)] <- "Intermediate predator"

# Test
# unique(ins_bio[ins_bio$group == "Intermediate predator", ]$morphotype)

# Abundnaces
herbAbundance <- tapply(ins_bio[ins_bio$group == "Herbivore", ]$amount, 
                        ins_bio[ins_bio$group == "Herbivore", ]$plot, 
                        sum, na.rm=T)
ipsAbundance <- tapply(ins_bio[ins_bio$group == "Intermediate predator", ]$amount, 
                    ins_bio[ins_bio$group == "Intermediate predator", ]$plot, 
                    sum, na.rm=T)

# Biomass

#Why there are so many NAs in herbivore biomass?
ins_bio[is.na(ins_bio$totbio), ]$totbio <- 0
# because there are no estimates of biomass for these species :/

herbBio <- tapply(ins_bio[ins_bio$group == "Herbivore", ]$totbio, 
                  ins_bio[ins_bio$group == "Herbivore", ]$plot, 
                  sum, na.rm=T)

ipsBio <- tapply(ins_bio[ins_bio$group == "Intermediate predator", ]$totbio, 
                  ins_bio[ins_bio$group == "Intermediate predator", ]$plot, 
                  sum, na.rm=T)

# Diversity
herbDivBio <- tapply(ins_bio[ins_bio$group == "Herbivore", ]$totbio, 
                  ins_bio[ins_bio$group == "Herbivore", ]$plot, 
                  vegan::diversity, index = "invsimpson")

ipsDivBio <- tapply(ins_bio[ins_bio$group == "Intermediate predator", ]$totbio, 
                 ins_bio[ins_bio$group == "Intermediate predator", ]$plot, 
                 vegan::diversity,index = "invsimpson")

# Richness [test this!!!]
herbRich <- tapply(ins_bio[ins_bio$group == "Herbivore", ]$morphotype, 
                        ins_bio[ins_bio$group == "Herbivore", ]$plot, 
                        length)
ipsRich <- tapply(ins_bio[ins_bio$group == "Intermediate predator", ]$morphotype, 
                       ins_bio[ins_bio$group == "Intermediate predator", ]$plot, 
                       length)

descriptorDf <- data.frame(plot = names(herbAbundance), 
                           habund = herbAbundance,
                           ipabund = ipsAbundance,
                           hbio = herbBio,
                           ipsbio = ipsBio,
                           hrich = herbRich,
                           iprich = ipsRich,
                           hdivbio = herbDivBio,
                           ipsdivbio = ipsDivBio)
rownames(treats) <- tolower(rownames(treats))
descriptorDf$treat <- treats[descriptorDf$plot, ]$treat

# ggplot dataset

panelData <- data.frame(value = c(herbAbundance,
                                  ipsAbundance,
                                  herbBio,
                                  ipsBio,
                                  herbRich,
                                  ipsRich,
                                  herbDivBio,
                                  ipsDivBio),
                        group = rep(c("Herbivore","AP",
                                    "Herbivore","AP",
                                    "Herbivore","AP",
                                    "Herbivore","AP"),each = 36),
                        descriptor = rep(c("Abundance",
                                           "Biomass",
                                           "Richness",
                                           "Diversity"), each = 2*36),
                        type = rep(c("Herbivore abundance",
                                 "AP abundance",
                                 "Herbivore biomass",
                                 "AP biomass",
                                 "Herbivore richness",
                                 "AP richness",
                                 "Herbivore diversity",
                                 "AP diversity"), each = 36),
                        treat = rep(treats[descriptorDf$plot, ]$treat, 8),
                        plots = rep(treats[descriptorDf$plot, ]$codes, 8))
panelData$block <- substr(panelData$plots,3,4)
panelData$tukey <- "" # set up column for labels from Tukey test

# panelData

# 1. Ready ----
library(lme4)
library(MASS)
library(lmerTest)
library(emmeans)
library(multcomp)
library(AER)


# glmer

### Biomass ----

# Group
herb <- panelData$group == "Herbivore"
ips <- panelData$group == "AP"

# Descriptors
bio <- panelData$descriptor == "Biomass"
abu <- panelData$descriptor == "Abundance"
rich <- panelData$descriptor == "Richness"
div <- panelData$descriptor == "Diversity"

# Selected treatments
trtsel <- panelData$treat %in% treats_to_plot

# Specific function that performs the analysis and appends the Tukey results to a dataset.
analyzeAndAppend2 <- function(Conditions, 
                             FUN,
                             dispTest = F,
                             ...){
  
  # Specific function that analyses the chunk of my data and
  # appends results of a Tukey post hoc test to panelData
  
  #Call the function
  model <- FUN(..., data = panelData[Conditions, ])
  
  print("model calculated")
  if(dispTest){
    print(dispersiontest(model,trafo=1))
  }
  
  # Get model summary and pairwise comparisons
  inter.test <- emmeans(model, "treat", 
                        data = panelData[Conditions, ])
  print("inter.test passed")
  pairwise <- cld(inter.test, Letter="abcdefghijklm")
  
  print("pairwise computed")
  
  # Extract letters and append to the dataset
  ltrs <- data.frame(pt = pairwise$treat,
                     pg = pairwise$.group)
  rownames(ltrs) <- ltrs$pt
  labs <- ltrs[as.character(panelData[Conditions, "treat"]),]$pg
  panelData[Conditions, ]$tukey <<- as.character(labs)
  
  print("data appended")
  
  # Print model results
  print(summary(model))
  
  return(model)
}


# Analyses

# library(merTools)

## 1A. Herbivore - biomass ----
analyzeAndAppend2(Conditions = (herb & bio & trtsel),
                 FUN = glmer,
                 formula = value~treat+(1|block),
                 family = gaussian(link="log"),
                 control=glmerControl(optimizer="bobyqa",
                                      optCtrl=list(maxfun=2e5)))

# Variance test
# data <- panelData[herb & bio & trtsel, ]
# var.test(value~treat, data = data)
# bartlett.test(value~treat, data = data)
# car::leveneTest(value~treat, data = data)

## 1B. IPs - biomass ----
analyzeAndAppend2(Conditions = (ips & bio & trtsel),
                  FUN = glmer,
                  formula = value~treat+(1|block),
                  family = gaussian(link="log"),
                  control=glmerControl(optimizer="bobyqa",
                                       optCtrl=list(maxfun=2e5)))

# data <- panelData[(ips & bio & trtsel), ]
# var.test(value~treat, data = data)
# bartlett.test(value~treat, data = data)
# car::leveneTest(value~treat, data = data)

## 2A. Herbivore - Abundance ----
# no random effects because I was getting singular var-cov matrix
analyzeAndAppend2(
  Conditions = (herb & abu & trtsel),
  FUN = glm.nb,
  formula = value~treat)

# data <- panelData[(herb & abu & trtsel), ]
# var.test(log(value)~treat, data = data)
# bartlett.test(log(value)~treat, data = data)
# car::leveneTest(value~treat, data = data)

## 2B. IPs - Abundnace ----
# no random effects because I was getting singular var-cov matrix
# Barr DJ, Levy R, Scheepers C, Tily HJ. Random effects structure for confirmatory hypothesis testing: Keep it maximal. Journal of Memory and Language, 68(3):255–278, April 2013.
analyzeAndAppend2(
  Conditions = (ips & abu & trtsel),
  FUN = glm.nb,
  formula = value~treat)

# data <- panelData[(ips & abu & trtsel), ]
# var.test(log(value)~treat, data = data)
# bartlett.test(log(value)~treat, data = data)
# car::leveneTest(value~treat, data = data)

## 3A. Herbivore - Richness ----
analyzeAndAppend2(
  Conditions = (herb & rich & trtsel),
  FUN = glmmPQL,
  family = quasipoisson(),
  random = ~ 1|block,
  fixed = value~treat,
  dispTest = F)

# data <- panelData[(herb & rich & trtsel), ]
# var.test(log(value)~treat, data = data)
# bartlett.test(log(value)~treat, data = data)
# car::leveneTest(value~treat, data = data)

## 3B. IPs - Richness ----
# overdispersed random effect models
analyzeAndAppend2(
  Conditions = (ips & rich & trtsel),
  FUN = glmmPQL,
  family = quasipoisson(),
  random = ~ 1|block,
  fixed = value~treat)

# data <- panelData[(ips & rich & trtsel), ]
# var.test(log(value)~treat, data = data)
# bartlett.test(log(value)~treat, data = data)
# car::leveneTest(value~treat, data = data)

## 4A. Herbivore - diversity SW ----
analyzeAndAppend2(Conditions = (herb & div & trtsel),
                  FUN = lmer,
                  formula = value~treat+(1|block))

library(nlme)
analyzeAndAppend2(Conditions = (herb & div & trtsel),
                  FUN = lme,
                  fixed = value~treat,
                  random= ~1|block)


# data <- panelData[(herb & div & trtsel), ]
# var.test(log(value)~treat, data = data)
# bartlett.test(log(value)~treat, data = data)
# car::leveneTest(value~treat, data = data)

## 4B. IPs - diversity SW ----
analyzeAndAppend2(Conditions = (ips & div & trtsel),
                  FUN = lm,
                  formula = value~treat)

# data <- panelData[(ips & div & trtsel), ]
# var.test(log(value)~treat, data = data)
# bartlett.test(log(value)~treat, data = data)
# car::leveneTest(value~treat, data = data)

# 2. Plots ----

# Add density to the panel data
dens <- read.table("datasets/densities.txt")
densPD <- data.frame(value = c(dens$hdens,dens$apdens),
                     group = rep(c("Herbivore", "AP"), 
                                 each = dim(dens)[1]),
                     descriptor = "Density",
                     type = rep(c("Herbivore density",
                                  "AP density"),each = dim(dens)[1]),
                     treat = rep(toupper(dens$trt), 
                                 2),
                     plots = rep(dens$plot, 2),
                     block = rep(dens$block,2),
                     tukey = "")

panelDataDens <- rbind(panelData, densPD)
panelDataDens[panelDataDens$treat == "EXCLOSURE", ]$treat <- "PREDATOR"

panelDataDens$treat
treats_to_plot
library(ggplot2)

# panelData$treat <- factor(panelData$treat,
#                           levels = treats_to_plot)

panelDataDens$treat <- factor(panelDataDens$treat,
                          levels = treats_to_plot)


# ann_text = 

# panelData[panelData$treat %in% treats_to_plot,] -> pD

panelDataDens[panelDataDens$treat %in% treats_to_plot,] -> pD

pD$type <- factor(pD$type, levels = levels(pD$type)[c(5,6,7,8,10,1,2,3,4,9)],
                  labels = c("H-A","H-B","H-D","H-R","H-Dens",
                             "AP-A","AP-B","AP-D","AP-R","AP-Dens"
                           ))

pD[grep("B", pD$type), ]$value <- pD[grep("B", pD$type), ]$value/1000

colvals <-  c("black", "black",
                              "black", "red",
                              "black", "gold",
                              "black", "black",
              "black","red",
                              "black", "black",
                              "black", "red",
                              "black", "black",
                              "black", "black",
              "black","gold")

pDplot <- ggplot(pD, 
       aes(x = treat, 
           y = value,
           group = group,
           label = tukey))+
  geom_jitter(width = 0.05, 
              col = rgb(10,10,10,50,
                        maxColorValue = 255),
              size = 3)+
  geom_line(aes(group  = block), lty = 2, 
            col = rgb(0,0,0,0.1))+
  facet_wrap(~type, ncol=5, scales = "free") + 
  stat_summary(fun.data=mean_cl_boot, 
               geom="pointrange",col = colvals,
               lwd=1) +
  theme_bw() +
  theme(strip.text = element_text(face = "bold", size = 12))+
  scale_x_discrete(labels = c("C", "Ex"))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6))+
  xlab("")+ylab("")

svg("ms1/draft_3/figures/fig1.svg", width = 6, height= 4)
pDplot
dev.off()
# Get real values
ggplot_build(pDplot)$data[[3]]
biotest <- pD[pD$type == "H-B", ]

# mean_cl_boot(biotest$value)
# mean_cl_boot(biotest[biotest$treat == "CONTROL", ]$value)
