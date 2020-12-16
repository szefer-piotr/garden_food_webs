# RDA

# Prepare data ----
rm(list=ls())

library(emmeans)
library(multcomp)
source("code/data_processing_code.R")

# Original datasets in case I need to go back to these
plantsOrig <- contingencyTable2(main_biomass, "CODE", "SP_CODE", "WEIGHT")
abumatOrig <- contingencyTable2(ins_bio,"plot","morphotype","amount")
biomatOrig <- contingencyTable2(ins_bio,"plot","morphotype","totbio")

# Hellinger transformed datasets
abumat <- decostand(abumatOrig, "hel")
plants <- decostand(plantsOrig, "hel")
biomat <- decostand(biomatOrig, "hel")

# Swith full/only herbivore dataset
# Remove intermediate predators from insect datasets
abumat_noip <- abumat[,-grep("aran|mant", colnames(abumat))]
biomat_noip <- biomat[,-grep("aran|mant", colnames(biomat))]


# Use only data for woody trees to remove unsampled herbaceous plants
tree_names <- main_biomass[main_biomass$LIFE.FORM %in% c("shrub","tree"),]$SP_CODE
plants_woody <- plants[,unique(tree_names)]

# Create a dummy variable dataset
rownames(treats) <- treats$codes
treats$block <- substr(treats$codes, 3,4)
library(psych)
treats_dummy <- dummy.code(treats$treat)
rownames(treats_dummy) <- treats$codes
treats_dummy <- as.data.frame(treats_dummy)
treats_dummy$block <- substr(rownames(treats_dummy), 3,4)

treats_to_plot <- as.character(unique(treats$treat))[c(3,4)]

# Predator vs control only
treats_trimmed <- treats[(treats$treat %in% treats_to_plot), ]
treats_trimmed$codes <- as.character(treats_trimmed$codes)
treats_dummy <- dummy.code(as.character(treats_trimmed$treat))
treats_trimmed <- cbind(treats_trimmed, treats_dummy)

# Datasets containing only selected treatments
treats_trimmed$sites <- rownames(treats_trimmed)
abumat_trimmed <- abumat_noip[rownames(abumat_noip) %in% rownames(treats_trimmed), ]
biomat_trimmed <- biomat_noip[rownames(biomat_noip) %in% rownames(treats_trimmed), ]
plants_trimmed <- plants_woody[rownames(plants_woody) %in% rownames(treats_trimmed), ]

# Intermediate predator subset
ipabu_trimmed <- abumat[treats_trimmed$sites,
                        grep("aran|mant", colnames(abumat))]
ipbio_trimmed <- biomat[treats_trimmed$sites,
                        grep("aran|mant", colnames(biomat))]

# treats_trimmed <- treats_trimmed[,!(colnames(treats_trimmed) %in% c("CONTROL"))]
treats_to_formula <- as.character(unique(treats_trimmed$treat))
treats_to_formula <- treats_to_formula[!(treats_to_formula %in% c("CONTROL"))]

# Add general descriptions for IPs to treats_trimmed
IPabundance <- rowSums(abumatOrig[,grep("aran|mant", 
                                        colnames(abumatOrig))])
IPbiomass <- rowSums(biomatOrig[,grep("aran|mant", 
                                      colnames(biomatOrig))])
IPdiversity <- vegan::diversity(abumatOrig[,grep("aran|mant", 
                                                 colnames(abumatOrig))], MARGIN = 1)
IPrichness <-  rowSums(abumatOrig[,grep("aran|mant", 
                                        colnames(abumatOrig))] > 0)

# Herbivore community descriptors
Habundance <- rowSums(abumatOrig[,-grep("aran|mant", 
                                        colnames(abumatOrig))])
Hbiomass <- rowSums(biomatOrig[,-grep("aran|mant", 
                                      colnames(biomatOrig))])
Hdiversity <- vegan::diversity(abumatOrig[,-grep("aran|mant", 
                                                 colnames(abumatOrig))], MARGIN = 1)
Hrichness <-  rowSums(abumatOrig[,-grep("aran|mant", 
                                        colnames(abumatOrig))] > 0)

treats_trimmed$habu <- Habundance[treats_trimmed$sites]
treats_trimmed$hbio <- Hbiomass[treats_trimmed$sites]
treats_trimmed$hdiv <- Hdiversity[treats_trimmed$sites]
treats_trimmed$hric <- Hrichness[treats_trimmed$sites]

treats_trimmed$ipabu <- IPabundance[treats_trimmed$sites]
treats_trimmed$ipbio <- IPbiomass[treats_trimmed$sites]
treats_trimmed$ipdiv <- IPdiversity[treats_trimmed$sites]
treats_trimmed$ipric <- IPrichness[treats_trimmed$sites]

pairs(treats_trimmed[,7:14], 
      bg=as.numeric(treats_trimmed$treat)+1,pch=21)

library(PerformanceAnalytics)
chart.Correlation(treats_trimmed[,7:14], 
                  bg=treats_trimmed$treat,
                  pch=21)


# Plant PCA
pldat <- as.data.frame(plants_trimmed[, colSums(plants_trimmed)!=0])
plPDC <- rda(pldat ~ Condition(block), data=treats_trimmed)
plant_ort_sites <- plPDC$CA$u
# Treatments and plant ortogonal 
treats_trimmedPCA <- cbind(treats_trimmed, plant_ort_sites)


# Test the significance of each variable on the community composition of herbioves and IAPS

# For herbivores withouth the effect of site and exclosure
prdaNulla <- rda(abumat_trimmed~1+Condition(block+treat), 
                data = treats_trimmedPCA)
prdaNullb <- rda(biomat_trimmed~1+Condition(block+treat), 
                data = treats_trimmedPCA)

# Test each PC axis separately
prdaScopea <- rda(abumat_trimmed~.+Condition(block+treat),
                 data = treats_trimmedPCA[, c("block","treat","habu","hbio","hdiv",
                                              "hric","ipabu","ipbio",
                                              "ipdiv","ipric","PC1",
                                              "PC2","PC3","PC4","PC5","PC6")])
prdaScopeb <- rda(biomat_trimmed~.+Condition(block+treat),
                  data = treats_trimmedPCA[, c("block","treat","habu","hbio","hdiv",
                                               "hric","ipabu","ipbio",
                                               "ipdiv","ipric","PC1",
                                               "PC2","PC3","PC4","PC5","PC6")])
ordistep(prdaNulla, prdaScopea)
ordistep(prdaNullb, prdaScopeb)

# For IAPs
prdaNulla <- rda(ipabu_trimmed~1+Condition(block+treat), 
                data = treats_trimmedPCA)
prdaNullb <- rda(ipbio_trimmed~1+Condition(block+treat), 
                data = treats_trimmedPCA)

prdaScopea <- rda(ipabu_trimmed~.+Condition(block+treat),
                  data = treats_trimmedPCA[, c("block","treat","habu","hbio","hdiv",
                                               "hric","ipabu","ipbio",
                                               "ipdiv","ipric","PC1",
                                               "PC2","PC3","PC4","PC5","PC6")])
# habu is significant

prdaScopeb <- rda(ipbio_trimmed~.+Condition(block+treat),
                  data = treats_trimmedPCA[, c("block","treat","habu","hbio","hdiv",
                                               "hric","ipabu","ipbio",
                                               "ipdiv","ipric","PC1",
                                               "PC2","PC3","PC4","PC5","PC6")])

aipform <- ordistep(prdaNulla, prdaScopea)
aipform$call
bipform <- ordistep(prdaNullb, prdaScopeb)
bipform$call
# Full RDA for herbivores and IAPs
abuHerbRDA <- rda(abumat_trimmed~CONTROL+Condition(block), treats_trimmedPCA)
bioHerbRDA <- rda(biomat_trimmed~CONTROL+Condition(block), treats_trimmedPCA)

abuIpRDA <- rda(ipabu_trimmed~CONTROL+Condition(block + habu), treats_trimmedPCA)
bioIpRDA <- rda(ipbio_trimmed~CONTROL+Condition(block), treats_trimmedPCA)

checkPredEffectStrength <- function(abuHerbRDA, treshold = 10){
  # Coord plot
  ordcompdf <- data.frame(rdacoords = (-1) * summary(abuHerbRDA)$species[,1])
  ordcompdf$rdaorders <- substr(rownames(ordcompdf),1,4)
  
  # Remove zeros as these are not informative
  ordcompdf_nz <- ordcompdf[ordcompdf$rdacoords != 0,]
  library(ggplot2)
  ggplot(ordcompdf_nz, aes(x = rdaorders, y=rdacoords))+
    geom_jitter(width = 0.1)
  
  # morphospecies with reasonable abundance
  # only C and P
  ibcp <- ins_bio[ins_bio$plot %in% treats_trimmed$sites, ]
  ibcp$morph <- as.character(ibcp$morph)
  abundances <- tapply(ibcp$amount, ibcp$morph, sum, na.rm = T)
  # treshold = 10
  tspec <- names(abundances[abundances >= treshold])
  rdact <- ordcompdf_nz[rownames(ordcompdf_nz) %in% tspec, ] # tresholded ds
  comp <- lm(rdacoords~rdaorders, rdact, weights = abundances[rownames(rdact)])
  print(summary(comp))
  
  tukeyComp <- emmeans(comp, "rdaorders")
  pairwise <- cld(tukeyComp, Letter="abcdefghijklm")
  
  pairwise <- as.data.frame(pairwise)
  pairwise$.group <- as.character(pairwise$.group)
  
  # For a boxplot as.data.frame(pairwise)
  
  p <- ggplot(rdact, aes(x = rdaorders, y=rdacoords))
  p +  geom_jitter(width = 0.1, col = rgb(10,10,10,80,maxColorValue = 255))+
    stat_summary(fun.data=mean_cl_boot, 
                 geom="pointrange", lwd=0.8) +
    stat_summary(fun=mean, geom="point",cex = 2)+
    geom_text(data = pairwise, aes(x = rdaorders, y = emmean,
                                   label = .group), nudge_y = 0.1)+
    geom_hline(yintercept = 0, lty = 2, 
               col = rgb(10,10,10,80,maxColorValue = 255))
}
library(ggpubr)

checkPredEffectStrength(abuHerbRDA)
checkPredEffectStrength(bioHerbRDA)
checkPredEffectStrength(abuIpRDA)
checkPredEffectStrength(bioIpRDA)

ggarrange(checkPredEffectStrength(abuHerbRDA),
          checkPredEffectStrength(bioHerbRDA),
          checkPredEffectStrength(abuIpRDA),
          checkPredEffectStrength(bioIpRDA),
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)



efvars <- treats_trimmedPCA[, c("block","treat","habu","hbio","hdiv",
                                "hric","ipabu","ipbio",
                                "ipdiv","ipric")]
# Indiviidual species changes with RDA
ahsef <- envfit(abuHerbRDA, abumat_trimmed) # species
ahvef <- envfit(abuHerbRDA, efvars) # variables - nothing

bhsef <- envfit(bioHerbRDA, biomat_trimmed) # species
bhvef <- envfit(bioHerbRDA, efvars) # variables - nothing

anova(abuHerbRDA)
anova(bioHerbRDA)
anova(abuIpRDA)
anova(bioIpRDA)
