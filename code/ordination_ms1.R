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

abufull.trimmed <- abumat[rownames(abumat) %in% rownames(treats_trimmed),]
abufull.trimmed <- abufull.trimmed[, colSums(abufull.trimmed) != 0]
biofull.trimmed <- biomat[rownames(biomat) %in% rownames(treats_trimmed),]
biofull.trimmed <- biofull.trimmed[, colSums(biofull.trimmed) != 0]

# text(fullRDA, display = "species")


# Intermediate predator subset
ipabu_trimmed <- abumat[treats_trimmed$sites,
                        grep("aran|mant", colnames(abumat))]
ipbio_trimmed <- biomat[treats_trimmed$sites,
                        grep("aran|mant", colnames(biomat))]

# treats_trimmed <- treats_trimmed[,!(colnames(treats_trimmed) %in% c("CONTROL"))]
treats_to_formula <- as.character(unique(treats_trimmed$treat))
treats_to_formula <- treats_to_formula[!(treats_to_formula %in% c("CONTROL"))]

# Separate analyses for AP and HERB with plant PCA and The opposite community PCA.




# Add general descriptions for IPs to treats_trimmed ----


# IPabundance <- rowSums(abumatOrig[,grep("aran|mant", 
#                                         colnames(abumatOrig))])
# IPbiomass <- rowSums(biomatOrig[,grep("aran|mant", 
#                                       colnames(biomatOrig))])
# IPdiversity <- vegan::diversity(abumatOrig[,grep("aran|mant", 
#                                                  colnames(abumatOrig))], MARGIN = 1)
# IPrichness <-  rowSums(abumatOrig[,grep("aran|mant", 
#                                         colnames(abumatOrig))] > 0)
# 
# # Herbivore community descriptors
# Habundance <- rowSums(abumatOrig[,-grep("aran|mant", 
#                                         colnames(abumatOrig))])
# Hbiomass <- rowSums(biomatOrig[,-grep("aran|mant", 
#                                       colnames(biomatOrig))])
# Hdiversity <- vegan::diversity(abumatOrig[,-grep("aran|mant", 
#                                                  colnames(abumatOrig))], MARGIN = 1)
# Hrichness <-  rowSums(abumatOrig[,-grep("aran|mant", 
#                                         colnames(abumatOrig))] > 0)
# 
# treats_trimmed$habu <- Habundance[treats_trimmed$sites]
# treats_trimmed$hbio <- Hbiomass[treats_trimmed$sites]
# treats_trimmed$hdiv <- Hdiversity[treats_trimmed$sites]
# treats_trimmed$hric <- Hrichness[treats_trimmed$sites]
# 
# treats_trimmed$ipabu <- IPabundance[treats_trimmed$sites]
# treats_trimmed$ipbio <- IPbiomass[treats_trimmed$sites]
# treats_trimmed$ipdiv <- IPdiversity[treats_trimmed$sites]
# treats_trimmed$ipric <- IPrichness[treats_trimmed$sites]






# Pairwise comparisons of the descriptors ----

# pairs(treats_trimmed[,7:14], 
#       bg=as.numeric(treats_trimmed$treat)+1,pch=21)
# 
# library(PerformanceAnalytics)
# chart.Correlation(treats_trimmed[,7:14], 
#                   pch=treats_trimmed$treat)

# Descriptor correlation plot
# library(ggcorrplot)
# inddat <- treats_trimmed[,7:14]
# inddatcor <- cor(inddat,
#                   method = "pearson", 
#                   use = "pairwise.complete.obs")
# inddatpval <- cor_pmat(inddat)
# 
# ggcorrplot(inddatcor,
#            hc.order = F,
#            type = "upper",
#            p.mat = inddatpval,
#            outline.color = "white",
#            ggtheme = ggplot2::theme_gray)


# Plant PCA ----

pldat <- as.data.frame(plants_trimmed[, colSums(plants_trimmed)!=0])
plPDC <- rda(pldat ~ Condition(block), data=treats_trimmed)
plant_ort_sites <- plPDC$CA$u


# Herbivore PCA main axes of herbivore community change ----

# hdat <- as.data.frame(plants_trimmed[, colSums(plants_trimmed)!=0])
hPCA <- rda(abumat_trimmed ~ Condition(block), data=treats_trimmed)
h_ort_sites <- hPCA$CA$u
colnames(h_ort_sites) <- paste(colnames(h_ort_sites), "h", sep = "")


# AP PCA main axes of AP community change ----

apPCA <- rda(ipabu_trimmed ~ Condition(block), data=treats_trimmed)
ap_ort_sites <- apPCA$CA$u
colnames(ap_ort_sites) <- paste(colnames(ap_ort_sites), "ap", sep = "")

t.t.all.PCA <- cbind(treats_trimmed,
                     plant_ort_sites,
                     h_ort_sites,
                     ap_ort_sites)


# Forward selection for H and AP with plant and opposite community PCA----

# Herbivore abundance
# prda.h.null <- rda(abumat_trimmed~1+Condition(block+treat),
#                    data = t.t.all.PCA)
# prda.h.scope <- rda(abumat_trimmed~ . +Condition(block+treat),
#                     data = t.t.all.PCA[, c("block","treat",
#                                            "PC1","PC2",
#                                            "PC3","PC4",
#                                            "PC5","PC6",
#                                            "PC1ap","PC2ap",
#                                            "PC3ap","PC4ap",
#                                            "PC5ap","PC6ap")])
# ordistep(prda.h.null, prda.h.scope) # no significant variables

# Herbivore biomass
# prda.h.null <- rda(biomat_trimmed~1+Condition(block+treat),
#                    data = t.t.all.PCA)
# prda.h.scope <- rda(biomat_trimmed~ . +Condition(block+treat),
#                     data = t.t.all.PCA[, c("block","treat",
#                                            "PC1","PC2",
#                                            "PC3","PC4",
#                                            "PC5","PC6",
#                                            "PC1ap","PC2ap",
#                                            "PC3ap","PC4ap",
#                                            "PC5ap","PC6ap")])
# ordistep(prda.h.null, prda.h.scope) # no significant variables

# AP abundance
# prda.ap.null <- rda(ipabu_trimmed ~ 1 + Condition(block+treat),
#                    data = t.t.all.PCA)
# prda.ap.scope <- rda(ipabu_trimmed ~ . + Condition(block+treat),
#                     data = t.t.all.PCA[, c("block","treat",
#                                            "PC1","PC2",
#                                            "PC3","PC4",
#                                            "PC5","PC6",
#                                            "PC1h","PC2h",
#                                            "PC3h","PC4h",
#                                            "PC5h","PC6h")])
# ordistep(prda.ap.null, prda.ap.scope) # no significant variables

# AP biomass
# prda.ap.null <- rda(ipbio_trimmed~1+Condition(block+treat),
#                    data = t.t.all.PCA)
# prda.ap.scope <- rda(ipbio_trimmed ~ . +Condition(block+treat),
#                     data = t.t.all.PCA[, c("block","treat",
#                                            "PC1","PC2",
#                                            "PC3","PC4",
#                                            "PC5","PC6",
#                                            "PC1h","PC2h",
#                                            "PC3h","PC4h",
#                                            "PC5h","PC6h")])
# ordistep(prda.ap.null, prda.ap.scope) # no significant variables




# Treatments and plant ortogonal 
# treats_trimmedPCA <- cbind(treats_trimmed, plant_ort_sites)

# Forward selectin of variables for the complete community.
# prdaNulla <- rda(abufull.trimmed~1+Condition(block+treat), 
#                  data = treats_trimmedPCA)
# prdaNullb <- rda(biofull.trimmed~1+Condition(block+treat), 
#                  data = treats_trimmedPCA)

# Test each PC axis separately
# prdaScopea <- rda(abufull.trimmed ~ . + Condition(block+treat),
#                   data = treats_trimmedPCA[, c("block","treat",
#                                                "PC1","PC2",
#                                                "PC3","PC4",
#                                                "PC5","PC6")])
# prdaScopeb <- rda(biofull.trimmed ~ . + Condition(block+treat),
#                   data = treats_trimmedPCA[, c("block","treat",
#                                                "PC1","PC2",
#                                                "PC3","PC4",
#                                                "PC5","PC6")])
# ordistep(prdaNulla, prdaScopea)
# ordistep(prdaNullb, prdaScopeb)



# Full RDA ----
fullRDA <- rda(abufull.trimmed ~ CONTROL + Condition(block), 
               data = treats_trimmed)
summary(fullRDA)
anova(fullRDA, by = "terms", permutations = 999)


fullRDA.PCA <- rda(biofull.trimmed ~ CONTROL + Condition(block+PC1+PC2+PC3+PC4),
    data = t.t.all.PCA[, c("CONTROL", "block","treat",
                           "PC1","PC2",
                           "PC3","PC4")])
plot(fullRDA.PCA)


#Design based permutation
h <- how(blocks = t.t.all.PCA$block, nperm = 9999)



# Design based permutation
anova(fullRDA.PCA, by = "terms", permutations = 999)

# PERMANOVA results show also no effect
# anosim(abufull.trimmed, treats_trimmed$CONTROL)
colors <- ifelse(grepl("aran|mant",
                       colnames(abufull.trimmed)), 
                 rgb(0.8,0.4,0,0.5), 
                 rgb(0,0.6,0.5,0.5))

plot(fullRDA.PCA, display="species", type = "n",
     xlab = "RDA (39.8%)", ylab = "PC (60.2%)")
points(fullRDA.PCA, display = "species", 
       col = colors, pch = 19, cex = 1.2)
shape::Arrows(0,0,-0.5363*0.1,0,lwd = 2, arr.type = "triangle",
              col = "grey20")
text(-0.5363*0.1,0.01,labels = "Effect of predators",
              col = "grey20")

summary(fullRDA.PCA)
fit <- envfit(fullRDA.PCA, t.t.all.PCA[, c("CONTROL")])

plot(fit)
# More complicated models are studied below
# Within order analysis - no significance
ordnm <- "cole"
for (fm in unique(ins_bio$family)){
  dset <- abufull.trimmed[,grepl(fm, colnames(abufull.trimmed))]
  ordrds <- rda(dset~CONTROL+Condition(block), data = treats_trimmed)
  print(anova(ordrds))
  points
}


# For herbivores withouth the effect of site and exclosure
prdaNulla <- rda(abumat_trimmed~1+Condition(block+treat), 
                data = treats_trimmedPCA)
prdaNullb <- rda(biomat_trimmed~1+Condition(block+treat), 
                data = treats_trimmedPCA)

# Test each PC axis separately
prdaScopea <- rda(abumat_trimmed ~ . + Condition(block+treat),
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

prdaScopeb <- rda(ipbio_trimmed~.+Condition(block+treat),
                  data = treats_trimmedPCA[, c("block","treat","habu","hbio","hdiv",
                                               "hric","ipabu","ipbio",
                                               "ipdiv","ipric","PC1",
                                               "PC2","PC3","PC4","PC5","PC6")])

aipform <- ordistep(prdaNulla, prdaScopea)
aipform$call # ipbio now is significant
bipform <- ordistep(prdaNullb, prdaScopeb)
bipform$call # 

# Full RDA for herbivores and IAPs
abuHerbRDA <- rda(abumat_trimmed~CONTROL+Condition(block), treats_trimmedPCA)
bioHerbRDA <- rda(biomat_trimmed~CONTROL+Condition(block), treats_trimmedPCA)

abuIpRDA <- rda(ipabu_trimmed~CONTROL+Condition(block+ipbio), treats_trimmedPCA)
bioIpRDA <- rda(ipbio_trimmed~CONTROL+Condition(block), treats_trimmedPCA)

anova(abuHerbRDA)
anova(bioHerbRDA)
anova(abuIpRDA) # This is marginally significant
anova(bioIpRDA)

# Plot of the effects
library(ggplot2)


checkPredEffectStrength <- function(abuHerbRDA, treshold = 10, ...){
  # Coord plot
  ordcompdf <- data.frame(rdacoords = (-1) * summary(abuHerbRDA)$species[,1])
  ordcompdf$rdaorders <- substr(rownames(ordcompdf),1,4)
  
  # Remove zeros as these are not informative
  ordcompdf_nz <- ordcompdf[ordcompdf$rdacoords != 0,]
  
  # morphospecies with reasonable abundance
  # only C and P
  ibcp <- ins_bio[ins_bio$plot %in% treats_trimmed$sites, ]
  ibcp$morph <- as.character(ibcp$morph)
  abundances <- tapply(ibcp$amount, ibcp$morph, sum, na.rm = T)
  # treshold = 10
  tspec <- names(abundances[abundances >= treshold])
  rdact <- ordcompdf_nz[rownames(ordcompdf_nz) %in% tspec, ] # tresholded ds
  comp <- lm(rdacoords~0+rdaorders, rdact, weights = abundances[rownames(rdact)])
  print(summary(comp))
  
  tukeyComp <- emmeans(comp, "rdaorders")
  pairwise <- cld(tukeyComp, Letter="abcdefghijklm")
  
  pairwise <- as.data.frame(pairwise)
  pairwise$.group <- as.character(pairwise$.group)
  
  # For a boxplot as.data.frame(pairwise)
  
  p <- ggplot(rdact, aes(x = rdaorders, y=rdacoords))
  p +  geom_jitter(width = 0.1, col = rgb(10,10,10,80,maxColorValue = 255))+
    stat_summary(fun.data=mean_cl_boot, 
                 geom="pointrange", lwd=0.8, ...) +
    stat_summary(fun=mean, geom="point",cex = 2, ...)+
    geom_text(data = pairwise, aes(x = rdaorders, y = emmean,
                                   label = .group), nudge_y = 0.1)+
    geom_hline(yintercept = 0, lty = 2, 
               col = rgb(10,10,10,80,maxColorValue = 255))
}
library(ggpubr)

# cols <- c(rgb(),
#           rgb())

p1 <- checkPredEffectStrength(abuHerbRDA, colour = c("red","red","black","red","black"))
p2 <- checkPredEffectStrength(bioHerbRDA, colour = c("red","red","black","red","black"))
p3 <- checkPredEffectStrength(abuIpRDA, colour = c("black","red"))
p4 <- checkPredEffectStrength(bioIpRDA, colour = c("black","gold"))

ggarrange(p1,p2,p3,p4,
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)



efvars <- treats_trimmedPCA[, c("block","treat","habu","hbio","hdiv",
                                "hric","ipabu","ipbio",
                                "ipdiv","ipric")]

# Individual species changes with RDA
ahsef <- envfit(abuHerbRDA, abumat_trimmed) # species
ahvef <- envfit(abuHerbRDA, efvars) # variables - nothing

bhsef <- envfit(bioHerbRDA, biomat_trimmed) # species
bhvef <- envfit(bioHerbRDA, efvars) # variables - nothing

anova(abuHerbRDA)
anova(bioHerbRDA)
anova(abuIpRDA)
anova(bioIpRDA)