# RDA

# Prepare data ----
rm(list=ls())
source("code/data_processing_code.R")

# Original datasets in case I need to go back to these
plantsOrig <- contingencyTable2(main_biomass, "CODE", "SP_CODE", "WEIGHT")
abumatOrig <- contingencyTable2(ins_bio,"plot","morphotype","amount")
biomatOrig <- contingencyTable2(ins_bio,"plot","morphotype","totbio")

# ins_bio[ins_bio$plot == "w1g6p5" & ins_bio$morphotype == "lepi013",]
# ins_bio[ins_bio$morphotype == "lepi013",]

# Hellinger transformed datasets
abumat <- decostand(abumatOrig, "hel")
plants <- decostand(plantsOrig, "hel")
biomat <- decostand(biomatOrig, "hel")

# Swith full/only herbivore dataset
# Remove intermediate predators from insect datasets
abumat_noip <- abumat[,-grep("aran|mant", colnames(abumat))]
biomat_noip <- biomat[,-grep("aran|mant", colnames(biomat))]
# # abumat_noip <- abumat
# biomat_noip <- biomat

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

# Select only focal treatments.
# Trimm the data such that it consist only of CONTROL, PREDATOR
# treats_trimmed <- treats_dummy[, c("CONTROL", 
#                                    "PREDATOR",
#                                    "WEEVIL25",
#                                    "WEEVIL125",
#                                    "INSECTICIDE",
#                                    "block")]
# 
# treats_trimmed <- treats_dummy[, c("CONTROL", 
#                                    "PREDATOR",
#                                    "block")]


treats_to_plot <- as.character(unique(treats$treat))[c(3,4,2)]
treats_to_plot <- as.character(unique(treats$treat))[c(3,4,5,2)]

# Predator vs control only
treats_trimmed <- treats[(treats$treat %in% treats_to_plot), ]
treats_trimmed$codes <- as.character(treats_trimmed$codes)
treats_dummy <- dummy.code(as.character(treats_trimmed$treat))
treats_trimmed <- cbind(treats_trimmed, treats_dummy)

# Add predator effect main effect (bird effect analysis)
# treats_to_plot <- as.character(unique(treats$treat))[c(3,4,2,5)]
# treats_trimmed <- treats[(treats$treat %in% treats_to_plot), ]
# treats_trimmed$codes <- as.character(treats_trimmed$codes)
# treats_trimmed$birdef <- "no_bird"
# treats_trimmed[treats_trimmed$treat == "CONTROL",]$birdef <- "bird"

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

treats_trimmed$ipabu <- IPabundance[treats_trimmed$sites]
treats_trimmed$ipbio <- IPbiomass[treats_trimmed$sites]
treats_trimmed$ipdiv <- IPdiversity[treats_trimmed$sites]
treats_trimmed$ipric <- IPrichness[treats_trimmed$sites]

# 1.Plants ordination ----
formula_string <- paste("plants_trimmed~", 
                    paste(paste(treats_to_formula,
                          collapse = "+"), 
                          "Condition(block)", 
                          sep = "+"), 
                    sep="")
rda_formula <- formula(formula_string)
plantTreat <- rda(rda_formula, data = treats_trimmed)

anova(plantTreat, by="terms")
# plot(plantTreat)
plantFit <- envfit(plantTreat, plants_trimmed)

summary(plantTreat)


# This maybe represents the variability well... for meany of these plants, standard deviation is zero.

# 2. Herbivore community ordination ----

# Same formula for herbivores
formula_string <- paste("biomat_trimmed~", 
                        paste(paste(treats_to_formula,
                                    collapse = "+"), 
                              "Condition(block)", 
                              sep = "+"), 
                        sep="")
rda_formula <- formula(formula_string)

herbTreat <- rda(rda_formula, data = treats_trimmed)
anova(herbTreat, by="terms")
# plot(herbTreat)

# See which herbivores respond to treatments
herbFit <- envfit(herbTreat, biomat_trimmed)

# Abundance based ordination for herbivores
formula_string <- paste("abumat_trimmed~", 
                        paste(paste(treats_to_formula,
                                    collapse = "+"), 
                              "Condition(block)", 
                              sep = "+"), 
                        sep="")
rda_formula <- formula(formula_string)
herbTreat <- rda(rda_formula, data = treats_trimmed)
anova(herbTreat, by="terms")
# plot(herbTreat)

# 3a. Removing effect of plant community PLANT pPCA ----

# Bosc suggested conditioning only on block (Appendix S1). 
pPCA <- rda(plants_trimmed ~ Condition(block), data=treats_trimmed)

# Site axes
ort_sites <- pPCA$CA$u
ort_sites <- as.data.frame(ort_sites)

selectOrthogonalVars <- function(maxa, 
                                 ort_sites, 
                                 additional_variables = NULL){
  
  # Orthogonal axes selection on herbivore dataset (no ip)
  # maxa <- 14 # 14 is maximal... why?
  pcs <- paste("PC", seq(1:maxa), sep="")
  pcseq <- paste(pcs, collapse="+")
  a_vars <- paste(additional_variables, collapse="+")
  form <- paste("biomat_trimmed", "~", pcseq,"+",a_vars, sep="")
  rdaform <- with(ort_sites, {as.formula(form)})

  # Combine treatments with PC axes
  treats_pc <- cbind(treats_trimmed, ort_sites[,1:maxa])

  # See which axes of plant variability influence herbivore community based
  # of abundance

  nullRDA <- rda(abumat_trimmed ~ 1 +
                   Condition(block+treat),
                 data = treats_pc)

  scopeRDA <- rda(update.formula(rdaform, . ~ . + Condition(block+treat)),
                  data = treats_pc)

  # Forward secelction
  fselbioins <- ordiR2step(nullRDA, scope = scopeRDA,
                           direction = "forward")

  fselbioins$anova
  pca_plants_formula <- fselbioins$call$formula
  print(pca_plants_formula)
  return(list(treatments = treats_pc,
              selected = fselbioins))
}

fselbioins <- selectOrthogonalVars(maxa = 14, ort_sites)

# With 14 PCs 2 are significant: PC2,PC1 for abundance
# PC1 and PC2 are significant for abundance

# Only PC4 is significant for biomass
# But only if I use first 4 PC's

# See which species were responsible for PC7 and PC10 axis
par(mfrow=c(1,2))
# plot(pPCA, choices = c(1, 2), display = "species")
# plot(pPCA, choices = c(7, 10), display = "species")
par(mfrow=c(1,1))

# 3b. Intermediate Predator PCA ----
# ipppPCA <- rda(ipbio_trimmed ~ Condition(block), data=treats_trimmed)
ipppPCA <- rda(ipabu_trimmed ~ Condition(block), data=treats_trimmed)

# Quick check wether there is bird effect on intermediate preds
# ipRDA <- rda(ipbio_trimmed ~ treat + Condition(block), data=treats_trimmed)
# ipRDA <- rda(ipbio_trimmed ~ birdef + Condition(block), data=treats_trimmed)
# anova(ipRDA)
# anova(ipRDA, by="terms")
# anova(ipRDA, by="axis")

# Site axes
ipPCsites <- ipppPCA$CA$u
ipPCsites <- as.data.frame(ipPCsites)

selectOrthogonalVars(10, ort_sites = ipPCsites,
                     additional_variables = c("ipabu","ipbio","ipdiv","ipric"))

# No effect of IPs PCs, nor general descriptors on herbivore community

# Partial RDA for herbivores abundance ----
form <- fselbioins$selected$call$formula
formula_string <- as.character(form)[3]
splitted_formula <- strsplit(formula_string, split=" ")
orthogonal_axes <- splitted_formula[[1]][grep("PC", splitted_formula[[1]])]
final_herg_abu_formula <- paste("abumat_trimmed", "~", 
                                 "PREDATOR+WEEVIL125+WEEVIL25", "+Condition(block+", paste(orthogonal_axes, collapse = "+"), ")", sep = "")

# Treatment effect on herbivorous communities conditioned on plant composition.
# There was also no effect of IPs community on herbivore communities.

herbAbuConditioned <-rda(formula(final_herg_abu_formula),
                         data=fselbioins$treatments)

anova(herbAbuConditioned, by="terms", permutations = 999)

herbivoresResponses <- envfit(herbAbuConditioned, abumat_trimmed)

sigHerbs <- names(herbivoresResponses$vectors$pvals)[herbivoresResponses$vectors$pvals <= 0.05]
herbivoresResponses$vectors$control


# Plot the graph
plot(herbAbuConditioned, type="n")
points(herbAbuConditioned, display = "sites",
       scaling = 2, pch=19, cex = 1.5,
       col = fselbioins$treatments$treat)

text(herbAbuConditioned, display = "species",
       scaling = 1, select = sigHerbs)

points(herbAbuConditioned, display="cn", scaling = 3,
       col = "red",
       pch = 25)

text(herbAbuConditioned, display="cn", scaling = 3,
     col = "red")


# Pairwise comparisons

# Predator as base
# final_herg_abu_formula <- paste("abumat_trimmed", "~", 
#                                 "CONTROL+WEEVIL125+WEEVIL25", "+Condition(block+", paste(orthogonal_axes, collapse = "+"), ")", sep = "")
# herbAbuConditioned <-rda(formula(final_herg_abu_formula),
#                          data=fselbioins$treatments)
# 
# anova(herbAbuConditioned, by="terms", permutations = 999)
# 
# 
# # W25 as base
# final_herg_abu_formula <- paste("abumat_trimmed", "~", 
#                                 "PREDATOR+WEEVIL125+CONTROL", "+Condition(block+", paste(orthogonal_axes, collapse = "+"), ")", sep = "")
# herbAbuConditioned <-rda(formula(final_herg_abu_formula),
#                          data=fselbioins$treatments)
# 
# anova(herbAbuConditioned, by="terms", permutations = 999)
# 
# # W125 as base
# final_herg_abu_formula <- paste("abumat_trimmed", "~", 
#                                 "CONTROL+PREDATOR+WEEVIL25", "+Condition(block+", paste(orthogonal_axes, collapse = "+"), ")", sep = "")
# herbAbuConditioned <-rda(formula(final_herg_abu_formula),
#                          data=fselbioins$treatments)
# 
# anova(herbAbuConditioned, by="terms", permutations = 9999)
# 
# # All included
# final_herg_abu_formula <- paste("abumat_trimmed", "~", 
#                                 "CONTROL+PREDATOR+WEEVIL25+WEEVIL125", "+Condition(block+", paste(orthogonal_axes, collapse = "+"), ")", sep = "")
# herbAbuConditioned <-rda(formula(final_herg_abu_formula),
#                          data=fselbioins$treatments)
# 
# anova(herbAbuConditioned, by="terms", permutations = 9999)

# BIRD EFFECT ----
biomat_bef <- biomat[treats_trimmed$sites, ]
abumat_bef <- abumat[treats_trimmed$sites, ]


rda_bio <-rda(biomat_bef ~ treat+Condition(block+PC4),data=treats_pc) 
rda_abu <-rda(abumat_bef ~ treat+Condition(block+PC4),data=treats_pc) 

# Bird effect only
rda_biob <-rda(biomat_bef ~ birdef+Condition(block+PC4),data=treats_pc)
rda_abub <-rda(abumat_bef ~ birdef+Condition(block+PC4),data=treats_pc)

rdatoplot <- rda_bio

anova(rdatoplot)
anova(rdatoplot, by="axis")
anova(rdatoplot, by="terms")

RsquareAdj(rdatoplot)

herbnames <- rownames(rdatoplot$CCA$v)
colors <- rep("forestgreen", length(herbnames))
colors[grep("aran|mant", herbnames)] <- "red"

pchar <- colors

plot(rdatoplot, scaling=3, display="species", type="n")
points(rdatoplot, scaling=3, display="species", pch=4, cex = 0.8,
       col = colors)

# text(rda_abu, display = "bp")
# text(rdatoplot, display = "bp", scaling = 3)
# text(rdatoplot, display = "species", scaling = 3)
points(rdatoplot, display = "cn", scaling = 3, col="red")
text(rdatoplot, display = "cn", scaling = 3, col="red")

cntvals <- rdatoplot$CCA$centroids

x1 <- cntvals[2,1]
y1 <- cntvals[2,2]
x2 <- cntvals[3,1]
y2 <- cntvals[3,2]
x3 <- cntvals[4,1]
y3 <- cntvals[4,2]

centroid <- c((x1+x2+x3)/3,(y1+y2+y3)/3)
c_control <- c(cntvals[1,1], cntvals[1,2])
arrows(c_control[1], c_control[2], centroid[1],centroid[2],
       length = 0.1)

# Which insects best fit the rda model (no factors)
colnames(abumat_trimmed)
sp_fitsa <- envfit(rdatoplot, abumat_trimmed)
sig_arth <- sp_fitsa$vectors$pvals
r2_arth <- sp_fitsa$vectors$r
# plot(sort(r2_arth, decreasing = T))
# there are only few significanlty affected herbivores
# r2 of these might be used in weighted regression as these bring more information 
# than nonsignificant ones
names(r2_arth[names(sig_arth[sig_arth <= 0.05] )])
names(r2_arth[r2_arth >= 0.2363755])
# Find projection of points on predator effect axis

# Calculate projection on bird effect axis ----
library(LearnGeom)
predEff <- CreateLinePoints(c_control, centroid)

vec <- c_control
# standardized u
theta <- atan(vec[2]/vec[1]) # in radians

par(mfrow=c(1,2))
herbpoints <- rdatoplot$CCA$v[,c(1,2)]
# plot(herbpoints)
# abline(0,0.5200527)
# abline(h=0, lty=2)
# abline(v=0, lty=2)

Rmat <- matrix(c(cos(theta), -sin(theta), 
                 sin(theta),  cos(theta)), 
               byrow=T, nrow=2)

hp_rot <- herbpoints %*%  Rmat

# plot(hp_rot)
# abline(0,0.5200527)
# abline(h=0, lty=2)
# abline(v=0, lty=2)

projection_on_predef <- hp_rot[,1]

# Predator effect axis vs PDI----
vals <- projection_on_predef
summary(vals)
qts <- quantile(vals, c(0.375,0.625))
# hist(vals,breaks = 50)
# abline(v=qts[1], lty= 2, col=2)
# abline(v=qts[2], lty= 2, col=2)
names(sp_fitsa$vectors$pvals)[sp_fitsa$vectors$pvals <= 0.05]
neutral <- vals[vals<=qts[2] & vals>=qts[1]]
positive <- vals[vals>qts[2]]
negative <- vals[vals<qts[1]]

length(vals)
length(neutral)
length(positive)
length(negative)

source("code/pdi.R")

# highest for cole009

pdimat[, "cole009"]
pdimat[, "lepi073"]
db <- diet_breadth
db <- db[db > 0]
db[db == min(db)]

# only 14 herbivore abundances are significantly  

# Why is the length of PDI different from length of herbivore progjections 
# from pRDA
names(diet_breadth)
dim(pdimat) #data is the same
names(vals)

# subsample from pdi db values
diet_breadth_df <- data.frame(pdi = diet_breadth,
                              bef = vals[names(diet_breadth)])
diet_breadth_df$herbivores <- rownames(diet_breadth_df)
diet_breadth_df$birdef <- "neutral"
diet_breadth_df[diet_breadth_df$herbivores %in% names(positive), ]$birdef <- "positive"
diet_breadth_df[diet_breadth_df$herbivores %in% names(negative), ]$birdef <- "negative"

library(ggplot2)
library(AER)
library(blmeco)
library(emmeans)
library(multcomp)

# Insted of total ammount of food plants o could calculate PDI!!!
dbdf <- diet_breadth_df

dbdf$lpdi <- log(dbdf$pdi)
dbdf <- dbdf[dbdf$pdi != 0, ]
lm1 <- lm(lpdi~birdef, dbdf)
plot(lm1) # pretty good model

inter.test1 <- emmeans(lm1, "birdef")
plot(inter.test1)
phabu <- cld(inter.test1, Letter="abcdefghijklm")

# I am going to log them anyway to better visualize
# ggplot(dbdf, aes(x = birdef, y = log(pdi)))+
#   geom_jitter(width = 0.1, col = rgb(0,0,0,50,maxColorValue = 255))+
#   stat_summary(fun.y = mean, geom = "point", col= "red")+
#   stat_summary(fun.data = "mean_cl_boot", 
#                geom = "errorbar",
#                width=0.05, col="red", lwd=1.1)+
#   ggtitle("H")

# # DIET SWITCHING ----
# # Herbivores present in P and C treatments
# predsites <- treats[treats$treat == "PREDATOR",]$codes
# contsites <- treats[treats$treat == "CONTROL",]$codes
# 
# ips <- grep("aran|mant", ins_bio$morphotype)
# 
# ins_bioOrig <- ins_bioins_bio <- ins_bio[-ips, ]
# 
# pabumat <- contingencyTable2(ins_bio[(ins_bio$plot %in% predsites), ],
#                              "plot","morphotype","amount")
# pbiomat <- contingencyTable2(ins_bio[ins_bio$plot %in% predsites, ],
#                              "plot","morphotype","totbio")
# 
# cabumat <- contingencyTable2(ins_bio[ins_bio$plot %in% contsites, ],
#                              "plot","morphotype","amount")
# cbiomat <- contingencyTable2(ins_bio[ins_bio$plot %in% contsites, ],
#                              "plot","morphotype","totbio")
# 
# comparable <- colnames(cbiomat)[colnames(cbiomat) %in% colnames(pbiomat)]
# 
# # Remove intermediate predators
# comparable
# 
# # Food plants for comparable herbivores
# cp_treats <- treats_trimmed[treats_trimmed$treat %in% c("CONTROL","PREDATOR"),]$sites
# csites <- treats_trimmed[treats_trimmed$treat %in% c("CONTROL"),]$sites
# psites <- treats_trimmed[treats_trimmed$treat %in% c("PREDATOR"),]$sites
# 
# ins_bio_cp <- ins_bio[ins_bio$plot %in% cp_treats, ]
# ins_bio_cp_comparable <- ins_bio_cp[ins_bio_cp$morphotype %in% comparable,]
# ibc <- ins_bio_cp_comparable
# ibc <- ibc[complete.cases(ibc),]
# 
# ccompFood <- contingencyTable2(ibc[ibc$plot %in% csites, ],
#                                "tree",
#                                "morphotype",
#                                "totbio")
# 
# pcompFood <- contingencyTable2(ibc[ibc$plot %in% psites, ],
#                                "tree",
#                                "morphotype",
#                                "totbio")
# 
# dim(pcompFood)
# dim(ccompFood)
# 
# # Combine dataset 
# rownames(pcompFood) <- paste("p", rownames(pcompFood), sep="_")
# rownames(ccompFood) <- paste("c", rownames(ccompFood), sep="_")
# 
# compFood <- rbind(pcompFood,ccompFood)

#### UNHASH IN CASE OF EMERGENCY

envdat <- data.frame(treat = rep(c("predator", "control"), 
              c(dim(pcompFood)[1],
                dim(ccompFood)[1])))
rownames(envdat) <- rownames(compFood)

dietrda <- rda(compFood~treat, data=envdat)
anova(dietrda, by="axis")
plot(dietrda)

# DIet shift initial code
# compare_row <- function(row){
#   compvec <- compFood[,row]
#   pred <- grep("p_", names(compvec))
#   cont <- grep("c_", names(compvec))
#   # Predator treatment will always have more woody plants
#   # However, I need to fix the names.
#   predvec <- compvec[pred]
#   contvec <- compvec[cont]
#   # Remove indicator
#   names(predvec) <- substr(names(predvec), 3, 8)
#   names(contvec) <- substr(names(contvec), 3, 8)
#   # Merge two vectors
#   compDat <- data.frame(pred = predvec)
#   rownames(compDat) <- names(predvec)
#   compDat$cont <- contvec[rownames(compDat)]
#   # Here is another problem. What about palnts thet were, not
#   # present in the control plot but were in the predator exclosure.
#   # If they were not present, but favored would that count as a shift?
#   # I will initially make these values zeros.
#   compDat$cont[is.na(compDat$cont)] <- 0
#   # Caluculate proportions for the diets
#   compDatStand <- as.data.frame(apply(compDat, 2, function(x){x/sum(x)}))
#   # I calculate differences between proportions in utilizing different
#   # plant species. I use absolute values. This way values would vary between
#   # 0 (no shift) and 2 - compelete shift. Therefore I could divide by 2 to 
#   # standardize
#   compDatStand$diff <- abs(compDatStand$pred - compDatStand$cont)
#   shift <- sum(compDatStand$diff)/2
#   return(shift)
# }
# 
# shiftDf <- data.frame()
# for (row in 1:dim(compFood)[2]){
#   print(row)
#   shiftSpec <- data.frame(species = colnames(compFood)[row],
#                      shift = compare_row(row))
#   shiftDf <- rbind(shiftDf, shiftSpec)
# }
# 
# # See if it does what it says
# compFood[, 'orth051'] # yes it does

# inspect individual species
sp <- "orth052"
spibc <- ibc[ibc$morphotype == sp,]
spibc$treat <- treats_trimmed[as.character(spibc$plot),]$treat

# Compare predator effect on food web ----
# Role of species in the food web ----
par(mfrow=c(1,1))
library(bipartite)
main$SP_CODE <- as.character(main$SP_CODE)
biomassC <- tapply(main[main$TREAT == "CONTROL", c("WEIGHT")],
                  main[main$TREAT == "CONTROL",]$SP_CODE,sum)
names(biomassC) <- tolower(names(biomassC))

biomassP <- tapply(main[main$TREAT == "PREDATOR", c("WEIGHT")],
                   main[main$TREAT == "PREDATOR",]$SP_CODE,sum)
names(biomassP) <- tolower(names(biomassP))

# Plot food webs P vs C ----
plotweb(ccompFood, method = "normal", 
        low.abun = biomassC, labsize = 1.1,
        text.rot = 90)
plotweb(pcompFood, method = "normal", low.abun = biomassP)

# C and Z vals - roles of specis ----
mcfw <- computeModules(ccompFood)
czvalsc <- czvalues(mcfw)
# plot(czvalsc$c, czvalsc$z)

mpfw <- computeModules(pcompFood)
czvalsp <- czvalues(mpfw)
grayAlpha <- rgb(0,0,0,10,maxColorValue = 255)

colorsC <- rep(grayAlpha,length(colnames(ccompFood)))
nm <- names(czvalsc$z)[czvalsc$z == max(czvalsc$z)]
nm <- "cole131"
colorsC[colnames(ccompFood) == nm] <- "red"

colorsP <- rep(grayAlpha,length(colnames(pcompFood)))
colorsP[colnames(pcompFood) == nm] <- "red"

par(mfrow=c(1,2))
# plot(czvalsc$c, czvalsc$z, col = colorsC, pch = 19)
abline(h=2.5)
abline(v=0.62)

# plot(czvalsp$c, czvalsp$z, col = colorsP, pch = 19)
abline(h=2.5)
abline(v=0.62)

# Procrustes plot ----
proc  <- procrustes(cbind(czvalsp$c, czvalsp$z),
                    cbind(czvalsc$c, czvalsc$z), 
                    scale = T,
                    symmetric = T)
par(mfrow=c(1,1))
plot(proc)

# Biggest change experienced by:
max(residuals(proc))

# Smallest change
min(residuals(proc))
sort(residuals(proc))

# Diet dissimilarity plot ----
calcDietDissimilarity <- function(name, ...){
  before <- pcompFood[,name]
  after <- ccompFood[,name]
  
  bdf <- data.frame(spc = names(before),
                    val = before,
                    trt = "before")
  adf <- data.frame(spc = names(after),
                    val = after,
                    trt = "after")
  bamat <- contingencyTable2(rbind(bdf,adf), "spc", "trt", "val")
  res <- round(c(vegdist(t(bamat), ...)),3)
  names(res) <- name
  return(res)
}

# create diet dissimilarity database
dietdiss <- data.frame()
for(nms in colnames(ccompFood)){
  entry <- calcDietDissimilarity(nms, "bray")
  dsrow <- data.frame(species = names(entry),
                      braydiss = entry)
  dietdiss <- rbind(dietdiss, dsrow)
}

# Normalised PDI
source("code/pdi.R")
seldb <- diet_breadth[dietdiss$species]

dbdd <- cbind(dietdiss, seldb)
dbdd <- dbdd[dbdd$seldb != 0 , ]
# Hypothesis here is: higher generality higher the change
# plot(seldb~braydiss, data=dbdd, log="y")

# BC PDI PLOT ----
plot(braydiss~seldb, data=dbdd)

# 
# ggplot(dbdd, aes(x=braydiss, y=seldb) ) +
#   stat_density_2d(aes(fill = ..level..), geom = "raster")

ggplot(dbdd, aes(x=braydiss, y=seldb) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE)+
  scale_fill_distiller(palette= "Spectral", direction=1) + 
  geom_point(col = "white")

dbdd$lsdb <- log(dbdd$seldb)
dbdd <- dbdd[dbdd$lsdb != -Inf, ]

summary(glm(lsdb~braydiss, data=dbdd)) # not significant
proc$scale

# par(mfrow=c(1,1))
# plot(proc)
# text(proc)

# How dissimilar are these sites in case of plants?
Cdf <- data.frame(species=names(biomassC),
                  bio = biomassC,
                  site = "control")
Pdf <- data.frame(species=names(biomassC),
                  bio = biomassC,
                  site = "exclosure")

dfQuantity <- function(sel_name, pabumat, cabumat){
    ptrt <- pabumat[,sel_name]
    ptrt <- data.frame(qty = ptrt, site = names(ptrt))
    ptrt$trt <- "exclosure"
    ctrt <- cabumat[,sel_name]
    ctrt <- data.frame(qty = ctrt, site = names(ctrt))
    ctrt$trt <- "control"
    dftrt <- rbind(ptrt,ctrt)
    dftrt$block <- substr(dftrt$site,3,4)
    return(dftrt)
}

library(ggplot2)

# plotComparison <- function(dfdata){
#   ggplot(dfdata, aes(x = trt, y = qty, group = block))+
#     geom_jitter(width=0.0025, cex=4)+
#     geom_line()
# }

dfdata <- dfQuantity(comparable[14],pabumat, cabumat)
ggplot(dfdata, aes(x = trt, y = qty, group = block))+
  geom_jitter(width=0.025, cex=4)+
  geom_line(lty = 2)
library(lme4)

testR <- glmer.nb(qty~trt+(1|block), dfdata)
testN <- glm.nb(qty~trt, dfdata)

# Notes
#Different than other suggestions, I guess you would like a similarity function or index between networks.
# In this case, I can suggest employing a network or graph matching approach, which seeks for solving the isomorphism problem (Np-complete https://en.wikipedia.org/wiki/Graph_isomorphism_problem).
# There are two family methods: exact or error-tolerance network matching.
# In your case, maybe the tolerance to \epsilon error methods is suitable for your problem.
# One of the most simple and easy to use measure for graph similarity is the hamming distance, for instances.
# (remembering that a distance can be approached as a similarity function, where closer distances, more similar graphs)
# https://stackoverflow.com/questions/40919612/why-find-the-hamming-distance-in-dynamical-networks
# And, the advantage is that Hamming or Levenshtein measures are available packages and ready to use in many programming languages, like python and R.
# For more information, some interesting works in the area (not mine) are:
#   Quantification of network structural dissimilarities. (https://www.nature.com/articles/ncomms13928)
# Graph distance for complex networks. (https://www.nature.com/articles/srep34944)
# Graph matching and learning in pattern recognition in the last 10 years. (https://www.worldscientific.com/doi/abs/10.1142/S0218001414500013)
# Bye and success in your work. 