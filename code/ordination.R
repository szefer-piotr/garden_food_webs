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


treats_to_plot <- as.character(unique(treats$treat))[c(3,4)]
# treats_to_plot <- as.character(unique(treats$treat))[c(6,3,4,5,2)] # insecticide as well

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

# 1. Plants ordination ----
formula_string <- paste("plants_trimmed~", 
                    paste(paste(treats_to_formula,
                          collapse = "+"), 
                          "Condition(block)", 
                          sep = "+"), 
                    sep="")
rda_formula <- formula(formula_string)
plantTreat <- rda(rda_formula, data = treats_trimmed)

# Pairwise
# plantTreat <- rda(plants_trimmed ~ CONTROL + PREDATOR + WEEVIL25 + Condition(block)
# , data = treats_trimmed)
# anova(plantTreat, by="terms")
# 
# plantTreat <- rda(plants_trimmed ~ WEEVIL125 + CONTROL + WEEVIL25 + Condition(block), data = treats_trimmed)
# anova(plantTreat, by="terms")
# 
# plantTreat <- rda(plants_trimmed ~ WEEVIL125 + PREDATOR + CONTROL + Condition(block), data = treats_trimmed)
# anova(plantTreat, by="terms")

anova(plantTreat, by="terms") # as previously we have significant insecticide
# plot(plantTreat)

# plantFit <- envfit(plantTreat, plants_trimmed)

summary(plantTreat)
# This maybe represents the variability well... for meany of these plants, standard deviation is zero.

#### MANUSCRIPT I ----
# 1a. pRDA on CvsP only, no significant effect of plant composition nor IAPs ----
invert_abu <- cbind(ipabu_trimmed,abumat_trimmed)

# Check wether plant abundance explains some variability
pldat <- as.data.frame(plants_trimmed[, colSums(plants_trimmed)!=0])
plPDC <- rda(pldat ~ Condition(block), data=treats_trimmed)

plant_ort_sites <- plPDC$CA$u
treats_trimmedPCA <- cbind(treats_trimmed, plant_ort_sites)

abumat_trimmed_no_cole <- abumat_trimmed[,
                                         -(which(colnames(abumat_trimmed) == "cole001"))]

# For herbivores
prdaNull <- rda(abumat_trimmed_no_cole~1+Condition(block+treat), 
                data = treats_trimmedPCA)
prdaScope <- rda(abumat_trimmed~PC1+PC2+PC3+PC4+Condition(block+treat),
                 data = treats_trimmedPCA)
prdaScope <- rda(abumat_trimmed_no_cole~PC5+PC6+Condition(block+treat), 
                 data = treats_trimmedPCA)
prdaScope <- rda(abumat_trimmed_no_cole~ipabu+ipbio+ipdiv+ipric+Condition(block+treat), 
                 data = treats_trimmedPCA)

# For IAPs
prdaNull <- rda(ipabu_trimmed~1+Condition(block+treat), 
                data = treats_trimmedPCA)

prdaScope <- rda(ipabu_trimmed~PC1+PC2+PC3+PC4+Condition(block+treat),
                 data = treats_trimmedPCA)
prdaScope <- rda(ipabu_trimmed~PC5+PC6+Condition(block+treat), 
                 data = treats_trimmedPCA)
prdaScope <- rda(ipabu_trimmed~habu+hbio+hdiv+hric+Condition(block+treat), 
                 data = treats_trimmedPCA)

ordistep(prdaNull, prdaScope, direction = "forward")
ordiR2step(prdaNull, prdaScope, direction = "forward")

# No PC axis for herbivores was and no characteristics of the IAP community was significant.
# No PC axis for IAPs was and no characteristics of the IAP community was significant.

prdaHerb <- rda(invert_abu~PREDATOR+Condition(block), data = treats_trimmed)

pointcols <- rep(c(rgb(255,0,0,100,maxColorValue = 255),
                   rgb(0,255,0,100,maxColorValue = 255)), 
                 c(dim(ipabu_trimmed)[2], dim(abumat_trimmed)[2]))

# Species significantly responding to the treatment
envHerb <- envfit(prdaHerb, invert_abu)
names(envHerb$vectors$r)[envHerb$vectors$pvals <= 0.05]

# only melanolepis multigmlandulosa!

plot(prdaHerb, display = "species")
plot(prdaHerb, type = "n",display="species")
points(prdaHerb, display = "species", col = pointcols, pch = 19)
plot(prdaHerb, display = "bp")
anova(prdaHerb, by = "terms")

# plot(envHerb[])
# And the treatment was also not significant.


# EXAMPLE ----
library(vegan)
library(ggfortify)

data(varespec)
data(varechem)

#CCA
cca_model<-cca(varespec ~ .,data=varechem)
plot(cca_model,choices=c(1,2), display=c('sp','bp'), scaling=2)

#Get CCA scores
df_species  <- data.frame(summary(cca_model)$species[,1:2])# get the species CC1 and CC2 scores
df_environ  <- scores(cca_model, display = 'bp') #get the environment vars CC1 and CC2 scores

cca1_varex<-round(summary(cca_model)$cont$importance[2,1]*100,2) #Get percentage of variance explained by first axis
cca2_varex<-round(summary(cca_model)$cont$importance[2,2]*100,2) #Get percentage of variance explained by second axis

#Set a scaling variable to multiply the CCA values, in order to get a very similar plot to the the one generated by plot(cca_model). You can adjust it according to your data
scaling_factor <- 2

ggplot(df_species, 
       aes(x=CCA1, y=CCA2)) + 
  #Draw lines on x = 0 and y = 0
  geom_hline(yintercept=0, 
             linetype="dashed") +
  geom_vline(xintercept=0, 
             linetype="dashed") +
  coord_fixed()+
  #Add species text
  geom_text(data=df_species, 
            aes(x=CCA1,#Score in CCA1 to add species text
                y=CCA2,#Score in CCA2 to add species text
                label=rownames(df_species),
                hjust=0.5*(1-sign(CCA1)),#Set the text horizontal alignment according to its position in the CCA plot
                vjust=0.5*(1-sign(CCA2))),#Set the text vertical alignment according to its position in the CCA plot
            color = "forestgreen")+
  #Add environmental vars arrows
  geom_segment(data=df_environ, 
               aes(x=0, #Starting coordinate in CCA1 = 0 
                   xend=CCA1*scaling_factor,#Ending coordinate in CCA1  
                   y=0, #Start in CCA2 = 0
                   yend=CCA2*scaling_factor), #Ending coordinate in CCA2 
               color="firebrick1", #set color
               arrow=arrow(length=unit(0.01,"npc"))#Set the size of the lines that form the tip of the arrow
  )+
  #Add environmental vars text
  geom_text(data=df_environ, 
            aes(x=CCA1*scaling_factor, 
                y=CCA2*scaling_factor,
                label=rownames(df_environ),
                hjust=0.5*(1-sign(CCA1)),#Add the text of each environmental var at the end of the arrow
                vjust=0.5*(1-sign(CCA2))),#Add the text of each environmental var at the end of the arrow 
            color="firebrick1")+
  #Set bw theme
  theme_bw()+
  #Set x and y axis titles
  labs(x=paste0("CCA1 (",cca1_varex," %)"),
       y=paste0("CCA2 (",cca2_varex," %)"))

# END OF THE EXAMPLE ----

# invert_abu[, "lepi008"]
# text(prdaHerb,display = "species")

insmetamds <- metaMDS(invert_abu)
env <- envfit(insmetamds, treats_trimmed)

dim(env$vectors$arrows)

# No differences in community composition
anosim(invert_abu, treats_trimmed$PREDATOR)

rownames(invert_abu)

plot(insmetamds)
plot(env)
arrows(x0 <- 0,
       y0 <- 0,
       x1 = env$vectors$arrows[env$vectors$pvals >= 0.05, 1],
       y1 = env$vectors$arrows[env$vectors$pvals >= 0.05, 2]
       )

# There were no differences in abundances of woody plants

for(col in 1:dim(pldat)[2]){
  
  x <- pldat[which(treats_trimmed$PREDATOR==1), col]
  y <- pldat[which(treats_trimmed$CONTROL==1), col]
  
  if(sum(x) == 0 | sum(y) == 0 ){
    # print("Zeros")
    next
  }
  
  if(sum(x > 0) <= 3 | sum(y > 0) <= 3){
    # print("Only one in each")
    next
  }
  print(colnames(pldat)[col])
  print(cbind(x,y))
  print(t.test(x,y, paired = T)$p.value)
  print("W T")
  print(wilcox.test(x, y, paired = TRUE, alternative = "two.sided")$p.value)
  
}

t.test(x,y, paired = T)$p.value

#### END


# 1b. Removing effect of herbivores on plants ----

# Below migh be unnecessary as we don't anticipate that herbivores controll plant community
abumat_trimmed_no_cole <- abumat_trimmed[, colnames(abumat_trimmed) != "cole001"]
herbpPCA <- rda(abumat_trimmed_no_cole ~ Condition(block), data=treats_trimmed)

# Site axes
herb_ort_sites <- herbpPCA$CA$u
herb_ort_sites <- as.data.frame(herb_ort_sites)

#NOT tested yet
# focal_dataset <- abumat_trimmed
# maxa <- 14
# additional_variables = NULL

selectOrthogonalVars <- function(maxa,
                                 focal_dataset,
                                 ort_sites, 
                                 additional_variables = NULL){
  
  # focal data set <- set to which ort_sites will be fitted
  
  # Orthogonal axes selection on herbivore dataset (no ip)
  # maxa <- 14 # 14 is maximal... why?
  pcs <- paste("PC", seq(1:maxa), sep="")
  pcseq <- paste(pcs, collapse="+")
  a_vars <- paste(additional_variables, collapse="+")
  
  if(is.null(additional_variables)){
    form <- paste(deparse(quote(focal_dataset)), 
                  "~", pcseq,
                  sep="")
  }else{
    form <- paste(deparse(quote(focal_dataset)), 
                  "~", pcseq,
                  "+",a_vars, 
                  sep="")
  }
  
  rdaform <- with(ort_sites, {as.formula(form)})
  
  # Combine treatments with PC axes
  treats_pc <- cbind(treats_trimmed, ort_sites[,1:maxa])
  
  # See which axes of plant variability influence herbivore community based
  # of abundance
  
  nullRDA <- rda(focal_dataset ~ 1 +
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

plants_sel_form <- selectOrthogonalVars(maxa = 1,
                                   focal_dataset =  plants_trimmed,
                                   ort_sites = herb_ort_sites,
                                   additional_variables = c("hbio",
                                                            "hdiv",
                                                            "habu",
                                                            "hric"))
# Herbivores PC2, PC3, PC1 and PC8
par(mfrow=c(1,2))
plot(herbpPCA, choices = c(1,2), display="species")
text(herbpPCA, choices = c(1,2), display="species")

plot(herbpPCA, choices = c(3,8), display="species")
text(herbpPCA, choices = c(3,8), display="species")

plot(herbpPCA, display="species")
text(herbpPCA, display="species")

hform <- plants_sel_form$selected$call$formula
hformula_string <- as.character(hform)[3]
hsplitted_formula <- strsplit(hformula_string, split=" ")
horthogonal_axes <- hsplitted_formula[[1]][grep("PC", hsplitted_formula[[1]])]
hfinal_herg_abu_formula <- paste("plants_trimmed", "~",
                                "PREDATOR+WEEVIL125+WEEVIL25+INSECTICIDE", 
                                "+Condition(block+", 
                                paste(horthogonal_axes, 
                                      collapse = "+"), ")", sep = "")

# # Treatment effect on herbivorous communities conditioned on plant composition.
# # There was also no effect of IPs community on herbivore communities.
# 
plantAbuConditioned <-rda(formula(hfinal_herg_abu_formula),
                         data=plants_sel_form$treatments)
anova(plantAbuConditioned, by="terms", permutations = 999)

# if we include Insecticide we obtain significance of w125

# 2. Herbivore community ordination ----
# Same formula for herbivores

biomat_trimmed_no_cole <- biomat_trimmed[, colnames(biomat_trimmed) != "cole001"]

formula_string <- paste("biomat_trimmed_no_cole~", 
                        paste(paste(treats_to_formula,
                                    collapse = "+"), 
                              "Condition(block)", 
                              sep = "+"), 
                        sep="")
rda_formula <- formula(formula_string)

herbTreatBio <- rda(rda_formula, data = treats_trimmed)
anova(herbTreatBio, by="terms")
# withouth conditioning, and withouth cole001
# we obtain significant w125 on herbivore biomass.

# See which herbivores respond to treatments
# herbFit <- envfit(herbTreat, biomat_trimmed)

# Abundance based ordination for herbivores

# Mannual select
formula_string <- formula(abumat_trimmed_no_cole~PREDATOR+Condition(block))

formula_string <- paste("abumat_trimmed_no_cole~", 
                        paste(paste(treats_to_formula,
                                    collapse = "+"), 
                              "Condition(block)", 
                              sep = "+"), 
                        sep="")
rda_formula <- formula(formula_string)
herbTreatAbu <- rda(rda_formula, data = treats_trimmed)
anova(herbTreatAbu, by="terms")

addNamesAndPlot <- function(ipsResponses,ipsRDA,treats){
  sigIps <- names(ipsResponses$vectors$pvals)[ipsResponses$vectors$pvals <= 0.05]
  # Plot the graph
  plot(ipsRDA, type="n")
  points(ipsRDA, display = "sites",
         scaling = 2, pch=19, cex = 1.5,
         col = treats$treat)
  text(ipsRDA, display = "species",
       scaling = 1, select = sigIps)
  points(ipsRDA, display="cn", scaling = 3,
         col = "red",
         pch = 25)
  text(ipsRDA, display="cn", scaling = 3,
       col = "red")
} 

herbFitBio <- envfit(herbTreatBio, biomat_trimmed_no_cole)
herbFitAbu <- envfit(herbTreatAbu, abumat_trimmed_no_cole)

addNamesAndPlot(herbFitAbu, herbTreatAbu, treats_trimmed)
addNamesAndPlot(herbFitBio, herbTreatBio, treats_trimmed)

# ADD NAMES!!!
plot(herbTreatAbu)
text(herbTreatAbu, display="species", scaling=1)
plot(herbTreatBio)




# 2a. Removing effect of plant community PLANT pPCA ----
# Bosc suggested conditioning only on block (Appendix S1). 
pPCA <- rda(plants_trimmed ~ Condition(block), data=treats_trimmed)

# Site axes
ort_sites <- pPCA$CA$u
ort_sites <- as.data.frame(ort_sites)

fselbioins <- selectOrthogonalVars(maxa = 10, 
                                   focal_dataset =  biomat_trimmed_no_cole, 
                                   ort_sites = ort_sites,
                                   additional_variables = c("ipabu","ipbio","ipdiv","ipric"))
# nothing significant for biomass

fselbioins <- selectOrthogonalVars(maxa = 1, 
                                   focal_dataset =  abumat_trimmed_no_cole, 
                                   ort_sites = ort_sites,
                                   additional_variables = c("ipabu","ipbio","ipdiv","ipric"))

# With 10PCs and additional intermediate predator variables
# 4 are significant: PC2,PC1,PC3,PC7 for abundance

# Only PC4 is significant for biomass
# But only if I use first 4 PC's

# See which species were responsible for PC axes
par(mfrow=c(1,2))
plot(pPCA, choices = c(1, 2), display = "species", scaling =1)
plot(pPCA, choices = c(3, 7), display = "species")
par(mfrow=c(1,1))

# 2b. Intermediate Predator PCA ----
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

fselbioips <- selectOrthogonalVars(10, focal_dataset = abumat_trimmed_no_cole,
                     ort_sites = ipPCsites,
                     additional_variables = c("ipabu","ipbio","ipdiv","ipric"))

# Effect of IPs diversity on herbivore community, but withouth accounting for plant effect

# 2c. Partial RDA for herbivores abundance ----
form <- fselbioins$selected$call$formula
formula_string <- as.character(form)[3]
splitted_formula <- strsplit(formula_string, split=" ")
orthogonal_axes <- splitted_formula[[1]][grep("PC", splitted_formula[[1]])]
final_herb_abu_formula <- paste("abumat_trimmed_no_cole", "~", 
                                 "PREDATOR+WEEVIL125+WEEVIL25+INSECTICIDE", 
                                "+Condition(block+", paste(orthogonal_axes, collapse = "+"), ")", sep = "")

# Treatment effect on herbivorous communities conditioned on plant composition.
# There was also no effect of IPs community on herbivore communities.

# remove cole001 from the dataset
# abumat_trimmed <- abumat_trimmed[, -which(colnames(abumat_trimmed) == "cole001")]

herbAbuConditioned <-rda(formula(final_herb_abu_formula),
                         data=fselbioins$treatments)

# Because of small sample size I use the model based permutations
# https://fromthebottomoftheheap.net/2014/11/03/randomized-complete-block-designs-and-vegan/
anova(herbAbuConditioned, by="terms", permutations = 999)

herbivoresResponses <- envfit(herbAbuConditioned, abumat_trimmed)

# NAME ADDING LINE 
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

# 3. IPs ordiantion ----

# Removing effect of plants and herbivores on ips
selips <- selectOrthogonalVars(10, 
                     focal_dataset = ipabu_trimmed,
                     ort_sites = ort_sites,
                     additional_variables = c("habu","hbio","hdiv","hric"))

# GENERAL FUNCTION FOR OBTAINING THE FORMULA
getFormula <- function(focalDS = "abumat_trimmed_no_cole",
                       treats_to_formula,
                       selectionFormula = fselbioins){
  form <- fselbioins$selected$call$formula
  formula_string <- as.character(form)[3]
  splitted_formula <- strsplit(formula_string, split=" ")
  orthogonal_axes <- splitted_formula[[1]][grep("PC", splitted_formula[[1]])]
  final_herb_abu_formula <- paste(focalDS, "~", 
                                  paste(treats_to_formula, collapse = "+"), 
                                  "+Condition(block+", 
                                  paste(orthogonal_axes, 
                                        collapse = "+"), ")", sep = "")
  return(final_herb_abu_formula)
}

ipspRDAform <- getFormula(focalDS = "ipabu_trimmed",
           treats_to_formula = treats_to_formula,
           selectionFormula = selips)

ipsRDA <- rda(formula(ipspRDAform),
              data = selips$treatments)
# ipspca <- rda(ipabu_trimmed~CONTROL+WEEVIL125+WEEVIL25+Condition(block + PC6), 
#               data = treats_trimmed)
anova(ipsRDA, by="terms")
# There is an effect of w25

#####*************

ipsResponses <- envfit(ipsRDA, ipabu_trimmed)

# NAME ADDING LINE 
sigIps <- names(ipsResponses$vectors$pvals)[ipsResponses$vectors$pvals <= 0.05]

# Plot the graph
plot(ipsRDA, type="n")
points(ipsRDA, display = "sites",
       scaling = 2, pch=19, cex = 1.5,
       col = selips$treatments$treat)
text(ipsRDA, display = "species",
     scaling = 1, select = sigIps)
points(ipsRDA, display="cn", scaling = 3,
       col = "red",
       pch = 25)
text(ipsRDA, display="cn", scaling = 3,
     col = "red")


#####*************

plot(ipsRDA, display="species")
plot(ipsRDA)
anova(ipspca, by="terms")

# 4. Insect ordiantion: Pairwise comparisons ----
dim(abumat_trimmed)
dim(ipabu_trimmed)
abumat_arthropod <- cbind(abumat_trimmed, ipabu_trimmed)

arthropca <- rda(abumat_arthropod ~ Condition(block), data = treats_trimmed)
plot(arthropca, display = "species")
text(arthropca, display = "species")

# 5. Correlation plot for families ----
abuFamOrig <- contingencyTable2(ins_bio,"plot","family","amount")
bioFamOrig <- contingencyTable2(ins_bio,"plot","family","totbio")

abuFam <- decostand(abuFamOrig, method = "hel")
bioFam <- decostand(bioFamOrig, method = "hel")

abuFam_trimmed <- abuFam[treats_trimmed$codes, ]
bioFam_trimmed <- bioFam[treats_trimmed$codes, ]

artFampca <- rda(abuFam_trimmed ~ Condition(block+treat), data = treats_trimmed)
artFamRDA <- rda(abuFam_trimmed ~ PREDATOR+WEEVIL125+WEEVIL25 + Condition(block), data = treats_trimmed)
bioartFampca <- rda(bioFam_trimmed ~ Condition(block+treat), data = treats_trimmed)

par(mfrow=c(1,2))
biplot(artFampca, display = "species")
summary(artFampca)

anova(artFamRDA, by="terms")
plot(artFamRDA, display="species")
plot(artFamRDA)

biplot(bioartFampca,display = "species")
summary(bioartFampca)
par(mfrow=c(1,1))

panel.cor <- function(x, y){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(x, y), digits=2)
  txt <- paste0("R = ", r)
  cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt)
}
# Customize upper panel
upper.panel<-function(x, y){
  points(x,y, pch = 19)
}
# Create the plots
# pairs(abuFam_trimmed, 
#       lower.panel = panel.cor,
#       upper.panel = upper.panel)

# pairs(abuFam_trimmed)

abuFam_trim_paired <- abuFamOrig[treats_trimmed$codes, ]
bioFam_trim_paired <- bioFamOrig[treats_trimmed$codes, ]

library(psych)
# pairs.panels(log(abuFam_trim_paired+1), 
#              method = "pearson", # correlation method
#              hist.col = "#00AFBB",
#              density = TRUE,  # show density plots
#              ellipses = TRUE # show correlation ellipses
# )

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

# Correlations between species and families
abumat_trimmed

# 6. Herbivore family ordination ----
insfam_mat_abu <- contingencyTable2(ins_bio, "plot", "family", "amount")
herbfam_mat_abu <- insfam_mat_abu[, -grep("aran|mant", colnames(insfam_mat_abu))]

herbfam_mat_abu_trim <- herbfam_mat_abu[rownames(treats_trimmed), ]

famRda <- rda(herbfam_mat_abu_trim~PREDATOR+Condition(block),
              treats_trimmed)
anova(famRda, by = "terms")
plot(famRda)

# remove effect of plants on orders
conditionalPCs <- selectOrthogonalVars(10, herbfam_mat_abu_trim,
                     ort_sites = ort_sites,
                     additional_variables = c("ipabu","ipbio","ipdiv","ipric"))

famRda <- rda(herbfam_mat_abu_trim~PREDATOR+WEEVIL125+WEEVIL25+ Condition(block+ipric + PC2), conditionalPCs$treatments)
anova(famRda, by = "terms")
plot(famRda)

# 6b. Biomass based herbivore family ordination ----
insfam_mat_abu <- contingencyTable2(ins_bio, "plot", "family", "totbio")
herbfam_mat_abu <- insfam_mat_abu[, -grep("aran|mant", colnames(insfam_mat_abu))]

herbfam_mat_abu_trim <- herbfam_mat_abu[rownames(treats_trimmed), ]

famRda <- rda(herbfam_mat_abu_trim~PREDATOR+WEEVIL125+WEEVIL25+Condition(block),
              treats_trimmed)
anova(famRda, by = "terms")
plot(famRda)

# remove effect of plants on orders
conditionalPCs <- selectOrthogonalVars(10, herbfam_mat_abu_trim,
                                       ort_sites = ort_sites,
                                       additional_variables = c("ipabu","ipbio","ipdiv","ipric"))

famRda <- rda(herbfam_mat_abu_trim~PREDATOR+WEEVIL125+WEEVIL25+ Condition(block+ipric), conditionalPCs$treatments)
anova(famRda, by = "terms")
plot(famRda, display="species")
text(famRda, display="cn")

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
# Herbivores present in P and C treatments
source("code/pdi.R")
at <- 5 #abundance_treshold
ins_bio_at <- ins_bio[ins_bio$amount >= at, ] 
ins_bio_at <- ins_bio_at[-grep("aran|mant", ins_bio_at$morphotype), ]

predsites <- treats[treats$treat == "PREDATOR",]$codes
contsites <- treats[treats$treat == "CONTROL",]$codes

ins_bioOrig <- ins_bioins_bio <- ins_bio[-ips, ]

pabumat <- contingencyTable2(ins_bio_at[(ins_bio_at$plot %in% predsites), ],
                             "plot","morphotype","amount")
pbiomat <- contingencyTable2(ins_bio_at[ins_bio_at$plot %in% predsites, ],
                             "plot","morphotype","totbio")

cabumat <- contingencyTable2(ins_bio_at[ins_bio_at$plot %in% contsites, ],
                             "plot","morphotype","amount")
cbiomat <- contingencyTable2(ins_bio_at[ins_bio_at$plot %in% contsites, ],
                             "plot","morphotype","totbio")



comparable <- colnames(cbiomat)[colnames(cbiomat) %in% colnames(pbiomat)]

# 
# # Remove intermediate predators
# comparable
# 
# # Food plants for comparable herbivores
cp_treats <- treats_trimmed[treats_trimmed$treat %in% c("CONTROL","PREDATOR"),]$sites
csites <- treats_trimmed[treats_trimmed$treat %in% c("CONTROL"),]$sites
psites <- treats_trimmed[treats_trimmed$treat %in% c("PREDATOR"),]$sites

ins_bio_cp <- ins_bio_at[ins_bio_at$plot %in% cp_treats, ]
ins_bio_cp_comparable <- ins_bio_cp[ins_bio_cp$morphotype %in% comparable,]
ibc <- ins_bio_cp_comparable
ibc <- ibc[complete.cases(ibc),]

# Species in the exclosure treatment get _P note ant the end
ibc$morphotype <- as.character(ibc$morphotype)
ibc[ibc$plot %in% psites, ]$morphotype <- paste(ibc[ibc$plot %in% psites,]$morphotype, 
                                              "P", sep = "_")

# Only one contingency table for all woody species
compFood <- contingencyTable2(ibc,
                              "tree",
                              "morphotype",
                              "totbio")

dim(compFood)

# envdat <- data.frame(treat = rep(c("predator", "control"), 
#               c(dim(compFood)[1],
#                 dim(compFood)[1])))
# rownames(envdat) <- rownames(compFood)

species_point_color <- rep(rgb(255,194,10,150,
                               maxColorValue = 255), 
                           dim(compFood)[2])
species_point_color[grep("_P", colnames(compFood))] <- rgb(12,123,220,150,
                             maxColorValue = 255)

par(mfrow=c(1,1))
dietrda <- metaMDS(compFood)
foodDist <- vegdist(t(compFood))

plot(dietrda, type = "n", display = "species")
points(dietrda, display = "species", 
       col = species_point_color, pch=19, cex = 1.5)

# Shift vs pdi
#diet breadths of comparable species 
distspec <- as.matrix(foodDist)
comparable <- comparable[comparable %in% colnames(distspec)]
dbspec <- diet_breadth_ab[comparable]
shiftvals <- diag(distspec[colnames(distspec) %in% comparable,
         colnames(distspec) %in% paste(comparable,"P",sep ="_")])
shiftvals[shiftvals == 1] <- 0.999

# Linear regression may be weighted by abundance
cds <- ins_bio_cp[ins_bio_cp$morphotype %in% comparable, ]
cds$morphotype <- as.character(cds$morphotype)
coll_abu <- tapply(cds$amount, cds$morphotype, sum)

betareg_mod <- betareg::betareg(shiftvals~dbspec, 
                                weights = log(coll_abu),
                                type = "ML")
summary(betareg_mod)
# plot(betareg_mod)

plotdf <- as.data.frame(cbind(shiftvals, dbspec, coll_abu))

preddat <- predict(betareg_mod, 
                 newdata = plotdf,
                 type = "quantile", 
                 at = c(0.025, 0.975))

plotdf <- cbind(plotdf, preddat)

ggplot(plotdf, aes(x = dbspec, y = shiftvals)) +
  geom_point(size = log(coll_abu)) +
  geom_line(aes(y = predict(betareg_mod, plotdf))) +
  geom_ribbon(aes(ymin= q_0.025, ymax = q_0.975), alpha=0.2)


ggplot(plotdf, aes(y=shiftvals, x=dbspec)) +
  geom_point(size = log(coll_abu), shape = 21) +
  geom_line(aes(y = predict(betareg_mod, plotdf))) +
  theme_bw()
# 
# ggplot()+
#   geom_point(size = log(coll_abu))+
#   stat_smooth(method="lm")

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
plot(czvalsc$c, czvalsc$z, col = colorsC, pch = 19)
abline(h=2.5)
abline(v=0.62)

plot(czvalsp$c, czvalsp$z, col = colorsP, pch = 19)
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