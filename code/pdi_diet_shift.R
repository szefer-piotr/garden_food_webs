# RDA

# Prepare data ----
rm(list=ls())
source("code/data_processing_code.R")
library(psych)

treats_to_plot <- as.character(unique(treats$treat))[c(3,4)]

treats_trimmed <- treats[(treats$treat %in% treats_to_plot), ]
treats_trimmed$codes <- as.character(treats_trimmed$codes)
treats_dummy <- dummy.code(as.character(treats_trimmed$treat))
treats_trimmed <- cbind(treats_trimmed, treats_dummy)
treats_trimmed$sites <- rownames(treats_trimmed)

cp_treats <- treats_trimmed[treats_trimmed$treat %in% c("CONTROL","PREDATOR"),]$sites
csites <- treats_trimmed[treats_trimmed$treat %in% c("CONTROL"),]$sites
psites <- treats_trimmed[treats_trimmed$treat %in% c("PREDATOR"),]$sites

source("code/pdi.R")
at <- 5 #abundance_treshold
ins_bio_at <- ins_bio[ins_bio$amount >= at, ] 
ins_bio_at <- ins_bio_at[-grep("aran|mant", ins_bio_at$morphotype), ]

contsites <- treats[treats$treat == "CONTROL",]$codes
predsites <- treats[treats$treat == "PREDATOR",]$codes
cbiomat <- contingencyTable2(ins_bio_at[ins_bio_at$plot %in% contsites, ],
                             "plot","morphotype","totbio")
pbiomat <- contingencyTable2(ins_bio_at[ins_bio_at$plot %in% predsites, ],
                             "plot","morphotype","totbio")
comparable <- colnames(cbiomat)[colnames(cbiomat) %in% colnames(pbiomat)]

cp_tr <- c(as.character(contsites),
           as.character(predsites))

ins_bio_cp <- ins_bio_at[ins_bio_at$plot %in% cp_tr, ]
ins_bio_cp_comparable <- ins_bio_cp[ins_bio_cp$morphotype %in% comparable,]
ibc <- ins_bio_cp_comparable
ibc <- ibc[complete.cases(ibc),]

# Species in the exclosure treatment get _P note ant the end
ibc$morphotype <- as.character(ibc$morphotype)
ibc[ibc$plot %in% predsites, ]$morphotype <- paste(ibc[ibc$plot %in% predsites,]$morphotype,"P", sep = "_")

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

shiftvals
dbspec
coll_abu

w_sv <- rep(shiftvals, log(coll_abu))
w_db <- rep(dbspec, log(coll_abu))

betareg_mod <- betareg::betareg(shiftvals~dbspec, 
                                weights = log(coll_abu),
                                type = "ML")


b_m_test_w <- betareg::betareg(w_sv~w_db,
                                type = "ML")

betareg_mod_inv <- betareg::betareg(shiftvals~dbspec, 
                                weights = coll_abu,
                                type = "ML")

summary(betareg_mod)
summary(betareg_mod_inv)
summary(b_m_test_w)

# plot(betareg_mod)

plotdf <- as.data.frame(cbind(shiftvals, dbspec, coll_abu))

preddat <- predict(betareg_mod, 
                   newdata = plotdf,
                   type = "quantile", 
                   at = c(0.025, 0.975))

plotdf <- cbind(plotdf, preddat)

library(ggplot2)
ggplot(plotdf, aes(x = dbspec, y = shiftvals)) +
  geom_point(size = log(coll_abu), alpha = 0.5) +
  geom_line(aes(y = predict(betareg_mod, plotdf)), col = "grey20",lwd = 2) +
  geom_ribbon(aes(ymin= q_0.025, ymax = q_0.975), alpha=0.2) + 
  theme_bw()+
  xlab("Specialization (Paired Differences Index)")+
  ylab("Bray-Curtis diet dissimilarity")
