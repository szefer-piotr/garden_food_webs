# source("code/data_processing_code.R")

# Plants from C and P that are in both sites.
ins_bio_noip <- ins_bio[-grep("aran|mant", ins_bio$morphotype), ]
fullbio <- tapply(ins_bio$totbio , ins_bio$tree, sum, na.rm = T)
fullabu <- tapply(ins_bio$amount , ins_bio$tree, sum, na.rm = T)

no_herb <- table(ins_bio_noip$tree, ins_bio_noip$morphotype)
attractivness <- rowSums(no_herb)

specperbio <- attractivness/fullbio

sort(specperbio)

# More abundant species gather more species, this is almost a perfect relationship
cbind(attractivness,fullbio)
lm1  <- lm(log(attractivness)~log(fullbio))
plot(log(attractivness)~log(fullbio))
abline(lm1)

cbind(attractivness,fullabu)
lm1  <- lm(log(attractivness)~log(fullabu))
plot(log(attractivness)~log(fullabu))
abline(lm1)
# We may not have any specialists at all, but only undersampled species.

main_biomass
main

specqua <- main[ ,c("SP_CODE", "LDMC", "HERB", "SLA", "WATER", "CODE")]
specqua$plco <- paste(specqua$SP_CODE, specqua$CODE, sep="")
specqua <- specqua[!duplicated(specqua$plco),]
rownames(specqua) <- specqua$plco

plplca <- specqua[ ,c("LDMC", "HERB", "SLA", "WATER")]
plplca <- plplca[complete.cases(plplca),]
library(vegan)
pca1 <- rda(plplca, na.rm=T)
plot(pca1, display="species", choices = c(1,2))
colorsf <- as.factor(substr(rownames(plplca),1,6))
colors <- as.numeric(colorsf)
plot(pca1, display="sites", choices = c(1,2), 
     col = colors, pch=19)
ax1 <- pca1$CA$u[,1]
ax2 <- pca1$CA$u[,2]

total_quality <- data.frame(quality = ax1, palat = ax2)
total_quality$spec <- substr(rownames(plplca),1,6)
mq <-  tapply(total_quality$quality, total_quality$spec, mean, na.rm=TRUE)
mp <-  tapply(total_quality$palat, total_quality$spec, mean, na.rm=TRUE)
plant_quality <- data.frame(quality = mq, palat = mp)
plot(quality~palat, plant_quality)
summary(lm(quality~palat, plant_quality))
abline(lm(quality~palat, plant_quality))
rownames(plant_quality) <- tolower(rownames(plant_quality))
