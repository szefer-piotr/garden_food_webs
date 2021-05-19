source("code/data_processing_code.R")

ins_bio_noip <- ins_bio[-grep("aran|mant", ins_bio$morphotype), ]
fullbio <- tapply(ins_bio$totbio , ins_bio$tree, sum, na.rm = T)
fullabu <- tapply(ins_bio$amount , ins_bio$tree, sum, na.rm = T)

no_herb <- table(ins_bio_noip$tree, ins_bio_noip$morphotype)
attractivness <- rowSums(no_herb)

specperbio <- attractivness/fullbio

sort(specperbio)
sort(fullabu)

# More abundant species gather more species, this is almost a perfect relationship
# cbind(attractivness,fullbio)
# lm1  <- lm(log(attractivness)~log(fullbio))
# plot(log(attractivness)~log(fullbio)) # more biomass, more species
# 
# abline(lm1)
# 
# cbind(attractivness,fullabu)
# lm1  <- lm(log(attractivness)~log(fullabu)) # sampling curve for insects
# plot(log(attractivness)~log(fullabu))
# abline(lm1)

# We may not have any specialists at all, but only undersampled species.
# main_biomass
# sel_sp <- unique(main[main$LIFE.FORM %in% c("tree", "shrub"), ]$SP_CODE)
# sel_sp <- as.character(sel_sp)[sel_sp %in% c("PIPTAR", "MELAMU")]
sel_sp <- c("PIPTAR", "MELAMU")
sel_sp <- toupper(c(rev(names(sort(fullabu)))))[1:10]

specqua <- main[ ,c("SP_CODE", "LDMC", "HERB", "SLA", "WATER", "CODE")]
specqua_filtered <- specqua[specqua$SP_CODE %in% sel_sp, ]
specqua_filtered$SP_CODE <- as.character(specqua_filtered$SP_CODE)

# Stack dataset for a given set of species for a facet plot
specqua_stacked <- stack(specqua_filtered[, c("LDMC", "HERB", "SLA", "WATER")])
specqua_stacked$sp_code <- rep(specqua_filtered[, "SP_CODE"], 4)
specqua_stacked$code <- rep(specqua_filtered[, "CODE"], 4)

# Visualisation
# library(ggplot2)
# ggplot(specqua_stacked, aes(x = sp_code, y = values))+
#   # geom_jitter(width = 0.1)+
#   geom_boxplot(varwidth = TRUE)+
#   facet_wrap(vars(ind), scales = "free")+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

specqua$plco <- paste(specqua$SP_CODE, specqua$CODE, sep="")
specqua <- specqua[!duplicated(specqua$plco),]
rownames(specqua) <- specqua$plco

plplca <- specqua[ ,c("LDMC", "HERB", "SLA", "WATER")]
plplca <- plplca[complete.cases(plplca),]
library(vegan)
pca1 <- rda(plplca, na.rm=T)
# plot(pca1, display="species", choices = c(1,2))
colorsf <- as.factor(substr(rownames(plplca),1,6))
colors <- as.numeric(colorsf)
# plot(pca1, display="sites", choices = c(1,2), 
#      col = colors, pch=19)
ax1 <- pca1$CA$u[,1]
ax2 <- pca1$CA$u[,2]

total_quality <- data.frame(quality = ax1, palat = ax2)
total_quality$spec <- substr(rownames(plplca),1,6)
mq <-  tapply(total_quality$quality, total_quality$spec, mean, na.rm=TRUE)
mp <-  tapply(total_quality$palat, total_quality$spec, mean, na.rm=TRUE)
plant_quality <- data.frame(quality = mq, palat = mp)
# plot(quality~palat, plant_quality)
# summary(lm(quality~palat, plant_quality))
# abline(lm(quality~palat, plant_quality))
rownames(plant_quality) <- tolower(rownames(plant_quality))
