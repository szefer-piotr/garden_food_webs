# PDI

# Calculate PDI for herbivores and correlate with these coordinates
# source("code/rank_abundance.R")

# diet_breadth <- data.frame(herbivores = colnames(cmx), 
#                            db = rankabu)

# source("code/data_processing_code.R")

pdidata <- ins_bio
pdidata$trt <- as.character(treats[pdidata$plot, ]$treat)
pdidata <- pdidata[!(pdidata$trt %in% c("FUNGICIDE", "INSECTICIDE")),]
pdimat <- contingencyTable2(pdidata, "tree", "morphotype", "totbio")

# Abundance data instead of biomass data
# pdimat <- contingencyTable2(pdidata, "tree", "morphotype", "amount")

par(mfrow=c(1,1))
library(bipartite)
# plotweb(pdimat)

# pdimat[pdimat ==0 ] <- NA

# colSums(pdimat>0)
# nm <- "orth024"

rawvals <- pdimat[, "orth024"]
toremna <- pdimat[, "orth024"]
toremna[toremna == 0] <- NA
nona_vals <- toremna[!is.na(toremna)]

# There is a difference between having the whole 
PDI(rawvals,normalise = T, log = T)
PDI(nona_vals,normalise = T, log = T)


diet_breadth <- PDI(pdimat, normalise = T, log = T)

# dbnorm <- PDI(pdimat, normalise = T, log = T)
# dbrel <- PDI(pdimat, normalise = F, log = T)
# 
# # plot(dbnorm~dbrel, log="x")
# 
# min(diet_breadth)
# 
# db <- diet_breadth
# db <- db[db != 0]
# db[db == min(db)]

# I need to remove zeros
# pdimat[, "cole004"]

diet_breadth <- diet_breadth[!(diet_breadth==0)]

diet_breadth[diet_breadth==1]

# I need to remove zeros
pdimat[, "cole014"]

diet_breadth["orth037"]
pdimat[, "orth037"]
# pdimat[, "orth092"]
