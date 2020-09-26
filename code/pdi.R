# PDI

# Calculate PDI for herbivores and correlate with these coordinates
# source("code/rank_abundance.R")

# diet_breadth <- data.frame(herbivores = colnames(cmx), 
#                            db = rankabu)

# source("code/data_processing_code.R")

pdidata <- ins_bio
pdidata$trt <- as.character(treats[pdidata$plot, ]$treat)

# ON/OFF All plots or all except for fungicide and insecticide treatments used
# pdidata <- pdidata[!(pdidata$trt %in% c("FUNGICIDE", "INSECTICIDE")),]

pdimat <- contingencyTable2(pdidata, "tree", "morphotype", "totbio")
pdissmat <- contingencyTable2(pdidata, "tree", "morphotype", "amount")
pdiss <- colSums(pdissmat)

# Abundance data instead of biomass data
# pdimat <- contingencyTable2(pdidata, "tree", "morphotype", "amount")

# par(mfrow=c(1,1))
library(bipartite)
# plotweb(pdimat)

# pdimat[pdimat ==0 ] <- NA
# colSums(pdimat>0)
# nm <- "orth024"

# Test of the index performance in case of present but unused resources
# rawvals <- pdimat[, "orth024"]
# toremna <- pdimat[, "orth024"]
# toremna[toremna == 0] <- NA
# nona_vals <- toremna[!is.na(toremna)]

# There is a difference between having NA or 0. But we assume that all presented plants were there for taking
# PDI(rawvals,normalise = T, log = T)
# PDI(nona_vals,normalise = T, log = T)

# Biomass based
diet_breadth <- PDI(pdimat, normalise = T, log = T) # biomass based
diet_breadth_ab <- PDI(pdissmat, normalise = T, log = T) # abundance based

# dbnorm <- PDI(pdimat, normalise = T, log = T)
# dbrel <- PDI(pdimat, normalise = F, log = T)
# # plot(dbnorm~dbrel, log="x")
# min(diet_breadth)
# db <- diet_breadth
# db <- db[db != 0]
# db[db == min(db)]

# Remove zeros whenever no resources were found used
diet_breadth <- diet_breadth[!(diet_breadth==0)]
diet_breadth_ab <- diet_breadth_ab[!(diet_breadth_ab==0)]

# Relationship between PDI based on abundance and biomass.
# comp_nms_pdi <- names(diet_breadth)[names(diet_breadth) %in% names(diet_breadth_ab)]
# 
# plot(diet_breadth[comp_nms_pdi], 
#      diet_breadth_ab[comp_nms_pdi])
# abline(0, 1)
