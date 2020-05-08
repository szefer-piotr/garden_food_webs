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
# pdimat <- contingencyTable2(pdidata, "tree", "morphotype", "amount")

par(mfrow=c(1,1))
library(bipartite)
# plotweb(pdimat)

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
# 
# pdimat[, "orth037"]
# pdimat[, "orth092"]
