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
plotweb(pdimat)

diet_breadth <- PDI(pdimat, normalise = F, log = T)
