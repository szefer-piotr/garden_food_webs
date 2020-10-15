# Insect

insects <- read.table("datasets/arthropods_clean.txt")
treats  <- read.table("datasets/treatments_clean.txt")
plants  <- read.table("datasets/plants_clean.txt")
size_dat <-read.table("datasets/size_dat_bio.txt")

library(bipartite)
library(lme4)
library(lmerTest)
# source("code/bio_log_ratio.R")
source("code/data_processing_code.R")
source("code/contingencyTable.R")
source("code/weighted-modularity-LPAwbPLUS/code/R/LPA_wb_plus.R")
source("code/weighted-modularity-LPAwbPLUS/code/R/MODULARPLOT.R")
source("code/weighted-modularity-LPAwbPLUS/code/R/convert2moduleWeb.R")
source("code/weighted-modularity-LPAwbPLUS/code/R/GetModularInformation.R")

ctrpl <- treats[treats$treat == "CONTROL",]$codes
ctrpl <- treats[!(treats$treat %in% c("FUNGICIDE", "INSECTICIDE")),] $codes
ctrlgardnets <- gardnets[ctrpl]
cg1 <- ctrlgardnets[[1]]
cg2 <- ctrlgardnets[[2]]
cg3 <- ctrlgardnets[[3]]
cg4 <- ctrlgardnets[[4]]
cg5 <- ctrlgardnets[[5]]
cg6 <- ctrlgardnets[[6]]

#
matrix_to_data <- function(mat){
  df <- data.frame()
  colnm <- colnames(mat)
  rownm <- rownames(mat)
  for(i in 1:nrow(mat)){
    for(j in 1:ncol(mat)){
      dfrow <- data.frame(rowname = rownm[i],
                          colname = colnm[j],
                          val = mat[i,j])
      df <- rbind(df, dfrow)
    }
  }
  return(df)
}

dfcg <- data.frame()
for(nm in names(ctrlgardnets)){
  cgx <- matrix_to_data(ctrlgardnets[[nm]])
  dfcg <- rbind(dfcg, cgx)
}

library(reshape)
cmx <- cast(dfcg, rowname~colname, fun.aggregate = sum)
rownames(cmx) <- cmx[,1]
cmx <- cmx[,-1]
cmx[cmx>0] <- 1
cmx
rankabu <- colSums(cmx)
hist(rankabu, breaks = 50)
table(rankabu)
treshold <- 2

# How often specialists vs generalists swith their resources under the treatment (predator removal)
sum(rankabu <= treshold) # number of specialists
sum(rankabu > treshold)  # rest

dim(cmx)

#co-inertia
# library(cocorresp)
# 
# data(beetles, plants)
# coin <- coinertia(beetles, plants)
# coin
# anova(coin)
