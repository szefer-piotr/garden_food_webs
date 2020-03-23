# Insect

insects <- read.table("datasets/arthropods_clean.txt")
treats  <- read.table("datasets/treatments_clean.txt")
plants  <- read.table("datasets/plants_clean.txt")
size_dat <-read.table("datasets/size_dat_bio.txt")

library(bipartite)
library(lme4)
library(lmerTest)
source("code/bio_log_ratio.R")
source("code/contingencyTable.R")
source("code/weighted-modularity-LPAwbPLUS/code/R/LPA_wb_plus.R")
source("code/weighted-modularity-LPAwbPLUS/code/R/MODULARPLOT.R")
source("code/weighted-modularity-LPAwbPLUS/code/R/convert2moduleWeb.R")
source("code/weighted-modularity-LPAwbPLUS/code/R/GetModularInformation.R")

ctrpl <- treats[treats$treat == "CONTROL",]$codes
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

dfcg1 <- matrix_to_data(cg1)
dfcg2 <- matrix_to_data(cg2)
dfcg3 <- matrix_to_data(cg3)
dfcg4 <- matrix_to_data(cg4)
dfcg5 <- matrix_to_data(cg5)
dfcg6 <- matrix_to_data(cg6)

dfcg <- rbind(dfcg1,
      dfcg2,
      dfcg3,
      dfcg4,
      dfcg5,
      dfcg6)

library(reshape)
cmx <- cast(dfcg, rowname~colname, fun.aggregate = sum)
rownames(cmx) <- cmx[,1]
cmx <- cmx[,-1]
cmx[cmx>0] <- 1
cmx
rankabu <- colSums(cmx)
hist(rankabu, breaks = 50)
table(rankabu)
