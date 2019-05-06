insects <- read.table("datasets/arthropods_clean.txt")
treats  <- read.table("datasets/treatments_clean.txt")
plants  <- read.table("datasets/plants_clean.txt")
size_dat <-read.table("datasets/size_dat_bio.txt")

library("bipartite")
source("code/bio_log_ratio.R")
# gardnets - networks in individual plots
# biofulldf biomass of all insects and trees

# Is it possible to study average number of plants per herbivore species using raw data?
# That could be done for arthropod species present in both predator and control plots. That would be problematic if there are differences in plant communities between blocks

# Simple vulnerability and generality patterns in networks
gnames <- as.character(treats$codes)
genvuldf <- data.frame()
for(gname in gnames){
  print(gname)
  nlres <- networklevel(gardnets[[gname]], index = "vulnerability")
  tt <- tryCatch(networklevel(gardnets[[gnames[gname]]], index = "vulnerability"),error=function(e) e, warning=function(w) w)
  ifelse(is(tt,"warning"),next,"OK")
  subgvdf <- data.frame(plot = gname, gen = nlres[1], vul = nlres[2])
  genvuldf <- rbind(genvuldf, subgvdf)
}

rownames(treats) <- treats$codes
genvuldf$trt <- treats[genvuldf$plot, ]$treat
genvuldf$block <- substr(genvuldf$plot, 3,4)

# But generality would also change if the number of plant at the plot is smaller!
library(ggplot2)
p <- ggplot(genvuldf, aes(x = trt, y = gen, col = block))
p + geom_point()

# What can I say using this data?