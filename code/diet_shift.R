# diet shift ordination

insects <- read.table("datasets/arthropods_clean.txt")
treats  <- read.table("datasets/treatments_clean.txt")
plants  <- read.table("datasets/plants_clean.txt")
size_dat <-read.table("datasets/size_dat_bio.txt")

library("bipartite")
library("vegan")
source("code/bio_log_ratio.R")

# This shift maybe should be beased on abundance rather than like in the plot
ins_bio$morphplot <- paste(ins_bio$family, ins_bio$plot, sep="")
ins_bio$tree <- as.character(ins_bio$tree)
abuinsmat <- contingencyTable2(ins_bio, "tree", "morphplot", "amount" , FALSE)
#dist(x, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)

# For the whole dataset
abunmds <- metaMDS(abuinsmat, distance = "bray")
plot(abunmds, display = "species")
text(abunmds, display = "species")
abuvd <- as.matrix(vegdist(t(abuinsmat), method = "bray"))

# Maybe this is not the best representation
tp <- as.character(treats[treats$treat %in% c("PREDATOR"), ]$codes)
tc <- as.character(treats[treats$treat %in% c("CONTROL") , ]$codes)

abuvd[substr(rownames(abuvd), 5,10) %in% tp, 
      substr(colnames(abuvd), 5,10) %in% tc]

# I just need to 
