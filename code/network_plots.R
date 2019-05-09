insects <- read.table("datasets/arthropods_clean.txt")
treats  <- read.table("datasets/treatments_clean.txt")
plants  <- read.table("datasets/plants_clean.txt")
size_dat <-read.table("datasets/size_dat_bio.txt")

library("bipartite")
source("code/bio_log_ratio.R")

gardnetsfam

# Example of the control - predator comparison
nc <- "w1g1p3"
np <- "w1g1p4"
netc <- gardnetsfam[[nc]]
netp <- gardnetsfam[[np]]
plntabudfc <- plants[(plants$CODE == nc), c("SP_CODE", "WEIGHT")]
plntabudfp <- plants[(plants$CODE == np), c("SP_CODE", "WEIGHT")]
plntabuc <- plntabudfc$WEIGHT
plntabup <- plntabudfp$WEIGHT
names(plntabuc) <- plntabudfc$SP_CODE
names(plntabup) <- plntabudfp$SP_CODE
par(mfrow=c(2,1))
plotweb(netc[, c("orth", "lepi","cole","hemi","homo")], low.abun = plntabuc)
plotweb(netp[, c("orth", "lepi","hemi","homo")], low.abun = plntabup)
# Seems like at least in the case of orthoptera majority of the interactions have shifted to a different resource
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
colc <- rep(cbPalette[4], length(colnames(netc)))
colc[!(colnames(netc) %in% c("orth", "lepi","cole","hemi","homo"))] <- cbPalette[7]
colp <- rep(cbPalette[4], length(colnames(netp)))
colp[!(colnames(netp) %in% c("orth", "lepi","hemi","homo"))] <- cbPalette[7]
plotweb(netc, low.abun = plntabuc, col.high = colc)
plotweb(netp, low.abun = plntabup, col.high = colp)

nc <- "w1g6p6"
np <- "w1g6p2"
netc <- gardnetsfam[[nc]]
netp <- gardnetsfam[[np]]
plntabudfc <- plants[(plants$CODE == nc), c("SP_CODE", "WEIGHT")]
plntabudfp <- plants[(plants$CODE == np), c("SP_CODE", "WEIGHT")]
plntabuc <- plntabudfc$WEIGHT
plntabup <- plntabudfp$WEIGHT
names(plntabuc) <- plntabudfc$SP_CODE
names(plntabup) <- plntabudfp$SP_CODE
par(mfrow=c(2,1))
plotweb(netc[, colnames(netc) %in% c("orth", "lepi","cole","hemi","homo")], low.abun = plntabuc)
plotweb(netp[, colnames(netp) %in% c("orth", "lepi","cole","hemi","homo")], low.abun = plntabup)
x11(4,4)
par(mfrow=c(2,1))
plotweb(netc, low.abun = plntabuc, col.high = colc)
plotweb(netp, low.abun = plntabup, col.high = colp)
par(mfrow=c(1,1))

# Example size distribution of individuals for Orthopterans at garden 1
nc <- "w1g1p3"
np <- "w1g1p4"

ic <- insects[insects$plot == nc & insects$family == "orth", 
              c("morphotype","amount")]
ip <- insects[insects$plot == np & insects$family == "orth",
              c("morphotype","amount")]

abumassc <- size_dat[ic$morphotype, "bio"]
abumassp <- size_dat[ip$morphotype, "bio"]

abuindc <- rep(abumassc, ic$amount)
abuindp <- rep(abumassp, ip$amount)

aidf <- data.frame(bio = c(abuindc,abuindp),
                   trt = rep(c("C","P"), c(length(abuindc),length(abuindp))))

t.test(log(abuindc), log(abuindp))

ggplot(aidf, aes(x=log(bio), fill = trt)) + 
  geom_histogram(alpha = 0.5, aes(y = ..density..), position = 'identity')

hist(log(abuindc), breaks = 10, xlim=c(-7, 0), ylim = c(0,50))
par(new=T)
hist(log(abuindp), breaks = 10, xlim=c(-7, 0), ylim = c(0,50))
