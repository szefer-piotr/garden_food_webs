rm(list = ls())

insects <- read.table("datasets/arthropods_clean.txt")
treats  <- read.table("datasets/treatments_clean.txt")
plants  <- read.table("datasets/plants_clean.txt")
size_dat <-read.table("datasets/size_dat_bio.txt")

library("bipartite")
library("igraph")
# source("code/bio_log_ratio.R")
source("code/data_processing_code.R")


# ********************************************************************
# Example of the control - predator comparison for garden 1
nc <- "w1g1p3" #control
np <- "w1g1p4" #predator

code <- "w1g1p3"

networkPlot <- function(garden_list, 
                        code,
                        hmltp = 1.5, 
                        pmltp=1.5,
                        emltp=1,
                        rem.ips = TRUE){
  net <- garden_list[[code]]
  if(rem.ips){
    net <- net[,-grep("aran|mant", colnames(net))]
  }
  print(net)
  
  netgraph <- graph_from_incidence_matrix(net,
                                           weighted = T)
  
  print("Net constructed")
  
  # 1. Customize the network
  
  
  # Weights
  # Plants
  pl_weight <- plants[(plants$CODE == code), 
                      c("SP_CODE", "WEIGHT")]
  print(pl_weight)
  
  print("Weights assigned")
  
  plw <- pl_weight$WEIGHT
  names(plw) <- pl_weight$SP_CODE
  
  # Herbivores/interactions
  hwc <- colSums(net)
  plws <- plw[names(plw) %in% names(V(netgraph))]
  hwcs <- hwc[names(hwc) %in% names(V(netgraph))]
  
  V(netgraph)$size <- c(plws*pmltp, hwcs*hmltp)
  
  E(netgraph)$width <- E(netgraph)$weight*emltp
  
  plts <- names(V(netgraph)) %in% names(plw)
  herbivores <- !plts
  
  print("bool for p and h")
  
  # Colors
  colors <- rep("green", length(V(netgraph)))
  colors[herbivores] <- "gold"
  V(netgraph)$color <- colors
  
  print("Colors assigned")
  print("Weights assigned")
  
  # 2. Plot weighted graph
  l <- layout_in_circle(netgraph)
  plot(netgraph, layout = l, 
       vertex.label = NA, 
       main = code)
  
}

networkPlot(abugardnets, "w1g3p4", hmltp = 0.6, emltp = 0.6,
            pmltp = 1)

# ****************************************************************
bipartieNetworkPlot <- function(garden_list, code, rem.ips=T, pmltp = 1){
  net <- garden_list[[code]]
  
  if(rem.ips){
    net <- net[,-grep("aran|mant", colnames(net))]
  }
  
  # Weights
  # Weights for vertices
  plntabudf <- plants[(plants$CODE == code), c("SP_CODE", "WEIGHT")]
  plntabu <- plntabudf$WEIGHT
  names(plntabu) <- plntabudf$SP_CODE

  plotweb(net, 
          low.abun = plntabu*pmltp)
  
}

bipartieNetworkPlot(abugardnets, pmltp = 5, "w1g1p3")

filterCodes <- function(trt){
  return(as.character(treats[treats$treat %in% trt, ]$codes))
}
x11(8,6)

# jpeg("PREDATOR_networks.jpg", width = 1200, height = 1000,
     # quality = 100, pointsize = 20)
par(mfrow = c(3,2))
for(codes in filterCodes("PREDATOR")){
  print(codes)
  bipartieNetworkPlot(abugardnets, pmltp = 5, codes)
}
# dev.off()





##### notes ******************************************

netc <- gardnetsfam[[nc]]
netp <- gardnetsfam[[np]]

# # igrpah graph
# netcgraph <- graph_from_incidence_matrix(netc,
#                                          weighted = T)
# class(netcgraph)

# Weights for vertices
plntabudfc <- plants[(plants$CODE == nc), c("SP_CODE", "WEIGHT")]
plntabudfp <- plants[(plants$CODE == np), c("SP_CODE", "WEIGHT")]
plntabuc <- plntabudfc$WEIGHT
plntabup <- plntabudfp$WEIGHT
names(plntabuc) <- plntabudfc$SP_CODE
names(plntabup) <- plntabudfp$SP_CODE
hwc <- colSums(netc)
hwp <- colSums(netp)

pnc <- rownames(netc)
pnp <- rownames(netp)
hnc <- colnames(netc)
hnp <- colnames(netp)

l <- layout_in_circle(netcgraph)

plot(netcgraph, layout = l) 

# par(mfrow=c(1,1))
plotweb(netc[, c("orth", "lepi","cole","hemi","homo")], 
        low.abun = plntabuc)
plotweb(netp[, c("orth", "lepi","hemi","homo")], 
        low.abun = plntabup)

# Seems like at least in the case of orthoptera majority of the interactions have shifted to a different resource
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
colc <- rep(cbPalette[4], length(colnames(netc)))
colc[!(colnames(netc) %in% c("orth", "lepi","cole","hemi","homo"))] <- cbPalette[7]
colp <- rep(cbPalette[4], length(colnames(netp)))
colp[!(colnames(netp) %in% c("orth", "lepi","hemi","homo"))] <- cbPalette[7]
x11(6,6)
par(mfrow=c(2,1))
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
x11(6,6)
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
