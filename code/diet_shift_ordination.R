# Herbivorius insects as sites

source("code/data_processing_code.R")

treats$block <- substr(treats$codes, 3,4)

# How dataset should be structured?
psites <- treats[treats$treat %in% "PREDATOR", ]$codes
csites <- treats[treats$treat %in% "CONTROL", ]$codes

sites <- data.frame(pred=psites, cont=csites)
sites$block <- unique(treats$block)

bl = "g1"

# herb before and herb after within block
maindat <- data.frame()
for(bl in unique(treats$block)){
  psite <- as.character(sites[sites$block == bl,"pred"])
  csite <- as.character(sites[sites$block == bl,"cont"])
  
  pgard <- gardnets[[psite]]
  cgard <- gardnets[[csite]]
  
  comp <- colnames(cgard)[colnames(cgard) %in% colnames(pgard)]
  
  sub_ins_bio <- ins_bio[ins_bio$plot %in% c(psite,csite), ]
  
  comp_ins_bio <- sub_ins_bio[sub_ins_bio$morphotype %in% comp,  ]
  
  desc <- "totbio"
  # desc <- "amount"
  
  pdat <- comp_ins_bio[comp_ins_bio$plot %in% psite, c("morphotype",
                                                       "tree", 
                                                       desc)]
  pdat$morphotype <- paste(pdat$morphotype, "p", sep = "_")
  pdat$trt <- "predator"
  
  cdat <- comp_ins_bio[comp_ins_bio$plot %in% csite, c("morphotype",
                                                       "tree", 
                                                       desc)]
  cdat$morphotype <- paste(cdat$morphotype, "c", sep = "_") 
  cdat$trt <- "control"
  
  cpdat <- rbind(pdat,cdat)
  cpdat$block <- bl
  
  maindat <- rbind(maindat, cpdat)
}

maindat$morphotype <- paste(maindat$morphotype, maindat$block, sep="_")
# rownames(maindat) <- maindat$morphotype
mainorddat <- contingencyTable2(maindat, "morphotype", "tree", "totbio")

# make trt data

maintrt <- data.frame()
for(mpt in unique(maindat$morphotype)){
  submd <- maindat[maindat$morphotype == mpt,]
  rowsubmd <- submd[1,c("morphotype","trt", "block")]
  maintrt <- rbind(maintrt, rowsubmd)
}

rownames(maintrt) <- maintrt$morphotype
mdstand <- decostand(mainorddat, "hel")
dsrda <- rda(mdstand~trt+Condition(block), maintrt)
anova(dsrda, by="terms")
plenv <- envfit(dsrda, mdstand, choices = c(1,2))
plot(dsrda, display="sites")
text(dsrda, display="species")

# The same but for plants as sites
bl = "g1"

# herb before and herb after within block
maindat <- data.frame()
for(bl in unique(treats$block)){
  psite <- as.character(sites[sites$block == bl,"pred"])
  csite <- as.character(sites[sites$block == bl,"cont"])
  
  pgard <- gardnets[[psite]]
  cgard <- gardnets[[csite]]
  
  comp <- rownames(cgard)[rownames(cgard) %in% rownames(pgard)]
  
  sub_ins_bio <- ins_bio[ins_bio$plot %in% c(psite,csite), ]
  
  comp_ins_bio <- sub_ins_bio[sub_ins_bio$tree %in% comp,  ]
  
  desc <- "totbio"
  # desc <- "amount"
  
  pdat <- comp_ins_bio[comp_ins_bio$plot %in% psite, c("tree",
                                                       "morphotype", 
                                                       desc)]
  pdat$tree <- paste(pdat$tree, "p", sep = "_")
  pdat$trt <- "predator"
  
  cdat <- comp_ins_bio[comp_ins_bio$plot %in% csite, c("tree",
                                                       "morphotype", 
                                                       desc)]
  cdat$tree <- paste(cdat$tree, "c", sep = "_") 
  cdat$trt <- "control"
  
  cpdat <- rbind(pdat,cdat)
  cpdat$block <- bl
  
  maindat <- rbind(maindat, cpdat)
}

maindat$tree <- paste(maindat$tree, maindat$block, sep="_")
# rownames(maindat) <- maindat$morphotype
mainorddat <- contingencyTable2(maindat, "tree", "morphotype", "totbio")

# make trt data
maintrt <- data.frame()
for(mpt in unique(maindat$tree)){
  submd <- maindat[maindat$tree == mpt,]
  rowsubmd <- submd[1,c("tree","trt", "block")]
  maintrt <- rbind(maintrt, rowsubmd)
}

rownames(maintrt) <- maintrt$tree
mdstand <- decostand(mainorddat, "hel")
dsrda <- rda(mdstand~trt+Condition(block), maintrt)
plot(dsrda)
anova(dsrda, by="terms")
herbresp <- envfit(dsrda, mdstand)
selherb <- names(herbresp$vectors$pvals)[herbresp$vectors$pvals <= 0.05]
plot(herbresp)
plot(dsrda, type = "n")
text(dsrda, display="species", select =selherb)
