# Script estimates biomasses of insects for each garden based on 
# available lenght-body mass equation for insects.

# 1. Load datasets (built using data_processing_code.R) ----
insects <- read.table("datasets/wng_arthro_clean.txt")
treats  <- read.table("datasets/treats_clean.txt")
plants  <- read.table("datasets/wng_main_bio.txt")
sizes <-  read.csv("datasets/corrected_measurements.csv", header=T)

plants$SP_CODE <- tolower(plants$SP_CODE)
insects$tree <- as.character(insects$tree)

# 2. Estimate insect sizes ----
# Size dataset, each morphotype's average size
size <- tapply(sizes$r_cm_size, sizes$morph, mean)
mft <- names(size)
group <- substr(names(size),1,4)
size_dat <- data.frame(morph = mft,
                       size = as.numeric(size),
                       group = group)

# Models used to estimate biomass
# Ganihar 1997
# Araneae: power;b0=-3.2105 (0.1075);b1=2.4681(0.0756)
# Orthoptera: power;b0=-3.5338(0.2668);b1=2.4619(0.1002)
# Hemiptera: power;b0=-3.8893(0.3387);b1=2.7642(0.3113)
# Homoptera: power;b0=-3.1984(0.1174);b1=2.3487(0.0779)
# Coleoptera: power;b0=-3.2689(0.0659);b1=2.4625(0.0415)

# Wardhough 
# (power model ln(weight) = ln(a) + b * length )
# weight <- exp(a) * size^b
# Mantodea:    a=-6.34(0.72);b=3.01(0.27)
# Araneae:     a=-2.13(0.15);b=2.23(0.11)
# Orthoptera:  a=-3.17(0.19);b=2.61(0.09)
# Hemiptera:   a=-3.01(0.17);b=2.59(0.09)
# Homoptera:   a=-3.20(0.12);b=2.35(0.08) (Ganihar 1997)
# Coleoptera:  a=-3.2(0.14); b=2.56(0.08)
# Lepidoptera: a=-5.44(); b=2.55

# Power model parameters for each family
allo_params <- data.frame(group = unique(size_dat$group),
                          a = c(-2.13,-3.2,-3.01,-3.20,-5.44,-6.34,-3.17),
                          b = c( 2.23,2.56, 2.59, 2.35, 2.55, 3.01, 2.61))

# Size of an average individual
size_dat$a <- 0
size_dat$b <- 0
for (grp in unique(allo_params$group)){
  aval <- allo_params[allo_params$group == grp, ]$a
  bval <- allo_params[allo_params$group == grp, ]$b
  size_dat[size_dat$group == grp, ]$a <- aval 
  size_dat[size_dat$group == grp, ]$b <- bval
}

size_dat$bio <- exp(size_dat$a) * (size_dat$size*10)^size_dat$b

size_dat[size_dat$morph == "lepi061",]

# exp(-3.58) * 31.95^2.42


write.table(size_dat, "datasets/size_dat_bio_corr.txt")

# Size distributions
library(ggplot2)
bio_dist <- ggplot(size_dat, aes(x = log(bio), fill=group))
bio_dist + geom_histogram(binwidth=0.5)


# Biomass based networks

# 3. Plant species names corrections ----

examine <- function(plt){
  #example plot
  sdp <- plants[plants$CODE == plt,]
  # see which plants are in the same 
  sdi <- insects[insects$plot == plt, ]
  sdp$SP_CODE
  unique(sdi$tree)
  # sdp[, c("WEIGHT","SP_CODE")]
  bool <- sum(unique(sdp$SP_CODE) %in% unique(sdi$tree)) == length(unique(sdi$tree))
  return(list(bool, sdp$SP_CODE, unique(sdi$tree)))
}

pcodes <- unique(plants$CODE)

plotnames <- c()
for(int in 1:length(pcodes)){
  ee <- examine(as.character(pcodes[int]))
  print(as.character(pcodes[int]))
  print(ee[[1]])
  if (ee[[1]] == FALSE){plotnames <- c(plotnames, as.character(pcodes[int]))}
}

# Names of plots where insect tree labels dont match ones from the plant dataset
plotnames

# w1g1p1
# what could this be?
examine("w1g1p1")
# insects[(insects$tree == "solatu" & insects$plot == "w1g1p1"), ]$tree # <- "melamu" 
# it has to be melamu the most abundant on the plot
insects[(insects$tree == "trems1" & insects$plot == "w1g1p1"), ]$tree  <- "tremor" 
#RESOLVED

# w1g1p2
# insects[insects$tree == "mimodi", ] #this doesn't appear in the plant dataset
# plants[plants$CODE == "w1g1p2", ]
# unique(plants$SPEC)
# 
# # Check if there are othher plots with this plant
# plants[plants$SP_CODE == "breyce", ]
# insects[insects$tree == "breyce", ]
# insects[insects$tree == "tricpl", ]

# mimodi is either breyce or tricpl, I can also remove insects and the plant,
# or i can assign these insects to a second most abundant species which is 
# tricpl (I will assign to tricpl)
examine("w1g1p2")
insects[(insects$tree == "mimodi" & insects$plot == "w1g1p2"), ]$tree <- "tricpl"
insects[(insects$tree == "trems1" & insects$plot == "w1g1p2"), ]$tree <- "tremor"
# RESOLVED

# w1g1p3
examine("w1g1p3")
insects[(insects$tree == "trems1" & insects$plot == "w1g1p3"), ]$tree <- "tremor"
# RESOLVED

# "w1g1p4"
examine("w1g1p4")
insects[(insects$tree == "trems1" & insects$plot == "w1g1p4"),  ]$tree <- "tremor"
# RESOLVED

# "w1g1p5" 
examine("w1g1p5")[[3]][!(examine("w1g1p5")[[3]] %in% examine("w1g1p5")[[2]])]
insects[(insects$tree == "trems1" & insects$plot == "w1g1p5"),  ]$tree <- "tremor"
# RESOLVED

# "w1g3p2" 
cd <- "w1g3p2"
examine(cd)
examine(cd)[[3]][!(examine(cd)[[3]] %in% examine(cd)[[2]])]
plants[plants$CODE == cd, ]
insects[insects$plot == cd, ]
insects[(insects$tree == "ficuco" & insects$plot == cd),  ]$tree <- "ficucp"
# RESOLVED

# "w1g3p3" 
cd <- "w1g3p3"
examine(cd)
examine(cd)[[3]][!(examine(cd)[[3]] %in% examine(cd)[[2]])]
insects[(insects$tree == "ficuco" & insects$plot == cd),  ]$tree <- "ficucp"
# RESOLVED

# "w1g3p4" 
cd <- "w1g3p4"
examine(cd)
examine(cd)[[3]][!(examine(cd)[[3]] %in% examine(cd)[[2]])]
insects[(insects$tree == "ficucg" & insects$plot == cd),  ]$tree <- "ficuco"
# RESOLVED

# "w1g3p6"
cd <- "w1g3p6"
examine(cd)
examine(cd)[[3]][!(examine(cd)[[3]] %in% examine(cd)[[2]])]
insects[(insects$tree == "ficuco" & insects$plot == cd),  ]$tree <- "ficupa"
# RESOLVED

# "w1g4p6"
cd <- "w1g4p6"
examine(cd)
examine(cd)[[3]][!(examine(cd)[[3]] %in% examine(cd)[[2]])]
insects <- insects[!(insects$tree == "karipa"),  ]
# RESOLVED

# "w1g5p2" 
cd <- "w1g5p2"
examine(cd)
examine(cd)[[3]][!(examine(cd)[[3]] %in% examine(cd)[[2]])]
insects[(insects$tree == "solatu" & insects$plot == cd),  ]$tree <- "solas1"
insects[(insects$tree == "trems1" & insects$plot == cd),  ]$tree <- "tremor"
# RESOLVED

# "w1g5p3"
cd <- "w1g5p3"
examine(cd)
examine(cd)[[3]][!(examine(cd)[[3]] %in% examine(cd)[[2]])]
insects[(insects$tree == "premob" & insects$plot == cd),  ]$tree <- "prems1"
# RESOLVED

# "w1g5p5" 
cd <- "w1g5p5"
examine(cd)
examine(cd)[[3]][!(examine(cd)[[3]] %in% examine(cd)[[2]])]
insects[(insects$tree == "ficucg" & insects$plot == cd),  ]$tree <- "ficuco"
# RESOLVED

# "w1g2p2"
cd <- "w1g2p2"
examine(cd)[[3]][!(examine(cd)[[3]] %in% examine(cd)[[2]])]
plants[plants$CODE == cd, ]
plants[plants$SP_CODE == "viteco",]
insects[insects$plot == cd, ]
# Figure out where these plants are
# plants[plants$SP_CODE ==  "viteco", ]$CODE # ok... there is only one place w1g2p6 
# plants[plants$SP_CODE ==  "macata", ]$CODE 
# # it is everywhere :/ but within the same garden it would be in w1g2p6
# 
# plants[plants$SP_CODE ==  "melamu", ]$CODE # might go there? -> w1g2p6
insects[(insects$tree %in% c("macata","viteco","melamu") & insects$plot == cd),  ]$plot <- "w1g2p6"
cd <- "w1g2p6"
examine(cd)
# RESOLVED

# "w1g6p4"
cd <- "w1g6p4"
examine(cd)
examine(cd)[[3]][!(examine(cd)[[3]] %in% examine(cd)[[2]])]
# how big of a deal is it?
insects <- insects[!(insects$tree %in% c("breyce", "piptar") & insects$plot == cd),  ] #only few ind
# I could remove it
# ??? breyce piptar
# RESOLVED

# "w1g2p3" 
cd <- "w1g2p3"
examine(cd)[[3]][!(examine(cd)[[3]] %in% examine(cd)[[2]])]
insects[(insects$tree == "trems1" & insects$plot == cd),  ]$tree <- "tremor"
plants[plants$CODE == cd, ]
insects[insects$plot == cd, ]
insects[(insects$tree %in% c("cordte") & insects$plot == cd),  ]
# definietly remove
insects <- insects[!(insects$tree == "cordte" & insects$plot == cd),  ]
# ??? cordte
# RESOLVED


write.table(insects, "datasets/arthropods_clean.txt")
write.table(treats,  "datasets/treatments_clean.txt")
write.table(plants,  "datasets/plants_clean.txt")
write.table(sizes,   "datasets/sizes_clean.txt")
