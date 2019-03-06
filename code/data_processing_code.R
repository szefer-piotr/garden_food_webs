# Contingency table function
contingencyTable2 <- function(dataset, ROW, COL, VALUE){
  # Get rid of the empty factors
  dataset[, colnames(dataset) == ROW] <- as.character(dataset[, colnames(dataset) == ROW])
  dataset[, colnames(dataset) == COL] <- as.character(dataset[, colnames(dataset) == COL])
  # Make a table, get rid of the empty rows and columns
  plants <- table(dataset[, colnames(dataset) == ROW], dataset[, colnames(dataset) == COL])
  plants <- plants[rowSums(plants) != 0, colSums(plants) != 0]
  # See where to insert values
  allSpecCodes <- colnames(plants)
  allPlotCodes <- rownames(plants)
  entries <- which(plants != 0, arr.ind = TRUE)
  # Loop through the entries and insert values
  for (entry in 1:dim(entries)[1]){
    plot <- entries[entry,1]
    plant <- entries[entry,2]
    specCode <- allSpecCodes[plant]
    plotCode <- allPlotCodes[plot]
    #res <- dataset[dataset$ROW == plotCode & dataset$COL == specCode,VALUE]
    res <- dataset[dataset[,ROW] == plotCode & dataset[,COL] == specCode,VALUE]
    # print(sum(res))
    plants[plot,plant] <- sum(res, na.rm = TRUE)
  }
  plants[is.na(plants)] <- 0
  
  # Change the table to a matrix or data.frame
  
  mat_a <- matrix(0, nrow = dim(plants)[1], ncol = dim(plants)[2])
  colnames(mat_a) <- colnames(plants)
  rownames(mat_a) <- rownames(plants)
  for (row in 1:dim(plants)[1]){
    for (col in 1:dim(plants)[2])
      mat_a[row,col] <- plants[row,col]
  }
  
  #return(plants)
  return(mat_a)
}

# herbivores <- read.csv("datasets/wng_insects/wng_bugs.csv",header=TRUE)
# 
# sort(unique(herbivores$code)) #3 plots are missing
# 
# dat13 <- herbivores[herbivores$code == "w1g1p3", ]
# 
# net13 <- contingencyTable2(dat13,"plant","morph","abu")
# 
# library(bipartite)
# plotweb2(net13)

measur <- read.csv("datasets/csv_measurments_all.csv")
arthro <- read.csv("datasets/csv_wng_all.csv")

arthro$family <- substr(arthro$morphotype, 1, 4)
measur$family <- substr(measur$morphotype, 1, 4)

table(arthro$family)
table(measur$family)

# Fix the some entries
measur$Size <- as.numeric(measur$Size)
names(measur) <- c("morphotype","no","size","scale","scl","notes")
measur$morphotype <- as.character(measur$morphotype)

# Cheeck the plat names
longnames <- which(sapply(as.character(measur$morphotype), nchar)>7)
chlist <- strsplit(as.character(measur[longnames, ]$morphotype), split=",")

for (i in 1:length(longnames)){
  row <- longnames[i]
  measur[row, 1] <- as.character(chlist[[i]][1])
  measur[row, 2] <- as.numeric(chlist[[i]][2])
  measur[row, 3] <- as.numeric(chlist[[i]][3])
  measur[row, 5] <- as.numeric(chlist[[i]][6])
}
measur <- measur[,c(1,2,3,5,6)]
measur$rsize <- measur$size * measur$scl
head(measur)

summary(arthro)
unique(arthro$morphotype) # 607 morphotypes
sort(unique(arthro$tree))

# Sida rhombifolia
arthro <- arthro[arthro$tree != "sida",]
dim(arthro)

# Change the names of trees
tochange <- data.frame(a = sort(unique(arthro$tree)),
           b=sort(unique(arthro$tree)))

# himibi??? "mimodi"
# arthro[arthro$tree == "premna",]
arthro[arthro$tree == "costsp",]
# arthro[arthro$plot == "w1g1p2",]
# 
# premna

tochange$b[1] <- "breyce"
tochange$b[5] <- "cordte"
tochange$b[10] <- "ficuco"
tochange$b[12] <- "mimodi"
tochange$b[14] <- "homano"
tochange$b[15] <- "homano"
tochange$b[25] <- "pipeum"
tochange$b[29] <- "prems1"
tochange$b[32] <- "solatu"
tochange$b[36] <- "tremor"
tochange$b[41] <- "viteco"

arthro$tree <- as.character(arthro$tree)

# Then use the corrected names to change names in the dataset
for(name in tochange$a){
  print(name)
  arthro[arthro$tree == name, ]$tree <- as.character(tochange[tochange$a == name, ]$b)
}

sort(unique(arthro$tree))

# Some notes:
# Check what Costus sp is.
# compare tree species with th eones observed in the previous study
# what ficuses are there?
# ARTOCO ARTOLA BARRS1 BREYCE CARIPA COMMBA DENDLO DENDS1 DRACLA ENDOLA
# FICUCO FICUCP FICUHA FICUHI FICUPA FICUVA FICUWA GUIOCO HOMANO LEEAIN
# MACAAL MACABI MACAFA MACAQU MACATA MANIES MELAMU MELOS1 MERRME MUSSCY
# PIPEAD PIPEUM PIPTAR PISOLO PREMOB PREMS1 TOURSA TREMOR TRICPL VITECO

# Plant biomass!!!
main <- read.table("datasets/wng_main_clean.txt", header = T)
# each plant in each garden has a biomass, so when i create my table, i could
# already use arthtopods biomass and supplement it with plant biomass at a given 
# plot

# Arthropod biomass
measur$fam <- substr(measur$morphotype, 1,4)
head(measur)

# Models used to estimate biomass
# Ganihar 1997
# Araneae: power;b0=-3.2105 (0.1075);b1=2.4681(0.0756)
# Orthoptera: power;b0=-3.5338(0.2668);b1=2.4619(0.1002)
# Hemiptera: power;b0=-3.8893(0.3387);b1=2.7642(0.3113)
# Homoptera: power;b0=-3.1984(0.1174);b1=2.3487(0.0779)
# Coleoptera: power;b0=-3.2689(0.0659);b1=2.4625(0.0415)

# Wardhough 
# (power model ln(weight) = ln(a) + b * length )
# Mantodea:   a=-6.34(0.72);b=3.01(0.27)
# Araneae:    a=-2.13(0.15);b=2.23(0.11)
# Orthoptera: a=-3.17(0.19);b=2.61(0.09)
# Hemiptera:  a=-3.01(0.17);b=2.59(0.09)
# Homoptera: 
# Coleoptera: a=-3.2(0.14); b=2.56(0.08)

params <- list(mant = c(a = -6.34, b = 3.01),
     aran = c(a = -2.13, b = 2.23))

est_bio <- function(fam, size){
  fam <- as.character(fam) #transform the input
  a <- params[[fam]]["a"]
  b <- params[[fam]]["b"]
  return(exp(a)*size^b) # weight=a*length^b
}

fam <- "aran"
size <- 40          # size has to be in [mm]
est_bio(fam, size) # what is the unit???? [mg]

# Average measurments for each morphotype,
# Or should I rather use median?
head(measur)
avmed <- tapply(measur$rsize, measur$morphotype, median)
estims <- data.frame(av_morp =avmed)
estims$fam <- substr(rownames(estims),1,4)
for (row in 1:dim(estims)[1]){
  estims$weight[row] <- est_bio(estims[row,]$fam, 
                                estims[row,]$av_morp)
}
