# Insects analysis
# source("C:\\Users\\Piotr Szefer\\Desktop\\Work\\garden experiment\\code\\Data_Processing_Script_v2.R")
# setwd("C:\\Users\\Piotr Szefer\\Desktop\\Work\\garden experiment\\datasets\\wng_insects")
# insects <- read.table("csv_wng_arthropods.csv", sep=",", header=T, skip=2)
# insects$group <- substr(insects$Morphotype,1,4)

insects <- read.table("datasets/wng_arthro_clean.txt")
sizes <- read.table("datasets/wng_measurements.txt", header=T)
sizes$length <- sizes$rsize
sizes$group <- substr(sizes$morphotype,1,4)

head(insects)
head(sizes)

# Size dataset
size <- tapply(sizes$length, sizes$morphotype, mean)
group <- substr(names(size),1,4)
size_dat <- cbind(size,group)

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

# Test the equations for Mantodea
mant <- sizes[sizes$group == "mant", ]
aran <- sizes[sizes$group == "aran", ]
homo <- sizes[sizes$group == "homo", ]
hemi <- sizes[sizes$group == "hemi", ]
cole <- sizes[sizes$group == "cole", ]

aran$morphotype <- as.character(aran$morphotype)
mant$morphotype <- as.character(mant$morphotype)
homo$morphotype <- as.character(homo$morphotype)
hemi$morphotype <- as.character(hemi$morphotype)
cole$morphotype <- as.character(cole$morphotype)

aran_size <- tapply(aran$length, aran$morphotype, mean)
mant_size <- tapply(mant$length, mant$morphotype, mean)
homo_size <- tapply(homo$length, homo$morphotype, mean)
hemi_size <- tapply(hemi$length, hemi$morphotype, mean)
cole_size <- tapply(cole$length, cole$morphotype, mean)

# T0 estimate the body size use equations on the individuals!
mant_ind <- insects[insects$group == "mant", ]
extra_row <- mant_ind[1,]
extra_row$Plot <- "w1g4p1"
extra_row$Amount <- 0
mant_ind <- rbind(mant_ind, extra_row)
mant_ind[which(mant_ind$Plot == "wg3p6"),]$Plot <- "w1g3p6"
mant_ind$Plot <- as.character(mant_ind$Plot)
mant_ind$bio <- 0
mant_bio <- exp(-6.34)*mant_size^3.01
for (morph in names(mant_bio)){
  mant_ind[mant_ind$Morphotype==morph,]$bio <- mant_bio[which(names(mant_bio) == morph)]
}
# calculate biomass
mant_ind$est_bio <- mant_ind$Amount * mant_ind$bio

# Prepare the dataset with gardens and plots
mant_plot <- tapply(mant_ind$est_bio,mant_ind$Plot,sum, na.rm=TRUE)
mant_gard <- names(mant_plot)
mant_data <- data.frame(code = mant_gard, bio = mant_plot)
mant_data$block <- substr(mant_data$code, 3,4)
mant_data <- cbind(mant_data, WNGtreat[-c(6),c(3,4)])
mant_data[18,]$bio <- 0.01

# Logged values of total biomass
mant_bio <- ggplot(mant_data, aes(x = TREATMENT, y=log(bio), group=block))+
  geom_line(aes(linetype = block), size = 0.5, alpha=0.5)+
  geom_jitter(width = 0.1, size=1.5)
mant_bio



# # weevil addition 125 "destabilized" the community of arachnids the most
# 

# Araneae biomass
aran_ind <- insects[insects$group == "aran", ]
aran_ind$Plot <- as.character(aran_ind$Plot)
aran_ind$bio <- 0
aran_bio <- exp(-2.13)*aran_size^2.23
for (morph in names(aran_bio)){
  aran_ind[aran_ind$Morphotype==morph,]$bio <- aran_bio[which(names(aran_bio) == morph)]
}

# calculate biomass
aran_ind$est_bio <- aran_ind$Amount * aran_ind$bio

# Prepare the dataset with gardens and plots
aran_plot <- tapply(aran_ind$est_bio,aran_ind$Plot,sum, na.rm=TRUE)
aran_gard <- names(aran_plot)
aran_data <- data.frame(code = aran_gard, bio = aran_plot)
aran_data$block <- substr(aran_data$code, 3,4)
aran_data <- cbind(aran_data, WNGtreat[-6,c(3,4)])

# Logged values of total biomass
aran_bio <- ggplot(aran_data, aes(x = TREATMENT, y=log(bio), group=block))+
  geom_line(aes(linetype = block), size = 0.5, alpha=0.5)+
  geom_jitter(width = 0.1, size=1.5)
aran_bio



# weevil addition 125 destabilized the community of arachnids the most

# Investigate the zero values 
# No mantoids in w1g4p1 plot biomass value is 0,replaced with 0.01

predator <- cbind(mant_data, aran_data$bio)
predator$sum <- predator$bio+predator$'aran_data$bio'

pred_bio <- ggplot(predator, aes(x = TREATMENT, y=log(sum), group=block))+
  geom_line(aes(linetype = block), size = 0.5, alpha=0.5)+
  geom_jitter(width = 0.1, size=1.5)
pred_bio


# Herbivores
# Homopterans biomass
homo_ind <- insects[insects$group == "homo", ]
homo_ind$Plot <- as.character(homo_ind$Plot)
homo_ind$bio <- 0
homo_bio <- exp(-3.1984)*homo_size^2.3487 #b0=-3.1984(0.1174);b1=2.3487(0.0779)
for (morph in names(homo_bio)){
  homo_ind[homo_ind$Morphotype==morph,]$bio <- homo_bio[which(names(homo_bio) == morph)]
}

# calculate biomass
homo_ind$est_bio <- homo_ind$Amount * homo_ind$bio

# Prepare the dataset with gardens and plots
homo_plot <- tapply(homo_ind$est_bio,homo_ind$Plot,sum, na.rm=TRUE)
homo_gard <- names(homo_plot)
homo_data <- data.frame(code = homo_gard, bio = homo_plot)
homo_data$block <- substr(homo_data$code, 3,4)
homo_data <- cbind(homo_data, 
                   WNGtreat[rownames(WNGtreat) %in% toupper(gsub("w1","w",
                                                                 homo_data$code)),
                            c(3,4)])

# Logged values of total biomass
homo_bio <- ggplot(homo_data, aes(x = TREATMENT, y=log(bio), group=block))+
  geom_line(aes(linetype = block), size = 0.5, alpha=0.5)+
  geom_jitter(width = 0.1, size=1.5)
homo_bio



# Hemiptera
hemi_ind <- insects[insects$group == "hemi", ]
hemi_ind$Plot <- as.character(hemi_ind$Plot)
hemi_ind$Morphotype <- as.character(hemi_ind$Morphotype)
hemi_ind[hemi_ind$Morphotype == "hemi29",]$Morphotype <- "hemi029" 
hemi_ind$bio <- 0
hemi_bio <- exp(-3.01)*hemi_size^2.59 #-3.01(0.17);b=2.59(0.09)
for (morph in names(hemi_bio)){
  if(dim(hemi_ind[hemi_ind$Morphotype == morph,])[1] != 0){
    hemi_ind[hemi_ind$Morphotype==morph,]$bio <- hemi_bio[which(names(hemi_bio) == morph)]
    print("done")
  }
}

# calculate biomass
hemi_ind$est_bio <- hemi_ind$Amount * hemi_ind$bio

# Prepare the dataset with gardens and plots
hemi_plot <- tapply(hemi_ind$est_bio,hemi_ind$Plot,sum, na.rm=TRUE)
hemi_gard <- names(hemi_plot)
hemi_data <- data.frame(code = hemi_gard, bio = hemi_plot)
hemi_data$block <- substr(hemi_data$code, 3,4)
hemi_data <- cbind(hemi_data, 
                   WNGtreat[rownames(WNGtreat) %in% toupper(gsub("w1","w",
                                                                 hemi_data$code)),
                            c(3,4)])

# Logged values of total biomass
hemi_bio <- ggplot(hemi_data, aes(x = TREATMENT, y=log(bio), group=block))+
  geom_line(aes(linetype = block), size = 0.5, alpha=0.5)+
  geom_jitter(width = 0.1, size=1.5)
hemi_bio

summary(lmer(log(bio)~TREATMENT+(1|block),data=hemi_data))

#################### Coleoptera ##########################################
# Coleoptera
cole_ind <- insects[insects$group == "cole", ]
cole_ind$Plot <- as.character(cole_ind$Plot)
cole_ind$Morphotype <- as.character(cole_ind$Morphotype)
cole_ind$bio <- 0
cole_bio <- exp(-3.2)*cole_size^2.56 #a=-3.2(0.14); b=2.56(0.08)
for (morph in names(cole_bio)){
  if(dim(cole_ind[cole_ind$Morphotype == morph,])[1] != 0){
    cole_ind[cole_ind$Morphotype==morph,]$bio <- cole_bio[which(names(cole_bio) == morph)]
    print("done")
  }
}

# calculate biomass
cole_ind$est_bio <- cole_ind$Amount * cole_ind$bio

# Prepare the dataset with gardens and plots
cole_plot <- tapply(cole_ind$Amount,cole_ind$Plot,sum, na.rm=TRUE)
cole_gard <- names(cole_plot)
cole_data <- data.frame(code = cole_gard, bio = cole_plot)
cole_data$block <- substr(cole_data$code, 3,4)
cole_data <- cbind(cole_data, 
                   WNGtreat[rownames(WNGtreat) %in% toupper(gsub("w1","w",
                                                                 cole_data$code)),
                            c(3,4)])

# Logged values of total biomass
cole_bio <- ggplot(cole_data, aes(x = TREATMENT, y=bio, group=block))+
  geom_line(aes(linetype = block), size = 0.5, alpha=0.5)+
  geom_jitter(width = 0.1, size=1.5)
cole_bio

summary(lmer(log(bio)~TREATMENT+(1|block),data=cole_data))



##########################################################################


# Cumulative (need to oknow which are zeros and which are missing. Seems like only g1p6 is missing)
# g2p6 and g6p6 exist, missing is g1p6 
dim(hemi_data)
dim(homo_data)

hemi_data$code <- as.character(hemi_data$code)
homo_data$code <- as.character(homo_data$code)

herbivores <- vector()
for(row in 1:dim(hemi_data)[1]){
  plot <- hemi_data[row,]$code
  print(plot)
  he <- as.numeric(hemi_data[hemi_data==plot,]$bio)
  ho <- as.numeric(homo_data[homo_data==plot,]$bio)
  print(c(he,ho))
  print(he+ho)
  if(dim(as.matrix(ho))[1] == 0){
    herbivores <- c(herbivores, he)
  }
  herbivores <- c(herbivores, (he+ho))
}

hemi_data$cumulative <- herbivores
herb_bio <- ggplot(hemi_data, aes(x = TREATMENT, y=log(cumulative), group=block))+
  geom_line(aes(linetype = block), size = 0.5, alpha=0.5)+
  geom_jitter(width = 0.1, size=1.5)
herb_bio
summary(lmer(log(cumulative)~TREATMENT+(1|block),data=hemi_data))

par(mfrow=c(1,1))

#windows(800,600)
#pdf("insects.pdf", height = 800,width=800)
require(cowplot)
plot_grid(aran_bio, mant_bio,
          pred_bio, hemi_bio,
          homo_bio,herb_bio,
          labels = c('Arachnids', 'Mantoids',
                     'All predators','Hemipterans',
                     'Homopterans','All herbivores'),las=2)
#dev.off()
predator$bio <- predator$sum
pred_data <- predator[,1:5]
herb_data <- hemi_data
herb_data$bio <- herbivores

herb_ggplot <- rbind(aran_data,mant_data,pred_data,
                     hemi_data[1:5],homo_data,herb_data[1:5])
herb_ggplot$plot <- rep(c("Aranea","Mantoidea",
                          "Predators","Hemiptera",
                          "Homoptera","Herbivores"), 
                        c(dim(aran_data)[1],dim(mant_data)[1],dim(pred_data)[1],
                          dim(hemi_data[1:5])[1],dim(homo_data)[1],dim(herb_data[1:5])[1]))

p <- ggplot(herb_ggplot, aes(x = TREATMENT, y = log(bio), group=block)) + 
  geom_line(aes(linetype = block), size = 0.5, alpha=0.5)+
  geom_jitter(width = 0.1, size=1.5)
p + facet_wrap(~plot, scales="free", drop=TRUE) +
  theme_bw() + 
  theme(axis.text.x=element_text(angle=90, hjust=1))

# Something is not ok in the garden 6 control plot. There are only some
# hemipterans there, not... i dont know... seems unprobable...
#herb_ggplot[herb_ggplot$plot == "Herbivores" & herb_ggplot$TREATMENT == "CONTROL",]
herb_ggplot <- herb_ggplot[!(herb_ggplot$plot == "Herbivores" & herb_ggplot$TREATMENT == "CONTROL"), ]

# Test for herbivores and predators (trophic levels)
Herb <- herb_ggplot[herb_ggplot$plot == "Herbivores",]

herb_plot <- ggplot(Herb, aes(x = TREATMENT, y = log(bio), group=block)) + 
  geom_line(aes(linetype = block), size = 0.5, alpha=0.5)+
  geom_jitter(width = 0.1, size=1.5)
herb_plot

lmer_herb <- lmer(log(bio)~TREATMENT + (1|block), data=Herb)
summary(lmer_herb)

Pred <- herb_ggplot[herb_ggplot$plot == "Predators",]

pred_plot <- ggplot(Pred, aes(x = TREATMENT, y = log(bio), group=block)) + 
  geom_line(aes(linetype = block), size = 0.5, alpha=0.5)+
  geom_jitter(width = 0.1, size=1.5)
pred_plot

lmer_pred <- lmer(log(bio)~TREATMENT + (1|block), data=Pred)
summary(lmer_pred)

kable(summary(lmer(log(bio)~TREATMENT+(1|block),data=aran_data))$coef)
summary(lmer(log(bio)~TREATMENT+(1|block),data=mant_data))
summary(lmer(log(sum)~TREATMENT+(1|block),data=predator))
summary(lmer(log(bio)~TREATMENT+(1|block),data=homo_data))
summary(lmer(log(bio)~TREATMENT+(1|block),data=hemi_data))
summary(lmer(log(cumulative)~TREATMENT+(1|block),data=hemi_data))

plot(lmer(log(bio)~TREATMENT+(1|block),data=mant_data))
plot(lmer(log(sum)~TREATMENT+(1|block),data=predator))
plot(lmer(log(bio)~TREATMENT+(1|block),data=homo_data))
plot(lmer(log(bio)~TREATMENT+(1|block),data=hemi_data))
plot(lmer(log(cumulative)~TREATMENT+(1|block),data=hemi_data))

# Is low herbivore load correlated with descriptors of plants?
bio_data <- fig1data[fig1data$TYPE == "Biomass",]

dim(herb_data)
x <- herb_data$bio

how <- order(bio_data$PLOT_CODE)
bio_data <- bio_data[how,]

dim(bio_data[-6,])
y <- bio_data[-6,]$VALUE

plot(log(herb_data$bio)~bio_data[-6,]$VALUE)
abline(lm(log(herb_data$bio)~bio_data[-6,]$VALUE))
summary(lm(log(herb_data$bio)~bio_data[-6,]$VALUE))

multi <- data.frame("herb" = log(herb_data$bio),
                    "bio" = bio_data[-6,]$VALUE)
multi$treat <- herb_data$TREATMENT
multi$gard <- herb_data$block

multi_reg <- ggplot(multi, aes(x =bio, y = herb, col=treat))
multi_reg + geom_point() + geom_smooth(aes(x = bio,y = herb, col=treat), 
                               method=lm, se=F)
glm_herb <- glm(herb~bio*treat, data=multi)
summary(glm_herb)
