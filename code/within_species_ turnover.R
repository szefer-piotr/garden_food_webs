rm(list=ls())
source("code/data_processing_code.R")
source("code/pdi.R")
source("code/diet_shift.R")

library(betapart)
library(reshape2)
library(glmmTMB)
library(emmeans)
library(multcomp)
library(ggplot2)

# Go through all gardens and all treatment comparisons
comparisons <- list(c("control", "predator"),
                    c("control", "insecticide"),
                    c("control", "weevil125"),
                    c("control", "weevil25"))

ins_bio$block <- substr(ins_bio$plot,  3, 4)

treats$block <- substr(treats$codes,3,4)

blocks <- unique(ins_bio$block)

within_data <- list()

# not use
comp = 1
trts = toupper(comparisons[[comp]])
bl = "g1"

for(comp in 1:length(comparisons)){
  trts <- toupper(comparisons[[comp]])
  print(trts)
  
  ref_sites <- treats[treats$treat == trts[1], ]$codes
  trt_sites <- treats[treats$treat == trts[2], ]$codes
  
  # create a subset for a given block and  treatment
  for(bl in blocks){
    print(bl)
    
    # make boolean filters
    ref_site <- as.character(ref_sites[grep(bl, ref_sites)])
    trt_site <- as.character(trt_sites[grep(bl, trt_sites)])
    
    ref_plant <- ins_bioOrig[ins_bioOrig$plot == ref_site, ]
    trt_plant <- ins_bioOrig[ins_bioOrig$plot == trt_site, ]
    
    ref_plant_names <- as.character(unique(ref_plant$tree))
    trt_plant_names <- as.character(unique(trt_plant$tree))
    
    sel_names <- ref_plant_names[ref_plant_names %in% trt_plant_names]
    
    if(length(sel_names) == 0){
      print(paste("In", bl, paste(trts,collapse=" "), "have no plants in common"))
      next
    }
    
    ref_plant_filtered <- ref_plant[ref_plant$tree %in% sel_names, ]
    trt_plant_filtered <- trt_plant[trt_plant$tree %in% sel_names, ]
    
    ref_plant_filtered$treetrt <- paste(ref_plant_filtered$tree,
                                          trts[1], sep = "_")
    trt_plant_filtered$treetrt <- paste(trt_plant_filtered$tree,
                                        trts[2], sep = "_")
    
    joined_data <- rbind(ref_plant_filtered,
                         trt_plant_filtered)
    
    joined_data$block <- bl
    
    within_data[[paste(paste(comparisons[[comp]],collapse="_"),
                       bl,sep="_")]] <- joined_data
  }
}

# Extract biomass of a given plant at given plot
plantBio <- function(species, plot){
  plbiomass <- main_biomass[main_biomass$CODE == plot, ]
  plsitebio <- plbiomass[plbiomass$SP_CODE == toupper(species), ]
  bio <- (plsitebio$WEIGHT)
  if(length(bio) != 0){
    print(paste("Weight of:",species, "at the plot", plot, "is", bio))
    return(bio)
  }else{
    print(paste("There is no",species, "at the plot", plot))
    return(0)
  }
}

# append bio data to a list of datasets
appendBioData <- function(within_data){
  for(i in 1:length(within_data)){
    print(within_data[[i]])
    for(row in 1:dim(within_data[[i]])[1]){
      print(row)
      within_data[[i]]$treebio[row] <- plantBio(within_data[[i]]$tree[row],
                                           within_data[[i]]$plot[row])
    }
  }
  return(within_data)
}

within_data <- appendBioData(within_data)

# for(ind in 1:length(within_data)){
#   print(within_data[[ind]]$morphotype)
# }


# Make it into df
makeDistanceIntoDf <- function(turnover_dat_ip, type = "beta.bray"){
  # Possible types: "beta.bray.bal" "beta.bray.gra" "beta.bray"
  df <- na.omit(melt(as.matrix(turnover_dat_ip[[type]])))
  df <- df[order(df$Var1), ]   # ordering
  colnames(df) <- c("col", "row", "value") # setting colnames
  df <- df[!(df$col == df$row), ]
  df <- df[grep("CONTROL", df$col), ]
  
  df$a_plant <- substr(df$col, 1,6)
  df$b_plant <- substr(df$row, 1,6)
  
  # Between treatment comparisons
  df_between <- df[df$a_plant == df$b_plant, ]
  
  # Given plant species vs all the rest
  df_within <- df[!(df$a_plant == df$b_plant), ]
  
  ##### THIS SHOULD COMPARE SPECIES IN THE CONTROL WITH
  ##### OTHER SPECIES IN THE TREATMENT PLOT
  
  return(list(between = df_between,
              within = df_within))
}


extractBetaDiversityData <- function(within_data, 
                                     comparison = comparisons[[1]],
                                     beta_type = "beta.bray"){
  
  indices <- grep(paste(comparison, collapse = "_"), 
       names(within_data))
  
  beta_list <- list()
  
  for(ind in indices){
    
    joined_data <- within_data[[ind]]
    
    bl <- unique(joined_data$block)
    
    comp_dat <- contingencyTable2(joined_data, 
                                  "morphotype",
                                  "treetrt",
                                  "amount")
    
    ips <- grep("aran|mant", rownames(comp_dat))
    
    comp_dat_ip <- comp_dat[ips, ]
    comp_dat_noip <- comp_dat[-ips, ]
    
    turnover_dat_ip <- beta.pair.abund(t(comp_dat_ip))
    turnover_dat_noip <- beta.pair.abund(t(comp_dat_noip))
    
    turn_dat_ip_df <- makeDistanceIntoDf(turnover_dat_ip,
                                         type = beta_type)
    turn_dat_noip_df <- makeDistanceIntoDf(turnover_dat_noip, 
                                           type = beta_type)
    indlist <- list(noip = turn_dat_noip_df,
                ip = turn_dat_ip_df)
    
    beta_list[[bl]] <- indlist
  }
  
  beta_list_nonull <- beta_list[!sapply(beta_list,is.null)]
  
  # Attach plant bio, check theh number of the comparison
  # noip between
  for(nm in names(beta_list_nonull)){
    print(nm)
    for(row in 1:dim(beta_list_nonull[[nm]]$noip$between)[1]){
      print(row)
      plant_name <- beta_list_nonull[[nm]]$noip$between[row, ]$a_plant
      plant_name
      
      ref_bool <- (treats$block == nm & treats$treat == toupper(comparison[1]))
      trt_bool <- (treats$block == nm & treats$treat == toupper(comparison[2]))
      
      ref_plot <- as.character(treats[ref_bool,]$codes)
      trt_plot <- as.character(treats[trt_bool,]$codes)
      
      ref_plant_bio <- plantBio(plant_name, ref_plot)
      trt_plant_bio <- plantBio(plant_name, trt_plot)
      
      beta_list_nonull[[nm]]$noip$between$ref_pl_bio <- ref_plant_bio
      beta_list_nonull[[nm]]$noip$between$trt_pl_bio <- trt_plant_bio
      
      beta_list_nonull[[nm]]$noip$between$block <- nm
    }
  }
  
  # ip between
  for(nm in names(beta_list_nonull)){
    print(nm)
    for(row in 1:dim(beta_list_nonull[[nm]]$ip$between)[1]){
      print(row)
      plant_name <- beta_list_nonull[[nm]]$ip$between[row, ]$a_plant
      plant_name
      
      ref_bool <- (treats$block == nm & treats$treat == toupper(comparison[1]))
      trt_bool <- (treats$block == nm & treats$treat == toupper(comparison[2]))
      
      ref_plot <- as.character(treats[ref_bool,]$codes)
      trt_plot <- as.character(treats[trt_bool,]$codes)
      
      ref_plant_bio <- plantBio(plant_name, ref_plot)
      trt_plant_bio <- plantBio(plant_name, trt_plot)
      
      beta_list_nonull[[nm]]$ip$between$ref_pl_bio <- ref_plant_bio
      beta_list_nonull[[nm]]$ip$between$trt_pl_bio <- trt_plant_bio
      
      beta_list_nonull[[nm]]$ip$between$block <- nm
    }
  }
  
  # noip within
  for(nm in names(beta_list_nonull)){
    print(nm)
    for(row in 1:dim(beta_list_nonull[[nm]]$noip$within)[1]){
      print("noip within")
      # plant_name <- beta_list_nonull[[nm]]$ip$within[row, ]$a_plant
      # plant_name
      # 
      # ref_bool <- (treats$block == nm & treats$treat == toupper(comparison[1]))
      # trt_bool <- (treats$block == nm & treats$treat == toupper(comparison[2]))
      # 
      # ref_plot <- as.character(treats[ref_bool,]$codes)
      # trt_plot <- as.character(treats[trt_bool,]$codes)
      # 
      # ref_plant_bio <- plantBio(plant_name, ref_plot)
      # trt_plant_bio <- plantBio(plant_name, trt_plot)
      # 
      # beta_list_nonull[[nm]]$ip$within$ref_pl_bio <- ref_plant_bio
      # beta_list_nonull[[nm]]$ip$within$trt_pl_bio <- trt_plant_bio
      
      if(dim(beta_list_nonull[[nm]]$noip$within)[1] == 0){
        next
      }
      
      beta_list_nonull[[nm]]$noip$within$block <- nm
    }
  }
  
  # ip within
  for(nm in names(beta_list_nonull)){
    print(nm)
    for(row in 1:dim(beta_list_nonull[[nm]]$ip$within)[1]){
      print("ip within")
      # plant_name <- beta_list_nonull[[nm]]$ip$within[row, ]$a_plant
      # plant_name
      # 
      # ref_bool <- (treats$block == nm & treats$treat == toupper(comparison[1]))
      # trt_bool <- (treats$block == nm & treats$treat == toupper(comparison[2]))
      # 
      # ref_plot <- as.character(treats[ref_bool,]$codes)
      # trt_plot <- as.character(treats[trt_bool,]$codes)
      # 
      # ref_plant_bio <- plantBio(plant_name, ref_plot)
      # trt_plant_bio <- plantBio(plant_name, trt_plot)
      # 
      # beta_list_nonull[[nm]]$ip$within$ref_pl_bio <- ref_plant_bio
      # beta_list_nonull[[nm]]$ip$within$trt_pl_bio <- trt_plant_bio
      
      if(dim(beta_list_nonull[[nm]]$ip$within)[1] == 0){
        next
      }
      
      beta_list_nonull[[nm]]$ip$within$block <- nm
    }
  }
  
  return(beta_list_nonull)
}


# Make dataset for a plot from a list
makeDataFrameFromDiversityList <- function(beta_data, 
                                           dataset = "noip", 
                                           type = "between"){
  
  data <- data.frame()
  
  for(nm in names(beta_data)){
    data_row <- beta_data[[nm]][dataset][[1]][type][[1]]
    data <- rbind(data,data_row)
  }
  
  return(data)
}


comparisons <- list(c("control", "predator"),
                    c("control", "insecticide"),
                    c("control", "weevil125"),
                    c("control", "weevil25"))

b_type <- "beta.bray"


beta_data <- beta_data_w125

# removeWithinControl <- function(beta_data){
#   
#   # within dfs
#   
#   # for both dataset = "noip" 
#   type = "within"
#   for (nm in names(beta_data)){
#     print(nm)
#     dset <- beta_data[[nm]][][[1]][type][[1]]
#     
#     dsetcontr <- dset[-grep("CONTROL", dset$row), ]
#     
#     beta_data[[nm]][][[1]][type][[1]] <- dsetcontr
#     
#   }
#   return(beta_data)
# }

beta_data_w125 <- extractBetaDiversityData(within_data, 
                                      comparison = comparisons[[3]],
                                      beta_type = b_type)
beta_data_w125_noctr <- removeWithinControl(beta_data_w125)

beta_data_pred <- extractBetaDiversityData(within_data, 
                                           comparison = comparisons[[1]],
                                           beta_type = b_type)
beta_data_ins <- extractBetaDiversityData(within_data, 
                                           comparison = comparisons[[2]],
                                          beta_type = b_type)
beta_data_w25 <- extractBetaDiversityData(within_data, 
                                           comparison = comparisons[[4]],
                                          beta_type = b_type)

# I should switch the names: between is within species trt and within is between

dataPrep <- function(beta_data){
  plotDatBetween <- makeDataFrameFromDiversityList(beta_data) 
  plotDatWithin <- makeDataFrameFromDiversityList(beta_data, type = "within")
  
  # I change it here but it should be changed in the code above for clarity
  plotDatBetween$type <- "Within"
  plotDatWithin$type <- "Between"
  joinedPlotDat <- rbind(plotDatBetween[,colnames(plotDatWithin)],
                         plotDatWithin)
  return(joinedPlotDat)
}

dw125 <- dataPrep(beta_data_w125)
dw125$comparison <- "Control vs. Weevil 125"
dw125 <- dw125[-grep("CONTROL", dw125$row), ]

dw25 <- dataPrep(beta_data_w25)
dw25$comparison <- "Control vs. Weevil 25"
dw25 <- dw25[-grep("CONTROL", dw25$row), ]

dins <- dataPrep(beta_data_ins)
dins$comparison <- "Control vs. Insecticide"
dins <- dins[-grep("CONTROL", dins$row), ]

dpred <- dataPrep(beta_data_pred)
dpred$comparison <- "Control vs. Exclosure"
dpred <- dpred[-grep("CONTROL", dpred$row), ]

fullPlotDat <- rbind(dw125,
                     dw25,
                     dins,
                     dpred)

# Test

# Beta approximation
approx_beta_1 <- 0.99999999

# W 25 - significant
w25dset <- fullPlotDat[fullPlotDat$comparison == "Control vs. Weevil 25", ]
w25dset
w25_brand <- glmmTMB(value ~ type + (1|block), 
                  data = w25dset, 
                  family= beta_family(link = "logit"))
summary(w25_brand)

# W 125 - marginally significant/ when approximation of 1 is increased it was gone
w125dset <- fullPlotDat[fullPlotDat$comparison == "Control vs. Weevil 125", ]
w125dset[w125dset$value == 1, ]$value <- approx_beta_1
w125_brand <- glmmTMB(value ~ type + (1|block), 
                     data = w125dset, 
                     family= beta_family(link = "logit"))
summary(w125_brand)

# Ins - not significant
winsdset <- fullPlotDat[fullPlotDat$comparison == "Control vs. Insecticide", ]
winsdset[winsdset$value == 1, ]$value <- approx_beta_1
wins_brand <- glmmTMB(value ~ type + (1|block), 
                     data = winsdset, 
                     family= beta_family(link = "logit"))

summary(wins_brand)

# Pred - not significant
wpreddset <- fullPlotDat[fullPlotDat$comparison == "Control vs. Exclosure", ]
wpreddset[wpreddset$value == 1, ]$value <- approx_beta_1
wpred_brand <- glmmTMB(value ~ type + (1|block), 
                     data = wpreddset, 
                     family= beta_family(link = "logit"))
summary(wpred_brand)


# Within vs between beta diversity

clrs <- rgb(0,0,0,alpha = 50,maxColorValue = 255)
ggplot(fullPlotDat, aes(x = type, y=value))+
  geom_jitter(width = 0.1, size = 4,
              fill = clrs,
              stroke = 0,
              shape = 16,
              alpha =0.1)+
  stat_summary(fun.data=mean_cl_boot, 
               geom="pointrange", lwd=0.8,
               col = c("grey20","grey20",
                       "grey20","grey20",
                       "grey20","grey20",
                       "red","red")) +
  stat_summary(fun=mean, geom="point",cex = 2,
               col = c("grey20","grey20",
                       "grey20","grey20",
                       "grey20","grey20",
                       "red","red"))+
  theme_bw()+
  facet_wrap(~comparison)

ggplot(fullPlotDat[fullPlotDat$type == "Between", ], aes(x = a_plant, y= value))+
  geom_jitter()+
  stat_summary(fun.data=mean_cl_boot, 
               geom="pointrange", lwd=0.8,
               col = "red") +
  stat_summary(fun=mean, geom="point",cex = 2,
               col = "red")

ggplot(joinedPlotDat, aes(x = type, y=value,
                           label = type,
                           colour = block))+
  geom_jitter(width = 0.1, size = 4)+
  stat_summary(fun.data=mean_cl_boot, 
               geom="pointrange", lwd=0.8) +
  stat_summary(fun=mean, geom="point",cex = 2)+
  stat_summary(fun=mean, geom="text", hjust = 2,
               vjust = 1)
  

ggplot(joinedPlotDat, aes(x = type, y=value,
                          label = type,
                          colour = block))+
  geom_jitter(width = 0.1)+
  stat_summary(fun.data=mean_cl_boot, 
               geom="pointrange", lwd=0.8) +
  stat_summary(fun=mean, geom="point",cex = 2)+
  stat_summary(fun=mean, geom="text", hjust = 2,
               vjust = 1)+
  facet_wrap(~a_plant)

# Log ration might affect the beta diversity measures
ggplot(plotDatBetween, aes(x=log(trt_pl_bio/ref_pl_bio), 
                           y=value, color = block))+
         geom_point()

# But the general biomass tends to influence turnover rates
ggplot(plotDatBetween, aes(x=log(trt_pl_bio+ref_pl_bio), 
                           y=value, color = block))+
  geom_point()
