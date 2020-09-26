rm(list=ls())
source("code/data_processing_code.R")
source("code/pdi.R")
source("code/diet_shift.R")
library(reshape2)

# Go through all gardens and all treatment comparisons
comparisons <- list(c("control", "predator"),
                    c("control", "insecticide"),
                    c("control", "weevil125"),
                    c("control", "weevil25"))

ins_bio$block <- substr(ins_bio$plot,  3, 4)

blocks <- unique(ins_bio$block)

within_data <- list()

# not use
# comp = 1
# trts = toupper(comparisons[[comp]])
# bl = "g1"

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
  
  return(list(between = df_between,
              within = df_within))
}


extractBetaDiversityData <- function(within_data, comparison = comparisons[[1]]){
  
  indices <- grep(paste(comparison, collapse = "_"), 
       names(within_data))
  
  beta_list <- list()
  
  for(ind in indices){
    joined_data <- within_data[[ind]]
    
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
                                         type = "beta.bray")
    turn_dat_noip_df <- makeDistanceIntoDf(turnover_dat_noip, 
                                           type = "beta.bray")
    indlist <- list(noip = turn_dat_noip_df,
                ip = turn_dat_ip_df)
    
    beta_list[[ind]] <- indlist
  }
  
  beta_list_nonull <- beta_list[!sapply(beta_list,is.null)]
  
  return(beta_list_nonull)
}

comparisons <- list(c("control", "predator"),
                    c("control", "insecticide"),
                    c("control", "weevil125"),
                    c("control", "weevil25"))

beta_data <- extractBetaDiversityData(within_data, 
                                      comparison = comparisons[[4]])


# Make dataset for a plot from a list
beta_data[[2]]$noip$between
