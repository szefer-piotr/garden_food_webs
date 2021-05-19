# Colonization success adn predator action

library(MASS)
library(ggplot2)
library(EnvStats)

ord.n <- 150
ord.names <- c(paste("ord", 1:ord.n, sep ="_"))
set.seed(126)

# How many insects flow to a plot
influx <- rnegbin(ord.n, mu = 5000, theta = 0.5)

plot.n <- 1

# From previous power analysi
mubio <- 0.09570716 * 1000 # in kg
sdbio <- 1.9
plant.biomass <- rlnorm(plot.n, 
                        log(mubio/2), 
                        log(sdbio))

# Abundance - biomass relationship
# Plant biomass on the plots
vals <- 10^(-0.4783 + log10(plant.biomass*1000)*0.6694)

# Based on the plabt biomass relationship this calculates how many insects a plot can suport
ins.capacity <- c()
for (i in 1:plot.n){
  ins.capacity <- c(ins.capacity, rnegbin(1,mu = vals[i], theta = 10))
}

# For a plot summarising this, plot.n should have been higher
pldf <- as.data.frame(cbind(ins.capacity, plant.biomass))

# ggplot(pldf,
#        aes(y = ins.capacity, x=plant.biomass))+
#   geom_point()+geom_smooth(method = "lm")
#   scale_x_continuous(trans = 'log10') + 
#   scale_y_continuous(trans = 'log10')

# Density dependent influx
influx.com.const <- rep(ord.names, influx)
influx.com <- list()
for (ic in 1:plot.n){
  influx.com[[ic]] <- influx.com.const
}
# for each plot

# 1. Community sample at full capacity
# com1 <- sample(influx.com, pldf$ins.capacity[1], replace = F)

# Density independent
# com1 <- sample(influx.com, pldf$ins.capacity[1], replace = T)
# com1 <- sample(ord.names, pldf$ins.capacity[1], replace = T)

# 2. Predation

# Density independent - fixed diet
predator.prop <- rpareto(ord.n, location = 2,1) 
predator.diet <- predator.prop/sum(predator.prop)
names(predator.diet) <- ord.names
# What theoretically a predator would eat

updateCommunity <- function (ord.names,
                            init.com,
                            predator.fr, 
                            predator.diet,
                            ins.capacity,
                            density.dep = FALSE,
                            ...){
  
  # Create a list taht would contain updated communities
  
  upd.com.list <- list()
  # print("list created")
  
  for (ic in 1:length(ins.capacity)){
    
    # print(paste("plot no", ic))
    
    # Sample diet of predaotr
    if (!density.dep){
      # Density independent - fixed diet
      remcom1.di <- sample(ord.names, predator.fr*ins.capacity[ic], 
                         replace = T, prob = predator.diet)
      # print("sampled predator diet density indep.")
    }
    
    if (density.dep){
      # Density dependent - opportunistic
      remcom1.di <- sample(init.com[[ic]], ins.capacity[1]*predator.fr)
      # print("density dep. predator diet sampled")
    }
    
    # FORAGING - What would be left in herbivore community 
    left.com1.di <- vecsets::vsetdiff(init.com[[ic]], remcom1.di)
    # print("remainig community obtained")
    
    # IMIGRATION - Refill the community the influx community
    com1.sam <- sample(influx.com[[ic]], 
                       size = (sum(table(init.com[[ic]])) - sum(table(left.com1.di))), 
                       replace = F)
    # print("immigration")
    
    new.com <- c(left.com1.di,com1.sam)
    # print("new community computed")
    
    upd.com.list[[ic]] <- new.com
    # print("list updated")
  
  }
  
  return(upd.com.list) 
  
}


# Assume that influx is the same for all sites
# new.coms <- updateCommunity(ord.names,
#                             init.com = influx.com,
#                             predator.fr = 0.1,
#                             predator.diet,
#                             ins.capacity)
# 
# lapply(influx.com, table) # Constant influx
# lapply(new.coms, table)

# Loop n times through this
iterateTimes <- function(init.com,
                         steps = 100, ...){
  
  ord.df <- data.frame(as.list(lapply(init.com, table)[[1]]))
  
  pre.iter.com <- init.com
  
  for (t in 1:steps){
    print(t)
    post.iter.com <- updateCommunity(ord.names,
                                init.com = pre.iter.com,
                                predator.diet,
                                ins.capacity, ...)
    pre.iter.com <- post.iter.com
    
    ord.df.row <- data.frame(as.list(lapply(pre.iter.com, table)[[1]]))
    ord.df <- dplyr::bind_rows(ord.df, ord.df.row)
    
    # Calculate characteristics
    # Richness
    # Diversity
  
  }
  
  return(ord.df)
  
}

# Strong predation density independent feeding
res.com.di <- iterateTimes(influx.com, steps = 200, 
                        predator.fr = 0.5)

# Strong predation density dependent feeding
res.com.dd <- iterateTimes(influx.com, steps = 200, 
                           predator.fr = 0.5, density.dep = TRUE)

par(mfrow = c(1,2))
matplot(log(res.com.di), type = "l", lwd = 2, lty = 1)
matplot(log(res.com.dd), type = "l", lwd = 2, lty = 1)
