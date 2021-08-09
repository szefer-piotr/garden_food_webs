woody_list <- toupper(c("melamu",
                "macaqu",
                "breyce",
                "premob",
                "ficuhi",
                "ficucp",
                "pipead",
                "homano",
                "macabi",
                "ficuwa",
                "macaal",
                "endola"))
woody_list

library(dplyr)

sum_dat <- data.frame()

ntd <- main %>%
  filter(SP_CODE %in% woody_list) %>%
  filter(TREAT %in% c("CONTROL", "INSECTICIDE"))


for(bl in unique(ntd$BLOCK)){
  print(bl)
  
  bl_dat <- ntd %>%
    filter(BLOCK == bl)
  
  i_dat <- bl_dat[bl_dat$TREAT == "INSECTICIDE",]
  c_dat <- bl_dat[bl_dat$TREAT == "CONTROL",]
  
  i_nm <- as.character(i_dat$SP_CODE)
  c_nm <- as.character(c_dat$SP_CODE)
  
  nmfiter <- i_nm[i_nm %in% c_nm]
  
  if(length(nmfiter)==0){
    next
  }

  fidat <- i_dat[i_dat$SP_CODE %in% nmfiter, ]
  fcdat <- c_dat[c_dat$SP_CODE %in% nmfiter, ]
  
  print(as.character(fidat$SP_CODE))
  print(as.character(fcdat$SP_CODE))
  
  data <- data.frame(block = bl, 
                     spec = nmfiter,
                     lratio = log(fcdat$WEIGHT/fidat$WEIGHT))
  
  print(data)  
  
  sum_dat <- rbind(sum_dat, data)
  
}

sum_dat
sum_dat$nc <- c(2.843, 4.147, 4.147,4.147,4.147, 2.947)

summary(lm(lratio~nc, data = sum_dat))

library(ggplot2)
library(ggpubr)

ggplot(sum_dat, aes(x = nc, y = lratio))+
  geom_point()+
  geom_smooth(method = "lm")+
  stat_cor(label.y = 4.5)
# Insects have more positive effect on plants with higher nitrogen content.


nitro <- read.csv("/home/piotrszefer/thesis/defence_presentation/data/sla_n.csv", header = T)


nitro
summary(nitro)
names(nitro) <- c("n","c","wd","ash","sla")

ggplot(nitro, aes(x = sla, y = n))+
  geom_point()+
  geom_smooth(method = "lm", lty = 2, se = F)+
  stat_cor(label.y = 4.5)
