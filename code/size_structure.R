# Size structure change

source("code/data_processing_code.R")

psites <- as.character(treats[treats$treat %in% c("PREDATOR"), ]$codes)
csites <- as.character(treats[treats$treat %in% c("CONTROL"), ]$codes)
ins_bio$block <- substr(ins_bio$plot, 3, 4)
ibcp <- ins_bio[ins_bio$plot %in% c(psites,csites), ]

ibcp$mode <- "herbivore"
ibcp[grep("aran|mant", ibcp$morphotype),]$mode <- "predator"

ssdf <- data.frame()
for(family in unique(ins_bio$family)){
  for(bl in unique(ins_bio$block)){
    print(bl)
    sub <- ibcp[ibcp$family == family & ibcp$block == bl, ]
    subp <- sub[sub$plot %in% psites, ]
    subc <- sub[sub$plot %in% csites, ]
    
    subpdf <- rep(subp$bio, subp$amount)
    subcdf <- rep(subc$bio, subc$amount)
    
    # hist(log(subpdf), breaks = 25)
    # hist(log(subcdf), breaks = 25)
    
    # dp <- density(log(subpdf), na.rm = T, adjust = 6, kernel = "gaussian")
    # dc <- density(log(subcdf), na.rm = T, adjust = 6, kernel = "gaussian")
    # 
    # yl <- pmax(max(dp$y),max(dc$y))
    # xl <- c(pmin(min(dp$x),min(dc$x)), pmax(max(dp$x),max(dc$x)))
    # 
    # plot(dc,
    #      lty=2,xlim = xl, ylim = c(0,yl),
    #      lwd=2, col="gray80", main = paste(family))
    # par(new=TRUE)
    # plot(dp,
    #      lty=1,xlim = xl, ylim = c(0,yl),
    #      lwd=2, col="red", main="", 
    #      xlab = "", ylab = "")
    
    print(family)
    dat <- data.frame(lbio = c(log(subpdf),log(subcdf)),
                      type = rep(c("predator", "control"),
                                 c(length(subpdf),
                                   length(subcdf))),
                      block = bl,
                      fam = family)
    ssdf <- rbind(ssdf, dat)
  }
}

library(lmer)
library(lmerTest)

# We have change in the size structure towards larger herbivores!
# fm = "aran"
for(fams in unique(ssdf$fam)){
  print(fams)
  mod1 <- summary(glmer(exp(lbio)~type+(1|block), 
                       data=ssdf[ssdf$fam == fams,],
                       family = gaussian(link="log"), start = 0))
  # print(summary(mod1))
  sm <- summary(mod1)
  c1 <- sm$coefficients[1,1]
  c2 <- sm$coefficients[2,1]
  pval <- sm$coefficients[2,4]
  
  print(paste(fams, "average mass: ", round(exp(c1), 5)))
  print(paste(fams, "in exclosure mass: ", round(exp(c1 + c2),5)))
  
  if(pval<=0.05){
    if(round(exp(c1), 5) < round(exp(c1 + c2),5)){
      print("INCREASE")
    } else {
      print("DECREASE")
    }
  } else {
    print("NS")
  }
  
}


library(ggplot2)
library(RColorBrewer)

crs <- alpha(RColorBrewer::brewer.pal(6,"PiYG"),0.1)

ggplot(ssdf, aes(y=lbio, x = type, color = block))+
  geom_jitter(width=0.1)+
  scale_color_manual(values=crs)+
  stat_summary(fun.y = mean, geom = "point", col= "red")+
  stat_summary(fun.data = mean_se, #"mean_cl_boot"
                 geom = "errorbar",
                 width=0.05, col="red", lwd=1.1)+
  facet_wrap(~fam)

ggplot(ssdf, aes(x=lbio, fill=type))+
  geom_density(position = "identity", 
               adjust=2,
               alpha=0.3,
               lwd=0.5)+
  # scale_fill_manual(crs[c(1,2)])+
  # scale_color_manual(values=crs)+
  # stat_summary(fun.y = mean, geom = "point", col= "red")+
  # stat_summary(fun.data = mean_se, #"mean_cl_boot"
  #              geom = "errorbar",
  #              width=0.05, col="red", lwd=1.1)+
  facet_wrap(~fam)
# summary(lmer(lbio~type+(1|block/fam), data=ssdf)) # family nested within the block
IC <- mod1$coefficients[1,1]
C <- mod1$coefficients[2,1]
exp(IC) # average size of a herbivore
exp(IC+ C)# <- average size in predator exclosure

# Calculate the same for each family, especially IP
