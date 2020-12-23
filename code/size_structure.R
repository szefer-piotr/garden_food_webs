# Size structure change

source("code/data_processing_code.R")

psites <- as.character(treats[treats$treat %in% c("PREDATOR"), ]$codes)
csites <- as.character(treats[treats$treat %in% c("CONTROL"), ]$codes)

# Only adults?
ins_bio <- ins_bio[ins_bio$adult.larvae != "L", ]
ins_bio$block <- substr(ins_bio$plot, 3, 4)
ibcp <- ins_bio[ins_bio$plot %in% c(psites,csites), ]

ibcp$mode <- "herbivore"
ibcp[grep("aran|mant", ibcp$morphotype),]$mode <- "predator"

ssdf <- data.frame()
# wtp <- c()
for(family in unique(ins_bio$family)){
  for(bl in unique(ins_bio$block)){

    print(bl)
    sub <- ibcp[ibcp$family == family & ibcp$block == bl, ]
    
    for(plsp in unique(sub$tree)){
      
      print(plsp)
      
      # DEBUG
      pdat <- sub[sub$tree == plsp, ]
      
      subp <- pdat[pdat$plot %in% psites, ]
      subc <- pdat[pdat$plot %in% csites, ]
      
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
      
      if(length(subpdf) != 0 & length(subpdf) != 0){
        print("enough data")
        dat <- data.frame(lbio = c(log(subpdf),log(subcdf)),
                          type = rep(c("predator", "control"),
                                     c(length(subpdf),
                                       length(subcdf))),
                          block = bl,
                          fam = family,
                          plant = plsp)
        ssdf <- rbind(ssdf, dat)
    }
    }
  }
}


# Species size dsitribution
ndssdf <- ssdf[!duplicated(ssdf),]


# W30  aran ficucp   434.5  77 0.001819995
# W31  aran costsp   434.5  77 0.001819995
# W32  aran pipead   434.5  77 0.001819995
# W33  aran macabi   434.5  77 0.001819995

ndssdf[ndssdf$family == "aran" & ndssdf$plant == "ficucp", ]

library(lmer)
library(lmerTest)

npar_res <- data.frame()
npar_tree_res <- data.frame()

# Ther is a change in the size structure towards larger herbivores!

for(fams in unique(ndssdf$fam)){
  print(fams)
  
  # predict data from the mod
  data <- ndssdf[ndssdf$fam == fams,]
  nddata <- data[!duplicated(data),]
  nddata$plant <- as.character(nddata$plant)

  for(sp in unique(nddata$plant)){
    print(sp)
    
    spnddat <- nddata[nddata$plant == sp, ]
    
    print(dim(spnddat)[1])
    print(head(spnddat))
    
    if(length(unique(spnddat$type)) == 1){
      next
    }else if(dim(spnddat)[1] <= 2){
      next
    }
    
    
    nptest <- wilcox.test(lbio~type, data = spnddat)
    kwtest <- kruskal.test(lbio~type, data = spnddat)
    pval = nptest$p.value
    
    npar_tree_res <- rbind(npar_tree_res,
                      data.frame(order = fams,
                                 plant = sp,
                                 W = nptest$statistic,
                                 N = dim(spnddat)[1],
                                 pval = nptest$p.value))
    
  }
  
  #Friedman Rrandomized block design test
  friedtest <- wilcox.test(lbio~type, data = nddata)
  kwtest <- kruskal.test(lbio~type, data = nddata)
  pval = friedtest$p.value
  
  npar_res <- rbind(npar_res,
                    data.frame(order = fams,
                               W = friedtest$statistic,
                               N = dim(nddata)[1],
                               pval = friedtest$p.value))
  
  # EVALUATE THE DIFFERENCE
  # cbio <- exp(mean(nddata[nddata$type == "control",]$lbio))
  # pbio <- exp(mean(nddata[nddata$type == "predator",]$lbio))
  # print(paste(fams, "average mass: ",round(cbio,5)))
  # print(paste(fams, "in exclosure mass: ", round(pbio,5)))
  # 
  # # Print some results out
  # if(pval<=0.05){
  #   if(cbio < pbio){
  #     print("INCREASE")
  #   } else {
  #     print("DECREASE")
  #   }
  # } else {
  #   print("NS")
  # }
}


npar_tree_res

names(npar_res) <- c("fam", "W","N","pval")
npar_res$label <- with(npar_res, paste("W = ", W,
                                       ", N = ", N,
                                       ", P-value = ", round(pval, 3), 
                                       sep = ""))
npar_res$type = "predator"




library(ggplot2)
library(RColorBrewer)

crs <- alpha(RColorBrewer::brewer.pal(6,"PiYG"),0.1)

# ndividuals - I think yes
# ggplot(ssdf, aes(y=lbio, x = type, color = block))+
#   geom_jitter(width=0.1)+
#   scale_color_manual(values=crs)+
#   stat_summary(fun.y = mean, geom = "point", col= "red")+
#   stat_summary(fun.data = mean_se, #"mean_cl_boot"
#                  geom = "errorbar",
#                  width=0.05, col="red", lwd=1.1)+
#   facet_wrap(~fam)

ggplot(ndssdf, aes(x=lbio, fill=type))+
  geom_density(position = "identity", 
               adjust=1,
               alpha=0.3,
               lwd=0.5)+
  facet_grid(vars(fam))+
  geom_text(data = npar_res, 
            mapping = aes(x=-Inf, y=-Inf, label = label),
            hjust = -1.3,
            vjust = -8)
  # scale_fill_manual(crs[c(1,2)])+
  # scale_color_manual(values=crs)+
  # stat_summary(fun.y = mean, geom = "point", col= "red")+
  # stat_summary(fun.data = mean_se, #"mean_cl_boot"
  #              geom = "errorbar",
  #              width=0.05, col="red", lwd=1.1)+
  
# summary(lmer(lbio~type+(1|block/fam), data=ssdf)) # family nested within the block
# IC <- mod1$coefficients[1,1]
# C <- mod1$coefficients[2,1]
# exp(IC) # average size of a herbivore
# exp(IC+ C)# <- average size in predator exclosure

# Check whether number of individualt match the data
with(ssdf, {
  nind <- dim(ssdf[fam == "orth" & type == "predator" & block == "g5", ])[1]
  print(nind)
})

with(ins_bio, {
  sum(ins_bio[family == "orth" & (plot %in% psites) & block == "g5", ]$amount)
})

# Violoin plots

ggplot(ssdf, aes(x=type, y=lbio))+
  geom_violin(scale="width")+
  # scale_fill_manual(crs[c(1,2)])+
  # scale_color_manual(values=crs)+
  # stat_summary(fun.y = mean, geom = "point", col= "red")+
  # stat_summary(fun.data = mean_se, #"mean_cl_boot"
  #              geom = "errorbar",
  #              width=0.05, col="red", lwd=1.1)+
  facet_wrap(~fam)

# Geom ridges
library(ggridges)
ggplot(ndssdf, aes(x = lbio, y = type, fill=type)) + 
  geom_density_ridges(quantile_lines = T, lty = 1,
                      quantiles = c(0.5), alpha=0.5, scale =3)+
  # facet_grid(vars(fam), vars(block))+
  facet_grid(vars(fam))+
  geom_text(data = npar_res, 
            mapping = aes(x=-Inf, y=-Inf, label = label),
            hjust = -1,
            vjust = -0.25)+
  theme_bw()+
  theme(legend.position = "none")



npar_res




ann_text <- data.frame(mpg = 15,wt = 5,lab = "Text",
                       cyl = factor(8,levels = c("4","6","8")))
