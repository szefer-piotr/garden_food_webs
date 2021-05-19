# Size structure change
rm(list=ls())
source("code/data_processing_code.R")

library(ggplot2)

psites <- as.character(treats[treats$treat %in% c("PREDATOR"), ]$codes)
csites <- as.character(treats[treats$treat %in% c("CONTROL"), ]$codes)

# Only adults?
# ins_bio <- ins_bio[ins_bio$adult.larvae != "L", ]
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

# Check whether number of individualt match the data

with(ssdf, {
  nind <- dim(ssdf[fam == "orth" & type == "predator" & block == "g5", ])[1]
  print(nind)
})

with(ins_bio, {
  sum(ins_bio[family == "orth" & (plot %in% psites) & block == "g5", ]$amount)
})

ann_text <- data.frame(mpg = 15,wt = 5,lab = "Text",
                       cyl = factor(8,levels = c("4","6","8")))

# Individual size distribution ----
ibcp 
ibcp.expanded <- ibcp[rep(row.names(ibcp), ibcp$amount), 
                      c(3,4,7,9,11)]
library(dplyr)

ibcp.expanded <- ibcp.expanded %>%
  mutate(Treatment = ifelse(plot %in% csites,
                            "C", "Ex"))

ndssdf <- ndssdf %>%
  mutate(Treatment = ifelse(type == "predator",
                            "Ex", "C"))

ord.labs <- c("Orthoptera", "Aranea","Homoptera","Heteroptera",
              "Mantodea", "Coleoptera","Lepidoptera")

levels(ibcp.expanded$family)
levels(ibcp.expanded$family) <- ord.labs[c(2,6,4,3,7,5,1)]

ibcp.expanded$family <- factor(ibcp.expanded$family, 
                               levels = ord.labs)
ibcp.expanded <- ibcp.expanded %>%
  mutate(lbio = log(bio))

sigcol <- "red"
nsigcol <- "black"
msigcol <- "gold"
ind.b.cols <- c(sigcol,sigcol,
                sigcol,sigcol,
                sigcol,sigcol,
                nsigcol,nsigcol,
                sigcol,sigcol,
                sigcol,sigcol,
                sigcol,sigcol)

sb.b.cols <- c(nsigcol,nsigcol,
                sigcol,sigcol,
                nsigcol,nsigcol,
                nsigcol,nsigcol,
                nsigcol,nsigcol,
                nsigcol,nsigcol,
                nsigcol,nsigcol)

# Individual based plots
mp1 <- ggplot(data = ibcp.expanded) +
  geom_pointrange(mapping = aes(y = log(bio), x = Treatment),
                  stat = "summary",
                  fun.min = function(z) {quantile(z,0.05)},
                  fun.max = function(z) {quantile(z,0.95)},
                  fun = median,
                  col = ind.b.cols)+
  facet_wrap(vars(family),ncol = 7)+ xlab("")+
  ylab("Log[individual body size]")
# 
# levels(ndssdf$fam) <- ord.labs
# 
# # Species based
# mp2 <- ggplot(data = ndssdf) +
#   geom_pointrange(mapping = aes(y=lbio, x=Treatment),
#                   stat = "summary",
#                   fun.min = function(z) {quantile(z,0.05)},
#                   fun.max = function(z) {quantile(z,0.95)},
#                   fun = median,
#                   col = sb.b.cols)+
#   facet_wrap(vars(fam),ncol = 7) + xlab("")+
#   ylab("Log[species body size]")
# 
# ggpubr::ggarrange(mp1, mp2, labels = c("A","B"), nrow = 2)

# U-Mann whitney test for families
# Hodges Lemann centrality instead of median
# Based on individuals and raw values

generalUW <-data.frame()
siteUW <- data.frame()
plantUW <- data.frame()
for (fam in unique(ibcp.expanded$family)){
  print(fam)
  subdat <- ibcp.expanded[ibcp.expanded$family == fam,]
  wt <- wilcox.test(subdat$bio~subdat$Treatment, conf.int = T)
  generalUW <- rbind(generalUW, data.frame(
    Order = fam,
    NC = nrow(subdat[subdat$Treatment == "C", ]),
    NEx = nrow(subdat[subdat$Treatment == "Ex", ]),
    HL = wt$estimate,
    MC = median(subdat[subdat$Treatment == "C", ]$lbio,na.rm=T),
    MEx = median(subdat[subdat$Treatment == "Ex", ]$lbio,na.rm=T),
    # HLLCL = wt$conf.int[[1]],
    # HLUCL = wt$conf.int[[2]],
    P = round(wt$p.value, 3)
  )) 
  
  for(st in unique(subdat$block)){
    print(st)
    sdst <-subdat[subdat$block == st, ]
    if(length(unique(sdst$Treatment)) < 2){
      print("no values in the other treatment")
      next
    }
    wt <- wilcox.test(sdst$bio~sdst$Treatment, conf.int = T)
    siteUW <- rbind(siteUW, data.frame(
      Order = fam,
      Site = st,
      NC = nrow(sdst[sdst$Treatment == "C", ]),
      NEx = nrow(sdst[sdst$Treatment == "Ex", ]),
      HL = wt$estimate,
      HLLCL = wt$conf.int[[1]],
      HLUCL = wt$conf.int[[2]],
      P = round(wt$p.value, 3)
    ))
  }
  
  for(pltnm in unique(subdat$tree)){
    print(pltnm)
    sdpl <-subdat[subdat$tree == pltnm, ]
    
    TT <- table(sdpl$Treatment)
    
    if(dim(TT) == 1 | TT[1]<2 | TT[2]<2 ){
      print("no species present in both treatment")
      next
    }
    wt <- wilcox.test(sdpl$bio~sdpl$Treatment, conf.int = T)
    plantUW <- rbind(plantUW, data.frame(
      Order = fam,
      Plant = pltnm,
      NC = nrow(sdpl[sdpl$Treatment == "C", ]),
      NEx = nrow(sdpl[sdpl$Treatment == "Ex", ]),
      HL = wt$estimate,
      HLLCL = wt$conf.int[[1]],
      HLUCL = wt$conf.int[[2]],
      P = round(wt$p.value, 3)
    ))
  }
}

generalUW <- generalUW %>% mutate(pred.eff = ifelse(P <= 0.05, 
                                       ifelse(HL < 0, 
                                              "decrease",
                                              "increase"), "ns"))
siteUW <- siteUW %>% mutate(pred.eff = ifelse(P <= 0.05, 
                                    ifelse(HL < 0, 
                                           "decrease",
                                           "increase"), "ns"))
plantUW <- plantUW %>% mutate(pred.eff = ifelse(P <= 0.05, 
                                     ifelse(HL < 0, 
                                            "decrease",
                                            "increase"), "ns"))


# Based on species

ndssdf$bio <- exp(ndssdf$lbio)

generalUWsb <-data.frame()
siteUWsb <- data.frame()
plantUWsb <- data.frame()

for (fm in unique(ndssdf$fam)){
  print(fm)
  subdat <- ndssdf[ndssdf$fam == fm,]
  wt <- wilcox.test(subdat$bio~subdat$Treatment, conf.int = T)
  generalUWsb <- rbind(generalUWsb, data.frame(
    Order = fm,
    NC = nrow(subdat[subdat$Treatment == "C", ]),
    NEx = nrow(subdat[subdat$Treatment == "Ex", ]),
    HL = wt$estimate,
    HLLCL = wt$conf.int[[1]],
    HLUCL = wt$conf.int[[2]],
    P = round(wt$p.value, 3)
  )) 
  
  for(st in unique(subdat$block)){
    print(st)
    sdst <-subdat[subdat$block == st, ]
    if(length(unique(sdst$Treatment)) < 2){
      print("no values in the other treatment")
      next
    }
    wt <- wilcox.test(sdst$bio~sdst$Treatment, conf.int = T)
    siteUWsb <- rbind(siteUWsb, data.frame(
      Order = fm,
      Site = st,
      NC = nrow(sdst[sdst$Treatment == "C", ]),
      NEx = nrow(sdst[sdst$Treatment == "Ex", ]),
      HL = wt$estimate,
      HLLCL = wt$conf.int[[1]],
      HLUCL = wt$conf.int[[2]],
      P = round(wt$p.value, 3)
    ))
  }
  
  for(pltnm in unique(subdat$plant)){
    print(pltnm)
    sdpl <-subdat[subdat$plant == pltnm, ]
    if(length(unique(sdpl$Treatment)) < 2){
      print("no species present in both treatment")
      next
    }
    wt <- wilcox.test(sdpl$bio~sdpl$Treatment, conf.int = T)
    plantUWsb <- rbind(plantUWsb, data.frame(
      Order = fm,
      Plant = pltnm,
      NC = nrow(sdpl[sdpl$Treatment == "C", ]),
      NEx = nrow(sdpl[sdpl$Treatment == "Ex", ]),
      HL = wt$estimate,
      HLLCL = wt$conf.int[[1]],
      HLUCL = wt$conf.int[[2]],
      P = round(wt$p.value, 3)
    ))
  }
}

generalUWsb <- generalUWsb %>% 
  mutate(pred.eff = ifelse(P <= 0.05,
                           ifelse(HL < 0,
                                  "decrease",
                                  "increase"), "ns"))
siteUWsb <- siteUWsb %>% 
  mutate(pred.eff = ifelse(P <= 0.05,
                           ifelse(HL < 0,
                                  "decrease",
                                  "increase"), "ns"))
plantUWsb <- plantUWsb %>% 
  mutate(pred.eff = ifelse(P <= 0.05,
                           ifelse(HL < 0,
                                  "decrease",
                                  "increase"), "ns"))

write.table(generalUW, "ms1/draft_3/tables/general_uw.txt")
write.table(siteUW, "ms1/draft_3/tables/site_uw.txt")
write.table(plantUW, "ms1/draft_3/tables/plant_uw.txt")

write.table(generalUWsb, "ms1/draft_3/tables/general_uw_sb.txt")
write.table(siteUWsb, "ms1/draft_3/tables/site_uw_sb.txt")
write.table(plantUWsb, "ms1/draft_3/tables/plant_uw_sb.txt")

