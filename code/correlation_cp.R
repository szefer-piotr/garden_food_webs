rm(list=ls())
source("code/data_processing_code.R")

ccodes <- as.character(treats[treats$treat %in% "CONTROL", ]$codes)
pcodes <- as.character(treats[treats$treat %in% "PREDATOR", ]$codes)

cpib <- ins_bio[ins_bio$plot %in% c(ccodes,pcodes), ]
hib <- cpib[-grep("aran|mant", cpib$morphotype), ]
iip <- cpib[grep("aran|mant", cpib$morphotype), ]


# 141 samples from which we collected insects
vals <- tapply(hib$amount, hib$plot, sum, na.rm =T)
valsnona <- vals[!is.na(vals)]

ivals <- tapply(iip$amount, iip$plot, sum, na.rm =T)
ivalsnona <- ivals[!is.na(ivals)]


pvals <- vals[names(vals) %in% pcodes]
cvals <- vals[names(vals) %in% ccodes]
ipvals <- ivals[names(ivals) %in% pcodes]
icvals <- ivals[names(ivals) %in% ccodes]

cordf <- cbind(pvals,cvals,ipvals,icvals)

pairs(cordf)


library(psych)
pairs.panels(cordf, 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
)
library(PerformanceAnalytics)
chart.Correlation(cordf, pch="+")

chart.Correlation(iris[-5], bg=iris$Species, pch=21)

# Both not significant
cor.test(cordf[,1], cordf[,3])
cor.test(cordf[,2], cordf[,4])

# Broken into oorders
matamo <- contingencyTable2(cpib, "plot", "family", "amount")
matbio <- contingencyTable2(cpib, "plot", "family", "bio")

rownames(treats) <- treats$codes

cptrt <- as.character(treats[rownames(matamo), ]$treat)

matamodf<-as.data.frame(matamo) 
matamodf$trt <- cptrt

matbiodf <- as.data.frame(matbio)
matbiodf$trt <- cptrt

library(corrplot)
corcont <- cor(matamodf[matamodf$trt == "CONTROL",-8], 
                    method = "pearson", 
                    use = "pairwise.complete.obs")
corcont_sig <- cor.mtest(matamodf[matamodf$trt == "CONTROL",-8], 
                          conf.level = .95,
                          use = "pairwise.complete.obs")
corrplot(corcont,
         diag=F,
         p.mat = corcont_sig$p, 
         insig = "label_sig",
         sig.level = c(.001, .01, .05), 
         pch.cex = .9, 
         pch.col = "white",
         type = "upper")
corpred <- cor(matamodf[matamodf$trt == "PREDATOR",-8], 
               method = "pearson", 
               use = "pairwise.complete.obs")
corpred_sig <- cor.mtest(matamodf[matamodf$trt == "PREDATOR",-8], 
                         conf.level = .95,
                         use = "pairwise.complete.obs")
corrplot(corpred,
         p.mat = corpred_sig$p, 
         insig = "label_sig",
         sig.level = c(.001, .01, .05), 
         pch.cex = .9, 
         pch.col = "white",
         type = "upper",
         diag = FALSE)

par(mfrow=c(1,2))



corcont <- cor(matbiodf[matbiodf$trt == "CONTROL",-8], 
               method = "pearson", 
               use = "pairwise.complete.obs")
corcont_sig <- cor.mtest(matbiodf[matbiodf$trt == "CONTROL",-8], 
                         conf.level = .95,
                         use = "pairwise.complete.obs")
corrplot(corcont,
         diag=F,
         p.mat = corcont_sig$p, 
         insig = "label_sig",
         sig.level = c(.001, .01, .05), 
         pch.cex = .9, 
         pch.col = "white",
         type = "upper")
corpred <- cor(matbiodf[matbiodf$trt == "PREDATOR",-8], 
               method = "pearson", 
               use = "pairwise.complete.obs")
corpred_sig <- cor.mtest(matbiodf[matbiodf$trt == "PREDATOR",-8], 
                         conf.level = .95,
                         use = "pairwise.complete.obs")
corrplot(corpred,
         p.mat = corpred_sig$p, 
         insig = "label_sig",
         sig.level = c(.001, .01, .05), 
         pch.cex = .9, 
         pch.col = "white",
         type = "upper",
         diag = FALSE)
