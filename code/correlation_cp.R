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
# pairs.panels(cordf, 
#              method = "pearson", # correlation method
#              hist.col = "#00AFBB",
#              density = TRUE,  # show density plots
#              ellipses = TRUE # show correlation ellipses
# )
library(PerformanceAnalytics)
# chart.Correlation(cordf, pch="+")

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
corcontabu <- cor(matamodf[matamodf$trt == "CONTROL",-8], 
                    method = "pearson", 
                    use = "pairwise.complete.obs")
corcontabu_sig <- cor.mtest(matamodf[matamodf$trt == "CONTROL",-8], 
                          conf.level = .95,
                          use = "pairwise.complete.obs")
# corrplot(corcont,
#          diag=F,
#          p.mat = corcont_sig$p, 
#          insig = "label_sig",
#          sig.level = c(.001, .01, .05), 
#          pch.cex = .9, 
#          pch.col = "white",
#          type = "upper")

corpredabu <- cor(matamodf[matamodf$trt == "PREDATOR",-8], 
               method = "pearson", 
               use = "pairwise.complete.obs")
corpredabu_sig <- cor.mtest(matamodf[matamodf$trt == "PREDATOR",-8], 
                         conf.level = .95,
                         use = "pairwise.complete.obs")
# corrplot(corpred,
#          p.mat = corpred_sig$p, 
#          insig = "label_sig",
#          sig.level = c(.001, .01, .05), 
#          pch.cex = .9, 
#          pch.col = "white",
#          type = "upper",
#          diag = FALSE)

# par(mfrow=c(1,1))



corcontbio <- cor(matbiodf[matbiodf$trt == "CONTROL",-8], 
               method = "pearson", 
               use = "pairwise.complete.obs")
corcontbio_sig <- cor.mtest(matbiodf[matbiodf$trt == "CONTROL",-8], 
                         conf.level = .95,
                         use = "pairwise.complete.obs")
# p3 <- corrplot(corcont,
#          diag=F,
#          p.mat = corcont_sig$p, 
#          insig = "label_sig",
#          sig.level = c(.001, .01, .05), 
#          pch.cex = .9, 
#          pch.col = "white",
#          type = "upper")

corpredbio <- cor(matbiodf[matbiodf$trt == "PREDATOR",-8], 
               method = "pearson", 
               use = "pairwise.complete.obs")
corpredbio_sig <- cor.mtest(matbiodf[matbiodf$trt == "PREDATOR",-8], 
                         conf.level = .95,
                         use = "pairwise.complete.obs")
# p4 <-corrplot(corpred,
#          p.mat = corpred_sig$p, 
#          insig = "label_sig",
#          sig.level = c(.001, .01, .05), 
#          pch.cex = .9, 
#          pch.col = "white",
#          type = "upper",
#          diag = FALSE)

library(ggcorrplot)

nms <- c("Aranea",
  "Coleoptera",
  "Heteroptera",
  "Homoptera",
  "Lepidoptera",
  "Mantodea",
  "Orthoptera",
  "trt")

colnames(matamodf) <- nms
colnames(matbiodf) <- nms

corcontabu <- cor(matamodf[matamodf$trt == "CONTROL",-8], 
                  method = "pearson", 
                  use = "pairwise.complete.obs")
corcontabu_sig <- cor_pmat(matamodf[matamodf$trt == "CONTROL",-8])
corpredabu <- cor(matamodf[matamodf$trt == "PREDATOR",-8], 
                  method = "pearson", 
                  use = "pairwise.complete.obs")
corpredabu_sig <- cor_pmat(matamodf[matamodf$trt == "PREDATOR",-8])
corcontbio <- cor(matbiodf[matbiodf$trt == "CONTROL",-8], 
                  method = "pearson", 
                  use = "pairwise.complete.obs")
corcontbio_sig <- cor_pmat(matbiodf[matbiodf$trt == "CONTROL",-8])
corpredbio <- cor(matbiodf[matbiodf$trt == "PREDATOR",-8], 
                  method = "pearson", 
                  use = "pairwise.complete.obs")
corpredbio_sig <- cor_pmat(matbiodf[matbiodf$trt == "PREDATOR",-8])



# outline.color = "white",
# ggtheme = ggplot2::theme_gray,
# colors = c("#F0E442", "#E69F00", "#56B4E9")

COLVEC <- c("#D55E00","#F0E442", "#56B4E9")
p1 <- ggcorrplot(corcontabu,
                 hc.order = F,
                 type = "upper",
                 p.mat = corcontabu_sig,
                 outline.color = "white",
                 ggtheme = ggplot2::theme_gray,
                 colors = COLVEC)
p2 <- ggcorrplot(corpredabu,
                 hc.order = F,
                 type = "upper",
                 p.mat = corpredabu_sig,
                 outline.color = "white",
                 ggtheme = ggplot2::theme_gray,
                 colors =  COLVEC)
p3 <- ggcorrplot(corcontbio,
                 hc.order = F,
                 type = "upper",
                 p.mat = corcontbio_sig,
                 outline.color = "white",
                 ggtheme = ggplot2::theme_gray,
                 colors =  COLVEC)
p4 <- ggcorrplot(corpredbio,
                 hc.order = F,
                 type = "upper",
                 p.mat = corpredbio_sig,
                 outline.color = "white",
                 ggtheme = ggplot2::theme_gray,
                 colors =  COLVEC)

library(ggpubr)
ggarrange(p1,p2,p3,p4,
          labels = c("A", "B", "C","D"),
          legend = "right", common.legend = T)
