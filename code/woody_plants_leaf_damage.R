rm(list=ls())
source("code/data_processing_code.R")

c1 <- main$LIFE.FORM %in% c("tree","shrub")

csites <- treats[treats$treat %in% c("CONTROL"), ]$codes
psites <- treats[treats$treat %in% c("PREDATOR"), ]$codes
csites <- as.character(csites)
psites <- as.character(psites)

c2 <- main$CODE %in% c(psites,csites)

main_dam <- main[c1 & c2, c("CODE",
                            "TREAT",
                            "BLOCK",
                            "SP_CODE",
                            "HERB",
                            "WATER",
                            "LDMC",
                            "SLA",
                            "WEIGHT")]

main_dam$SP_CODE <- as.character(main_dam$SP_CODE)

rowSums(table(main_dam$SP_CODE, main_dam$CODE))

library(ggplot2)

ggplot(main_dam, aes(x = TREAT, y=log(WEIGHT)))+
  geom_jitter(width = 0.01)+
  facet_wrap(~SP_CODE, scales = "free")

library(lme4)
library(lmerTest)

ggplot(main_dam, aes(x = TREAT, y=HERB))+
  geom_boxplot()

# This is not generally true.
lmer1 <- lmer(HERB~TREAT+(1|BLOCK), data = main_dam)
summary(lmer1)

# This is not generally true.
lmer1 <- lmer(HERB~TREAT+(1|BLOCK), data = main_dam)
summary(lmer1)

# For individual species

for(sp in unique(main_dam$SP_CODE)){
  print(sp)
  
  dset <- main_dam[main_dam$SP_CODE == sp, ]
  table(dset$TREAT) -> tbl
  
  if(dim(tbl) == 1){
    print("Single treatment")
    next
  }
  
  if(tbl[1] == 1 | tbl[2] == 1){
    print("Only one observation in one of the treatment")
    next
  }else{
    print("run model")
    # print(dset)
    lmer2 <- lm(HERB~TREAT, data = dset)
    print(summary(lmer2)$coefficients)
  }
}
