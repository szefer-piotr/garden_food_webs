# Summary for insects and plants

c1 <- main$LIFE.FORM %in% c("tree","shrub")

csites <- treats[treats$treat %in% c("CONTROL"), ]$codes
psites <- treats[treats$treat %in% c("PREDATOR"), ]$codes
csites <- as.character(csites)
psites <- as.character(psites)

c2 <- main$CODE %in% c(psites,csites)
sum(main[c1 & c2,]$WEIGHT)

main_la <- main %>%
  mutate(leaf_area = (LEAVES*1000 * SLA))

cmsq <- sum(main_la$leaf_area, na.rm = T)
cmsq/10000


sum(ins_bio[ins_bio$plot %in% c(psites,csites), ]$totbio, na.rm=T)
