source("code/data_processing_code.R")

# Plants from C and P that are in both sites.

ins_bio_noip <- ins_bio[-grep("aran|mant", ins_bio$morphotype), ]
fullbio <- tapply(ins_bio$totbio , ins_bio$tree, sum, na.rm = T)
no_herb <- table(ins_bio_noip$tree, ins_bio_noip$morphotype)
attractivness <- rowSums(no_herb)

specperbio <- attractivness/fullbio

sort(specperbio)

cbind(attractivness,fullbio)
plot(log(attractivness)~log(fullbio))

main_biomass
main
