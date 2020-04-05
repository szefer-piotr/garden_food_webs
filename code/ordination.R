# RDA

source("code/data_processing_code.R")

ins_bio

#gardnets
#gardnetsfam
#abugardnets
#abugardnetsfam


# All below is woong because it uses all plants fromo all treatments.
abumat <- contingencyTable2(ins_bio,"tree","morphotype","amount")
abumat <- decostand(abumat, "hel")

rda_abu <-rda(abumat) 
biplot(rda_abu)
text(rda_abu, display = "sites")
text(rda_abu, display = "species")

# Which insects best fit the rda model (no factors)
sp_fits <- envfit(rda_abu, abumat)

rda_bio <-rda()

# Co-inertia?

rdap # rda for plants
rdai # rda for insects