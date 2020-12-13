rm(list=ls())
source("code/data_processing_code.R")

# 141 samples from which we collected insects
sum(table(ins_bio$tree, ins_bio$plot)>0)

# 27 tree species
spsp <- toupper(unique(ins_bio$tree))

# Withouth fungicide
sites <- unique(main[main$TREAT != "FUNGICIDE", ]$CODE)
sum(table(ins_bio[ins_bio$plot %in% sites, ]$tree,
          ins_bio[ins_bio$plot %in% sites, ]$plot)>0)


onlyTree <- main[main$SP_CODE %in% spsp,]
onlyTree$SP_CODE <- as.character(onlyTree$SP_CODE)

onlyTree <- onlyTree[onlyTree$TREAT != "FUNGICIDE", ]

sampl_no <- sum(table(onlyTree$SP_CODE))

sumdryw <- tapply(onlyTree$DRY, onlyTree$SP_CODE, sum, na.rm =T)
mindryw <- tapply(onlyTree$DRY, onlyTree$SP_CODE, min, na.rm =T)
maxdryw <- tapply(onlyTree$DRY, onlyTree$SP_CODE, max, na.rm =T)

min(mindryw)
max(maxdryw)
# 28 species in approximately 141 samples (dryed leaves)

sample_data <- as.data.frame(cbind(sampl_no,sumdryw,mindryw,maxdryw))
sum(sample_data$sampl_no)
