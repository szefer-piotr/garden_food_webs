sum(table(ins_bio$tree, ins_bio$plot)>0)

spsp <- toupper(unique(ins_bio$tree))

onlyTree <- main[main$SP_CODE %in% spsp,]
onlyTree$SP_CODE <- as.character(onlyTree$SP_CODE)

sampl_no <- table(onlyTree$SP_CODE)

sumdryw <- tapply(onlyTree$DRY, onlyTree$SP_CODE, sum, na.rm =T)
mindryw <- tapply(onlyTree$DRY, onlyTree$SP_CODE, min, na.rm =T)
maxdryw <- tapply(onlyTree$DRY, onlyTree$SP_CODE, max, na.rm =T)

min(mindryw)
max(maxdryw)
# 28 species in approximately 141 samples (dryed leaves)

sample_data <- cbind(sampl_no,sumdryw,mindryw,maxdryw)
