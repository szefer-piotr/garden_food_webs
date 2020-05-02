source("code/data_processing_code.R")

psites <- as.character(treats[treats$treat %in% c("PREDATOR"), ]$codes)
csites <- as.character(treats[treats$treat %in% c("CONTROL"), ]$codes)
ibcp <- ins_bio[ins_bio$plot %in% c(psites, csites), ]
ipcp_noip <- ibcp[-grep("aran|mant", ibcp$morphotype), ]


gardnets_noip <- gardnets

for(gard in names(gardnets){
  
}