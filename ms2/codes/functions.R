
TukeyGroups <- function(model, data, ...){
  # Get model summary and pairwise comparisons
  inter.test <- emmeans(model, "treat", 
                        data = data)
  print("inter.test passed")
  pairwise <- cld(inter.test, Letter="abcdefghijklm", type = "response")
  yvalue <- aggregate(.~treat, 
                      data=data[, c("treat","val")], 
                      mean)
  tdat <- merge(yvalue, pairwise[,c(2,5,6,7)])
  return(pairwise)
}

AppendGroups <- function(pairwise, data){
  # Extract letters and append to the dataset
  panelData <- data
  panelData$tukey <- "none"
  ltrs <- data.frame(pt = pairwise[, "treat"],
                     pg = pairwise[, ".group"])
  rownames(ltrs) <- ltrs$pt
  panelData$tukey <- ltrs[panelData$treat, ]$pg
  print("Data appended")
  return(panelData)
}

# Functions

# AnalyzeAndAppend <- function(data,
#                              Conditions, 
#                              FUN,
#                              ...){
#   
#   # Conditions <- gdp$guild == "herbivore" & gdp$ind == "bio"
#   
#   # Log link gaussian.
#   mod <- glmer(val~treat+(1|block),family = gaussian(link = "log"),
#                   data = data[Conditions, ],
#                   control=glmerControl(optimizer="bobyqa",
#                                        optCtrl=list(maxfun=2e5)))
#   summary(mod)
#   
#   tg <- TukeyGroups(mod, data[Conditions, ], ...)
#   data[Conditions, ] <- AppendGroups(tg, data[Conditions, ], ...)
#   
#   return(data)
# }
