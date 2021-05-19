
#

library(dplyr)
library(ggplot2)


n = 10
a = 20

var.clinic <- 10
var.patient <- 15

mean.clinic <- 0
mean.patient <- 10

population.size <- n * a

clinic.effect <- rnorm(a, 0, sqrt(var.clinic))

mock.data <- data.frame(patient = rep(paste("patient", 
                                            1:10, 
                                            sep = ""),n),
                        clinic = rep(1:a, each = n),
                        clinic_eff = rep(clinic.effect, each = n))

mock.data <- mock.data %>%
  mutate(yij = rnorm(dim(mock.data)[1], 
                     mean.patient, 
                     sqrt(var.patient))+clinic_eff)

mock.data <- mock.data %>%
  mutate(seizures = rpois(dim(mock.data)[1], 
                     mean.patient+clinic_eff))


ggplot(mock.data, aes(x=clinic, y = seizures))+
  geom_jitter(width = 0.1)

# Poisson distribution
mu <- 10
linear_predictor <- rnorm()