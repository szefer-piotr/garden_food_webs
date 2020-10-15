# Create food web diagrams

library(diagram)

par(mar = c(1, 1, 1, 1), mfrow = c(1, 1))
names <- c("Top predators", "Herbivores", "Arthropod predators", "Plants",
           "Fungi+", "Fungi-")

M <- matrix(nrow = 6, ncol = 6, byrow = TRUE, data = 0)

M[1, 2] <- 1
M[1, 3] <- 1 
M[2, 4] <- 1
M[3, 2] <- 1 
M[4, 5] <- 1
M[6, 4] <- 1

x11(8,8)
plotmat(M, pos = c(1, 2,1,2), curve = 0, name = names, lwd = 1,
        box.lwd = 2, cex.txt = 0, box.type = "square", box.prop = 0.5)
