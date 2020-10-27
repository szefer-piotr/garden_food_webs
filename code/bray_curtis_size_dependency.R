
create_sp_names <- function(number){
  return(paste("sp", seq(1,number), sep=""))
}

number = 20

get_probabilities <- function(number, adjust = 5, noise = T){
  y <- rexp(100000,1/(number/adjust))
  hpr <- hist(y, breaks = number+5, freq = F)
  probs <- hpr$density[1:number]
  return(probs)
}

size = 20

commdf <- data.frame(spec=create_sp_names(size),
                     probs = get_probabilities(size))

commOfSize <- function(n, noise = 1){
  
  d=commdf$probs
  logit=function(x){log(x/(1-x))}
  ld <- logit(d)
  ns <- rnorm(length(d),mean=0,sd=noise)
  ldns <- ld + ns
  de=exp(ldns)/(1+exp(ldns))
  noise_probs <- de/sum(de)
  
  # plot(commdf$probs, pch = 20, col="blue")
  # par(new=T)
  # plot(noise_probs, pch=19, col = "red")
  
  sample(commdf$spec, n, prob = noise_probs,
         replace = T)
}

# Initiate community:
initial_comm <- commOfSize(1000, noise=0)

rand_no = 20

plotData <- data.frame()

for(noise in seq(0,1.25, by = 0.25)){
  for(csize in seq(1000, 10, by = -5)){
    
    # print each community size that will be compared with the 
    # initial community
    print(csize)
    
    # Set an initial communtity matrix
    sizeComm <- summary(initial_comm)
    # Get few estimations of a random community
    for(rand in 1:rand_no){
      sizeComm <- cbind(sizeComm, 
                        summary(commOfSize(csize, noise = noise)))
    }
    
    colnames(sizeComm) <- c("init", paste("rep", seq(1:rand_no), sep="_"))
    
    values <- as.matrix(beta.pair.abund(t(sizeComm))$beta.bray)[,1]
    
    row <- data.frame(size = csize, 
                      bc = values[-1],
                      noise = noise)
    
    plotData <- rbind(plotData, row)
    
  }
  
}

library(ggplot2)
ggplot(plotData, aes(x = size, y = bc))+
  geom_point(alpha = 0.1)+
  facet_wrap(~noise)
