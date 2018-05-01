### Notation ####
### status: censor indicator ###
### observed time: time ###
### Y: measurements ###
### I: # patients 
### J: max time ###
### Simulation Setup ###
I <- Nobs
J <- 12

#T <- rep(NA, I)
### True value for parameters###
sigma1_true <- 1
sigma2_true <- 0.02
sigma3_true <- 0.5

alpha_true <- c(0, 0)
Z <- cbind(rep(1, I), rep(1, I))

gamma21_true <- 1
gamma22_true <- 1
gamma3_true <- -0.13

gamma2_true <- c(gamma21_true, gamma22_true)

betaMean_true <- c(2.5, -0.2)
beta1_true <- rnorm(I, betaMean_true[1], sd = sigma1_true)
beta2_true <- rnorm(I, betaMean_true[2], sd = sigma2_true)
beta_true <- cbind(beta1_true, beta2_true)

### 
simulatedtime <- rweibull(I, shape = 1,
                 scale=exp((beta_true%*%gamma2_true))) 
status <- as.numeric(simulatedtime< J-1)
maxtime <- J-1
time <- ifelse(status==0, maxtime, simulatedtime)
index <- which(time < maxtime)
d <- length(unique(time))-1

Y <- list()
xi <- list()
for(i in 1:I){
  xi[[i]] <- rbind(rep(1, floor(time[i]+1)), 1:floor(time[i]+1))
  #xi[[i]] <- 1:floor(time[i]+1)
  Y[[i]] <- beta_true[i,]%*%xi[[i]] + rnorm(floor(time[i]+1), 0, sigma3_true)
}


