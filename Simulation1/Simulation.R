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
sigma1_true <- 14.91
sigma2_true <- 2.89
sigma3_true <- 0.85

#alpha_true <- c(0, 0)
alpha_true <- c(1, -0.17, 5.39, -3.69, -6.45, 7.56, 10.75, 16.52, -0.93, -2.35, -1.54, -0.27, -0.52, -0.67, -0.96, -1.15)
#Z <- cbind(rep(1, I), rep(1, I))
Z0 <- rep(1, I)
Z1 <- rnorm(I, mean = 40, sd = 8)
Z2 <- as.numeric(runif(I, 0, 1)>0.41)
Z3 <- as.numeric(runif(I, 0, 1)>0.30)
Z4 <- as.numeric(runif(I, 0, 1)>0.77)
## TX ###
Z5 <- as.numeric(runif(I, 0, 1)>0.9)
Z6 <- as.numeric(runif(I, 0, 1)>0.77)
Z7 <- as.numeric(runif(I, 0, 1)>0.6)
## PKPRA ###
Z8 <- as.numeric(runif(I, 0, 1)>0.29)
Z9 <- as.numeric(runif(I, 0, 1)>0.76)
## HLA ###
Z10 <- as.numeric(runif(I, 0, 1)>0.2)
## DIA ###
Z11 <- as.numeric(runif(I, 0, 1)>0.78)
Z12 <- as.numeric(runif(I, 0, 1)>0.62)
Z13 <- as.numeric(runif(I, 0, 1)>0.83)
Z14 <- as.numeric(runif(I, 0, 1)>0.90)
## txtype2 ###
Z15 <- as.numeric(runif(I, 0, 1)>0.61)
#Z5 <- rnorm(I, mean = 5, sd = 2.5)
Z <- cbind(Z0, Z1, Z2, Z3, Z4, Z5, Z6, Z7, Z8, Z9, Z10, Z11, Z12, Z13, Z14, Z15)

#gamma21_true <- 1
#gamma22_true <- 1
#gamma3_true <- -0.13
W <- rep(0, I)
indexW <- sample(1:100, 20)
W[indexW] <- 1
Z_gamma <- cbind(Z, W)
gamma21_true <- -0.07
gamma22_true <- -0.21
gamma3_true <- -0.13

gamma10_true <- 4
gamma11_true <- 0.02
gamma12_true <- -0.23
gamma13_true <- 0.17
gamma14_true <- 0.23
gamma15_true <- -0.29
gamma16_true <- -0.76
gamma17_true <- -0.95
gamma18_true <- 0.06
gamma19_true <- 0.27
gamma110_true <- 0.24
gamma111_true <- 0.07
gamma112_true <- 0.08
gamma113_true <- 0.33
gamma114_true <- 0.38
gamma115_true <- 0.14


gamma1_true <- c(gamma10_true, gamma11_true, gamma12_true, gamma13_true, gamma14_true, gamma15_true, gamma16_true, gamma17_true, gamma18_true, gamma19_true, gamma110_true, gamma111_true, gamma112_true, gamma113_true, gamma114_true, gamma115_true, gamma3_true)
gamma2_true <- c(gamma21_true, gamma22_true)


gamma2_true <- c(gamma21_true, gamma22_true)

betaMean_true <- c(48.94, -1.36)
beta1_true <- rnorm(I, betaMean_true[1], sd = sigma1_true)
beta2_true <- rnorm(I, betaMean_true[2], sd = sigma2_true)
beta_true <- cbind(beta1_true, beta2_true)

### 
simulatedtime <- rweibull(I, shape = 1,
                 scale=exp((Z_gamma%*%gamma1_true + beta_true%*%gamma2_true))) 
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
  Y[[i]] <- as.vector(alpha_true%*%Z[i,]) + beta_true[i,]%*%xi[[i]] + rnorm(floor(time[i]+1), 0, sigma3_true)
}


