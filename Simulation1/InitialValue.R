###Fit mixed effects model ###
niter <- 20
#sigmaBeta_store <- matrix(NA, nr = niter, nc = 2)
#sigma3_store <- rep(NA, niter)
#betaMean_store <- matrix(NA, nr = niter, nc = 2)
nSampleVec <- (1:niter)*200
sigma1_temp <- sigma1_true + rnorm(1, 0, 0.1)
sigma2_temp <- sigma2_true + rnorm(1, 0, 0.005)
sigma3_temp <- sigma3_true + rnorm(1, 0, 0.1)
sigma_beta_temp <- c(sigma1_temp, sigma2_temp)
alpha_temp <- alpha_true + rnorm(16, 0, c(0.1, 0.01, rep(0.2,6), rep(0.05,8)))
betaMean_temp <- betaMean_true + rnorm(2, 0, c(0.2, 0.1))
WeightsList <- list()

for(iter in 1:niter){
  ### E step ###
  sampledRF <- sampleRE(nSampleVec[iter], sigma3_temp, sigma_beta_temp, alpha_temp, betaMean_temp, Y, Z, xi, I)
  ### M step ###
  for(i in 1:I){
    WeightsList[[i]] <- rep(1/nSampleVec[iter], nSampleVec[iter])
  }
  betaEstimateTemp <- betaEstimate(sampledRF, WeightsList, I)
  #betaEstimateTemp <- betaEstimate(sampledRF, I)
  
  btemp <- bestimate(betaEstimateTemp)
  #print(btemp)
  #betaMean_store[iter,] <- btemp
  Btemp <- Bestimate(sampledRF, btemp,WeightsList, I)
  #Btemp <- Bestimate(sampledRF, btemp, I)
  
  sigma_beta_temp <- sqrt(Btemp)
  #sigmaBeta_store[iter, ] <- sigma_beta_temp
  #print(sigma_beta_temp)
  
  alpha_temp <- alpha_estimate(Z, Y, betaEstimateTemp, xi, I)
  #alpha_store[iter,] <- alpha_temp
  
  sigma3_temp <- sigma_estimate(Z, Y, betaEstimateTemp, alpha_temp, xi, I)
  #sigma3_store[iter] <- sigma3_temp
  
  #print(sigma3_temp)
}

betaEstimate_init <- betaEstimateTemp
btemp_init <- btemp
sigma_beta_init <- sigma_beta_temp
alpha_init <- alpha_temp
sigma3_init <- sigma3_temp

###Then fit survival model ###
###To deal with the multiple modes issue, we use multiple initial value ###
bestLL <- 1000000000000
#gamma_initial <- rep(NA, length(gamma2_true))
gamma_initial <- rep(NA, length(c(gamma1_true, gamma2_true)))
lambda0_initial <- rep(NA,d) 


Cj <- rep(NA, d)
index <- which(time < maxtime)
for(itry in 1:10){
  gamma1_temp <- gamma1_true + rnorm(17, 0, sd = c(0.05, 0.001, rep(0.01, 6), 0.003, rep(0.01, 2), rep(0.005, 6)))
  
  gamma2_temp <- gamma2_true + rnorm(2, 0, sd = c(0.001, 0.005))

  for(iter in 1:5){
    phi <- exp(Z_gamma%*%gamma1_temp + betaEstimate_init%*%gamma2_temp)
    PHI <- (phi*time)[index]
    ##u(k) ###
    ordered_PHI <- c(0, sort(PHI))
    mu <- phi*time
    
    leftInterval <- ordered_PHI[-length(ordered_PHI)]
    rightInterval <- ordered_PHI[-1]
    Interval <- cbind(leftInterval, rightInterval)
    
    Cj_numerator <- rep(1,d) 
    Cj_denominator <- (rightInterval - leftInterval)*unlist(lapply(rightInterval, function(x)(length(which(mu >= x)))))
    lambda0_temp <- Cj_numerator/Cj_denominator
    
    ###optimize gamma ###
    #par <- c(gamma2_temp)
    par <- c(gamma1_temp, gamma2_temp)
    rt <- nlminb(start = par, optimizeInitial, control=list(trace=FALSE, abs.tol = 0.0001, iter.max = 5), Z = Z_gamma, beta = betaEstimateTemp, time = time, status = status, lambda0 = lambda0_temp, leftInterval = leftInterval, rightInterval = rightInterval, index = index)
    gamma_temp <- rt$par
    gamma_temp1 <- gamma_temp[1:length(gamma1_true)]
    gamma_temp2 <- gamma_temp[(length(gamma1_true)+1):(length(gamma1_true)+length(gamma2_true))]
  }
  #gamma1_temp <- gamma1_true + rnorm(5, 0, sd = c(0.05, 0.001, 0.01,0.01, 0.05))

  if(bestLL > rt$objective){
    bestLL <- rt$objective
    gamma_initial <- rt$par
    lambda0_initial <- lambda0_temp
  }
}




