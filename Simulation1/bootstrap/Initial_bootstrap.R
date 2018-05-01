###bootstrap: initial value ###
Bootstrap_initialValue <- function(status, time, index, d, Y, xi, I){
  niter <- 10
  nSampleVec <- (1:niter)*200
  sigma1_temp <- 1 + rnorm(1, 0, 0.1)
  sigma2_temp <- 0.02 + rnorm(1, 0, 0.005)
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
    
    btemp <- bestimate(betaEstimateTemp)
    
    Btemp <- Bestimate(sampledRF, btemp,WeightsList, I)
    
    sigma_beta_temp <- sqrt(Btemp)
    
    alpha_temp <- alpha_estimate(Z, Y, betaEstimateTemp, xi, I)
    
    sigma3_temp <- sigma_estimate(Z, Y, betaEstimateTemp, alpha_temp, xi, I)
    
  }
  
  betaEstimate_init <- betaEstimateTemp
  btemp_init <- btemp
  sigma_beta_init <- sigma_beta_temp
  sigma3_init <- sigma3_temp
  alpha_init <- alpha_temp
  
  ###Then fit survival model ###
  ###To deal with the multiple modes issue, we use multiple initial value ###
  bestLL <- 1000000000000
  gamma_initial <- rep(NA, length(c(gamma1_true, gamma2_true)))
  lambda0_initial <- rep(NA,d) 
  
  Cj <- rep(NA, d)
  #index <- which(time < maxtime)
  for(itry in 1:10){
    gamma1_temp <- gamma1_true + rnorm(17, 0, sd = c(0.05, 0.001, rep(0.01, 6), 0.003, rep(0.01, 2), rep(0.005, 6)))
    gamma2_temp <- gamma2_true + rnorm(2, 0, sd = c(0.05, 0.05))
    
    for(iter in 1:6){
      phi <- exp(Z_gamma%*%gamma1_temp + betaEstimate_init%*%gamma2_temp)
      PHI <- (phi*time)[index]
      ##u(k) ###
      ordered_PHI <- c(0, sort(PHI))
      mu <- phi*time
      
      leftInterval <- ordered_PHI[-length(ordered_PHI)]
      rightInterval <- ordered_PHI[-1]
      Interval <- cbind(leftInterval, rightInterval)
      
      #Cj_numerator <- rep(1,d)
      ### check it ###
      Cj_numerator <- unlist(lapply(rightInterval, function(x)(length(which(mu == x)))))
      Cj_denominator <- (rightInterval - leftInterval)*unlist(lapply(rightInterval, function(x)(length(which(mu >= x)))))
      #print(length = Cj_denominator)
      lambda0_temp <- Cj_numerator/Cj_denominator
      
      ###optimize gamma ###
      #par <- c(gamma2_temp)
      par <- c(gamma1_temp, gamma2_temp)
      rt <- nlminb(start = par, optimizeInitial, control=list(trace=FALSE, abs.tol = 0.0001, iter.max = 10), Z = Z_gamma, beta = betaEstimateTemp, time = time, status = status, lambda0 = lambda0_temp, leftInterval = leftInterval, rightInterval = rightInterval, index = index)
      gamma_temp <- rt$par
      gamma_temp1 <- gamma_temp[1:length(gamma1_true)]
      gamma_temp2 <- gamma_temp[(length(gamma1_true)+1):(length(gamma1_true)+length(gamma2_true))]
      #gamma2_temp <- rt$par
    }
    
    if(bestLL > rt$objective){
      bestLL <- rt$objective
      gamma_initial <- rt$par
      lambda0_initial <- lambda0_temp
      rightInterval_initial <- rightInterval
      leftInterval_initial <- leftInterval 
    }
  }
  return(list(betaEstimate_init = betaEstimate_init, btemp_init = btemp_init, sigma_beta_init = sigma_beta_init, sigma3_init = sigma3_init, gamma_initial = gamma_initial, lambda0_initial = lambda0_initial, rightInterval = rightInterval, leftInterval = leftInterval ))
}



