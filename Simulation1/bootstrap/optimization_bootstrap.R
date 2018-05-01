Bootstrap_optimization <- function(betaEstimate_init, btemp_init,sigma_beta_init, sigma3_init, gamma_initial, lambda0_initial, alpha_init, status, time, index, d, Y, xi, I, leftInterval, rightInterval){
  betaEstimateTemp <- betaEstimate_init 
  btemp <- btemp_init
  btemp1 <- btemp
  sigma_beta_temp <- sigma_beta_init
  sigma_beta_temp1 <- sigma_beta_temp
  alpha_temp <- alpha_init
  sigma3_temp <- sigma3_init
  sigma3_temp1 <- sigma3_temp
  lambda0_temp <- lambda0_initial
  lambda0_temp1 <- lambda0_temp
  gamma_temp <- gamma_initial
  #gamma2_temp <- gamma_initial
  #gamma2_temp1 <- gamma2_temp
  eps = 1
  iter <- 0
  maxIterperSim <- 30
  while((eps > 0.002)&(iter <= maxIterperSim)){
    iter <- iter + 1
    print(iter)
    gamma1_temp <- gamma_temp[1:length(gamma1_true)]
    gamma2_temp <- gamma_temp[(length(gamma1_true)+1):(length(gamma1_true)+length(gamma2_true))]
    ### E step ###
    sampledRF <- sampleRE(iter*200, sigma3_temp, sigma_beta_temp, alpha_temp, btemp, Y, Z, xi, I)
    WeightsList <- WeightsMatrix(sampledRF, Z_gamma, gamma1_temp, gamma2_temp, time, status, lambda0_temp, leftInterval, rightInterval, index, I)
    
    betaEstimateTemp <- betaEstimate(sampledRF, WeightsList, I)
    
    btemp <- bestimate(betaEstimateTemp)
    
    Btemp <- Bestimate(sampledRF, btemp,WeightsList, I)
    
    sigma_beta_temp <- sqrt(Btemp)
    
    #print(sigma_beta_temp)
    alpha_temp <- alpha_estimate(Z, Y, betaEstimateTemp, xi, I)
    
    sigma3_temp <- sigma_estimate(Z, Y, betaEstimateTemp, alpha_temp, xi, I)
    
    #d <- length(unique(time))-1
    Cj <- rep(NA, d)
    #index <- which(time < maxtime)
    
    phi <- exp(Z_gamma%*%gamma1_temp + betaEstimate_init%*%gamma2_temp)
    PHI <- (phi*time)[index]

    ordered_PHI <- c(0, sort(PHI))
    mu <- phi*time
    
    leftInterval <- ordered_PHI[-length(ordered_PHI)]
    rightInterval <- ordered_PHI[-1]
    Interval <- cbind(leftInterval, rightInterval)
    
    #Cj_numerator <- rep(1,d) 
    Cj_numerator <- unlist(lapply(rightInterval, function(x)(length(which(mu == x)))))
    Cj_denominator <- (rightInterval - leftInterval)*unlist(lapply(rightInterval, function(x)(length(which(mu >= x)))))
    lambda0_temp <- Cj_numerator/Cj_denominator
    
    #par <- c(gamma2_temp)
    par <- c(gamma_temp)
    rt <- nlminb(start = par, optimizeInitial, control=list(trace=FALSE, abs.tol = 0.0001, iter.max = 10), Z = Z_gamma, beta = betaEstimateTemp, time = time, status = status, lambda0 = lambda0_temp, leftInterval = leftInterval, rightInterval = rightInterval, index = index)
    #gamma2_temp <- rt$par
    gamma_temp <- rt$par
    
    dif <- c(btemp1 - btemp, sigma_beta_temp1- sigma_beta_temp, sigma3_temp1 - sigma3_temp, gamma_temp1 - gamma_temp)
    dif2 <- c(btemp1, sigma_beta_temp1, sigma3_temp1, gamma_temp1)
    eps <- max(abs(dif)/(abs(dif2)+0.001))
    #print(eps)
    
    btemp1 <- btemp
    
    sigma_beta_temp1 <- sigma_beta_temp
    
    sigma3_temp1 <- sigma3_temp
    
    lambda0_temp1 <- lambda0_temp
    
    gamma_temp1 <- gamma_temp
    
    alpha_temp1 <- alpha_temp
  }
  return(list(btemp = btemp, sigma_beta_temp = sigma_beta_temp, sigma3_temp = sigma3_temp, gamma_temp = gamma_temp, lambda0_temp = lambda0_temp, alpha_temp = alpha_temp))
}


