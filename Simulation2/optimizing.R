### sample random effects ###
### Ni step at ith iteration ###
### converge criteria eps = 0.0001 ###
#Niter <- 30
#NSampleVec <- floor((1:Niter)*200)
betaEstimateTemp <- betaEstimate_init 
btemp <- btemp_init
btemp1 <- btemp
sigma_beta_temp <- sigma_beta_init
sigma_beta_temp1 <- sigma_beta_temp
alpha_temp <- c(0, 0)
sigma3_temp <- sigma3_init
sigma3_temp1 <- sigma3_temp
lambda0_temp <- lambda0_initial
lambda0_temp1 <- lambda0_temp
#gamma_temp <- gamma_initial
#gamma1_temp <- gamma_initial[1:length(gamma1_true)]
gamma2_temp <- gamma_initial
gamma2_temp1 <- gamma2_temp
eps = 1
iter <- 0
while((eps > 0.002)&(iter <= maxIterperSim)){
  iter <- iter + 1
  print(iter)
  ### E step ###
  sampledRF <- sampleRE(iter*200, sigma3_temp, sigma_beta_temp, alpha_temp, btemp, Y, Z, xi, I)
  WeightsList <- WeightsMatrix(sampledRF, Z, c(0,0), gamma2_temp, time, status, lambda0_temp, leftInterval, rightInterval, index, I)
  
  betaEstimateTemp <- betaEstimate(sampledRF, WeightsList, I)
  
  btemp <- bestimate(betaEstimateTemp)
  
  Btemp <- Bestimate(sampledRF, btemp,WeightsList, I)
  
  sigma_beta_temp <- sqrt(Btemp)
  
  print(sigma_beta_temp)
  
  #alpha_temp <- alpha_estimate(Z, Y, betaEstimateTemp, xi, I)
  
  sigma3_temp <- sigma_estimate(Z, Y, betaEstimateTemp, alpha_temp, xi, I)
  
  #intervals <- uIntervals(gamma2_temp, betaEstimateTemp, time, index)
  #leftInterval <- intervals$leftInterval
  #rightInterval <- intervals$rightInterval
  
  #lambda0_temp <- lambda0Fun(time, index, gamma2_temp, leftInterval, rightInterval, sampledRF, WeightsList, I)
  #rt <- nlminb(start = c(gamma2_temp), optimizeFun, control=list(trace=TRUE, abs.tol = 0.0001, iter.max = 1), Z = Z, betaSample = sampledRF, time = time, status = status, lambda0 = lambda0_temp, leftInterval = leftInterval, rightInterval = rightInterval, index = index, Weights = WeightsList, I = I, Nsample = NSampleVec[iter]) 
  #gamma_temp <- rt$par
  #gamma2_temp <- gamma_temp
  d <- length(unique(time))-1
  Cj <- rep(NA, d)
  index <- which(time < maxtime)
  
  phi <- exp(betaEstimateTemp%*%gamma2_temp)
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
  
  par <- c(gamma2_temp)
  rt <- nlminb(start = par, optimizeInitial, control=list(trace=FALSE, abs.tol = 0.0001, iter.max = 10), Z = Z, beta = betaEstimateTemp, time = time, status = status, lambda0 = lambda0_temp, leftInterval = leftInterval, rightInterval = rightInterval, index = index)
  gamma2_temp <- rt$par
  
  dif <- c(btemp1 - btemp, sigma_beta_temp1- sigma_beta_temp, sigma3_temp1 - sigma3_temp, gamma2_temp1 - gamma2_temp)
  dif2 <- c(btemp1, sigma_beta_temp1, sigma3_temp1, gamma2_temp1)
  eps <- max(abs(dif)/abs(dif2+0.001))
  print(eps)
  
  btemp1 <- btemp
 
  sigma_beta_temp1 <- sigma_beta_temp
 
  sigma3_temp1 <- sigma3_temp

  lambda0_temp1 <- lambda0_temp

  gamma2_temp1 <- gamma2_temp
  
}
#for(iter in 1:Niter){

#}

