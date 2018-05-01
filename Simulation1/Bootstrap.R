###Bootstrap resampling###
ndata <- Nobs 
Nbts <- Nbootstrap 
indexMatrix <- matrix(NA, nr = Nbts, nc = ndata)
for(i in 1:Nbts){
  indexMatrix[i,] <- sample.int(ndata, size = ndata, replace = TRUE)
}

sigmaBeta_bootstrap <- matrix(NA, nr = Nbts, nc = 2)
sigma3_bootstrap <- rep(NA, Nbts)
betaMean_bootstrap <- matrix(NA, nr = Nbts, nc = 2)
gamma_bootstrap <- matrix(NA, nr = Nbts, nc = 19)
alpha_bootstrap <- matrix(NA, nr = Nbts, nc = 16)

###
for(iterboots in 1:Nbts){
  iMatrixboots <- indexMatrix[iterboots,]
  statusboots <- status[iMatrixboots]
  timeboots <- time[iMatrixboots]
  indexboots <- unique(which(timeboots < maxtime))
  #dboots <- length(unique(timeboots))-1
  dboots <- length(indexboots)
  
  Y_boots <- list()
  xi_boots <- list()
  for(i in 1:ndata){
    xi_boots[[i]] <- xi[[iMatrixboots[i]]]
    Y_boots[[i]] <- Y[[iMatrixboots[i]]]
  }
  
  InitialValue <- Bootstrap_initialValue(statusboots, timeboots, indexboots, dboots, Y_boots, xi_boots,I)
  betaEstimate_initBootstrap <- InitialValue$betaEstimate_init
  btemp_initBootstrap <- InitialValue$btemp_init
  sigma_beta_initBootstrap <- InitialValue$sigma_beta_init
  sigma3_initBootstrap <- InitialValue$sigma3_init
  gamma_initialBootstrap <- InitialValue$gamma_initial
  lambda0_initialBootstrap <- InitialValue$lambda0_initial
  leftInterval_temp <- InitialValue$leftInterval
  rightInterval_temp <- InitialValue$rightInterval
  ###Finish optimization step ###
  BootstrapEM <- Bootstrap_optimization(betaEstimate_initBootstrap, btemp_initBootstrap, sigma_beta_initBootstrap, sigma3_initBootstrap, gamma_initialBootstrap, lambda0_initialBootstrap, alpha_init, statusboots, timeboots, indexboots, dboots, Y_boots, xi_boots,I, leftInterval_temp, rightInterval_temp)
  
  sigmaBeta_bootstrap[iterboots,] <- BootstrapEM$sigma_beta_temp
  sigma3_bootstrap[iterboots] <- BootstrapEM$sigma3_temp
  betaMean_bootstrap[iterboots,] <- BootstrapEM$btemp
  gamma_bootstrap[iterboots,] <- BootstrapEM$gamma_temp
  alpha_bootstrap[iterboots,] <- BootstrapEM$alpha_temp
  print(iterboots)
}

sd_sigmaBeta <- apply(sigmaBeta_bootstrap, 2, sd)
sd_sigma3 <- sd(sigma3_bootstrap)
sd_betaMean <- apply(betaMean_bootstrap, 2, sd)
sd_gamma <- apply(gamma_bootstrap, 2, sd) 
sd_alpha <- apply(alpha_bootstrap, 2, sd)

