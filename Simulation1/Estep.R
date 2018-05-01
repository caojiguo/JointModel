### Monte Carlo EM ####
### E-step: we sample from p(beta|Y) and evaluate the weights ###
### The posterior of beta is Multi-variate normal ###
betapos_cov <- function(sigma, sigma_beta, xi){
  sig_inv <- xi%*%t(xi)/sigma^2 + diag(1/sigma_beta^2)
  return(solve(sig_inv))
}

betapos_mean <- function(betapos_cov, alpha, b, sigma_beta,sigma, y, z, xi){
  t_mu <- (t(y-as.vector(z%*%alpha))%*%t(xi)/sigma^2 + t(b)%*%diag(1/sigma_beta^2))%*%betapos_cov
  return(as.vector(t_mu))
}

### sample beta_{i} ###
sampleRE <- function(Nsample, sigma, sigma_beta, alpha, b, Y, Z, xi, I){
  nobs <- I
  nc <- length(sigma_beta)
  #betaSample <- matrix(NA, nr = nobs*Nsample, nc = nc)
  betaiMean <- matrix(NA, nr = nobs, nc = nc)
  betaSample <- list()
  for(i in 1:nobs){
   # print(xi[[i]])
    cov <- betapos_cov(sigma, sigma_beta, xi[[i]])
    mean <- betapos_mean(cov, alpha, b, sigma_beta,sigma, as.vector(Y[[i]]), as.vector(Z[i,]), xi[[i]])
    betaSample_temp <- mvrnorm(Nsample ,mean, cov)
    betaSample[[i]] <- betaSample_temp
    #betaSample[(((i-1)*Nsample+1):(i*Nsample)),] <- betaSample_temp
    betaiMean[i,] <- apply(betaSample_temp, 2, mean)
  }
  return(betaSample = betaSample)
}

### evaluate weights for beta_{i} ###
### evaluate the survival loglikelihood, store logWeights ###
logWeightsFun <- function(betai, Zi, gamma1, gamma2, timei, statusi, lambda0, leftInterval, rightInterval, indexi){
  ZiRep <- matrix(rep(Zi, nrow(betai)), nr = nrow(betai), byrow = TRUE)
  par <- c(gamma1, gamma2)
  phi <- exp(cbind(ZiRep, betai)%*%par)
  mu <- phi*timei
  Interval <- cbind(leftInterval, rightInterval)
  if(statusi == 1){
    rt <- log(unlist(lapply(mu, function(x)(sum(as.numeric(x>leftInterval & x<= rightInterval)*lambda0))))+exp(-100)) + log(phi)-unlist(lapply(mu, function(x)(sum(as.numeric(x>= rightInterval)*lambda0*(rightInterval-leftInterval)))))
  }else{
    rt <- -unlist(lapply(mu, function(x)(sum(as.numeric(x>= rightInterval)*lambda0*(rightInterval-leftInterval)))))
  }
  return(rt)
}

normalizeLogWeights <- function(logWeights){
  maxLogWeights <- max(logWeights)
  stLogWeights <- logWeights - maxLogWeights
  rt <- exp(stLogWeights - log(sum(exp(stLogWeights))))
  return(rt)
}

WeightsMatrix <- function(beta, Z, gamma1, gamma2, time, status, lambda0, leftInterval, rightInterval, index, I){
  WeightsRF <- list()
  for(i in 1:I){
    betai <- beta[[i]]
    Zi <- Z[i,]
    timei <- time[i]
    statusi <- status[i]
    indexi <- index[i]
    logWeightsi <- logWeightsFun(betai, Zi, gamma1, gamma2, timei, statusi, lambda0, leftInterval, rightInterval, indexi)
    WeightsRF[[i]] <- normalizeLogWeights(logWeightsi)
  }
  return(WeightsRF)
}
