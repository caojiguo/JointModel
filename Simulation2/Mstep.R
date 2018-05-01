### M-step ###
### optimize alpha, beta, sigma ###
###incorporate weights in this function ###
betaEstimate <- function(betaSample, WeightsList, I){
  E_beta <- matrix(NA, nr = I, nc = ncol(betaSample[[1]]))
  for(i in 1:I){
    E_beta[i,] <- apply(betaSample[[i]]*cbind(WeightsList[[i]],WeightsList[[i]]), 2, sum) ### incorporate weights here
  }
  return(E_beta)
}


# betaEstimate <- function(betaSample, I){
#   E_beta <- matrix(NA, nr = I, nc = ncol(betaSample[[1]]))
#   for(i in 1:I){
#     E_beta[i,] <- apply(betaSample[[i]], 2, mean) ### incorporate weights here
#   }
#   return(E_beta)
# }

bestimate <- function(E_beta){
  return(apply(E_beta, 2, mean))
}

### output are variance instead of sd ###
### incorporate weights here ###
Bestimate <- function(betaSample, meanSample,WeightsList, I){
  Nsample <- nrow(betaSample[[1]])
  Var_beta <- matrix(NA, nr = I, nc = ncol(betaSample[[1]]))
  meanSampleMatrix <- matrix(rep(meanSample, Nsample), nr = Nsample, byrow = TRUE)
  for(i in 1:I){
    Var_beta[i,] <- apply((betaSample[[i]] - meanSampleMatrix)*(betaSample[[i]] - meanSampleMatrix)*cbind(WeightsList[[i]],WeightsList[[i]]), 2, sum)
  }
  return(apply(Var_beta, 2, mean))
}


alpha_estimate <- function(Z, Y, E_beta, xi, I){
  dim <- ncol(Z)
  # nc <- ncol(Y); every Y has different length.
  part1 <- matrix(0, nr = dim, nc = dim)
  part2 <- rep(0, dim)
  for(i in 1:I){
    Zi <- matrix(rep(Z[i,], length(unlist(Y[[i]]))), nc = dim, byrow = TRUE)
    #Zi <- Zrep[((i-1)*nc+1):(i*nc),]
    part1 <- part1 + t(Zi)%*%Zi
    part2 <- part2 + t(Zi)%*%as.vector(Y[[i]]-as.vector(t(unlist(xi[[i]]))%*%E_beta[i,]))
  }
  alpha <- solve(part1)%*%part2
  return(alpha)
}

sigma_estimate <- function(Z, Y, E_beta, alpha, xi, I){
  dim <- ncol(Z)
  part1 <- 0
  for(i in 1:I){
    Zi <- matrix(rep(Z[i,], length(unlist(Y[[i]]))), nc = dim, byrow = TRUE)
    yy <- (as.vector(Y[[i]])-t(xi[[i]])%*%E_beta[i,]- Zi%*%alpha)
    part1 <- part1 + t(yy)%*%yy
  }
  denominator <- 0
  for(i in 1:I){
    denominator <- denominator + length(Y[[i]])
  }
  sigma <- part1/denominator
  return(as.vector(sqrt(sigma)))
}

### \lambda_{0} ###
### First estimate the left and right intervals ###
# meanBetai <- function(betaEstimate, Weights, I){
#   rt <- matrix(NA, nr = I, nc = ncol(betaEstimate[[1]]))
#   for(i in 1:I){
#     rt[i,] <- apply(betaEstimate[[i]], 2, sum)
#   }
#   return(rt)
# }
uIntervals <- function(gamma2, E_beta, time, index){
  phi <- exp(E_beta%*%gamma2)
  PHI <- (phi*time)[index]
  ordered_PHI <- c(0, sort(PHI))
  mu <- phi*time
  leftInterval <- ordered_PHI[-length(ordered_PHI)]
  rightInterval <- ordered_PHI[-1]
  return(list(leftInterval = leftInterval, rightInterval = rightInterval))
}

lambda0Fun <- function(time, index, gamma2, leftInterval, rightInterval, betaSample, Weights, I){
  d <- length(unique(time))-1
  ### use phi to denote u_{k} ###
  phi <- list()
  for(i in 1:I){
    phi[[i]] <- exp(betaSample[[i]]%*%gamma2)*time[i]
  }
  Cj_numerator <- matrix(NA, nr = I, nc = d)
  Cj_denominator <- matrix(NA, nr = I, nc = d)
  for(i in 1:I){
    if(status[i] == 0){
      Cj_numerator[i,] <- 0
    }else{
      for(k in 1:d){
        Cj_numerator[i,k] <- sum(unlist(lapply(phi[[i]], function(x)(as.numeric((leftInterval[k]<x)&(x<=rightInterval[k])))))*Weights[[i]])
      }
    }
    for(k in 1:d){
      Cj_denominator[i,k] <- (rightInterval[k] - leftInterval[k])*sum(unlist(lapply(phi[[i]], function(x)(as.numeric(x>=rightInterval[k]))))*Weights[[i]]) 
    }
  }
  lambda0 <- apply(Cj_numerator, 2, sum)/apply(Cj_denominator, 2, sum)
  return(lambda0)
}


###obtain initial value of lambda0 ###
lambda0_init <- function(par, Z, E_beta, time, status, index){
  phi <- exp(E_beta%*%gamma2_temp)
  PHI <- (phi*time)[index]
  ##u(k) ###
  ordered_PHI <- c(0, sort(PHI))
  mu <- phi*time
  
  leftInterval <- ordered_PHI[-length(ordered_PHI)]
  rightInterval <- ordered_PHI[-1]
  Interval <- cbind(leftInterval, rightInterval)
  
  Cj_numerator <- rep(1,d) 
  Cj_denominator <- (rightInterval - leftInterval)*unlist(lapply(rightInterval, function(x)(length(which(mu >= x)))))
  lambda0 <- Cj_numerator/Cj_denominator
  return(lambda0)
}


###obtain initial value of gamma ###
optimizeInitial <- function(par, Z, beta, time, status, lambda0, leftInterval, rightInterval, index){
  phi <- exp(beta%*%par)
  mu <- phi*time
  status1 <- status[index]
  mu1 <- mu[index]
  #print(length(lambda0)-length(rightInterval))
  Interval <- cbind(leftInterval, rightInterval)
  rt <- -(sum(log(unlist(lapply(mu1, function(x)(sum(as.numeric(x>leftInterval & x<= rightInterval)*lambda0))))+exp(-100)  )) + sum(status*(log(phi))-unlist(lapply(mu, function(x)(sum(as.numeric(x>= rightInterval)*lambda0*(rightInterval-leftInterval)))))))
  return(rt)
}

###optimize gamma ###
optimizeFun <- function(par, Z, betaSample, time, status, lambda0, leftInterval, rightInterval, index, Weights, I, Nsample){
  d <- length(unique(time))-1
  phi <- list()
  uk <- list()
  #gamma1 <- par[1:ncol(Z)]
  #gamma2 <- par[-(1:ncol(Z))]
  for(i in 1:I){
    phi[[i]] <- exp(betaSample[[i]]%*%par)
    uk[[i]] <- phi[[i]]*time[i]
  }
  ## loglikelihood is the summation of three parts
  ## we use part1, part2, part3 to denote
  rt <- 0
  for(i in 1:I){
    if(status[i] == 0){
      rt <- rt - sum(lambda0*(rightInterval-leftInterval)*unlist(lapply(rightInterval, function(x)(sum(as.numeric(x <= phi[[i]])*Weights[[i]])))))
        #sum(unlist(lapply(uk[[i]], function(x)(sum(as.numeric(x>= rightInterval)*lambda0*(rightInterval-leftInterval)))))*Weights[[i]])
    }else{
      part1 <- 0
      for(ii in 1:Nsample){
        part1 <- part1 + Weights[[i]][ii]*log(sum(lambda0*as.numeric(apply(cbind(leftInterval, rightInterval), 1, function(x)((uk[[i]][ii]>x[1])&(uk[[i]][ii]<=x[2])))))+exp(-100))
      }
      rt <- rt + sum(part1) + sum(log(phi[[i]])*Weights[[i]]) - sum(lambda0*(rightInterval-leftInterval)*unlist(lapply(rightInterval, function(x)(sum(as.numeric(x <= phi[[i]])*Weights[[i]])))))
    }
    
  }
  #status1 <- status[index]
  #mu1 <- mu[index]
  #Interval <- cbind(leftInterval, rightInterval)
  #rt <- -(sum(log(unlist(lapply(mu1, function(x)(sum(as.numeric(x>leftInterval & x<= rightInterval)*lambda0))))+exp(-100)  )) + sum(status*(log(phi))-unlist(lapply(mu, function(x)(sum(as.numeric(x>= rightInterval)*lambda0*(rightInterval-leftInterval)))))))
  return(-rt)
}
