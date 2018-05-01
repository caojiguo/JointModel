rm(list=ls())
### main function ###
library(MASS)
set.seed(3421)
source('Estep.R')
source('Mstep.R')
source('./bootstrap/Initial_bootstrap.R')
source('./bootstrap/optimization_bootstrap.R')

filename_sigmaBeta = './results/sigmaBeta'
filename_sigma3 = './results/sigma3'
filename_betaMean = './results/betaMean'
filename_gamma2 = './results/gamma2'
filename_Niters = './results/Niters'

filename_sigmaBetasd = './results/sigmaBetasd'
filename_sigma3sd = './results/sigma3sd'
filename_betaMeansd = './results/betaMeansd'
filename_gamma2sd = './results/gamma2sd'



Nsimulations = 100
Nobs <- 100
Nbootstrap <- 100
maxIterperSim <- 100
### we do bootstrap if 'DoBootstrap = 1' ###
DoBootstrap = 1
#alpha_store <- matrix(NA, nr = Nsimulations, nc = ncol(Z))
sigmaBeta_store <- matrix(NA, nr = Nsimulations, nc = 2)
sigma3_store <- rep(NA, Nsimulations)
betaMean_store <- matrix(NA, nr = Nsimulations, nc = 2)
#gamma1_store <- matrix(NA, nr = Nsimulations, nc = ncol(Z))
gamma2_store <- matrix(NA, nr = Nsimulations, nc = 2)
Niters_store <- rep(NA, Nsimulations)

sigmaBetasd_store <- matrix(NA, nr = Nsimulations, nc = 2)
sigma3sd_store <- rep(NA, Nsimulations)
betaMeansd_store <- matrix(NA, nr = Nsimulations, nc = 2)
#gamma1_store <- matrix(NA, nr = Nsimulations, nc = ncol(Z))
gamma2sd_store <- matrix(NA, nr = Nsimulations, nc = 2)
ifAppend=FALSE;
for(isim in 1:Nsimulations){
  if(isim == 1){
    ifAppend=FALSE
  }else{
    ifAppend=TRUE
  }
  
  source('Simulation.R')
  source('InitialValue.R')
  source('optimizing.R')
  #alpha_store[isim,] <- alpha_temp
  sigmaBeta_store[isim,] <- sigma_beta_temp
  sigma3_store[isim] <- sigma3_temp
  betaMean_store[isim,] <- btemp
  #gamma1_store[isim,] <- gamma1_temp
  gamma2_store[isim,] <- gamma2_temp
  Niters_store[isim] <- iter
  
  if(DoBootstrap == 1){
    source('Bootstrap.R')
    sigmaBetasd_store[isim,] <- sd_sigmaBeta
    sigma3sd_store[isim] <- sd_sigma3
    betaMeansd_store[isim,] <- sd_betaMean
    gamma2sd_store[isim,] <- sd_gamma2
    
    write.table(matrix(round(sd_sigmaBeta, digits = 6),nr=1), file = paste(filename_sigmaBetasd,".txt",sep=""), append = ifAppend, row.names = FALSE,
                col.names = FALSE);
    
    write.table(matrix(round(sd_sigma3, digits = 6),nr=1), file = paste(filename_sigma3sd,".txt",sep=""), append = ifAppend, row.names = FALSE,
                col.names = FALSE);
    
    write.table(matrix(round(sd_betaMean, digits = 6),nr=1), file = paste(filename_betaMeansd,".txt",sep=""), append = ifAppend, row.names = FALSE,
                col.names = FALSE);
    
    write.table(matrix(round(sd_gamma2, digits = 6),nr=1), file = paste(filename_gamma2sd,".txt",sep=""), append = ifAppend, row.names = FALSE,
                col.names = FALSE);
  }
  
  write.table(matrix(round(sigma_beta_temp, digits = 6),nr=1), file = paste(filename_sigmaBeta,".txt",sep=""), append = ifAppend, row.names = FALSE,
              col.names = FALSE);
  
  write.table(matrix(round(sigma3_temp, digits = 6),nr=1), file = paste(filename_sigma3,".txt",sep=""), append = ifAppend, row.names = FALSE,
              col.names = FALSE);

  write.table(matrix(round(btemp, digits = 6),nr=1), file = paste(filename_betaMean,".txt",sep=""), append = ifAppend, row.names = FALSE,
              col.names = FALSE);
  
  write.table(matrix(round(gamma2_temp, digits = 6),nr=1), file = paste(filename_gamma2,".txt",sep=""), append = ifAppend, row.names = FALSE,
              col.names = FALSE);
  
  write.table(matrix(round(iter, digits = 6),nr=1), file = paste(filename_Niters,".txt",sep=""), append = ifAppend, row.names = FALSE,
              col.names = FALSE);
  
}


