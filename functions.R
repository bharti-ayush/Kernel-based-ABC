# This file contains all the functions required to run the main file. Load these functions in the environment before running the main file.


############-----------Saleh-Valenzuela Model-----------#############
# Function Parameters:-
#       lambda_cluster:   cluster arrival rate  (\Lambda)
#       lambda_ray:       ray arrival rate  (\lambda)
#       gamma_cluster:    cluster decay rate (\Gamma)
#       gamma_ray:        ray decay rate  (\gamma)    
#       betaSquare_00:    average power of the first arriving component (Q)
#       sigma2_N:         noise variance (\sigma_W^2)
#       B:                bandwidth used for measuring the data
#       N:                number of channel realizations for a given parameter value
#       Ns:               number of frequency (or time) points in the data
#       tau_0:            time delay of the first arriving multipath component

# Output: simulated channel transfer function data, dimensions NxNs
SalehValenzuela <- function(lambda_cluster , lambda_ray , gamma_cluster , gamma_ray , betaSquare_00 ,sigma2_N, B = 4e9, N = 100, Ns = 801, tau_0 = 0){
  
  delta_f <- B / (Ns-1)
  t_end <- 1 / delta_f
  
  mu_clusterPoisson <- lambda_cluster * t_end
  
  H <- array(0, dim = c(N, Ns))
  
  for (n in 1:N) {
    
    repeat{
      nClusters <- rpois(1, mu_clusterPoisson)
      if (nClusters > 0) {
        break
      }
    }
    
    delays_cluster <- runif(nClusters - 1, tau_0, t_end)
    
    delays_cluster <- c(tau_0, delays_cluster[order(delays_cluster)])
    
    betaSquare_cluster <- betaSquare_00 * exp(- delays_cluster / gamma_cluster)
    
    nRays <- array(0, dim = c(nClusters, 1))
    mu_rayPoisson <- array(0, dim = c(nClusters, 1))
    
    total_delays <- c()
    total_gains <- c()
    
    for (i in 1:nClusters) {
      
      mu_rayPoisson[i] <- lambda_ray * (t_end - delays_cluster[i])  #Take this out of this loop
      
      nRays[i] <- rpois(1, mu_rayPoisson[i])  #Check if this can be taken outside the loop
      
      tau_kl <- c(0, runif(nRays[i], 0, t_end - delays_cluster[i])) #+ delays_cluster[i]
      
      betaSquare_kl <- betaSquare_cluster[i] * exp(- tau_kl / gamma_ray)
      
      alpha_kl <- array(0, dim = c(nRays[i] +1 ,1))
      
      for (l in 1:length(tau_kl)) {
        
        alpha_kl[l] <- (rnorm(1,mean = 0,sd = sqrt(betaSquare_kl[l]/2)) + rnorm(1,mean = 0,sd = sqrt(betaSquare_kl[l]/2))*1i)
        
      }
      
      total_delays <- append(total_delays, tau_kl + delays_cluster[i])
      total_gains <- append(total_gains, alpha_kl)
    }
    
    delays_ray_full <- total_delays[order(total_delays)]
    
    H[n,] <- exp(-1i*2*pi*delta_f * (0:(Ns-1))%*%t(total_delays)) %*% total_gains
    
  }
  
  Noise <- array(data = 0, dim = c(N, Ns))
  
  for (i in 1:N) {
    
    Noise[i,] <- rnorm(Ns, mean = 0, sd= sqrt(sigma2_N / 2)) + rnorm(Ns, mean = 0,  sd= sqrt(sigma2_N / 2)) * 1i
    
  }
  
  # Received signal in frequency domain
  
  Y <- H + Noise
  
  return(Y)
  # return(list(Y = Y,delays = total_delays))
}


##########--------Function to compute power delay profile of data Y-----------#########
pdp <- function(Y){
  library(pracma)
  
  y <- array(0, dim = dim(Y))
  p <- array(0, dim = dim(Y))
  
  for (i in 1:length(Y[,1])) {
    y[i,] <- ifft(Y[i,])
    p[i,] <- abs(y[i,])^2
  }
  
  return(10*log10(apply(p, 2, mean) ))
}

##########--------Function to compute K temporal moments of the data Y-----------#########
temporalMomentsGeneral <- function(Y, K=3, B = 4e9){
  library(pracma)
  
  N <- dim(Y)[1]
  Ns <- dim(Y)[2]
  
  delta_f <- B / (Ns-1)
  t_max <- 1 / delta_f
  
  tau <- seq(from = 0, to = t_max, length.out = Ns)
  
  m <- array(0, dim = c(N, K))
  
  for (k in 1:K) {
    for (i in 1:N) {
      y <- ifft(Y[i,])
      m[i,k] <- trapz(tau, tau^(k-1) * abs(y)^2)
    }
  }
  
  return(m)
}

##########--------Function to compute the Gram matrix of Z-----------#########
kernel_matrix <- function(Z,bw){
  
  return( exp(-dist(Z,method="euclidean", diag = TRUE, upper = TRUE)^2 / bw^2))
  
}


##########--------Function to compute the MMD between data-sets X and Y with lengthscale bw-----------#########
MMD_fast <- function(X,Y, bw){
  
  if (!is.matrix(X)) {
    X <- array(X, dim = c(length(X), 1))
    Y <- array(Y, dim = c(length(Y), 1))
  }
  
  m <- dim(Y)[1]
  n <- dim(X)[1]
  
  Z <- rbind(X,Y)
  
  B <- as.matrix(kernel_matrix(Z, bw))
  diag(B) <- 1
  
  K_XX <- B[1:n, 1:n]
  # diag(K_XX) <- 0
  K_YY <- B[(n+1) : (n+m), (n+1) : (n+m)]
  K_XY <- B[1:n, (n+1):(n+m) ]
  # diag(K_YY) <- 0
  
  
  one_n <- array(1, dim = c(n,1))
  one_m <- array(1, dim = c(m,1))
  
  estimated_MMD <- 1 / (n * (n-1)) * (t(one_n) %*% K_XX %*% one_n) -
    2 / (n * m) * t(one_n) %*% K_XY %*% one_m +
    1 / (m * (m-1))* t(one_m) %*% K_YY %*% one_m
  
  return(estimated_MMD)
}


plot_estimates <- function(param, param_true, weights = NULL){
  M <- dim(param)[1]
  if (is.null(weights)) {
    w <- rep(1 / M, M)
  } else {
    w <- weights
  }
  
  beta_start <- 1e-5
  beta_end <- 5e-4
  # 
  lambda_cluster_start <- 1e6
  lambda_cluster_end <- 1e7
  
  lambda_ray_start <- 1e7
  lambda_ray_end <- 2e9
  
  gamma_cluster_start <- 4e-9
  gamma_cluster_end <- 4e-8
  
  gamma_ray_start <- 4e-10
  gamma_ray_end <- 4e-8
  
  sigma2_N_start <- 1e-7
  sigma2_N_end <- 1e-6
  
  
  par(mfrow = c(2,3))
  plot(density(param[,1], weights = w), xlim = c(beta_start, beta_end), xlab = "Q", main = "")
  abline(v = param_true[1], col = 'green', lwd = 2)
  plot(density(param[,2], weights = w), xlim = c(lambda_cluster_start, lambda_cluster_end), xlab = expression(Lambda), main = "")
  abline(v = param_true[2], col = 'green', lwd = 2)
  plot(density(param[,3], weights = w), xlim = c(lambda_ray_start, lambda_ray_end), xlab = expression(lambda), main = "")
  abline(v = param_true[3], col = 'green', lwd = 2)
  plot(density(param[,4], weights = w), xlim = c(gamma_cluster_start, gamma_cluster_end), xlab = expression(Gamma), main = "")
  abline(v = param_true[4], col = 'green', lwd = 2)
  plot(density(param[,5], weights = w), xlim = c(gamma_ray_start, gamma_ray_end), xlab = expression(gamma), main = "")
  abline(v = param_true[5], col = 'green', lwd = 2)
  plot(density(param[,6], weights = w), xlim = c(sigma2_N_start, sigma2_N_end), xlab = expression(sigma[N]^2), main = "")
  abline(v = param_true[6], col = 'green', lwd = 2)
  par(mfrow = c(1,1))
  
}


##########--------Function to compute the means and covariances of the temporal moments-----------#########
summaryMoments <- function(x){
  library(matrixcalc)
  z <-  c(apply(x, 2, mean), upper.triangle(cov(x)))
  z <- z[z!=0]
  y <- array(z, dim = c(1, length(z)))
  # colnames(y) <- c("Mean(m0)", "Mean(m1)", "Mean(m2)", "Var(m0)", "Cov(m0,m1)", "Var(m1)", "Cov(m0,m2)", "Cov(m1,m2)", "Var(m2)")
  return(y)
}

##########--------Function running the ABC algorithm-----------#########
# Function Parameters:
# param:      set of parameter values sampled from the prior
# S:          set of simulated temporal moments corresponding to the parameter values
# m_obs:      observed temporal moments
# M_epsilon:  number of accepted parameter samples
# transf:     whether to logit transform the parameter values before regression adjustment or not (used for bounded priors)
# mismatch:   whether the model is misspecied or not
# B:          bandwidth
# Ns:         number of frequency or time points in the data

# Output a list containing:
# rej: parameter samples after rejection ABC
# reg: rejection ABC parameter samples after regression adjustment

abc_MMD <- function(param, S, m_obs, M_epsilon = 100, transf = TRUE, mismatch = TRUE, B = 4e9, Ns = 801){
  
  M <- dim(param)[1]  # number of parameter values sampled from the prior
  nParam <- dim(param)[2] # total number of parameters to be estimated
  
  acceptedSamples <- array(0, dim = c(M_epsilon, dim(param)[2]))
  distance <- array(0, dim = c(M,1))
  
  # Computing lengthscale of Gaussian Kernel
  a <- fields::rdist(m_obs, compact = TRUE)
  lengthscale <- sqrt(median(a^2 / 2))
  
  # Computing MMD
  for (i in 1:M) {
    distance[i] <- MMD_fast(m_obs, S[,,i], lengthscale)
  }
  
  max_accepted_distance <- sort(distance, decreasing = F)[M_epsilon]
  
  acceptedSamples <- param[which(distance <= max_accepted_distance),] # Accepting the parameter samples with lowest MMD values
  acceptedDistances <- distance[which(distance <= max_accepted_distance)]
  rej <- acceptedSamples
  
  
  # Summarising temporal moments into means and covariances
  sumStats <- array(0, dim = c(M,9))
  for (i in 1:M) {
    sumStats[i,] <- summaryMoments(S[,,i])
  }
  
  # Summarising observed temporal moments
  if (mismatch == TRUE) {
    #####------Adjusting regression adjustment---####
    theta_hat <- array(0, dim = c(nParam,1))
    for (i in 1:nParam) {
      d <- density(acceptedSamples[,i])
      theta_hat[i] <- d$x[which.max(d$y)]
    }
    
    Y_new <- SalehValenzuela(theta_hat[2] , theta_hat[3] , theta_hat[4],
                             theta_hat[5], theta_hat[1], theta_hat[6], N=1000, B = B, Ns = Ns)
    
    m_obs_new <- log(temporalMomentsGeneral(Y_new, 3))
    sumStats_obs <- summaryMoments(m_obs_new)
  } else {
    sumStats_obs <- summaryMoments(m_obs) 
  }
  ######---- Performing regression adjustment------#####
  # Logit transform
  if (transf == TRUE) {
    for (k in 1:6) {
      acceptedSamples[,k] <- (acceptedSamples[,k] - prior_start[k]) / (prior_end[k] - prior_start[k])
      acceptedSamples[,k] <- log(acceptedSamples[,k] / (1 - acceptedSamples[,k]))
    }
  }
  
  norm_sumStats <- array(0,dim = dim(sumStats))
  
  norm_sumStats_obs <- array(0, dim = c(length(sumStats_obs), 1))
  
  for (i in 1:length(sumStats_obs)) {
    norm_sumStats[,i] <- sumStats[,i] / mad(sumStats[,i]) 
    norm_sumStats_obs[i] <- sumStats_obs[i] / mad(sumStats[,i]) 
  }
  
  norm_sumStats_accepted <- norm_sumStats[which(distance <= max_accepted_distance),]
  
  s_obs_norm <- matrix(data = rep(norm_sumStats_obs, M), nrow =  M_epsilon, ncol = length(norm_sumStats_obs), byrow = TRUE)
  
  X <- cbind(array(1, dim = c(M_epsilon,1)), norm_sumStats_accepted - s_obs_norm)
  
  weights <- 1 - (acceptedDistances / max_accepted_distance)^2
  
  W <- diag(weights)
  
  solution <- ginv(t(X) %*% W %*% X) %*% (t(X) %*% W %*% acceptedSamples)
  
  beta <- solution[-1,]
  
  theta_adjusted <- acceptedSamples - (norm_sumStats_accepted - s_obs_norm) %*% beta
  
  if (transf == TRUE) {
    for (k in 1:6) {
      theta_adjusted[,k] <- exp(theta_adjusted[,k]) / (1 + exp(theta_adjusted[,k]))
      theta_adjusted[,k] <- theta_adjusted[,k] * (prior_end[k] - prior_start[k]) + prior_start[k]
      
    }
  }
  
  return(list("rej" = rej, "reg" = theta_adjusted))  
}
