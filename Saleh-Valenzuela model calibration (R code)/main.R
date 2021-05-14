library(rmatio)
library(pracma)
library(abc)
library(MASS)
library(mvtnorm)

# Loading Data
param <- readRDS(file = "param.rds")  # parameter values sampled from the prior for the first iteration of ABC algorithm
S <- readRDS(file = "S.rds")           # temporal moments simulated based on the sampled parameter values
# S <- S[1:100,,] # Taking the number of realizations of temporal moments to be N_sim = 100

#######-------- Parameter Settings -------########

B = 4e9 # Measurement Bandwidth
Ns = 801 # Number of frequency/time points 

N <- 100  # Number of realizations of simulated data
M <- 2000  # Number of simulated samples
M_epsilon <- 100  # Number of accepted samples
nParameters <- 6  # Number of parameters
T <- 10  # Number of iterations
tolerance <- 0.05
nSimulations1 <- 2000 # Number of parameter values sampled in the subsequent iterations of PMC-ABC

# Prior ranges
beta_start <- 1e-9
beta_end <- 1e-7

lambda_cluster_start <- 5e6
lambda_cluster_end <- 1e8

lambda_ray_start <- 5e6
lambda_ray_end <- 3e9

gamma_cluster_start <- 5e-9
gamma_cluster_end <- 5e-8

gamma_ray_start <- 5e-10
gamma_ray_end <- 5e-9

sigma2_N_start <- 2e-10
sigma2_N_end <- 2e-9

prior_start <- c(beta_start, lambda_cluster_start, lambda_ray_start, gamma_cluster_start, gamma_ray_start, sigma2_N_start)
prior_end <- c(beta_end, lambda_cluster_end, lambda_ray_end, gamma_cluster_end, gamma_ray_end, sigma2_N_end)


############----------GENERATING OBSERVED DATA----##########

# True values of the parameters
betaSquare_00_true <- 5e-8
lambda_cluster_true <- 2e7
lambda_ray_true <- 1e9
gamma_cluster_true <- 1e-8
gamma_ray_true <- 2e-9
sigma2_N_true <- 1e-9
# 
param_true <- c(betaSquare_00_true, lambda_cluster_true, lambda_ray_true, gamma_cluster_true, gamma_ray_true, sigma2_N_true)

# Generating pseudo-observed data from the S-V model
Y_obs <- SalehValenzuela(lambda_cluster_true , lambda_ray_true , gamma_cluster_true,
                         gamma_ray_true, betaSquare_00_true, sigma2_N_true, N=1000, B = B, Ns = Ns)

t <- 1 # Iteration counter
K <- 4 # Number of temporal moments to be used  

m_obs <- log(temporalMomentsGeneral(Y_obs, K)) # Computing observed temporal moments
saveRDS(m_obs, "m_obs_sim")

# Running the rejection ABC with regression adjustment
samples <- abc_MMD(param, S, m_obs, M_epsilon, transf = TRUE, mismatch = FALSE, B = B, Ns = Ns)
acceptedSamples <- samples$reg # Taking the adjusted parameter samples as the output of the function

# Saving parameter samples
saveRDS(samples$rej, file = paste0("samples_lund_rejection", t))
saveRDS(samples$reg, file = paste0("samples_regression_lund", t))

# Saving plots

pdf("plot_rejection_lund1.pdf", paper="USr", width = 10)
plot_estimates(samples$rej, param_true) # Plotting rejection ABC samples
dev.off()

pdf("plot_regression_lund1.pdf", paper="USr", width = 10)
plot_estimates(samples$reg, param_true) # Plotting regression adjusted samples
dev.off()


###########----------Population Monte Carlo ABC---------##########

sd_theta <- apply(acceptedSamples, 2, sd) * sqrt(2)
Sigma <- diag(sd_theta^2)
old_w <- array(1/M_epsilon, dim = c(M_epsilon,1)) # Equal weights for the 2nd iteration

ind <- sample(x = 1:length(acceptedSamples[,1]), size = M, replace = TRUE)

theta_accepted <- acceptedSamples

# Sampling new population of parameter values from the output of previous iteration
theta_star <- array(0, dim = c(M, nParameters))
theta_perturbed <- array(0, dim = c(M, nParameters))

for (i in 1:nParameters) {
  
  for (j in 1:M) {
    
    theta_star[j,i] <- theta_accepted[,i][ind[j]]
    
    repeat{
      theta_perturbed[j,i] <- rnorm(1, mean = theta_star[j,i], sd = sd_theta[i])
      if (theta_perturbed[j,i] >0 & theta_perturbed[j,i] > prior_start[i] & theta_perturbed[j,i] < prior_end[i]) {
        # if (theta_perturbed[j,i] >0){
        break
      }
    }
  }
}

old_samples <- samples$reg

# Saving new population of parameters
pdf(paste0("theta_perturbed_lund", t, ".pdf"), paper="USr", width = 10)
plot_estimates(theta_perturbed, param_true)
dev.off()
# 
saveRDS(theta_perturbed, file = paste0("theta_perturbed_lund", t))


ptm <- proc.time()

for (t in 2:T) {
  
  # Generating channel transfer functions for the new population of parameters
  channel <- array(0, dim = c(M, N, Ns))
  
  for (j in 1:M) {
    
    channel[j, , ] <- SalehValenzuela(lambda_cluster = theta_perturbed[j,2], lambda_ray =  theta_perturbed[j,3],
                                      gamma_cluster = theta_perturbed[j,4], gamma_ray = theta_perturbed[j,5],
                                      betaSquare_00 = theta_perturbed[j,1], sigma2_N = theta_perturbed[j,6], N = N, B = B, Ns = Ns)
    
    print(paste0(j, ":", t)) # Counter to keep track of completed simulations
    
  }
  
  # Computing temporal moments
  S <- array(0, dim = c(N, K ,M))
  
  for (i in 1:M) {
    S[,,i] <- log(temporalMomentsGeneral(channel[i,,], K))
  }
  
  saveRDS(S, file = paste0("S_lund", t))
  
  # rejection ABC with regression adjustment
  samples <- abc_MMD(theta_perturbed, S, m_obs, M_epsilon, transf = TRUE, mismatch = FALSE, B = B, Ns = Ns)
  
  # Saving parameter samples
  saveRDS(samples$rej, file = paste0("samples_rejection_lund", t))
  saveRDS(samples$reg, file = paste0("samples_regression_lund", t))
  
  # Saving plots
  
  pdf(paste0("plot_rejection_lund", t, ".pdf"), paper="USr", width = 10)
  plot_estimates(samples$rej, param_true)
  dev.off()
  
  pdf(paste0("plot_regression_lund", t, ".pdf"), paper="USr", width = 10)
  plot_estimates(samples$reg, param_true)
  dev.off()
  
  acceptedSamples <- samples$reg
  
  # Computing the weights
  w <- array(0, dim = c(M_epsilon, 1))
  
  for (i in 1:M_epsilon) {
    temp <- 0
    for (j in 1:M_epsilon) {
      temp <- old_w[j] *  dmvnorm(acceptedSamples[i,], mean = old_samples[j,], Sigma) + temp
    }
    w[i] <- 1 / temp
  }
  
  w_norm <- w / sum(w) # Normalizing the weights to sum to 1
  
  saveRDS(w_norm, file = paste0("weights", t))
  
  sd_theta <- apply(acceptedSamples, 2, sd) * sqrt(2)
  Sigma <- diag(sd_theta^2)
  
  indices <- sample.int(M_epsilon, size = M, replace = TRUE, prob = w_norm)
  
  theta_star <- acceptedSamples[indices,]
  
  for (i in 1:nParameters) {
    
    for (j in 1:M) {
      
      repeat{
        theta_perturbed[j,i] <- rnorm(1, mean = theta_star[j,i], sd =  sd_theta[i])
        
        if (theta_perturbed[j,i] >0 & theta_perturbed[j,i] > prior_start[i] & theta_perturbed[j,i] < prior_end[i]) {
          break
        }
      }
    }
  }
  saveRDS(theta_perturbed, file = paste0("theta_perturbed_lund", t))
  # Open a pdf file
  pdf(paste0("theta_perturbed_lund", t, ".pdf"), paper="USr", width = 10)
  plot_estimates(theta_perturbed, param_true)
  # Close the pdf file
  dev.off()
  
  old_w = w_norm
  old_samples = acceptedSamples
  
}
proc.time() - ptm
