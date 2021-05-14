% This is the main file for running the kernel-based ABC calibration of
% propagation graph model from A. Bharti, F.-X. Briol, and T. Pedersen, 
% “A general method for calibrating stochastic radio channel models with kernels,”
% IEEE Trans. on Antennas and Propag., 2021.(available at https://arxiv.org/abs/2012.09612).

close all;
clc;
clear all;
tic
%% Loading data
load('param1.mat'); 
load('data_moments.mat');
nMoments = 4; % number of temporal moments to be used for calibration
param = param(:,1:2000); % parameter values sampled from prior
S = Moments(:,1:nMoments,1:2000); 
S = log(S); % Corresponding log temporal moments simulated from the graph model
%% Parameters
B = 4e9; % Bandwidth
numR = 625; % number of realizations of the channel for 5x5 virtual array at both Tx and Rx
Ns = 801; % number of frequency points in the data
K = 2000; % number of parameter values sampled from the prior

delta_t = 1/B; % time resolution

% True parameter values
g_true = 0.6; % reflection coefficient
N_true = 15; % number of scatterers
Pvis_true = 0.6; % Probability of visibility
sigma_N_true = 1e-9; % noise variance

numR_obs = 625; % number of realizations in observed data
nacc = 100; % number of accepted parameter samples in rejection ABC
T = 10; % total number o iterations of PMC-ABC 
nParam = 4; % number of parameters to estimate
prob_factor = 1e4; % used to generate new population in PMC-ABC

% Priors
g_min = 0; g_max = 1;
N_min = 5; N_max = 35;
P_min = 0;  P_max = 1;
sigma_N_min = 2e-10; sigma_N_max = 2e-9;

prior_min = [g_min; N_min; P_min; sigma_N_min];
prior_max = [g_max; N_max; P_max; sigma_N_max];

bounds = [g_min g_max; N_min-1 N_max+1; P_min P_max; sigma_N_min sigma_N_max];

param_true = [g_true; N_true; Pvis_true; sigma_N_true];

Taxis = (0:Ns-1)*delta_t; % time axis

xaxis_names = {'Reflection coefficient, g', 'No. of scatterers, N_{s}', 'Probability of visibility, P_{vis}', 'Noise variance, \sigma_N^2'};

xlimits = [g_min g_max; N_min N_max; P_min P_max; sigma_N_min sigma_N_max];

%% Generating observed statistics (uncomment to generate pseudo-observed data again)
% m_obs1 = log(generateSummaries(param_true, 1, Taxis, nMoments));
% m_obs2 = log(generateSummaries(param_true, 1, Taxis, nMoments));
% m_obs3 = log(generateSummaries(param_true, 1, Taxis, nMoments));
% m_obs4 = log(generateSummaries(param_true, 1, Taxis, nMoments));
% 
% m_obs = [m_obs1; m_obs2; m_obs3; m_obs4];
% save('m_obs_simulated', 'm_obs');

load m_obs_simulated;
%% Running ABC
[abc_samples_rejection, abc_samples] = abc_MMD(param, S, m_obs, nacc, 'logit', bounds); 

save(['samples_rejection', num2str(1)], 'abc_samples_rejection');
save(['samples_regression', num2str(1)], 'abc_samples');
abc_plot(abc_samples, xaxis_names, xlimits, param_true, 1);
print(['plot_simulated', num2str(1)],'-dpdf')

%% Generating next population of samples
sd_theta = sqrt(2) * std(abc_samples);
Sigma = diag(sd_theta.^2);

ind = datasample(1:length(abc_samples(:,1)), K);

theta_accepted = abc_samples;
theta_star = zeros(K, nParam);
theta_perturbed = zeros(K, nParam);

for i = 1:nParam
   for j = 1:K
       
      theta_star(j,i) = abc_samples(ind(j), i);
      
      while 1
          theta_perturbed(j,i) = theta_star(j,i) + sd_theta(i) * randn(1);
          if theta_perturbed(j,i) > prior_min(i) && theta_perturbed(j,i) < prior_max(i)
              break
          end
      end
       
   end
end
theta_perturbed(:,2) = round(theta_perturbed(:,2));
save(['thetaPerturbed', num2str(1)], 'theta_perturbed');

old_w = ones(nacc, 1) .* 1 / (nacc);
old_samples = abc_samples;

%% Running sequential ABC
for t = 2:T
    %% Generating new channel responses
    S = generateSummaries(theta_perturbed', K, Taxis, nMoments);
    S = log(S);
    
    save(['S', num2str(t)], 'S');
     
    %% Perform abc
    [abc_samples_rejection, abc_samples] = abc_MMD(theta_perturbed, S, m_obs, nacc, 'logit', bounds); 
    
    save(['samples_regression', num2str(t)], 'abc_samples');
    save(['samples_rejection', num2str(t)], 'abc_samples_rejection');
    
    % Plot results
    abc_plot(abc_samples, xaxis_names, xlimits, param_true, t);
    print(['plot_simulated', num2str(t)],'-dpdf')

    %% Compute weights for next population
    
    w = zeros(nacc, 1);
    
    for l = 1:nacc
        w(l) = 1 / (sum(old_w .* mvnpdf(abc_samples(l,:), old_samples(1:nacc,:), Sigma)));
    end
    
    w_norm = w./sum(w);   % Normalizing weights to be probabilities
    
    save(['weights', num2str(t)], 'w_norm');
    
    sd_theta = sqrt(2) * std(abc_samples);    % Computing standard deviation
    Sigma = diag(sd_theta.^2);
    
    theta_star = abc_samples(randsample(nacc, K, true, w_norm), :);
      
    for i = 1:nParam
        for j = 1:K
            while 1
                theta_perturbed(j,i) = theta_star(j,i) + sd_theta(i) * randn(1);
                if theta_perturbed(j,i) > prior_min(i) && theta_perturbed(j,i) < prior_max(i)
                    break
                end
            end
        end
    end
    theta_perturbed(:,2) = round(theta_perturbed(:,2));
    save(['thetaPerturbed', num2str(t)], 'theta_perturbed');
    
    old_w = w_norm;
    old_samples = abc_samples;
    
end
toc