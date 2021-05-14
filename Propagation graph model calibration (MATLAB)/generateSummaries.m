%% This function generates channel realizations for the given set of parameter values
% param:        set of parameter values
% K:            number of parameters sampled from the prior
% t:            time-axis
% nMoments:     number of temporal moments to compute (set as 4 in the paper)

function [Dat, rho] = generateSummaries(param, K, t, nMoments)

%%
Dat = zeros(625, nMoments,K);
% finalStats = zeros(10,K);

for i = 1:K
    % generate transfer functions from the PG model
       [H,~, rho] =  graphModelStatistics(param(1:4,i),  1);
       % Add noise to the transfer functions
       Y = H + randn(size(H)) * sqrt(param(4,i) / 2) + 1i * randn(size(H)) * sqrt(param(4,i) / 2);
       disp(i)
       
       y = inverseFourier(Y); % Transform frequency-domain data to time-domain
       
       h_vv = reshape(y, 625,801);
       % Compute the temporal moments
       Dat(:,:,i) = computemoments(t',h_vv, nMoments);
              
end
    

