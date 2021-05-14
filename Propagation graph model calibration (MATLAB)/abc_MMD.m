function [acceptedSamples, theta_adjusted] = abc_MMD(param, S, m_obs, nacc, transf, bounds)

target_size = size(m_obs);
if target_size(1) < target_size(2)
   m_obs = m_obs'; 
end

param_size = size(param);
if param_size(1) < param_size(2)
   param = param'; 
end

% Rejection ABC using MMD
[M, nParam] = size(param);
distance = zeros(M,1);

a = pdist(m_obs);
lengthscale = sqrt(median(a.^2 / 2));

for i=1:M
   distance(i) = MMD(m_obs, squeeze(S(:,:,i)), lengthscale); 
   disp(i)
end

[sorted_distance, ~] = sort(distance);
sorted_distance = sorted_distance(1:nacc);
max_accepted_distance = sorted_distance(nacc);
samples_rej = param(distance <= max_accepted_distance,:);
acceptedSamples = samples_rej;

% Regression adjustment
obs_stats = summaryMoments(m_obs');

I = target_size(2);
nStatistics = (I^2 + 3*I) / 2;
sumStats = zeros(M,nStatistics);
for i = 1:M
   sumStats(i,:) = summaryMoments(squeeze(S(:,:,i))' );
end

if transf == 'logit'
    for k = 1:nParam
        samples_rej(:,k) = ( samples_rej(:,k) - bounds(k,1) ) / (bounds(k,2) - bounds(k,1));
        samples_rej(:,k) = log(samples_rej(:,k) ./ (1 - samples_rej(:,k)));
    end
end

    norm_factor = 1.4826 * mad(sumStats, 1);

    norm_sumstats = zeros(size(sumStats));
    
    for i=1:nStatistics
        norm_sumstats(:,i) = sumStats(:,i) / norm_factor(i);
        obs_stats_norm(i) = obs_stats(i) / norm_factor(i) ;
    end
    
%     sum1 = 0;
%     
%     for j = 1:nStatistics
%         sum1 = sum1 + (norm_sumstats(:,j) - obs_stats_norm(j)).^2;
%     end
%     
%     distance = sqrt(sum1);
% 
%     [sorted_distance, ~] = sort(distance);
%     sorted_distance = sorted_distance(1:nacc);
%     max_accepted_distance = sorted_distance(nacc);
    theta_star = samples_rej;
%     theta_star = param(distance <= max_accepted_distance,:);

    norm_sumstats_star = norm_sumstats(distance <= max_accepted_distance,:);

    s_obs_norm = repmat(obs_stats_norm, nacc, 1);

    X = [ones(nacc, 1) , norm_sumstats_star - s_obs_norm];

    weights = 1 - (distance(distance <= max_accepted_distance) / sorted_distance(end)).^2; %Epanechnikov kernel

    W = diag(weights);

    solution = (X' * W * X) \ (X' * W * theta_star);
    beta = solution(2:end,:);

     theta_adjusted = theta_star - (norm_sumstats_star - s_obs_norm) * beta;

if transf == 'logit'
    for k = 1:nParam
        theta_adjusted(:,k) = exp(theta_adjusted(:,k)) ./ (1 + exp(theta_adjusted(:,k))) ;
        theta_adjusted(:,k) = theta_adjusted(:,k) * (bounds(k,2) - bounds(k,1) ) + bounds(k,1);
    end
end

end