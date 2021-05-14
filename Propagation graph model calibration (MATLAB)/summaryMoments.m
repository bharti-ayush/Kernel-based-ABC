function y = summaryMoments(x)

target_size = size(x);
if target_size(1) < target_size(2)
   x = x'; 
end

mu = mean(x,1);

C = triu(cov(x));

y = C(:);
y(y==0) = [];

y = [mu';y];


