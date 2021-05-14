
function res = kernel_matrix(Z, bw)

D = squareform(pdist(Z));

res = exp(- D.^2 ./ bw^2);
