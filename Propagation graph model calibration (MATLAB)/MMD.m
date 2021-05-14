% This function computes the MMD between X and Y for a given lengthscale bw
function estimated_MMD = MMD(X, Y, bw)

[m,~] = size(Y);
[n,~] = size(X);

Z = [X;Y];

B = kernel_matrix(Z, bw);

K_XX = B(1:n, 1:n);
K_YY = B((n+1) : (n+m), (n+1) : (n+m));
K_XY = B(1:n, (n+1) : (n+m));

one_n = ones(n,1);
one_m = ones(m,1);

estimated_MMD = 1 / (n * (n-1)) * ( one_n' * K_XX * one_n ) - 2 / (n * m) * (one_n' * K_XY * one_m) + 1 / (m * (m-1)) * (one_m' * K_YY * one_m);
