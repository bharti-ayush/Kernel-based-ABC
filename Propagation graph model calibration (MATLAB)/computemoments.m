
function M = computemoments(tt,H, nMoments)
 [N1,N2]=size(H);
if N1 < N2
   H =  H';
end

[~,N]=size(H);
   
M = zeros(N, nMoments);

for k = 1:nMoments
    for ii = 1:N
         P = abs(H(:,ii)).^2;
         M(ii, k) = trapz(tt,tt.^(k-1) .* P);
    end
end

end
