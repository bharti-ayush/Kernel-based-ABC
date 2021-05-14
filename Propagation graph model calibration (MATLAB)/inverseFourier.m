%% Converting transfer function to impulse response 
function y_noisy = inverseFourier(Y)
y_noisy = zeros(size(Y));
dataDim = size(Y);
% 
% if length(dataDim) == 4
%     dataDim(5) = 1;
% end

% for ii=1:dataDim(end)
    for uu = 1:dataDim(1)
            for vv = 1:dataDim(2)
%                 for jj = 1:dataDim(1)
                
                    y_noisy(uu,vv,:) = ifft((squeeze(Y(uu,vv,:))));
                    
%                 end
            end
    end

% end