function y = Dy(M)
% Input:  M --- [m, n]
% Output: y --- [m-1, n]
    y = M(2:end, :) - M(1:end-1, :);
end