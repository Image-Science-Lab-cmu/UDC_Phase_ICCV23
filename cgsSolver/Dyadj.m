function x = Dyadj(y)
% Input:  y --- [m-1, n]
% Output: x --- [m, n]
    x = -y;
    x(end+1, :) = 0;
    x(2:end,:) = x(2:end,:) + y;
end