function x = Dxadj(y)
% Input:  y --- [m, n-1]
% Output: x --- [m, n]
    x = -y;
    x(:, end+1) = 0;
    x(:,2:end) = x(:,2:end) +y;
end