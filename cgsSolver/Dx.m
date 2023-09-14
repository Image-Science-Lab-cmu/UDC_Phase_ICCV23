function y = Dx(M)

% Input:  M --- [m, n]
% Output: y --- [m, n-1]
    y = M(:, 2:end) - M(:, 1:end-1);
end