function dn_vec = downSample1D(vec, ratio)

N = size(vec, 2);
N = N - mod(N, ratio);
vec = vec(:, 1: N);

% dn_vec = resample(vec, 1, ratio);
vec = imfilter(vec, ones(1, ratio)/ ratio);
% mid = floor(ratio / 2) + 1;
% Anqi @ Oct 27 to keep consistent with python version
if mod(ratio, 2) == 0
    mid = floor(ratio / 2);
else
    mid = floor(ratio / 2) + 1;
end 
dn_vec = vec(:, mid:ratio:end);

% Note: python conv2d samples the center point of a 
% convolution filter, while we take the first sample
% in Matlab previously.