function recons = LS_L2Gradient_Deblurring_3Ch(y, omega, init, lambda, operator, noise_variance)

[H, W, ~] = size(y);        %shape of image
recons = zeros(H, W, 3);
if strcmp(operator, 'matmul')
    [~, nH, ~] = size(omega);
    recons = zeros(nH, W, 3);
end
for cc = 1:3
    if isempty(noise_variance)
        recons(:,:,cc) = leastSquareDeblurring(y(:,:,cc,:), omega(:,:,cc,:), init, lambda, operator, []);
    else
        recons(:,:,cc) = leastSquareDeblurring(y(:,:,cc,:), omega(:,:,cc,:), init, lambda, operator, noise_variance(:,:,cc,:));
    end
end

end



function recons = leastSquareDeblurring(y, omega, init, lambda, operator, noise_covariance)
cgs_iters = 150;  % Number of conjugate gradient descent iterations
cgs_tol = 1e-6;  % cgs tolerance
% lambda = 1e-3;   % Regularization parameter

% todo: read in blurry image
% img = imread('cameraman.tif');
% img = double(img)/255;
% 
[H, W] = size(y); %shape of image
[kH, kW] = size(omega);
% y = AOmega(img, omega);

% If Nvidia GPU available, send y and omega
% to GPU and all relavent computation will
% be conducted on GPU.
if gpuDeviceCount > 0 
    omega = gpuArray(omega);
    y = gpuArray(y);
end

vec = @(x) x(:);

if strcmp(operator, 'matmul')
    
    % y = Omega * x
    % y:    [H, W]
    % Omega [H, kH]
    % x:    [kH, W]
    
    [H, W] = size(y); %shape of image
    [~, kH] = size(omega);
    
    AOmega = @(x) omega*x;
    AOmega_adj = @(y) omega'*y;

    Aforward = @(x) vec(AOmega(reshape(x,kH,W)));
    Aadj = @(y) vec(AOmega_adj(reshape(y,H,W)));
    
    Dxforward = @(x) vec(Dx(reshape(x,kH,W)));
    DxAdj = @(x) vec(Dxadj(reshape(x,kH,W-1)));

    Dyforward = @(x) vec(Dy(reshape(x,kH,W)));
    DyAdj = @(x) vec(Dyadj(reshape(x,kH-1,W)));

elseif strcmp(operator, 'conv')
    AOmega = @(x) myConv2(x, omega, 'valid');
    AOmega_adj = @(x) myConv2(x, omega(end:-1:1, end:-1:1), 'full');
    
    Aforward = @(x) vec(AOmega(reshape(x,H+kH-1,W+kW-1)));
    Aadj = @(x) vec(AOmega_adj(reshape(x,H,W)));
    
    Dxforward = @(x) vec(Dx(reshape(x,H+kH-1,W+kW-1)));
    DxAdj = @(x) vec(Dxadj(reshape(x,H+kH-1,W+kW-1-1)));
    
    Dyforward = @(x) vec(Dy(reshape(x,H+kH-1,W+kW-1)));
    DyAdj = @(x) vec(Dyadj(reshape(x,H+kH-1-1,W+kW-1)));


else
    error('Error: undefined operator');
end


A = @(x) Aadj(Aforward(x)) + lambda * ((DxAdj(Dxforward(x))) + (DyAdj(Dyforward(x))));
b = Aadj(y(:));


if ~isempty(init)
    recons = cgs(A, b, cgs_tol, cgs_iters, [], [], vec(init));
else
    [recons, ~, ~, ~, resvec] = cgs(A, b, cgs_tol, cgs_iters);
end

switch operator
    case 'matmul'
        recons = reshape(recons, kH, W);
    case 'conv'
        recons = reshape(recons, H+kH-1, W+kW-1);
        recons = recons(floor(kH/2)+1: floor(kH/2)+H, floor(kW/2)+1:floor(kW/2)+W);
    otherwise
        error('Unknown operator');
end

if gpuDeviceCount > 0
    recons = gather(recons);
end
end