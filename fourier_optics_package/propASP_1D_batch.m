function [x2, Uout] = propASP_1D_batch(Uin, wvl, d1, d2, Dz)
% function [x2 y2 Uout] ...
% = ang_spec_prop(Uin, wvl, d1, d2, Dz)
ft = @(x, delta) fftshift(fft(fftshift(x, 2), [], 2), 2) * delta; % apply ft to each row of x
ift = @(x, delta_f) ifftshift(ifft(ifftshift(x, 2), [], 2), 2) * length(x) *(delta_f); % apply ift to each row of x

[~, N] = size(Uin); % assume square grid
k = 2*pi/wvl; % optical wavevector
% source-plane coordinates
x1 = (-N/2 : 1 : N/2 - 1) * d1;
if gpuDeviceCount > 0; x1 = gpuArray(x1); end
r1sq = x1.^2;

% spatial frequencies (of source plane)
df1 = 1 / (N*d1);
fX = (-N/2 : 1 : N/2 - 1) * df1;
if gpuDeviceCount > 0; fX = gpuArray(fX); end
fsq = fX.^2;

% scaling parameter
m = d2/d1;
% observation-plane coordinates
x2 = (-N/2 : 1 : N/2 - 1) * d2;
if gpuDeviceCount > 0; x2 = gpuArray(x2); end
r2sq = x2.^2;

% quadratic phase factors
Q1 = exp(1i*k/2*(1-m)/Dz*r1sq);
Q2 = exp(-1i*pi^2*2*Dz/m/k*fsq);
Q3 = exp(1i*k/2*(m-1)/(m*Dz)*r2sq);

% compute the propagated field
Uout = Q3.* ift(Q2 .* ft(Q1 .* Uin / m, d1), df1);
% Uout = ift(Q2 .* ft(Q1 .* Uin / m, d1), df1);