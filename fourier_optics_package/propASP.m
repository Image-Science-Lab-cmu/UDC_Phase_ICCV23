function [x2, y2, Uout] = propASP(Uin, wvl, d1, d2, Dz)
% function [x2 y2 Uout] ...
% = ang_spec_prop(Uin, wvl, d1, d2, Dz)
ft2 = @(x, delta) fftshift(fft2(fftshift(x))) * delta^2;
ift2 = @(x, delta_f) ifftshift(ifft2(ifftshift(x))) * (size(x,1) * delta_f)^2;

N = size(Uin,1); % assume square grid
k = 2*pi/wvl; % optical wavevector
% source-plane coordinates
[x1, y1] = meshgrid((-N/2 : 1 : N/2 - 1) * d1);
% if gpuDeviceCount > 0; x1 = gpuArray(x1); y1=gpuArray(y1); end
r1sq = x1.^2 + y1.^2;
% spatial frequencies (of source plane)
df1 = 1 / (N*d1);
[fX, fY] = meshgrid((-N/2 : 1 : N/2 - 1) * df1);
% if gpuDeviceCount > 0; fX = gpuArray(fX); fY=gpuArray(fY); end
fsq = fX.^2 + fY.^2;
% scaling parameter
m = d2/d1;
% observation-plane coordinates
[x2, y2] = meshgrid((-N/2 : 1 : N/2 - 1) * d2);
% if gpuDeviceCount > 0; x2 = gpuArray(x2); y2=gpuArray(y2); end
r2sq = x2.^2 + y2.^2;
% quadratic phase factors
Q1 = exp(1i*k/2*(1-m)/Dz*r1sq);
Q2 = exp(-1i*pi^2*2*Dz/m/k*fsq);
Q3 = exp(1i*k/2*(m-1)/(m*Dz)*r2sq);
% compute the propagated field
Uout = Q3.* ift2(Q2 .* ft2(Q1 .* Uin / m, d1), df1);