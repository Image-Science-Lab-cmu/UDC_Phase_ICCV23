function[u2,L2,x2]=propFF_1D_batch(u1,L1,lambda,f,d)
% propagation - Fraunhofer pattern
% assumes uniform sampling
% u1 - source plane field
% L1 - source plane side length
% lambda - wavelength
% f - propagation distance
% L2 - observation plane side length
% u2 - observation plane field
%
[~, M]=size(u1); %get input field array size
dx1=L1/M; %source sample interval
% k=2*pi/lambda; %wavenumber
% %
% L2=lambda*f/dx1; %obs sidelength
% dx2=lambda*f/L1; %obs sample interval
% % x2=-L2/2:dx2:L2/2-dx2; %obs coords
% x2 = (0:1:M-1)*dx2 - L2/2; %obs coords
% if gpuDeviceCount > 0; x2 = gpuArray(x2); end

% todo: since we only measure intensity, omit pure phase change.
% c=1/(1j*lambda*f)*exp(1j*k*(1-d/f)/(2*f)*(x2.^2));
c=1/(1j*lambda*f);

u2=ifftshift(fft(fftshift(u1, 2), [], 2), 2)*dx1;
u2=c.*u2;

% placeholder
L2 = 0;
x2 = 0;
end