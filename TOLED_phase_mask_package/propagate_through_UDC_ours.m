function [I3] = propagate_through_UDC_ours(params, thetas, lambda)

[root, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(root, '../utils'));
addpath(fullfile(root, '../fourier_optics_package'));

%% Prepare the setup
% optical system parameters
batchSize = length(thetas); 
k = 2*pi/lambda;                % wavenumber
M = floor(lambda*params.f/params.delta2/params.delta3/2)*2;

%  Load display aperture and display pattern
x1 = (-M/2: 1: M/2-1) * params.delta1;
displayAperture = rect(x1/params.displayApertureSize);
displayApertureM = sum(displayAperture);
display = load_display(displayApertureM, params.delta1, params.pitch, params.openRatio);

% Load phase mask profiles
microLensArray1 = load_microLensArray(k, displayApertureM, params.delta1, params.pitch, params.f1, params.DOE_m, params.DOE_lambda0);
microLensArray2 = load_microLensArray(k, displayApertureM, params.delta1, params.pitch, params.f2, params.DOE_m, params.DOE_lambda0);

% Pad display and phase mask to required number of samples M
padM = ceil((M-displayApertureM)/2);
display = padarray(display, [0,padM], 1,'both');
microLensArray1 = padarray(microLensArray1, [0,padM], 1,'both');
microLensArray2 = padarray(microLensArray2, [0,padM], 1,'both');
display = display(1:M);microLensArray1=microLensArray1(1:M); microLensArray2=microLensArray2(1:M);

% Prepare input wavefront
u1 = ones(1, M);
u1 = u1 .* displayAperture;
u1 = repmat(u1, [batchSize, 1]);
% tilted wavefront from angle thetas
thetas = 90-thetas;         % propagating direction (degree)
alphas = cos(thetas / 180 * pi);
tilted_phases = exp(1i*k*alphas'*x1);
u1 = u1 .* tilted_phases;
clear tilted_phases;

%% LensArray1: phase mask before display plane
u1 = u1 .* microLensArray1;
clear microLensArray1;

%% Fresnel Propagation: from phase mask to display plane by f1 (m)
[~,u2] = propASP_1D_batch(u1, lambda, params.delta1, params.delta2, params.f1);
clear u1;

%% Modulated by display plane 
u2 = u2 .* display .* displayAperture;

%% Fresnel Propagation: from display plane to the second phase mask by f2 (m)
[~,u2] = propASP_1D_batch(u2, lambda, params.delta1, params.delta2, params.f2);

%% LensArray2: phase mask behind display plane
mainLensAperture = rect(x1/params.lensApertureSize);
u2 = u2 .* microLensArray2 .* mainLensAperture;
clear microLensArray2 mainLensAperture;
   
%% Focus by main lens
%  propagate from main lens to sensor.
[u3,~,~] = propFF_1D_batch(u2, M*params.delta2, lambda, params.f, 0);
I3=abs(u3).^2;
clear u2 u3;

%% crop PSF to sensor region and downsample to sensor resolutions
%  Crop to sensor region of interest [m]
x3 = (-M/2: 1: M/2-1);
sensorWindow = rect(x3 / (params.sensorRes * params.scale));
I3 = I3(:, sensorWindow);
%  Downsample to sensor resolution
I3 = downSample1D(I3, params.scale);
I3 = gather(I3);