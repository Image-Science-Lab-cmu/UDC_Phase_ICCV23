function [I3_lensOnly] = propagate_through_UDC_TOLED(params, thetas, lambda)
[root, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(root, '../utils'));
addpath(fullfile(root, '../fourier_optics_package'));

%% optical system parameters
batchSize = length(thetas);     % number of samples
k = 2*pi/lambda;                % wavenumber
M = floor(lambda*params.f/params.delta2/params.delta3/2)*2;

%% Display plane: initialize aperture and display pattern
% display aperture
x1 = (-M/2: 1: M/2-1) * params.delta1;
displayAperture = rect(x1/params.displayApertureSize);
displayApertureM = sum(displayAperture);

% load OLED display
display = load_display(displayApertureM, params.delta1, params.pitch, params.openRatio);

% pad patterns to M
padM = ceil((M-displayApertureM)/2);
display = padarray(display, [0,padM], 1,'both');
display = display(1:M);


% main lens aperture
mainLensAperture = rect(x1/params.lensApertureSize);
u1_1 = display .* mainLensAperture;
u1_1 = repmat(u1_1, [batchSize, 1]);
clear display mainLensAperture;

% Tilted wavefront
thetas = 90-thetas;         % propagating direction (degree)
alphas = cos(thetas / 180 * pi);
tilted_phase = exp(1i*k*alphas'*x1);

u1_1 = u1_1 .* tilted_phase;
clear tilted_phase;

% lens propagation
[u3_lensOnly,~,~] = propFF_1D_batch(u1_1, M*params.delta1, lambda, params.f, 0);
I3_lensOnly = abs(u3_lensOnly).^2;
clear u3_lensOnly;

%% crop PSF to sensor region and downsample to sensor resolution
sensorRes = round(params.sensorSize / params.sensorPitch);
x3 = (-M/2: 1: M/2-1);
sensorWindow = rect(x3 / (sensorRes * params.scale));
I3_lensOnly = I3_lensOnly(:, sensorWindow);
I3_lensOnly = downSample1D(I3_lensOnly, params.scale);
I3_lensOnly = gather(I3_lensOnly);