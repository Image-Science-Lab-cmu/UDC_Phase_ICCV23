function [params] = init_system(f, f1, f2, openRatio, pitch, ...
    sensorPitch, magnification, DOE_m, DOE_lambda0)


% main lens aperture [m]
if ~exist('f', 'var')
    params.f = 4.67e-3;                          % main lens focal length (m)
else
    params.f = f;
end
params.lensApertureSize = params.f/2;            % FNumber 2

if ~exist('sensorPitch', 'var')
    params.sensorPitch = 2e-6;              % sensor pitch of a conventional camera 
else
    params.sensorPitch = sensorPitch;       % sensor pitch of customized camera
end

if ~exist('magnification', 'var')
    params.magnification = 1;              % no magnification
else
    params.magnification = magnification;       % specified magnification
end

params.n = 1.5262;                         % refractive index of propagating media

params.delta1 = 1.5e-8;                    % spacing on display plane [m]
params.delta2 = 1.5e-8;                    % spacing on lens plane [m]
params.delta3 = 1.45e-7;                   % spacing on sensor plane [m]
params.delta3 = params.delta3;             % keep M at 2^20

params.sensorSize = 1024*2e-6;              % diam of sensor [m] (1024*sensorPitch)
params.sensorRes = round(params.sensorSize / params.sensorPitch);                     % number of pixels on sensor
params.scale = round(params.sensorPitch / params.delta3);
                                           % downsample scale
params.pitch = pitch;                      % diam of display pixel [m]
params.openRatio = openRatio;              % open ratio of display pixel [m]

                             
params.f1 = f1;                             % first microlens focal length [m]
params.f2 = f2;                             % second microlens focal length [m]
params.DOE_m = DOE_m;                       % max height of phase plate = DOE_m * lambda_0
params.DOE_lambda0 = DOE_lambda0;           % lambda0 for each microlens

% Distance between one phase mask and display plane
params.z=f1;                                % Fresnel prop distance [m]
params.displayApertureSize = ...            % Display aperture size (slightly larger than lens aperture)
    params.lensApertureSize + 2*params.z*tan(50/180*pi);
end
