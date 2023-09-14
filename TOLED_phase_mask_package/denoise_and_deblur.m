function [lin_deblur] = denoise_and_deblur(lin_blur, omega, operator, ...
    noise_var, lambda)
addpath bm3d_matlab_package
addpath bm3d_matlab_package/bm3d
addpath cgsSolver

%% denoise
if noise_var ~= 0
    noiseProfile = BM3DProfile();
    noiseProfile.gamma = 0;
    noise_type =  'gw';
    seed = 0; % seed for pseudorandom noise realization
    [~, PSD, ~] = getExperimentNoise(noise_type, noise_var, seed, ...
        size(lin_blur));
    lin_blur = CBM3D(lin_blur, PSD, noiseProfile);
end
    
%% deblur
% Using cgs to solve both spatially varying and invariant blur.
lin_deblur = LS_L2Gradient_Deblurring_3Ch(lin_blur, omega, [], ...
    lambda, operator, []);