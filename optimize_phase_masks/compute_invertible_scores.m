% This script computes the invertibility score matrix
% for different DPI
addpath ../;
simulationType = 'densePSF';
SAVE = 1;

% dense wavelength
% Note: Use spectral response with the same energy for RGB
load('spectral_response_new.mat');
load('PhasePlate_Params.mat');

for id = [1, 2, 7: 12]

    f = 4.67e-3;            % main lens focal length
    opening = 1 / PhasePlate_Params(id).lensOpenRatio;
    magnification = PhasePlate_Params(id).magnification;
    f1 = PhasePlate_Params(id).f1;
    f2 = PhasePlate_Params(id).f2;
    pitch = PhasePlate_Params(id).lensPitch;
    DOE_m = PhasePlate_Params(id).DOE_m;
    DOE_lambda0 = PhasePlate_Params(id).DOE_lambda0;
    sensorPitch = 2e-6;     % Sensor pixel size
    

    % Image system for an under-display camera
    % f, f1, f2, openRatio, pitch, ...
    % sensorPitch, magnification, DOE_m, DOE_lambda0
    params = init_system(f, f1, f2, 1 / opening, pitch, ...
        sensorPitch, magnification, DOE_m, DOE_lambda0);
    sensor = get_sensor(params);
    
    invertible_scores = zeros(length(SR.spectral(1:10:end)));
    
    % Design phase plate at reference wavelength lambda0
    for lambda0 = SR.spectral(1:10:end)
        
        % Use lambda0 for all lenslets
        params.DOE_lambda0 = lambda0;
    
        for wvl = SR.spectral(1: 10: end)
            % get sensor information under these parameters
            I3_TOLED = propagate_through_UDC_TOLED(params, 0, wvl);
            I3_TOLED_energy = sum(I3_TOLED);
            
            % compute PSF for frontal parallel incident light
            pixelIds = [0];
            thetas_y = sensor.pixel2angle(pixelIds);
            I3_ours = propagate_through_UDC_ours(params, thetas_y, wvl);
            fprintf('PhaseMask focal length=%dum, folding with lambda0=%.0fnm, incident wavelength=%.0fnm\n', ...
                f1 * 1e6, lambda0 * 1e9, wvl * 1e9);
            
            % normalize energy by I3_TOLED
            I3_ours = I3_ours / I3_TOLED_energy;
            I3_TOLED = I3_TOLED / I3_TOLED_energy;

            % compute invertible score for using lambda0
            invertible_scores(lambda0 == SR.spectral(1: 10: end), ...
                wvl == SR.spectral(1: 10: end)) = compute_inv(I3_ours);
        end
    end

    if SAVE
        save(sprintf('invertible_scores_pitch%d_m_%d.mat', ...
        pitch * 1e6, DOE_m), 'invertible_scores');
    end

end
% end

function invertible_score = compute_inv(k)
% k is PSF
% - of a single wavelength
% - not normalized

k = k / sum(k);
mtf = abs(ifftshift(fft(fftshift(k))));

len = length(k);
mtf = mtf(floor(len/2) + 1: len);
invertible_score = sum(mtf);
end