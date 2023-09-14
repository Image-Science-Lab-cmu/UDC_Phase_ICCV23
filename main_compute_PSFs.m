% This script computes and saves spatially-varying 1D PSFs for
% the proposed setup TOLED + PhaseMask

function [] = main_compute_PSFs(ids)

addpath('TOLED_phase_mask_package');
% PhasePlate_Params saves several setup configurations used in the paper.
% You can specify configurations by setting 'ids'.
Params = load('TOLED_phase_mask_package/PhasePlate_Params.mat');
PhasePlate_Params = Params.PhasePlate_Params;

% If not specified, compute for the first configuration
if ~exist('ids', 'var'); ids = 1; end

simulationType = 'peakPSF';    % Choose from [densePSF | peakPSF]
save_rootdir = 'output';        % Output directory
mkdir(sprintf('%s/%s/', save_rootdir, simulationType));

% Note: spectral response with the same energy for RGB
switch simulationType
    case 'densePSF'
        SR = load('TOLED_phase_mask_package/spectral_response_new.mat');
        SR = SR.SR;
    case 'peakPSF'
        SR.spectral = [0.61e-6, 0.53e-6, 0.47e-6];
        SR.weights = [1, 0, 0;
            0, 1, 0;
            0, 0, 1];
    otherwise
        error('Choose simulation type between peakPSF and densePSF.');
end



for id =  ids

    f = 4.67e-3;                        % Focal length of the camera main lens

    opening = 1 / PhasePlate_Params(id).lensOpenRatio;
    magnification = PhasePlate_Params(id).magnification;
    f1 = PhasePlate_Params(id).f1;      % Focal length of the first phase mask
    f2 = PhasePlate_Params(id).f2;      % Focal length of the second phase mask
    pitch = PhasePlate_Params(id).lensPitch; % Display pixel pitch
    sensorPitch = 2e-6;                 % Sensor pixel size

    DOE_m = PhasePlate_Params(id).DOE_m;
    DOE_lambda0 = PhasePlate_Params(id).DOE_lambda0;

    save_dir = sprintf('%s/%s/TOLED_%.4fm_%.6fm_%.6fm_Opening_%d_Mag_%d_pitch_%d_phasePlate_M_%d_%s/', ...
        save_rootdir, simulationType, f, f1, f1, opening, magnification, ...
        pitch * 1e6, DOE_m, DOE_lambda0);
    mkdir(save_dir);

    % Image system for an under-display camera
    params = init_system(f, f1, f2, 1 / opening, pitch, ...
        sensorPitch, magnification, DOE_m, DOE_lambda0);
    sensor = get_sensor(params);

    PSFs_y = zeros(1024, sensor.M, 3);           % 1024 PSFs along the short edge of the sensor
    PSFs_y_TOLED = zeros(1024, sensor.M, 3);     % 1024 PSFs along the short edge of the sensor

    % Load lambda0 for each lenslet that controls its folding height
    switch DOE_lambda0
        % If 'fix', all lenslets share the same lambda0.
        case 'fix'
            params.DOE_lambda0 = [0.53e-6];
            % If 'optimized', load optimized lambda0s for all lenslets.
        case 'optimized'
            load(sprintf('optimize_phase_masks/lambda0s_pitch%d_m_%d.mat', pitch*1e6, params.DOE_m));
            params.DOE_lambda0 = lambda0s;
            % If 'uniform', uniformly sample lambda0s from 300nm to 700 nm.
        case 'uniform'
            repNum = ceil(params.displayApertureSize / pitch);
            params.DOE_lambda0 = rand(repNum, 1) * 0.3e-6 + 0.4e-6;
            params.DOE_lambda0 = sort(params.DOE_lambda0);
        otherwise
            fprintf('Unknown DOE_lambda type.')
    end

    % Traverse densely sampled wavelengths.
    for wvl = SR.spectral

        % Propagate wavefront throught UDC with TOLED display.
        I3_TOLED = propagate_through_UDC_TOLED(params, 0, wvl);
        I3_TOLED_energy = sum(I3_TOLED);

        % Compute spatially-varying PSFs.
        % pixelIds = -512: 1: 511; %todo!!
        pixelIds = 0; %todo!
        thetas_y = sensor.pixel2angle(pixelIds);
        numPixel = length(pixelIds);
        batchSize = 30;   % batch propagation for mulitple directions

        I3_ours = zeros(numPixel, sensor.M);

        format long;
        fprintf('Display pitch = %.0f um, Phase mask focal length = %d um, wavelength = %.0f nm...\n', ...
            pitch * 1e6, f1 * 1e6, wvl * 1e9);

        pre = 1;

        while pre <= numPixel

            post = min(pre + batchSize - 1, numPixel);
            I3_ours(pre: post, :) = propagate_through_UDC_ours(params, ...
                thetas_y(pre: post), wvl);
            pre = post + 1;

            % clear GPU RAM to avoid out-of-memory
            if exist('gpuId', 'var'); reset(gpu);end;
        end

        % Relative energy compared to that of UDC TOLED.
        I3_ours = I3_ours / I3_TOLED_energy;
        I3_TOLED = I3_TOLED / I3_TOLED_energy;

        % Save weighted averaged PSFs
        % Weights of the current wavelength to RGB measurements.
        weights = SR.weights(SR.spectral == wvl, :);
        PSFs_y(1: numPixel, :, :) = ...
            PSFs_y(1: numPixel, :, :) + ...
            cat(3, I3_ours * weights(1), ...
            I3_ours * weights(2), ...
            I3_ours * weights(3));
        PSFs_y_TOLED(1: numPixel, :, :) = ...
            PSFs_y_TOLED(1: numPixel, :, :) + ...
            cat(3, I3_TOLED * weights(1), ...
            I3_TOLED * weights(2), ...
            I3_TOLED * weights(3));

    end

    % Reduce light transmission by half due to polarization dependent
    % implementation of phase masks.
    PSFs_y = PSFs_y / 2;

    save([save_dir, 'PSFs.mat'], 'PSFs_y', 'PSFs_y_TOLED');
    fprintf('\n\n\n');

end
