function main_quantitative_evalution(ids, gpuId)

% Simulate and evaluate images captured under our setup.
%
% Input:
%   ids  : setup id that specifies configurations  
%   gpuId: (optional) GPU id if computed using a GPU.

addpath utils
addpath bm3d_matlab_package
addpath bm3d_matlab_package/bm3d
addpath cgsSolver

profile off;
profile on -history;

if exist('gpuId', 'var'); gpu = gpuDevice(gpuId); end
if ~exist('ids', 'var'); ids = 1; end


% Setup names
simulationNames = {
    'TOLED_0.0047m_0.000042m_0.000042m_Opening_4.200000e+00_Mag_4_pitch_42_phasePlate_M_5_fix', ...
    'TOLED_0.0047m_0.000042m_0.000042m_Opening_4.200000e+00_Mag_4_pitch_42_phasePlate_M_5_optimized', ...
    'TOLED_0.0047m_0.000042m_0.000042m_Opening_4.200000e+00_Mag_4_pitch_42_phasePlate_M_5_uniform', ...
    'TOLED_0.0047m_0.000042m_0.000042m_Opening_4.200000e+00_Mag_4_pitch_42_phasePlate_M_1_fix', ...
    'TOLED_0.0047m_0.000042m_0.000042m_Opening_4.200000e+00_Mag_4_pitch_42_phasePlate_M_1_optimized', ...
    'TOLED_0.0047m_0.000042m_0.000042m_Opening_4.200000e+00_Mag_4_pitch_42_phasePlate_M_1_uniform', ...
    'TOLED_0.0047m_0.000056m_0.000056m_Opening_4.200000e+00_Mag_4_pitch_56_phasePlate_M_5_optimized', ...
    'TOLED_0.0047m_0.000056m_0.000056m_Opening_4.200000e+00_Mag_4_pitch_56_phasePlate_M_1_optimized', ...
    'TOLED_0.0047m_0.000084m_0.000084m_Opening_4.200000e+00_Mag_4_pitch_84_phasePlate_M_5_optimized', ...
    'TOLED_0.0047m_0.000084m_0.000084m_Opening_4.200000e+00_Mag_4_pitch_84_phasePlate_M_1_optimized', ...
    'TOLED_0.0047m_0.000168m_0.000168m_Opening_4.200000e+00_Mag_4_pitch_168_phasePlate_M_5_optimized', ...
    'TOLED_0.0047m_0.000168m_0.000168m_Opening_4.200000e+00_Mag_4_pitch_168_phasePlate_M_1_optimized', ...
    };

operator = 'matmul';
LTR = 0.238*2;          % Light transmission rate of ours setup
refLTR = 0.238*2;       % Reference LTR (for comparison setups)
% Note: scale gain to gaurantee the same intensity for
% different display layouts.
% TOLED gain: 0.476/0.238 = 2
% Ours  gain: 0.476/0.476 = 1
% ZTE   gain: 0.476/0.750 = 0.63
gainRatio = refLTR / LTR;
% quantize captured image into #quantBit-bit linear image
quantBit = 12;

for id = ids
    simulationName = simulationNames{id};
    simulationType = 'densePSF';
    srcImgDir = '../under-display-camera/matlab/fig/HQ/test/';
    srcImgName = dir([srcImgDir, '*.png']);
    % save folder
    mkdir(sprintf('output/%s/%s/', ...
        simulationType, simulationName));
    % Load hyperparameter for iterative solver and denoiser
    load(sprintf('output/%s/%s/tuned_parameters.mat', ...
        simulationType, simulationName));

    %% generate PSF matrix
    load(sprintf('output/%s/%s/PSFs.mat', ...
        simulationType, simulationName));
   
    % Construct spatially-varying PSF operators
    % normalize largest energy to 1
    PSFs_y = PSFs_y / (LTR / 0.238); % PSF energy set to 1
    
    % compute valid field of view
    energy_y = sum(PSFs_y(:,:,2), 2); % each row is one PSF
    fov_y = (energy_y >= 0.5);        
    fov_x = (ones(2048, 1) == 1);
    
    % construct PSF matrix
    omega = construct_yPSFMatrix(PSFs_y);
    omega = omega(fov_y, fov_y, :);
    
    
    %% for SNR; for different images
    SNRs = 24:4:40;
    curr_ssims = zeros(length(SNRs), 1);
    curr_psnrs = zeros(length(SNRs), 1);
    
    for SNR = [40]
        
        mean_psnr = 0;
        mean_ssim = 0;
        imgIds = 1: 1: 30;

        % Load tuned parameters for image reconstruction
        lambda = best_lambdas(SNRs == SNR);
        noise_var = best_noise_vars(SNRs == SNR);
        
        for imgId = imgIds
            
            % Load latent sharp image
            img = im2double(imread([srcImgDir, '/', ...
                srcImgName(imgId).name]));
            img = img ./ (max(img(:)));
            img = img(fov_y, fov_x, :);   % crop valid field of view
            
            % Simulate blurry image and add noise
            % The intensity of imgBlur is equivalent to (img * LTR * gainRatio)
            imgBlurnoisy = capture(img * LTR, omega, operator, SNR, ...
                quantBit, gainRatio);
            
            % Denoise and deblur
            imgSharp = denoise_and_deblur(imgBlurnoisy, omega, ...
                operator, noise_var, lambda);
            imgSharp = imgSharp / LTR / gainRatio; % Intensity compensation

            % Compute PSNR and SSIM of the restored image
            psnrVal = psnr(imgSharp, img);
            [ssimVal, ssimMap] = ssim(imgSharp, img, 'Radius', 1.5);
            
            mean_psnr = mean_psnr + psnrVal;
            mean_ssim = mean_ssim + ssimVal;
            
            fprintf('%s, SNR=%d, lambda=%.5f, imgId=%d, PSNR=%.2f, SSIM=%.2f\n', ...
                simulationName, SNR, lambda, imgId, psnrVal, ssimVal);
            
            if mod(imgId-1, 5) == 0
                prefix = sprintf('output/%s/%s/%s', ...
                    simulationType, simulationName, srcImgName(imgId).name);
                imwrite(imgBlurnoisy, ...
                    sprintf('%s_blurImg_SNR_%d.png', prefix, SNR));
                imwrite(imgBlurnoisy / LTR / gainRatio, ...
                    sprintf('%s_blurImg_SNR_%d_compensated.png', ...
                    prefix, SNR));
                imwrite(imgSharp, ...
                    sprintf('%s_deblurImg_L2Gradient_lambda_%.5f_SNR_%d.png', ...
                    prefix, lambda, SNR));
                
                % visualize the ssimMap
                close all;
                figure('Renderer', 'painters', 'Position', ...
                    [10, 10, 1200, 500]);
                grid on; grid minor;
                set(gcf,'Color',[1 1 1], 'InvertHardCopy','off');
                imagesc(ssimMap(:,:,2)); colorbar; caxis([0, 1])
                saveas(gcf, sprintf('%s_ssimMap_SNR_%d.png', prefix, SNR));
            end

            % clear GPU RAM to avoid OOM
            if exist('gpuId', 'var'); reset(gpu);end;
            
        end
        
        mean_ssim = mean_ssim / length(imgIds);
        mean_psnr = mean_psnr / length(imgIds);
        
        curr_ssims(SNRs == SNR) = mean_ssim;
        curr_psnrs(SNRs == SNR) = mean_psnr;
        
        % 30 test images 
        save(sprintf('output/%s/%s/sweep_snr.mat', ...
            simulationType, simulationName), ...
            'SNR', 'curr_ssims', 'curr_psnrs');
    end
end
