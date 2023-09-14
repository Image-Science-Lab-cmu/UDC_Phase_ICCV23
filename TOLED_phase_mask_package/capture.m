function [imgBlurnoisy] = capture(img, PSF, operator, SNR, quantBit, gainRatio)

% blur the image
% July 12: We skip blur in y-direction for now
% since it multiplies quantum efficiency to each
% channel again.
% quantBit: Quantization of the captured image

switch operator
    case 'conv'
        [kH, kW, ~] = size(PSF);        % PSF size
        [H, W, ~] = size(img);          % image size
        
        % Convolve img and PSF with 'valid' boundary condition.
        img = padarray(img, [ceil(kH/2), ceil(kW/2)], 0, 'both');
        img = img(1:kH+H-1, 1: kW+W-1, :);
        imgBlur = zeros(H, W, 3);
        for cc = 1: 3
            imgBlur(:,:,cc) = conv2(img(:,:,cc), PSF(:,:,cc), 'valid');
        end
        
    case 'matmul'
        [nH, ~, ~] = size(PSF);         % PSF length
        [H, W, ~] = size(img);          % image size
        
        imgBlur = zeros(nH, W, 3);
        for cc = 1: 3
            imgBlur(:,:,cc) = PSF(:,:,cc) * img(:,:,cc);
        end
    otherwise
        msg('Unknown operator');
end


%% add noise
SNRs = 24:4:40;                         % SNR controls the noise level
Ls = [273, 654, 1608, 4005, 10024];     % Light level
sensor.capacity = 15506;
sensor.noise_std = 4.87;
if ~any(SNRs == SNR); error('Only SNR = [24, 28, 32, 36, 40] is valid.'); end
L = Ls(SNRs == SNR);
sensor.gain = 1/L * gainRatio;          % Scale gain by gainRatio

% Add noise
imgBlurnoisy = add_noise(imgBlur * L, 1, sensor);

% quantization noise
imgBlurnoisy = floor(imgBlurnoisy * (2^quantBit - 1)) / (2^quantBit - 1);