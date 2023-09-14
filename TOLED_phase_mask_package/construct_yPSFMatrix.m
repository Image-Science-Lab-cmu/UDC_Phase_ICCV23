function [yPSFMatrix] = construct_yPSFMatrix(PSFs_y)

% convert spatially varying PSF to a PSF matrix
% Input:
%   PSFs_y: number of PSF x PSF dim x 3 
%           each row is one PSF
%           each PSF is 1 x 1024
% Output:
%   yPSFMatrix: each column is one PSF, dimension 1024
%               number of column is number of spatially-varying PSF, ie. 1024

[numPSF, dimPSF, numChannels] = size(PSFs_y);
yPSFMatrix = zeros(dimPSF, numPSF, numChannels);

for cc = 1: numChannels
    yPSFMatrix(:, :, cc)=PSFs_y(:, :, cc)';
end
yPSFMatrix = yPSFMatrix(end: -1: 1, :, :);

end