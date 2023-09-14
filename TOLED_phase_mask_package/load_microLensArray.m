function [microLensArray, hmap] = load_microLensArray(k, M, delta, pitch, ...
    f, DOE_m, DOE_lambda0)

if ~exist('DOE_m', 'var')
    DOE_m = 0;
end

if ~exist('DOE_lambda0', 'var')
    DOE_lambda0 = '';
end

% Load phase plate profile
patternN = round(pitch / delta);
x1_2 = (-patternN/2: 1: patternN/2-1) * delta;
n = 1.515;  % Refractive index of phase masks
R = (n-1)*f;

% If use a fixed lambda0 for all lenslets
if length(DOE_lambda0) == 1
        
    hmax = sign(R) * DOE_lambda0 / (n-1) * DOE_m;
    
    hmap1 = mod(x1_2.^2/2/R, hmax);
    microLens1 = exp(-1i*k*(n-1)*hmap1);
    
    % create negative micro-lens array
    repNum = ceil(M * delta / pitch);
    microLensArray = repmat(microLens1, [1, repNum]);
    microLensArray = microLensArray(:, 1:M);
    hmap = repmat(hmap1, [1, repNum]);
    hmap = hmap(:, 1:M);

else
    % If use a different lambda0 for each lenslet    
    repNum = ceil(M * delta / pitch);
    if length(DOE_lambda0) < repNum; error('Not enougth DOE lambda0s.'); end
    
    hmax = sign(R) * DOE_lambda0 / (n-1) * DOE_m;
    hmap = [];
    for microLensId = 1: repNum
        hmap = [hmap, mod(x1_2.^2/2/R, hmax(microLensId))];
    end
    microLensArray = exp(-1i*k*(n-1)*hmap);
    microLensArray = microLensArray(:, 1:M);
    
end