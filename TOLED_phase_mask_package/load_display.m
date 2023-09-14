function [display] = load_display(M, delta1, pitch, openRatio)

% create unit OLED pattern
% patternSize = 168e-6;
patternN = round(pitch / delta1);

% slitSize = 34e-6;
slitSize = pitch * openRatio;
slitN = round(slitSize / delta1);

pattern = zeros(1, patternN);
N0 = round((patternN - slitN) / 2);
pattern(N0+1:N0+slitN) = 1;

% repeat OLED pattern to form display opennings
% x1 = (-M/2: 1: M/2-1) * delta1;
% displayAperture = rect(x1/D1);
repNum = ceil(M * delta1 / pitch);
u1_1 = repmat(pattern, [1, repNum]);
display = u1_1(:, 1:M);
% u1_1 = u1_1 .* displayAperture;