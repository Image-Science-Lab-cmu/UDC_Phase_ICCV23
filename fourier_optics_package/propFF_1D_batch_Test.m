% test input
addpath utils;
thetas = [0,0,0,0,0];
lambda = 0.5e-6;
batchSize = length(thetas); 
SINGLE = 0;
% D0 = 4.2e-3;                    % diam of the disp aperture [m]
D1 = 3.475e-3;                    % diam of the lens aperture [m]
% D1 = 3.277e-3;
k = 2*pi/lambda;                % wavenumber
f = 5.6e-3;                     % main lens focal length (m)

% delta1 = 0.05e-6;               % diam of the lens aperture [m]
% delta2 = 0.05e-6;               % diam of the lens aperture [m]
% M = 2^18;                       % number of samples

delta1 = 1.5e-8;               % diam of the lens aperture [m]
delta2 = 1.5e-8;               % diam of the lens aperture [m]
delta3 = 1.5e-7;               % todo: spacing on sensor plane
% M = 2^20;                       % number of samples
M = floor(lambda*f/delta2/delta3/2)*2;

D2 = M*delta2;                  % diam of sensor region of interest[m]
L1 = M*delta1;                  % diam of input plane [m]

f0 = -0.0002;                   % negative lens focal length (m)
f1 = 0.0004;                    % fresnel propagatoin distance (m) % todo:2x
opening = 5;
z=f1+f0;                        % positive lens focal length (m)
D0 = D1 + 2*z*tan(50/180*pi);   % diam of the disp aperture [m]

patternSize = 168e-6;
patternN = round(patternSize / delta1);

% slitSize = 34e-6;
slitSize = patternSize / opening;
slitN = round(slitSize / delta1);

pattern = zeros(1, patternN);
N0 = round((patternN - slitN) / 2);
pattern(N0+1:N0+slitN) = 1;

% repeat OLED pattern to form display opennings
x1 = (-M/2: 1: M/2-1) * delta1;
displayAperture = rect(x1/D0);
repNum = ceil(M * delta1 / patternSize);
u1_1 = repmat(pattern, [1, repNum]);
u1_1 = u1_1(:, 1:M);
u1_1 = u1_1 .* displayAperture;
u1_1 = repmat(u1_1, [batchSize, 1]);

Uin_box = displayAperture;  % test 1

% Tilted wavefront
thetas = 90-thetas;         % propagating direction (degree)
alphas = cos(thetas / 180 * pi);
tilted_phases = exp(1i*k*alphas'*x1);
Uin_squareWave_5_same = u1_1 .* tilted_phases; % test 3
Uin_squareWave_1 = Uin_squareWave_5_same(1,:); % test 2

thetas = [0, 5, 10, 15, 20];
thetas = 90-thetas;         % propagating direction (degree)
alphas = cos(thetas / 180 * pi);
tilted_phases = exp(1i*k*alphas'*x1);
Uin_squareWave_5_diff = u1_1 .* tilted_phases; % test 4

L1 = M * delta1;

%% Test 1: Box
[exp_u2, exp_L2, exp_x2] = propFF_1D(Uin_box, L1, lambda, f, 0);
[u2, L2, x2] = propFF_1D_batch(Uin_box, L1, lambda, f, 0);
assert(sum(abs(exp_u2 - u2)) == 0);
assert(sum(abs(exp_L2 - L2)) == 0);
assert(sum(abs(exp_x2 - x2)) == 0);

%% Test 2: square wave
[exp_u2, ~, ~] = propFF_1D(Uin_squareWave_1, L1, lambda, f, 0);
[u2, ~, ~] = propFF_1D_batch(Uin_squareWave_1, L1, lambda, f, 0);
assert(sum(abs(exp_u2 - u2)) == 0);

%% Test 3: square wave mini-batch same Uin
[exp_u2, ~, ~] = propFF_1D(Uin_squareWave_1, L1, lambda, f, 0);
[u2, ~, ~] = propFF_1D_batch(Uin_squareWave_5_same, L1, lambda, f, 0);
tol = 1e-10;
assert(max(Uin_squareWave_1 - Uin_squareWave_5_same(1,:)) < tol, 'Input not equal');
assert(max(abs(exp_u2 - u2(1, :))) < tol, 'The first output error');
assert(max(abs(exp_u2 - u2(2, :))) < tol, 'The second output error');
assert(max(abs(exp_u2 - u2(3, :))) < tol, 'The third output error');
assert(max(abs(exp_u2 - u2(4, :))) < tol, 'The fourth output error');
assert(max(abs(exp_u2 - u2(5, :))) < tol, 'The fifth output error');

%% Test 4: square wave mini-batch with different Uin
[u2, ~, ~] = propFF_1D_batch(Uin_squareWave_5_diff, L1, lambda, f, 0);

tol=1e-10;
[exp_u2, ~, ~] = propFF_1D(Uin_squareWave_5_diff(1,:), L1, lambda, f, 0);
assert(max(abs(exp_u2 - u2(1, :))) < tol, 'The first output error');

[exp_u2, ~, ~] = propFF_1D(Uin_squareWave_5_diff(2,:), L1, lambda, f, 0);
assert(max(abs(exp_u2 - u2(2, :))) < tol, 'The first output error');

[exp_u2, ~, ~] = propFF_1D(Uin_squareWave_5_diff(3,:), L1, lambda, f, 0);
assert(max(abs(exp_u2 - u2(3, :))) < tol, 'The first output error');

[exp_u2, ~, ~] = propFF_1D(Uin_squareWave_5_diff(4,:), L1, lambda, f, 0);
assert(max(abs(exp_u2 - u2(4, :))) < tol, 'The first output error');

[exp_u2, ~, ~] = propFF_1D(Uin_squareWave_5_diff(5,:), L1, lambda, f, 0);
assert(max(abs(exp_u2 - u2(5, :))) < tol, 'The first output error');

% conclusion: propFF_1D_batch functions the same as propFF_1D
