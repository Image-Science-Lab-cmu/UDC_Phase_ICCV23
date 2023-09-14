
clear;

% Number of lenslets / displayPixel within the aperture
numDispPixels = [58, 45, 31, 17];
% Pitch of each display pixel
pitches = [42e-6, 56e-6, 84e-6, 168e-6];

for pitch = pitches
    for DOE_m = [1, 5]

        numDispPixel = numDispPixels(pitches == pitch);

        load('spectral_response_new.mat');
        load(sprintf('invertible_scores_pitch%d_m_%d.mat', ...
            pitch * 1e6, DOE_m));

        SR_weights = SR.weights(1:10:end, :);
        SR_weights = SR_weights ./ sum(SR_weights,1);

        A = SR_weights' * invertible_scores';
        A = A ./ max(A(:));
        b = [1;1;1];

        % log file
        fid = fopen(sprintf('optimize_varying_lambda0s_pitch%d_m_%d.txt', ...
            pitch * 1e6, DOE_m), 'a');

        %% Use barrier methods
        % The goal is to optimize for vector x. Each element x_i represents
        % the propabability of folding lenslets using lambda_i.
        %
        % Optimization target:
        % min_x || Ax - b ||^2
        % s.t.  x >= 0
        %
        % Use log barrier approximation, the constrained problem can be written as
        % the following unconstrained problem,
        % min_x t || Ax - b ||^2 - log(x)
        %
        % Gradient g(x) = A'*(Ax - b) - 1./x
        % Hessian  h(x) = A'*A + (1./x^2 - 0);
        %
        % Algorithm:
        % Initialize with t_0, mu; compute x_0 = argmin l(t_0, x)
        % 1. Solve the barrier problem at t=t_k, using Newton initialized at x_(k-1)
        % 2. Stop if duality gap m/t < epsilon, else update t_(k+1) = mu*t.

        for q = [0.1] % Tuned hyperparameter for TOLED display


            mu = 2;
            epsilon = 1e4;

            h_ii = zeros(30, 30);
            for ii=1:30
                h_ii(ii, ii) = 1;
            end

            % gradient and hessian function
            f_func = @(x, t) t*(A*x-b)'*(A*x-b) + t*q * sum(x) - sum(my_log(h_ii'*x));
            g_func = @(x, t) t*(A'*(A*x-b)) + t*q - 1./x;
            h_func = @(x, t) t*A'*A + 0 + diag(1./x.^2);

            % initialization
            t(1) = 1e-2;
            x_init = ones(30, 1);
            x(:,1) = newton(x_init, t(1), f_func, g_func, h_func);

            for kk = 2: 1000

                t(kk) =  mu * t(kk-1);
                if t(kk) > epsilon; break; end;

                x(:, kk) = newton(x(:, kk-1), t(kk), f_func, g_func, h_func);

                loss = (A*x(:,kk)-b)'*(A*x(:,kk)-b);
                norm = sum(x(kk));
                fprintf('iter=%d loss=%.4f, norm=%.4f\n', kk, loss, norm);

            end

            x_optim = x(:, end);
            x_optim = x_optim / sum(x_optim);

            %% save
            % convert propabablity x_optim to display pixel counts.
            x_optim_quantized = round(x_optim * numDispPixel);
            % make sure x_optim_quantized sum to numDispPixel
            quantization_error = x_optim * numDispPixel - x_optim_quantized;
            [quantization_error, error_order] = sort(quantization_error, 'descend');
            for ii = 1: length(x_optim_quantized)
                if sum(x_optim_quantized) >= numDispPixel
                    break;
                else
                    x_optim_quantized(error_order(ii)) = ...
                        x_optim_quantized(error_order(ii)) + 1;
                end
            end
            for ii = length(x_optim_quantized): -1 : 1
                if sum(x_optim_quantized) <= numDispPixel
                    break;
                else
                    x_optim_quantized(error_order(ii)) = ...
                        x_optim_quantized(error_order(ii)) - 1;
                end
            end

            figure(1), stem(x_optim * numDispPixel, 'LineWidth', 2);
            figure(2), set(gcf,'Color',[1 1 1], 'InvertHardCopy','off');
            stem(x_optim_quantized, 'LineWidth', 2);

            % convert counts for each lambda0 into array of lambda0s
            wvls = SR.spectral(1: 10: end);
            lambda0s = [];
            for wvl_id = 1: length(wvls)
                count = x_optim_quantized(wvl_id);
                wvl = wvls(wvl_id);
                while count > 0
                    lambda0s = [lambda0s, wvl];
                    count = count - 1;
                end
            end

            %% Comparing optimized lambdas with 'fix' and 'uniform'
            % 'Uniform': Uniformly sample lambdas0 from 400nm to 700nm
            x_uniform = zeros(30,1);
            wvl_uniform = linspace(min(wvls), max(wvls), numDispPixel);
            for curr_wvl = wvl_uniform
                [~, ii] = min(abs(curr_wvl - wvls));
                x_uniform(ii) = x_uniform(ii) + 1;
            end
            % 'Fix': Sample all lambdas0 at 530nm
            x_fix = zeros(30,1); x_fix(14)=numDispPixel;

            RGB_optim = (A * x_optim_quantized)';
            RGB_uniform = (A * x_uniform)';
            RGB_fix = (A * x_fix)';

            fprintf('q=%.2f, loss=%.4f, l1_norm=%.4f, optimized RGB(%.2f, %.2f, %2f), uniform RGB(%.2f, %.2f, %2f), fix RGB(%.2f, %.2f, %2f)\n', ...
                q, loss, norm, ...
                RGB_optim(1), RGB_optim(2), RGB_optim(3), ...
                RGB_uniform(1), RGB_uniform(2), RGB_uniform(3), ...
                RGB_fix(1), RGB_fix(2), RGB_fix(3));
            fprintf(fid, 'q=%.2f, loss=%.4f, l1_norm=%.4f, optimized RGB(%.2f, %.2f, %2f), uniform RGB(%.2f, %.2f, %2f), fix RGB(%.2f, %.2f, %2f)\n', ...
                q, loss, norm, ...
                RGB_optim(1), RGB_optim(2), RGB_optim(3), ...
                RGB_uniform(1), RGB_uniform(2), RGB_uniform(3), ...
                RGB_fix(1), RGB_fix(2), RGB_fix(3));
            fclose(fid);
        end
        save(sprintf('lambda0s_pitch%d_m_%d.mat', pitch * 1e6, DOE_m), ...
            'lambda0s');

    end
end

%% Functions

function x_star = newton(x_init, t, f_func, g_func, h_func)

x_star = x_init;
for kk = 2: 1000
    x_pre = x_star;
    %     disp(f_func(x_pre, t));
    x_star = x_star - inv(h_func(x_star, t)) * g_func(x_star, t);
    if abs(f_func(x_pre, t) - f_func(x_star, t)) < 0.001
        break;
    end
    if f_func(x_star, t) == Inf
        x_star = x_pre;
        break;
    end
end
end

function y = my_log(x)

if x > 0
    y = reallog(x);
else
    y = -Inf;
end
end

