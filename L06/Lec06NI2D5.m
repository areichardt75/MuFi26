%% Lecture 5 - Numerical Integration 2D
% Two-Dimensional Numerical Integration Examples
% Methods: trapz (2D), integral2, and custom implementations
clearvars;
clc;

% 
%% Example 5: Monte Carlo Integration 
fprintf('=== EXAMPLE 5: Monte Carlo Integration ===\n');
fprintf('Using Monte Carlo to integrate f(x,y) = sin(x)*cos(y)\n');
fprintf('over [0,π] x [0,π]\n\n');

% Define the function
f5 = @(x, y) sin(x) .* cos(y);

% Integration limits
x_LO = 0; x_HI = pi;
y_LO = 0; y_HI = pi;

% Analytical: ∫∫sin(x)cos(y)dxdy = [-cos(x)]|_0^π * [sin(y)]|_0^π = 2 * 0 = 0
analytical5 = 0;
fprintf('Analytical Solution: %.8f\n', analytical5);

% Monte Carlo integration
n_samples = [1000, 10000, 100000, 1000000];
fprintf('\nMonte Carlo Results:\n');
for n = n_samples
    rng(now); % For reproducibility
    x_rand = x_LO + (x_HI - x_LO) * rand(n, 1);
    y_rand = y_LO + (y_HI - y_LO) * rand(n, 1);
    f_vals = f5(x_rand, y_rand);
    area = (x_HI - x_LO) * (y_HI - y_LO);
    result_mc = area * mean(f_vals);
    error_mc = abs(result_mc - analytical5);
    fprintf('n = %7d: %.8f (error: %.2e)\n', n, result_mc, error_mc);
end

% Compare with integral2
result5_int2 = integral2(f5, x_LO, x_HI, y_LO, y_HI);
fprintf('\nintegral2: %.8f (error: %.2e)\n\n', result5_int2, abs(result5_int2 - analytical5));

%% Visualizations

% Plot 6: Monte Carlo convergence
figure;
n_mc = logspace(2, 6, 50);
errors_mc = zeros(size(n_mc));
for i = 1:length(n_mc)
    rng(42);
    n_i = round(n_mc(i));
    x_r = x_LO + (x_HI - x_LO) * rand(n_i, 1);
    y_r = y_LO + (y_HI - y_LO) * rand(n_i, 1);
    area = (x_HI - x_LO) * (y_HI - y_LO);
    result_mc_i = area * mean(f5(x_r, y_r));
    errors_mc(i) = abs(result_mc_i - analytical5);
end
loglog(n_mc, errors_mc, 'b-', 'LineWidth', 2);
hold on;
loglog(n_mc, 1./sqrt(n_mc), 'r--', 'LineWidth', 2);
grid on;
xlabel('Number of Samples');
ylabel('Absolute Error');
title('Monte Carlo Convergence');
legend('Monte Carlo Error', '1/√n (theoretical)', 'Location', 'best');
hold off;
