%% Numerical Integration in one-dimension
% Example 2: Trigonometric Function
fprintf('=== EXAMPLE 2: Trigonometric Function ===\n');
fprintf('Integrating f(x) = sin(x) * exp(-x/5) from 0 to 2π\n\n');

% Define the function
f2 = @(x) sin(x) .* exp(-x/5);

% Integration limits
x_lower = 0;
x_higher = 2*pi;

% Using integral (most accurate)
result = integral(f2, x_lower, x_higher);
fprintf('Integral Result: %.8f\n\n', result);

%% Compare different numbers of points with trapz
n_values = [10, 50, 100, 500, 1000];
fprintf('Convergence study with Trapezoidal Rule:\n');
for n = n_values
    x2 = linspace(x_lower, x_higher, n);
    y2 = f2(x2);
    result_n = trapz(x2, y2);
    error_n = abs(result_n - result);
    fprintf('n = %4d: %.8f (error: %.2e)\n', n, result_n, error_n);
end

%% Plotting 
figure; 
% subplot(2,2,2);
x = linspace(x_lower, x_higher, 1000);
y = f2(x);
plot(x, y, 'b-', 'LineWidth', 2);
hold on;
% Color the area under graph (integrational value)
area(x, y, 'FaceColor', 'c', 'FaceAlpha', 0.3);
grid on;
xlabel('x');
ylabel('f(x)');
title('f(x) = sin(x) * exp(-x/5)');
hold off;

