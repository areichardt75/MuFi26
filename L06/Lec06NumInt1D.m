%% Numerical Integration in one-dimension
% Simple Polynomial Function
clearvars; 
clc;
fprintf('=== EXAMPLE 1: Polynomial Function ===\n');

% Define the function
% f1 = @(x) x.^2 + 2.*x + 1;
% f1 = @(x) x.^4+3*x.^3+x-1;
f1 = @(x) 2*exp(-x/2).*sin(2*pi/3*x);

% Integration limits
x_lower = 0;
x_higher = 3;

% Analytical solution for comparison
analytical1 = (3^3/3 + 3^2 + 3) - (0); % = 9 + 9 + 3 = 21
fprintf('Analytical Solution: %.6f\n\n', analytical1);

% Method 1: Trapezoidal Rule (trapz)
% Number of points - N1
N1 = 1000;
% generating points for x dimension
x1_trapezoid = linspace(x_lower, x_higher, N1);
% calculate function value at the points generated
y1_trapezoid = f1(x1_trapezoid);

result_trapezoid = trapz(x1_trapezoid, y1_trapezoid);
error_trapezoid = abs(result_trapezoid - analytical1);
fprintf('Trapezoidal Rule (n=%d): %.6f (error: %.2e)\n',...
  N1, result_trapezoid, error_trapezoid);

%% Method 2: Simpson's Rule (own implementation)
% 
result_simpson = simpson_rule(f1, x_lower, x_higher, N1);
error_simpson = abs(result_simpson - analytical1);
fprintf('Simpson''s Rule (n=%d): %.6f (error: %.2e)\n', N1, result_simpson, error_simpson);


%% Method 3: Adaptive Quadrature (integral)
result_integral = integral(f1, x_lower, x_higher);
error_integral = abs(result_integral - analytical1);
fprintf('Adaptive Quadrature: %.6f (error: %.2e)\n\n', result_integral, error_integral);

%% Visualization
% subplot(2,2,1);r
figure;
x_plot1 = linspace(x_lower, x_higher, 1000);
y_plot1 = f1(x_plot1);
plot(x_plot1, y_plot1, 'b-', 'LineWidth', 2);
hold on;
x_trap_viz = linspace(x_lower, x_higher, 5);
y_trap_viz = f1(x_trap_viz);
for i = 1:length(x_trap_viz)-1
    patch([x_trap_viz(i), x_trap_viz(i+1), x_trap_viz(i+1), x_trap_viz(i)], ...
          [0, 0, y_trap_viz(i+1), y_trap_viz(i)], 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'r');
end
grid on;
xlabel('x');
ylabel('f(x)');
title('Trapezoidal Rule: f(x) = x^2 + 2x + 1');
legend('Function', 'Trapezoidal Approximation');
hold off;


%% Custom functions of solution
function result = simpson_rule(f, a, b, n)
    % Ensure n is even
    if mod(n, 2) == 1
        n = n + 1;
    end
    
    h = (b - a) / n;
    x = linspace(a, b, n+1);
    y = f(x);
    
    result = y(1) + y(end);
    result = result + 4 * sum(y(2:2:end-1));
    result = result + 2 * sum(y(3:2:end-2));
    result = result * h / 3;
end