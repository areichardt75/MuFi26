%% Lecture 5 - Numerical Integration 2D
% Two-Dimensional Numerical Integration Examples
% Methods: trapz (2D), integral2, and custom implementations
clearvars;
clc;

%% Example 1: Simple 2D Polynomial
fprintf('=== EXAMPLE 1: 2D Polynomial Function ===\n');
fprintf('Integrating f(x,y) = x^2 + y^2 over [0,1] x [0,1]\n\n');

% Define the function
f1 = @(x, y) x.^2 + y.^2;

% Integration limits
x_low = 0; x_high = 1;
y_low = 0; y_high = 1;

% Analytical solution: ∫∫(x^2 + y^2)dxdy = [x^3/3 + xy^2]|_0^1 * dy = 
% = ∫(1/3 + y^2)dy = 1/3 + 1/3 = 2/3
analytical1 = 2/3;
fprintf('Analytical Solution: %.8f\n\n', analytical1);

% Method 1: Using integral2 (adaptive quadrature)
result_integral2 = integral2(f1, x_low, x_high, y_low, y_high);
error_integral2 = abs(result_integral2 - analytical1);
fprintf('integral2: %.8f (error: %.2e)\n', result_integral2, error_integral2);

% Method 2: Using trapz (2D trapezoidal rule)
N = 1000;
x1 = linspace(x_low, x_high, N);
y1 = linspace(y_low, y_high, N);
[X1, Y1] = meshgrid(x1, y1);
Z1 = f1(X1, Y1);
result_trapz2d = trapz(y1, trapz(x1, Z1, 2));
error_trapz2d = abs(result_trapz2d - analytical1);
fprintf('trapz (2D, n=%d): %.8f (error: %.2e)\n\n', N, result_trapz2d, error_trapz2d);

%% Visualization
% Plot 1: Simple 2D polynomial
% subplot(2,3,1);
figure;
[X_plot1, Y_plot1] = meshgrid(linspace(x_low, x_high, 50),...
  linspace(y_low, y_high, 50));
Z_plot1 = f1(X_plot1, Y_plot1);
surf(X_plot1, Y_plot1, Z_plot1);
colormap(jet);
shading interp;
xlabel('x'); ylabel('y'); zlabel('f(x,y)');
title('f(x,y) = x^2 + y^2');
colorbar;

%%
% Plot 4: Contour plot of polynomial
figure;
contourf(X_plot1, Y_plot1, Z_plot1, 20);
xlabel('x'); ylabel('y');
title('Contour: f(x,y) = x^2 + y^2');
colorbar;
