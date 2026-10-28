%% Lecture 5 - Numerical Integration
% Two-Dimensional Numerical Integration Examples
% Methods: trapz (2D), integral2, and custom implementations
clearvars;
clc;

%
%% Example 3: Non-Rectangular Domain (Circle)

fprintf('=== EXAMPLE 3: Integration over Circular Domain ===\n');
fprintf('Integrating f(x,y) = x*y over unit circle (x^2 + y^2 <= 1)\n\n');

% Define the function
f3 = @(x, y) x .* y;

% For a circular domain, we can use variable limits with integral2
% y ranges from -sqrt(1-x^2) to sqrt(1-x^2)
y_lower = @(x) -sqrt(1 - x.^2);
y_upper = @(x) sqrt(1 - x.^2);

% Integrate over circular domain
result3 = integral2(f3, -1, 1, y_lower, y_upper);
fprintf('integral2 (circular domain): %.8f\n', result3);
fprintf('(By symmetry, should be ≈ 0)\n\n');

%% Visualizations

% Plot 3: Circular domain
figure;
theta = linspace(0, 2*pi, 100);
r = linspace(0, 1, 30);
[THETA, R] = meshgrid(theta, r);
X_circ = R .* cos(THETA);
Y_circ = R .* sin(THETA);
Z_circ = f3(X_circ, Y_circ);
surf(X_circ, Y_circ, Z_circ);
shading interp;
xlabel('x'); ylabel('y'); zlabel('f(x,y)');
title('f(x,y) = xy (circular domain)');
colorbar;

