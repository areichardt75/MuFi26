%% Lecture 5 - Numerical Integration 2D
% Two-Dimensional Numerical Integration Examples
% Methods: trapz (2D), integral2, and custom implementations
clearvars;
clc;


%% Example 4: Variable Integration Limits
fprintf('=== EXAMPLE 4: Variable Integration Limits ===\n');
fprintf('Integrating f(x,y) = x + y over triangular region:\n');
fprintf('0 <= x <= 1, 0 <= y <= x\n\n');

% Define the function
f4 = @(x, y) x + y;

% Variable y limits
y_low = @(x) zeros(size(x));
y_high = @(x) x;

% Analytical: ∫[0,1]∫[0,x](x+y)dy dx = ∫[0,1][xy + y^2/2]|_0^x dx 
%           = ∫[0,1](x^2 + x^2/2)dx = ∫[0,1](3x^2/2)dx = [x^3/2]|_0^1 = 1/2
analytical4 = 1/2;
fprintf('Analytical Solution: %.8f\n', analytical4);

result4 = integral2(f4, 0, 1, y_low, y_high);
error4 = abs(result4 - analytical4);
fprintf('integral2: %.8f (error: %.2e)\n\n', result4, error4);

%% Visualizations


% Plot 5: Triangular integration region
figure;
x_tri = linspace(0, 1, 50);
y_tri = linspace(0, 1, 50);
[X_tri, Y_tri] = meshgrid(x_tri, y_tri);
Z_tri = f4(X_tri, Y_tri);
Z_tri(Y_tri > X_tri) = NaN; % Mask out region outside triangle
surf(X_tri, Y_tri, Z_tri);
shading interp;
xlabel('x'); ylabel('y'); zlabel('f(x,y)');
title('f(x,y) = x + y (triangular domain)');
colorbar;
