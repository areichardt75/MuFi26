%% Lecture 5. Numerical Integration
clearvars;
clc;

%% Example 2: Gaussian Function (2D Normal Distribution)
fprintf('=== EXAMPLE 2: 2D Gaussian Function ===\n');
fprintf('Integrating f(x,y) = exp(-(x^2 + y^2)) over [-2,2] x [-2,2]\n\n');

% Define the function
f2 = @(x, y) exp(-(x.^2 + y.^2));

% Integration limits
x_low = -2; x_high = 2;
y_low = -2; y_high = 2;

% Using integral2
result2 = integral2(f2, x_low, x_high, y_low, y_high);
fprintf('integral2 Result: %.8f\n', result2);
fprintf('(Theoretical over infinite domain: π = %.8f)\n\n', pi);

%% Visualizations
% Plot 2: Gaussian function
figure;
[X_plot2, Y_plot2] = meshgrid(linspace(x_low, x_high, 50),...
  linspace(y_low, y_high, 50));
Z_plot2 = f2(X_plot2, Y_plot2);
surf(X_plot2, Y_plot2, Z_plot2);
shading interp;
xlabel('x'); ylabel('y'); zlabel('f(x,y)');
title('f(x,y) = exp(-(x^2 + y^2))');
colorbar;
