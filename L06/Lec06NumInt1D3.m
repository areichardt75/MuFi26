%% Numerical Integration in one-dimension
% Example 3: Highly Oscillatory Function
fprintf('\n=== EXAMPLE 3: Oscillatory Function ===\n');
fprintf('Integrating f(x) = sin(20*x) / (1 + x^2) from 0 to 5\n\n');

% Define the function
f3 = @(x) sin(20*x) ./ (1 + x.^2);

% Integration limits
x_lower = 0;
x_higher = 5;

% Adaptive quadrature is best for oscillatory functions
result3 = integral(f3, x_lower, x_higher);
fprintf('Adaptive Quadrature: %.8f\n\n', result3);

%% Visualisation
% Plot 3: Oscillatory function
% subplot(2,2,3);
figure;
x_plot3 = linspace(x_lower, x_higher, 2000);
y_plot3 = f3(x_plot3);
plot(x_plot3, y_plot3, 'b-', 'LineWidth', 1);
hold on;
area(x_plot3, y_plot3, 'FaceColor', 'g', 'FaceAlpha', 0.3);
grid on;
xlabel('x');
ylabel('f(x)');
title('f(x) = sin(20x) / (1 + x^2)');
hold off;
