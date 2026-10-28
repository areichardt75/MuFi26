% Advanced Optimization: Portfolio Optimization Problem
% Minimize portfolio risk while achieving target return
% Subject to: sum of weights = 1, all weights >= 0, and minimum return constraint

clear; clc;

% Define 5 assets with historical returns (simulated data)
rng(42); % For reproducibility
rng(now)
n_assets = 5;
n_periods = 100;

% Generate random returns for 5 assets
returns = randn(n_periods, n_assets) .* [0.15, 0.20, 0.18, 0.12, 0.25] + ...
          [0.08, 0.12, 0.10, 0.06, 0.15];

% Calculate expected returns and covariance matrix
expected_returns = mean(returns)';
cov_matrix = cov(returns);

fprintf('Expected Annual Returns:\n');
for i = 1:n_assets
    fprintf('Asset %d: %.2f%%\n', i, expected_returns(i)*100);
end
%%
% Target return (10% annually)
target_return = 0.10;

% Objective: Minimize portfolio variance (risk)
% Portfolio variance = w' * Sigma * w
objective = @(w) w' * cov_matrix * w;

% Constraints:
% 1. Sum of weights = 1 (Aeq * w = beq)
Aeq = ones(1, n_assets);
beq = 1;

% 2. Portfolio return >= target_return (-A * w <= -b form)
A = -expected_returns';
b = -target_return;

% 3. All weights >= 0 (lower bound)
lb = zeros(n_assets, 1);
ub = ones(n_assets, 1);

% Initial guess (equal weights)
w0 = ones(n_assets, 1) / n_assets;

% Optimization options
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');

% Solve the optimization problem
[w_optimal, portfolio_variance] = fmincon(objective, w0, A, b, Aeq, beq, lb, ub, [], options);

% Calculate portfolio metrics
portfolio_return = expected_returns' * w_optimal;
portfolio_risk = sqrt(portfolio_variance);

%% Display results
fprintf('\n=== OPTIMAL PORTFOLIO ALLOCATION ===\n');
for i = 1:n_assets
    fprintf('Asset %d: %.2f%%\n', i, w_optimal(i)*100);
end
fprintf('\nPortfolio Expected Return: %.2f%%\n', portfolio_return*100);
fprintf('Portfolio Risk (Std Dev): %.2f%%\n', portfolio_risk*100);
fprintf('Sharpe Ratio (assuming rf=2%%): %.4f\n', (portfolio_return - 0.02)/portfolio_risk);

%% Visualize the efficient frontier
n_points = 50;
target_returns = linspace(min(expected_returns), max(expected_returns), n_points);
efficient_risks = zeros(n_points, 1);

for i = 1:n_points
    A_temp = -expected_returns';
    b_temp = -target_returns(i);
    
    try
        [~, var_temp] = fmincon(objective, w0, A_temp, b_temp, Aeq, beq, lb, ub, [], ...
                                optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp'));
        efficient_risks(i) = sqrt(var_temp);
    catch
        efficient_risks(i) = NaN;
    end
end

% Plot efficient frontier
figure;
subplot(1,2,1);
plot(efficient_risks*100, target_returns*100, 'b-', 'LineWidth', 2);
hold on;
plot(portfolio_risk*100, portfolio_return*100, 'ro', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
scatter(sqrt(diag(cov_matrix))*100, expected_returns*100, 100, 'filled', 'MarkerFaceColor', [0.5 0.5 0.5]);
grid on;
xlabel('Risk (Standard Deviation %)');
ylabel('Expected Return (%)');
title('Efficient Frontier');
legend('Efficient Frontier', 'Optimal Portfolio', 'Individual Assets', 'Location', 'best');
hold off;

% Plot portfolio weights
subplot(1,2,2);
bar(w_optimal*100);
xlabel('Asset Number');
ylabel('Allocation (%)');
title('Optimal Portfolio Weights');
grid on;

fprintf('\n=== Optimization Complete ===\n');