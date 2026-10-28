% Least squares polynomial fitting 
% More data, no linear fitting
% p(x) = c0+c1*x+x2*x^2+ ... + cn*x^n
% data points : (xi,yi), i=1...m
% 
% A. Find quadratic least squares fit to the data! 
%   number of 
% 
clearvars
xi = [0;1;2;3]
yi = [3;2;4;4]
A = [ones(size(xi)) xi xi.^2]
cp = (transpose(A)*A)\(transpose(A)*yi)

%% Assembling equations using built-in functions  
% Using vander() --> full matrix 
% last column is not used
Astar = fliplr(vander(xi));
Av = Astar(:,1:end-1);
cps = (transpose(Av)*Av)\(transpose(Av)*yi)

% Calculate sum of residual squares
rxv = yi-(cp(1)+cp(2)*xi+cp(3)*xi.^2);
rx = sum(norm(rxv))

%% Plot measuring data and fitted function
% 
x = linspace(min(xi), max(xi), 100);
y = cp(1)+cp(2)*x+cp(3)*x.^2;

figure;
  plot(xi,yi, 'ro','MarkerSize',8,'MarkerFaceColor','r');
  hold on;
  plot(x,y,'k--','LineWidth',2);
  xlabel('x'); 
  ylabel('y');
  title('Quadratic interpolation of measured data');
  legend('Measured data','Fitted curve','Location','southeast');

