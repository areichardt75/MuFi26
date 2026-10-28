% Polynom fitting on points 
% Use n+1 points for fitting an nth order polynom 
% All points are used and curve is on points
clearvars;
x = [-1; 1; 2; 3];
y = [4; 2; 1; 16];
A = [1 x(1) x(1)^2 x(1)^3;...
  1 x(2) x(2)^2 x(2)^3;...
  1 x(3) x(3)^2 x(3)^3;...
  1 x(4) x(4)^2 x(4)^3];
f = A\y
%% Graph solution
% alternative solution : use vander() for vandermont matrix creation


xfit = linspace(min(x), max(x), 1e3);
yfit = f(1)+f(2)*xfit+f(3)*xfit.^2+f(4)*xfit.^3;

figure; 
  plot(x,y, 'ro','MarkerFaceColor','r','MarkerSize',8);
  hold on;
  plot(xfit, yfit, 'k--');
  xlabel('x');
  ylabel('y');

%% Alternative solution : usage of vander() for Vandermont matrix
% 
V = vander(x);
Av = fliplr(V);
f2 = Av \ y;

%% Create plot different way
% Define interpolation function 
fx = @(x,fp) fp(1)+fp(2)*x+fp(3)*x.^2+fp(4)*x.^3;

yf = fx(x,f2);
figure; 
  plot(x,y, 'ro','MarkerFaceColor','r','MarkerSize',8);
  hold on;
  plot(xfit, yfit, 'k--');
  xlabel('x');
  ylabel('y');


