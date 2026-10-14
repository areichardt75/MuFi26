% Define potential function 
% Parameters
a = 2;
V0= 5;
% Potential function
fx = @(x) (abs(x)<=a).*(V0/a^2*(x-a).*(x+a))+V0;

%% Graph potential function
x = -2.5*a:a/100:2.5*a;
pot = fx(x);
figure; 
  plot(x, pot, 'k-','LineWidth',2);
  xlabel('distance');
  ylabel('Potential');
%% 
% p(x) = sum_{i=0}^4 p_i x^i
% p(-3a/2) = V0; p(-a/3) = 0; p(0) = V0/4;
% p(a/3) = 0; p(3a/2) = V0
a = 2; V0 = 5;
x = transpose([-3*a/2 -a/3 0 a/3 3*a/2]);
y = transpose([V0 0 V0/4 0 V0]);
A = [x.^[4 3 2 1 0]]
p = A \ y

%% Nonlinear potential function
potnl = @(x) 0.3774*x.^4-2.9802*x.^2+1.25;
x = -1.5*a:a/100:1.5*a;
Vx = potnl(x);
figure;
  plot(x, Vx, 'r-','LineWidth',2); 
  xlabel('x'); ylabel('V(x)');
  title('V(x) = A*x^4+B*x^2+C');
