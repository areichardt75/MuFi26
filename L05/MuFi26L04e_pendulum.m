%% Mufi 2025 - Lecture 8. - Solution of ODEs
% Simple pendulum 

% general form : Th'' = -g/ell*sin(Th)
% linearized form : Th'' = -g/ell*Th

ell = 1;
g = 9.81;

% odefun : x = [Th Th']
odefun = @(x,ell,g) [x(2);-g/ell*sin(x(1))];
odelin = @(x,ell,g) [x(2);-g/ell*x(1)];

% preparing simulation
T0 = 2*pi/sqrt(g/ell);
Tspan = 0:0.01:2.5*T0;
x0 = [pi/3;0];
opts = odeset('RelTol',1e-6);
% solving ODE for general case 
[T,Y] = ode45(@(t,x) odefun(x,ell,g), Tspan, x0,opts);
% and for linearized cases
[Tlin,Ylin] = ode45(@(t,x) odelin(x,ell,g), Tspan, x0,opts);

% taking variables out of solution
Th = Y(:,1); Thv = Y(:,2);
Thlin=Ylin(:,1); Thvlin=Y(:,2);

%% Visualization
figure; 
  plot(T,Th,'r-',Tlin,Thlin,'k-','LineWidth',2);
  xlabel('time'); 
  ylabel('Theta');
  legend('general','linearized','Location','northeast');

%% Phase space 
figure;
  plot(Th, Thv, 'r-o','LineWidth',2);
  xlabel('Theta'); ylabel('Theta velocity');
  title('Phase space of simple pendulum')