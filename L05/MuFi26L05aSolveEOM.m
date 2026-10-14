% Solving equations of motion
% 
% y = [x, v]
V0 = 5;
a = 2;
m = 1;
%% Solve linear oscillator problem
% pot = @(x,a,V0) (abs(x)<=a).*(V0/a^2*(x-a).*(x+a))+V0;
eom = @(x,a,V0,m) [x(2);-V0/(a^2*m)*2*x(1)];

Tspan = [0 10];
x0 = [a/3;0.1];
[T,Y] = ode45(@(t,x) eom(x,a,V0,m),Tspan, x0);

xlin = Y(:,1); 
vlin = Y(:,2); 

figure; 
  plot(T,xlin,'r-');

%% Solve a more advanced problem 
% V(x) = 0.3774*x^4-2.9802*x^2+1.25
% a) y0 = [a/3;0.1] 
% b) y0 = [a/3;-2]
% clearvars; % uncomment to get a fresh run
m=1; V0 = 5; a = 2;
fx = @(x) -0.3774*3*x.^3+2*2.9802*x;
eom2 = @(x) [x(2);fx(x(1))/m];

Tspan = [0 10];
% comment out x0= .. lines to get same intial state
% x0 = [a/3; -2];
x0 = [5*a/4; -4];
[T,Y] = ode45(@(t,x) eom2(x), Tspan, x0); 

x = Y(:,1);
v = Y(:,2);

figure;
  plot(T,x,'r-');
  xlabel('t'); ylabel('x');

%% Advanced visualization
figure;
  subplot(2,2,2); plot(T,x,'r-'); 
      xlabel('t'); ylabel('x');title('x(t)');
  subplot(2,2,4); plot(T,v,'b-'); 
      xlabel('t'); ylabel('v');title('v(t)');
  subplot(2,2,[1,3]); plot(x,v,'k-'); 
      xlabel('x'); ylabel('v');title('x(v)');


