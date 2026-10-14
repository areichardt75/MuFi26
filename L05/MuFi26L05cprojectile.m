% Motion of a projectile - 
% 
% A projectile starts at x0,y0 with vx0, vy0 
% Calculate the trajectory of the projectile

x0 = 0;
y0 = 0;
v0x = 10;
v0y = 20;
xinit = [x0;v0x;y0;v0y];
Tinit = [0 10];

eom = @(x,g) [x(2);0;x(4);-g];

Tspan = [0 10];
opts = odeset('AbsTol',1e-2);
[T,Y] = ode45(@(t,x) eom(x,9.81), Tspan, xinit, opts);

x = Y(:,1); y = Y(:,3); 
vx = Y(:,2); vy=Y(:,4);

%% Visualization
figure; plot(x,y,'r-');
  xlabel('x'); ylabel('y'); 


