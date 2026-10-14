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
% Tspan = [0 2];
opts = odeset('AbsTol',1e-2, 'Events',@eomevents);
[T,Y,te,ye,ie] = ode89(@(t,x) eom(x,9.81), Tspan, xinit, opts);
% Extract some information from 
x = Y(:,1);   y = Y(:,3); 
vx = Y(:,2);  vy=Y(:,4);

if ~isempty(ie)
  fprintf('Event at time : %8.5f !\n',te)
end
% 
%% Visualization
figure; plot(x,y,'r-');
  xlabel('x'); ylabel('y'); 
  


function [value,isterminal,direction] = eomevents(t,y)
% checks if the mass enters ground or not 
value = y(3);     % detect height = 0 (y(3)==0 )
isterminal = 1;   % stop the integration
direction = -1;   % negative direction

end

