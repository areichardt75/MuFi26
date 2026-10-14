%  In a region E and B is defined 
% E(x,y) = ex*Ex + ey*Ey
% B(x,y) = ez*B
% x(1)-->x, x(2)-->vx, x(3)-->y, x(4)-->vy
E0 = 0.1; alfa = 0;%alfa = pi/6;
qperm = 1;
Ex = E0*cos(alfa); Ey = E0*sin(alfa);
B = 2;
% force = @(v) qperm*[Ex+B*v(2);Ey-B*v(1)];
% eom = @(x) f = force([x(2),x(4)]); [x(2), f(1) ,x(4), f(2)];

Tspan = [0 100];
% x0 = [2 1 0 0];
x0 = [2 0 0 0];
[T,Y] = ode45( @(t,x) eom(x,Ex, Ey,B,qperm), Tspan, x0 );

x = Y(:,1); y = Y(:,3);
vx= Y(:,2); vy= Y(:,4);

%
% figure; plot(T,x,'r-');
figure; 
  subplot(221); plot(T,x,'r-o'); xlabel('t'); ylabel('x');
  subplot(222); plot(T,y,'r-o'); xlabel('t'); ylabel('y');
  subplot(223); plot(x,y,'r-o'); xlabel('x'); ylabel('y');
  subplot(224); plot(vx,vy,'r-o'); xlabel('vx'); ylabel('vy');
  



% provate functions of the script
function [dydt] = eom(x,Ex,Ey,B,qperm)
  dydt = zeros(size(x));
  v = [x(2) x(4)];
  f = qperm*[Ex+B*v(2);Ey-B*v(1)];
  % force(vv);
  dydt(1) = x(2);
  dydt(2) = f(1);
  dydt(3) = x(4);
  dydt(4) = f(2);
end
