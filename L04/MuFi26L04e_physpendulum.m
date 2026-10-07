%% Mufi 2025 - Lecture 8. - Solution of ODEs
% Physical pendulum 

fi0 = [pi/3;0];
ell=1;
% fi0 = [0;8+2];
Tspan = 0:0.01:2;
[~, x] = ode45(@(t,x) penfun(t, x, ell), Tspan, fi0);
fi = x(:,1); om = x(:,2);
[T,xlin] = ode45(@(t,x) penfunlin(t, x, ell), Tspan, fi0);
%%
figure;
  subplot(211);  plot(T,fi, 'r-'); xlabel('t'); ylabel('fi');
  subplot(212); plot(T, om, 'k-'); xlabel('t'); ylabel('om');
% figure; 
%   plot(T, fi, 'r-', T, xlin(:,1), 'k-' ); 
%     xlabel('t'); ylabel('fi');
%     legend('full','linearized','Location','best');

%% abrazolas
xpM = 2*ell/3*sin(fi); ypM = -2*ell/3*cos(fi);
xpm = -ell/3*sin(fi); ypm = ell/3*cos(fi);
minx = min([xpM;xpm]); maxx = max([xpM;xpm]);
miny = min([ypm;ypM]); maxy = max([ypm;ypM]);

%%
figure; 
  subplot(2,2,1); plot(T,fi,'k-');
    hold on;
    plot(T(1), fi(1), 'ro','MarkerSize',8,'MarkerFaceColor','r');
  subplot(2,2,3); plot(T,om,'k-');
    hold on;
    plot(T(1), om(1), 'ro','MarkerSize',8,'MarkerFaceColor','r');
  subplot(2,2,[2,4]); 
  plot(xpm(1), ypm(1), 'ko','MarkerFaceColor','k','MarkerSize',8);
  hold on;
  plot(xpM(1), ypM(1), 'bo','MarkerFaceColor','b','MarkerSize',8);
  line([xpm(1) xpM(1)],[ypm(1) ypM(1)],'Color','k');
  plot(0,0,'bo','MarkerSize',2);
  hold off;
  axis([minx maxx miny maxy]);
  axis equal;
  pause(1);

  for id=2:length(fi)
    subplot(2,2,1); plot(T,fi,'k-'); hold on;
    plot(T(id), fi(id), 'ro','MarkerSize',4,'MarkerFaceColor','r');
    hold off; title('fi'); xlabel('t'); ylabel('fi');
    subplot(2,2,3); plot(T,om,'k-'); hold on;
    plot(T(id), om(id), 'ro','MarkerSize',4,'MarkerFaceColor','r');
    hold off; title('om'); xlabel('t'); ylabel('om');
    subplot(2,2,[2,4]); 
    plot(xpM(id), ypM(id), 'bo','MarkerFaceColor','b','MarkerSize',8);
    hold on;
    plot(xpm(id), ypm(id), 'ko','MarkerFaceColor','k','MarkerSize',8);
    line([xpm(id) xpM(id)],[ypm(id) ypM(id)],'Color','k');
    plot(0,0,'bo','MarkerSize',2);
    hold off;
    axis equal;
    axis([minx maxx miny maxy]);
    pause(0.2);
  end


function [fiv] = penfun(t,fi,ell)
%PENDULUMFUN Differentia
%   fi'' =  -7*g/(3*1.76*ell)*cos fi 
% 
n=4;
g = 9.81;
fiv = zeros(size(fi));
fiv(1) = fi(2);
fiv(2) = -g/ell*(2*n-1)/1.88*sin(fi(1));
end

%% Functions used in script
function [fiv] = penfunlin(t,fi,ell)
%PENDULUMFUN Differentia
%   fi'' =  -7*g/(3*1.76*ell)*cos fi 
% 
n=4;
g = 9.81;
fiv = zeros(size(fi));
fiv(1) = fi(2);
fiv(2) = -g/ell*(2*n-1)/1.88*(fi(1));
end