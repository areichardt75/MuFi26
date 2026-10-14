%% Solution of van der Pol equation
% 
% y1'= y2
% y2'= 1000*(1-y1^2)y2-y1
% 
Tspan = [0 5000];
% Anonymous definition of the DE
vdP = @(x) [x(2);1000*(1-x(1)^2)*x(2)-x(1)];

%% Optional / solve using ode45 
% tic;
% [T45, Y45] = ode45(@(t,x) vdP(x), [0 1000],[2 0]);
% trun = toc; fprintf('Running time (ode45) : %8.4f\n',trun);
% figure; plot(T45, Y45(:,1), 'b-');
%   xlabel('t'); ylabel('y1, y2'); 

%% Solve using ode23t 
% Special solver for stiff problems
% 
tic; 
[T,Y] = ode23t(@(t,x) vdP(x), Tspan,[2 0]);
trun = toc; fprintf('Running time (ode23t) : %8.4f\n',trun);
figure; 
  plot(T,Y(:,1),'k-d','LineWidth',2,...
    'MarkerFaceColor','r',...
    'MarkerEdgeColor','r');
  xlabel('t'); ylabel('y1, y2'); 





