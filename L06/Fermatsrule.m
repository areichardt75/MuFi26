% Lecture 5. Optimization problem
% 
x0 = -10; y0 = 3;
x2 = 10; y2 = -3;
% level of interface
y1 = 0;

c1 = 1; c2 = 2;
t1sq = @(x,c) ((x-x0).^2+(y1-y0).^2)/(c^2);
t2sq = @(x,c) ((x2-x).^2+(y2-y1).^2)/(c^2);

fun = @(x,c1,c2) sqrt(t1sq(x,c1))+ sqrt(t2sq(x, c2));
opts = optimset('display','iter');

[x1] = fminunc(@(x) fun(x, c1, c2), 5, opts)
len1 = sqrt(t1sq(x1,c1)*(c1^2));
len2 = sqrt(t2sq(x1,c2)*(c2^2));
fprintf('Time1^2 : %.5f | time 2^2 : %.5f\n',t1sq(x1,c1),t2sq(x1,c2));
fprintf('Length1 : %.5f | Length2 : %.5f| Total length : %.5f\n',...
  len1, len2,len1+len2)
%%
figure; 
  plot([x0,x1,x2],[y0,y1,y2],'ro','MarkerSize',8,'MarkerFaceColor','r');
  hold on;
  l1 = line([x0,x1],[y0,y1],'color','b','LineWidth',2);
  l2 = line([x1, x2],[y1,y2], 'Color','b','LineWidth',2);
