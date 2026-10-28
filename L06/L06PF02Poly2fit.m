% Example of least squares fitting of a line
% Lecture 04. Example 3. Least Squares fitting
% 
clearvars
xi =transpose([-1 0 2.1 2.3 2.4 5.3 6 6.5 8]);
yi = transpose([-1.02 -0.52 0.55 0.7 0.7 2.13 2.52 2.82 3.54]);
figure;
plot(xi,yi,'ro','MarkerFaceColor','r');
%%
A = [ones(size(xi)) xi]
Anorm = transpose(A)*A
Bnorm = transpose(A)*yi
cpar = Anorm \ Bnorm
%%
x = linspace(min(xi), max(xi), 50);
y = cpar(1)+cpar(2)*x;
figure;
plot(xi,yi,'ro','MarkerFaceColor','r');
hold on;
ax = gca;
plot(x,y,'k--');

ax.XGrid = 'on'; ax.YGrid = 'off';
ax.XMinorGrid = 'on';
ax.YMinorGrid = 'on';
hl = line([x(1) x(end)],[0 0]);
hl.Color = 'k';hl.LineWidth = 1;
vl = line([0 0],[-2 4]);
vl.Color = 'k'; vl.LineWidth = 1;
xlabel('x');
ylabel('y');
title('Data points and fitted line');
legend('Data Points','Fitted Line','Location','north');
