% Simple Optimization
% Find minima of f(x)
% f(x) = x
f1 = @(x) (x-2).*(x-1).*(x+1).*(x-3);
x = -2:0.01:5;
x01 = 2;
opts = optimset('display','iter');
[x1,fval] = fminunc( f1, x01,opts);
figure; 
  plot(x, f1(x), 'k-', x01, f1(x01), 'bo',x1,f1(x1), 'ro');
  xlabel('x'); ylabel('y'); 
  title(sprintf('f1 - initial condition x0=%5.3f',x01));
fprintf('Solution : %6.3f | error = %6.3f\n',x1, abs(fval-f1(x1)) )
%% Use an other initial value
opts = optimset('display','iter');
x02 = 0;
[x2, fval] = fminunc(f1, x02, opts);
figure;
  plot(x,f1(x), 'k-', x02, f1(x02), 'bo', x2, f1(x2), 'ro');
  xlabel('x'); ylabel('y'); 
  title(sprintf('f1 - initial condition x=%5.3f',x02));

