% Find the point where the sin(x) function takes its minimum 
% in range 0<x<2*pi
% 

func = @sin;
xLO = 0;
xHI = 2*pi;
x = fminbnd(func, xLO, xHI)


%% Minimize with extra parameter
% 
a = 9/7;
fun = @(x) sin(x-a);
x = fminbnd(fun, 1, 2*pi)
